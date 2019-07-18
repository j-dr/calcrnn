#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>
#include <string.h>
#include <gsl/gsl_math.h>

#ifdef USEOPENMP
#include <omp.h>
#endif

#include "allheader.h"

/*
  Notes for how to do this.
  
  Queen - 
  1) read in a few hundred MB section of the halo file
  2) sort the halos so that they are ordered by which of the particle BBoxes they fall into - use the largest volume one
  3) send all halos in BBoxes from a single file to each worker - this should roughly make sure each worked does not exceed mem limits
        as long as it didn't when doing Rnn for the particles
  4) recv halo Rnn values from workers
  5) write Rnn values to disk
  
  Worker - 
  1) recv halos from queen
  2) compute halo bbox w/ buffers
  3) read all parts which intersect halo bbox 
  4) compute Rnn
  5) send Rnn back to queen
  
*/

//#define DEBUG_QWCOMM

static void do_rnn_halos_serial(BBox *boxes, int Nb);
static void do_rnn_halos_parallel(BBox *boxes, int Nb);
static void queen_code_rnn(BBox *boxes, int Nb);
static void worker_code_rnn(BBox *boxes, int Nb);
static void queen_code_bbox(BBox *boxes, int Nb);
static void worker_code_bbox(BBox *boxes, int Nb);
static long read_halo_chunk(Halo *halos, long Nh, FILE *fp, FILE *fpb);
static long sort_halo_chunk_into_bboxes(Halo *halos, long Nh, BBox *boxes, int Nb);
static void sort_halos_sortindex(Halo *halos, long Nh);
static int compHaloSortIndex(const void *a, const void *b);
static float get_mindist_bbox_point(BBox bx, float px, float py, float pz, float L);
static void reorder_halos_phcurve(Halo *halos, long Nh);
static void compute_bboxes_for_halos(Halo *halos, long Nh, BBox *boxes, int Nb);

static long read_halo_chunk(Halo *halos, long Nh, FILE *fp, FILE *fpb)
{
  long i,nr = 0;
  char fline[MAX_FILENAME];
  char flineb[MAX_FILENAME];
  long id;
  
  for(i=0;i<Nh;++i)
    {
      fgets(fline,MAX_FILENAME,fp);
      if(fpb != NULL)
	fgets(flineb,MAX_FILENAME,fpb);
      
      if(!feof(fp))
	{
	  if(strcmp(allData.HaloFileFormat,"SKELETON") == 0)
	    sscanf(fline,"%ld %*d %*g %*g %*g %*g %*g %*g %g %g %g %*g %*g %*g %*d\n",
		   &(halos[i].id),&(halos[i].px),&(halos[i].py),&(halos[i].pz));
	  else if(strcmp(allData.HaloFileFormat,"SUB") == 0)
	    {
	      // Scale Id Desc_scale Descid Num_prog Pid Upid Desc_pid Phantom Mvir Orig_Mvir Rvir Rs  Vrms Mmp Last_mm Vmax X  Y  Z  VX  VY  VZ  Macc Mpeak Vacc Vpeak
	      // %*g   %d %*g        %*g    %*g      %*g %*g  %*g      %*g     %*g  %*g       %*g  %*g %*g  %*g %*g     %*g  %g %g %g %*g %*g %*g %*g  %*g   %*g  %*g  \n
	      //            %*g %ld %*g %*g %*g %*g %*g %*g %*g %*g %*g %*g %*g %*g %*g %*g %*g %g %g %g %*g %*g %*g %*g %*g %*g %*g\n
	      sscanf(fline,"%*g %ld %*g %*g %*g %*g %*g %*g %*g %*g %*g %*g %*g %*g %*g %*g %*g %g %g %g %*g %*g %*g %*g %*g %*g %*g\n",
		     &(halos[i].id),&(halos[i].px),&(halos[i].py),&(halos[i].pz));
	    }
	  else 
	    {
	      fprintf(stderr,"Halo file format '%s' not supported!\n",allData.HaloFileFormat);
	      assert(0);
	    }
	  
	  if(fpb != NULL)
	    {
	      sscanf(flineb,"%ld %ld %ld\n\n",
		     &id,&(halos[i].bboxFileIndex),&(halos[i].bboxIndex));
	      assert(id == halos[i].id);
	    }
	  
	  ++nr;
	  
	  if(allData.BoxLength > 0.0)
	    {
	      halos[i].px = wrap_pos(halos[i].px,allData.BoxLength);
	      halos[i].py = wrap_pos(halos[i].py,allData.BoxLength);
	      halos[i].pz = wrap_pos(halos[i].pz,allData.BoxLength);
	    }
	}
      else
	break;
    }
  
  return nr;
}


void do_rnn_halos(BBox *boxes, int Nb)
{
  if(NTasks == 1)
    do_rnn_halos_serial(boxes,Nb);
  else
    do_rnn_halos_parallel(boxes,Nb);
}

static void do_rnn_halos_parallel(BBox *boxes, int Nb)
{
  if(ThisTask == 0)
    queen_code_bbox(boxes,Nb);
  else
    worker_code_bbox(boxes,Nb);
  
  ////////////////////////////
  MPI_Barrier(MPI_COMM_WORLD);
  ////////////////////////////
    
  if(ThisTask == 0)
    queen_code_rnn(boxes,Nb);
  else
    worker_code_rnn(boxes,Nb);
}

typedef struct {
  long bbind;
  long find;
  long index;
} BBoxData;

static void queen_code_bbox(BBox *boxes, int Nb)
{
  FILE *fph,*fp;
  long NhTot;
  long Nh,NhAlloc,NhRead;
  char fname[MAX_FILENAME],hline[MAX_FILENAME];
  long baseChunkSize,chunk,chunkStart,Nchunks;
  Halo *halos;
  long NumHaloWorkUnits,startWorkHalos,NworkHalos,workHaloBBoxFileIndex;
  long haloWorkUnit,i,j,k;
  long NumReadQueue,NumTasksWorking;
  int *readQueue;
  long NumDone;
  long NBBoxRecv,NBBoxDataTot = 0;
  BBoxData *bboxData = NULL,*tmpBBoxData;
  MPI_Status status;
    
  //alloc queues
  readQueue = (int*)malloc(sizeof(int)*NTasks);
  assert(readQueue != NULL);
  NumReadQueue = 0;

  //open halo file and get # of halos
  fph = fopen_retry(allData.HaloFile,"r");
  assert(fph != NULL);
  NhTot = fnumlines(fph) - 1;
  fgets(hline,MAX_FILENAME,fph); //read the comment line
  
  //get halo chunk size and mem
  NhAlloc = (allData.HaloChunkSizeMB*1024l*1024l)/sizeof(Halo);
  assert(NhAlloc > 0);
  
  halos = (Halo*)malloc(sizeof(Halo)*NhAlloc);
  assert(halos != NULL);
  
  //get # of halo chunks
  baseChunkSize = NhAlloc;
  Nchunks = NhTot/NhAlloc;
  if(NhTot != baseChunkSize*Nchunks)
    ++Nchunks;
  
  //open BBox file
  sprintf(hline,"%s",allData.HaloFile);
  extract_filename(hline);
  sprintf(fname,"%s/bbox_%s",allData.OutputPath,hline);
  fp = fopen_retry(fname,"w");
  assert(fp != NULL);
  fprintf(fp,"#ID FileInd BBoxInd\n");
  
  for(chunk=0;chunk<Nchunks;++chunk)
    {
      fprintf(stderr,"doing BBox for halo chunk %ld of %ld.\n",chunk+1,Nchunks);
      
      //////////////////////////////
      //get halo chunk
      //////////////////////////////
      
      //get parms of chunk
      Nh = baseChunkSize;
      chunkStart = baseChunkSize*chunk;
      if(chunkStart + Nh > NhTot)
	Nh = NhTot - 1 - chunkStart + 1;
      
      //read halos & assign file indexes to them
      NhRead = read_halo_chunk(halos,Nh,fph,NULL);
      assert(NhRead == Nh);
      for(i=0;i<Nh;++i)
	halos[i].haloFileIndex = chunkStart + i;
      
      //sort them on disk into work units
      NumHaloWorkUnits = NTasks-1;
      k = Nh/(NTasks-1);
      for(j=0;j<NumHaloWorkUnits;++j)
	{
	  if(j != NumHaloWorkUnits-1)
	    NhRead = j*k+k;
	  else
	    NhRead = Nh;
	  
	  for(i=j*k;i<NhRead;++i)
	    halos[i].bboxFileIndex = j;
	}
            
      //////////////////////////////
      //do rnn for each work unit
      //////////////////////////////
      
      //init for loop logic
      NumDone = 0;
      NumTasksWorking = 0;
      startWorkHalos = 0;
      haloWorkUnit = 0;
            
      while(haloWorkUnit < NumHaloWorkUnits || NumTasksWorking > 0)
	{
	  //see if anyone needs to read
	  while(NumReadQueue > 0 && haloWorkUnit < NumHaloWorkUnits)
	    {
	      //get all halos with the same work index
	      workHaloBBoxFileIndex = halos[startWorkHalos].bboxFileIndex;
	      NworkHalos = 0;
	      for(i=startWorkHalos;i<Nh;++i)
		{
		  if(halos[i].bboxFileIndex == workHaloBBoxFileIndex)
		    ++NworkHalos;
		  else
		    break;
		}
	      
	      //store current halo pos in sort index - this wasy can fill in index values when they are sent back
	      for(i=0;i<NworkHalos;++i)
		halos[i+startWorkHalos].sortIndex = i + startWorkHalos;
	      
#ifdef DEBUG_QWCOMM               
	      fprintf(stderr,"%d: sending halos to task %d (%ld tasks in read queue, %ld tasks working).\n",
		      ThisTask,readQueue[NumReadQueue-1],NumReadQueue,NumTasksWorking);
#endif

	      //send # of halos
	      MPI_Send(&NworkHalos,1,MPI_LONG,readQueue[NumReadQueue-1],TAG_WORKWORKER,MPI_COMM_WORLD);
	      
	      //send the halos
	      MPI_Send(halos+startWorkHalos,sizeof(Halo)*NworkHalos,MPI_BYTE,readQueue[NumReadQueue-1],TAG_HALO_DATA,MPI_COMM_WORLD);
	      
#ifdef DEBUG_QWCOMM               
	      fprintf(stderr,"%d: finished sending halos to task %d.\n",ThisTask,readQueue[NumReadQueue-1]);
#endif

	      --NumReadQueue;
	      ++NumTasksWorking;
	      ++haloWorkUnit;
	      startWorkHalos += NworkHalos;
	    }
              
	  //now recv messages from workers
	  MPI_Recv(&NBBoxRecv,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	  
	  if(status.MPI_TAG == TAG_READ_REQ)
	    {
	      readQueue[NumReadQueue] = status.MPI_SOURCE;
	      ++NumReadQueue;

#ifdef DEBUG_QWCOMM               
	      fprintf(stderr,"%d: added task %d to read queue (%ld tasks in read queue, %ld tasks working).\n",
		      ThisTask,status.MPI_SOURCE,NumReadQueue,NumTasksWorking);
#endif
	    }
	  else if(status.MPI_TAG == TAG_RNN_DONE)
	    {
	      //make sure have enough temp BBox memory
	      if(NBBoxRecv > NBBoxDataTot)
		{
		  NBBoxDataTot = NBBoxRecv + NBBoxRecv*1.2;
		  tmpBBoxData = (BBoxData*)realloc(bboxData,sizeof(BBoxData)*NBBoxDataTot);
		  assert(tmpBBoxData != NULL);
		  bboxData = tmpBBoxData;
		}
	      
	      //recv halo rnns from worker
	      MPI_Recv(bboxData,sizeof(BBoxData)*NBBoxRecv,MPI_BYTE,status.MPI_SOURCE,TAG_RNN_DATA,MPI_COMM_WORLD,&status);
	      
	      //put in global buffer of BBox values
	      for(i=0;i<NBBoxRecv;++i)
		{
		  halos[bboxData[i].index].bboxFileIndex = bboxData[i].find;
		  halos[bboxData[i].index].bboxIndex = bboxData[i].bbind;
		}
	      
	      --NumTasksWorking;
	      ++NumDone;
	      
#ifdef DEBUG_QWCOMM               
	      fprintf(stderr,"Halo BBox: finished %ld of %ld (%.2f percent) of halo work units.\n",
		      NumDone,NumHaloWorkUnits,((double) NumDone)/((double) NumHaloWorkUnits)*100.0);
#endif
	    }
	  else
	    {
	      fprintf(stderr,"recv weird tag at root! tag,src = %d|%d\n",status.MPI_TAG,status.MPI_SOURCE);
	      MPI_Abort(MPI_COMM_WORLD,123);
	    }
	}
          
      assert(NumDone == NumHaloWorkUnits);
      assert(haloWorkUnit == NumHaloWorkUnits);
          
      //make sure to get last of read requests
      while(NumReadQueue < NTasks-1)
	{
	  MPI_Recv(&NBBoxRecv,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
              
	  if(status.MPI_TAG == TAG_READ_REQ)
	    {
	      readQueue[NumReadQueue] = status.MPI_SOURCE;
	      ++NumReadQueue;

#ifdef DEBUG_QWCOMM               
	      fprintf(stderr,"%d: added task %d to read queue (%ld tasks in read queue, %ld tasks working).\n",
		      ThisTask,status.MPI_SOURCE,NumReadQueue,NumTasksWorking);
#endif
	    }
	  else
	    {
	      fprintf(stderr,"recv weird tag at root while cleaning up read requests! tag,src = %d|%d\n",status.MPI_TAG,status.MPI_SOURCE);
	      MPI_Abort(MPI_COMM_WORLD,123);
	    }
	}
          
      //////////////////////////////
      //write BBox to disk
      //////////////////////////////
      
      //write to disk
      for(i=0;i<Nh;++i)
        fprintf(fp,"%ld %ld %ld\n",halos[i].id,halos[i].bboxFileIndex,halos[i].bboxIndex);
      fflush(fp);
      
    }//for(chunk=0;chunk<Nchunks;++chunk)
  
  //now kill the worker tasks
  for(i=1;i<NTasks;++i)
    MPI_Send(&NworkHalos,1,MPI_LONG,i,TAG_KILLWORKER,MPI_COMM_WORLD);
    
  //clean up
  free(halos);
  fclose(fp);
  fclose(fph);
  free(readQueue);
  if(NBBoxDataTot > 0)
    free(bboxData);
}

static void worker_code_bbox(BBox *boxes, int Nb)
{
  MPI_Status status;
  long NworkHalosCurr = 0;
  long NworkHalos;
  long NBBoxSend,i;
  Halo *halos = NULL,*tmpHalos;
  long NhBuffMax = 0;
  BBoxData *bboxData = NULL,*tmpBBoxData;
  BBox *tmp_bboxes;

  //alloc boxes to contain modified halo boxes
  tmp_bboxes = (BBox*)malloc(sizeof(BBox)*allData.Ndiv*allData.Ndiv*allData.Ndiv);
  assert(tmp_bboxes != NULL);

  while(1)
    {
      if(NworkHalosCurr == 0)
	MPI_Send(&NBBoxSend,1,MPI_LONG,0,TAG_READ_REQ,MPI_COMM_WORLD);
              
      MPI_Recv(&NworkHalos,1,MPI_LONG,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
              
      if(status.MPI_TAG == TAG_KILLWORKER)
	break;
      else if(status.MPI_TAG == TAG_WORKWORKER)
	{

#ifdef DEBUG_QWCOMM               
	  fprintf(stderr,"%d: going to recv %ld halos.\n",ThisTask,NworkHalos);
#endif

	  NworkHalosCurr = NworkHalos;
	  NBBoxSend = NworkHalos;
	  
	  //make buffers for halos, halo positions and BBox data
	  if(NworkHalosCurr > NhBuffMax)
	    {
	      NhBuffMax = NworkHalosCurr;
	      	      
	      tmpBBoxData = (BBoxData*)realloc(bboxData,sizeof(BBoxData)*NhBuffMax);
	      assert(tmpBBoxData != NULL);
	      bboxData = tmpBBoxData;
	      
	      tmpHalos = (Halo*)realloc(halos,sizeof(Halo)*NhBuffMax);
	      assert(tmpHalos != NULL);
	      halos = tmpHalos;
	    }
	  
	  //recv halos
	  MPI_Recv(halos,sizeof(Halo)*NworkHalosCurr,MPI_BYTE,0,TAG_HALO_DATA,MPI_COMM_WORLD,&status);
	  
	  //fill halo pos & rnn buffers
	  for(i=0;i<NworkHalosCurr;++i)
	    bboxData[i].index = halos[i].sortIndex;
	  
#ifdef DEBUG_QWCOMM               
	  fprintf(stderr,"%d: read %ld halos for file %ld, going to read parts for halos.\n",ThisTask,NworkHalosCurr,halos[0].bboxFileIndex);
#endif
#ifdef DEBUG_QWCOMM               
	  fprintf(stderr,"%d: read parts for halos, doing bbox for halos.\n",ThisTask);
#endif
	  
	  //do bbox calc
	  compute_bboxes_for_halos(halos,NworkHalosCurr,boxes,Nb);
	    
	  //fill BBoxData to send back
	  for(i=0;i<NworkHalosCurr;++i)
	    {
	      bboxData[i].find = halos[i].bboxFileIndex;
	      bboxData[i].bbind = halos[i].bboxIndex;
	    }
	  
#ifdef DEBUG_QWCOMM               
	  fprintf(stderr,"%d: did bbox for halos, send bbox for halos.\n",ThisTask);
#endif
	  
	  //tell root done with BBox
	  MPI_Send(&NBBoxSend,1,MPI_LONG,0,TAG_RNN_DONE,MPI_COMM_WORLD); 
	  
	  //send root the data
	  MPI_Send(bboxData,sizeof(BBoxData)*NworkHalosCurr,MPI_BYTE,0,TAG_RNN_DATA,MPI_COMM_WORLD);
	  
	  NworkHalosCurr = 0;
	}
      else
	{
	  fprintf(stderr,"worker %d recv weird tag! tag,src = %d|%d\n",ThisTask,status.MPI_TAG,status.MPI_SOURCE);
	  MPI_Abort(MPI_COMM_WORLD,999);
	}
      
    } //end of while loop
  
  //free rnn data and halo pos buffers
  free(tmp_bboxes);
  if(NhBuffMax > 0)
    {
      free(halos);
      free(bboxData);
    }
}

typedef struct {
  float rnn;
  long index;
} RnnData;

static void queen_code_rnn(BBox *boxes, int Nb)
{
  FILE *fph,*fpb,*fp;
  long NhTot;
  long Nh,NhAlloc,NhRead;
  char fname[MAX_FILENAME],hline[MAX_FILENAME];
  long baseChunkSize,chunk,chunkStart,Nchunks;
  Halo *halos;
  long NumHaloWorkUnits,startWorkHalos,NworkHalos,workHaloBBoxFileIndex;
  long haloWorkUnit,i;
  long NumTasksDoingIO,NumReadQueue,NumTasksWorking;
  int *readQueue;
  long NumDone;
  long NRnnRecv,NRnnDataTot = 0;
  RnnData *rnnData = NULL,*tmpRnnData;
  MPI_Status status;
    
  //alloc queues
  readQueue = (int*)malloc(sizeof(int)*NTasks);
  assert(readQueue != NULL);
  NumReadQueue = 0;

  //open halo file and get # of halos
  fph = fopen_retry(allData.HaloFile,"r");
  assert(fph != NULL);
  NhTot = fnumlines(fph) - 1;
  fgets(hline,MAX_FILENAME,fph); //read the comment line
  
  sprintf(hline,"%s",allData.HaloFile);
  extract_filename(hline);
  sprintf(fname,"%s/bbox_%s",allData.OutputPath,hline);
  fpb = fopen_retry(fname,"r");
  assert(fpb != NULL);
  fgets(hline,MAX_FILENAME,fpb); //read the comment line
  
  //get halo chunk size and mem
  NhAlloc = (allData.HaloChunkSizeMB*1024l*1024l)/sizeof(Halo);
  assert(NhAlloc > 0);
  
  halos = (Halo*)malloc(sizeof(Halo)*NhAlloc);
  assert(halos != NULL);
  
  //get # of halo chunks
  baseChunkSize = NhAlloc;
  Nchunks = NhTot/NhAlloc;
  if(NhTot != baseChunkSize*Nchunks)
    ++Nchunks;
  
  //open Rnn file
  sprintf(hline,"%s",allData.HaloFile);
  extract_filename(hline);
  sprintf(fname,"%s/rnn_%s",allData.OutputPath,hline);
  fp = fopen_retry(fname,"w");
  assert(fp != NULL);
  fprintf(fp,"#ID Rnn [for N = %d]\n",allData.NRnn);
  
  //mak the I/O nice
  fprintf(stderr,"\n");
  
  for(chunk=0;chunk<Nchunks;++chunk)
    {
      fprintf(stderr,"doing Rnn for halo chunk %ld of %ld.\n",chunk+1,Nchunks);
      
      //////////////////////////////
      //get halo chunk
      //////////////////////////////
      
      //get parms of chunk
      Nh = baseChunkSize;
      chunkStart = baseChunkSize*chunk;
      if(chunkStart + Nh > NhTot)
	Nh = NhTot - 1 - chunkStart + 1;
      
      //read halos & assign file indexes to them
      NhRead = read_halo_chunk(halos,Nh,fph,fpb);
      assert(NhRead == Nh);
      for(i=0;i<Nh;++i)
	halos[i].haloFileIndex = chunkStart + i;
      
      //sort them on disk into work units
      NumHaloWorkUnits = sort_halo_chunk_into_bboxes(halos,Nh,boxes,Nb);
      
      fprintf(stderr,"finished sorting halo chunk %ld of %ld for Rnn.\n",chunk+1,Nchunks);
      
      //////////////////////////////
      //do rnn for each work unit
      //////////////////////////////
      
      //init for loop logic
      NumDone = 0;
      NumTasksDoingIO = 0;
      NumTasksWorking = 0;
      startWorkHalos = 0;
      haloWorkUnit = 0;
            
      while(haloWorkUnit < NumHaloWorkUnits || NumTasksWorking > 0)
	{
	  //see if anyone needs to read
	  while(NumReadQueue > 0 && NumTasksDoingIO < allData.NumTasksIOInParallel && haloWorkUnit < NumHaloWorkUnits)
	    {
	      //get all halos with the same work index
	      workHaloBBoxFileIndex = halos[startWorkHalos].bboxFileIndex;
	      NworkHalos = 0;
	      for(i=startWorkHalos;i<Nh;++i)
		{
		  if(halos[i].bboxFileIndex == workHaloBBoxFileIndex)
		    ++NworkHalos;
		  else
		    break;
		}
	      
	      //store current halo pos in sort index - this wasy can fill in Rnn values when they are sent back
	      for(i=0;i<NworkHalos;++i)
		halos[i+startWorkHalos].sortIndex = i + startWorkHalos;
	      
#ifdef DEBUG_QWCOMM               
	      fprintf(stderr,"%d: sending halos to task %d (%ld tasks in read queue, %ld tasks working).\n",
		      ThisTask,readQueue[NumReadQueue-1],NumReadQueue,NumTasksWorking);
#endif

	      //send # of halos
	      MPI_Send(&NworkHalos,1,MPI_LONG,readQueue[NumReadQueue-1],TAG_WORKWORKER,MPI_COMM_WORLD);
	      
	      //send the halos
	      MPI_Send(halos+startWorkHalos,sizeof(Halo)*NworkHalos,MPI_BYTE,readQueue[NumReadQueue-1],TAG_HALO_DATA,MPI_COMM_WORLD);
	      
#ifdef DEBUG_QWCOMM               
	      fprintf(stderr,"%d: finished sending halos to task %d.\n",ThisTask,readQueue[NumReadQueue-1]);
#endif

	      --NumReadQueue;
	      ++NumTasksDoingIO;
	      ++NumTasksWorking;
	      ++haloWorkUnit;
	      startWorkHalos += NworkHalos;
	    }
              
	  //now recv messages from workers
	  MPI_Recv(&NRnnRecv,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	  
	  if(status.MPI_TAG == TAG_READ_REQ)
	    {
	      readQueue[NumReadQueue] = status.MPI_SOURCE;
	      ++NumReadQueue;

#ifdef DEBUG_QWCOMM               
	      fprintf(stderr,"%d: added task %d to read queue (%ld tasks in read queue, %ld tasks working).\n",
		      ThisTask,status.MPI_SOURCE,NumReadQueue,NumTasksWorking);
#endif
	    }
	  else if(status.MPI_TAG == TAG_READ_DONE)
	    {
	      --NumTasksDoingIO;
	    }
	  else if(status.MPI_TAG == TAG_RNN_DONE)
	    {
	      //make sure have enough temp Rnn memory
	      if(NRnnRecv > NRnnDataTot)
		{
		  NRnnDataTot = NRnnRecv + NRnnRecv*1.2;
		  tmpRnnData = (RnnData*)realloc(rnnData,sizeof(RnnData)*NRnnDataTot);
		  assert(tmpRnnData != NULL);
		  rnnData = tmpRnnData;
		}
	      
	      //recv halo rnns from worker
	      MPI_Recv(rnnData,sizeof(RnnData)*NRnnRecv,MPI_BYTE,status.MPI_SOURCE,TAG_RNN_DATA,MPI_COMM_WORLD,&status);
	      
	      //put in global buffer of Rnn values
	      for(i=0;i<NRnnRecv;++i)
		halos[rnnData[i].index].rnn = rnnData[i].rnn;
	      
	      --NumTasksWorking;
	      ++NumDone;
	      
	      fprintf(stderr,"Rnn: finished %ld of %ld (%.2f percent) of halo work units.\n",
		      NumDone,NumHaloWorkUnits,((double) NumDone)/((double) NumHaloWorkUnits)*100.0);
	    }
	  else
	    {
	      fprintf(stderr,"recv weird tag at root! tag,src = %d|%d\n",status.MPI_TAG,status.MPI_SOURCE);
	      MPI_Abort(MPI_COMM_WORLD,123);
	    }
	}
          
      assert(NumDone == NumHaloWorkUnits);
      assert(haloWorkUnit == NumHaloWorkUnits);
          
      //make sure to get last of read requests
      while(NumReadQueue < NTasks-1)
	{
	  MPI_Recv(&NRnnRecv,1,MPI_LONG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
              
	  if(status.MPI_TAG == TAG_READ_REQ)
	    {
	      readQueue[NumReadQueue] = status.MPI_SOURCE;
	      ++NumReadQueue;

#ifdef DEBUG_QWCOMM               
	      fprintf(stderr,"%d: added task %d to read queue (%ld tasks in read queue, %ld tasks working).\n",
		      ThisTask,status.MPI_SOURCE,NumReadQueue,NumTasksWorking);
#endif
	    }
	  else
	    {
	      fprintf(stderr,"recv weird tag at root while cleaning up read requests! tag,src = %d|%d\n",status.MPI_TAG,status.MPI_SOURCE);
	      MPI_Abort(MPI_COMM_WORLD,123);
	    }
	}
          
      //////////////////////////////
      //write Rnn to disk
      //////////////////////////////
      
      //resort into file index order
      for(i=0;i<Nh;++i)
        halos[i].sortIndex = halos[i].haloFileIndex;
      sort_halos_sortindex(halos,Nh);

      //write to disk
      for(i=0;i<Nh;++i)
        fprintf(fp,"%ld %f\n",halos[i].id,halos[i].rnn);
      fflush(fp);
      
    }//for(chunk=0;chunk<Nchunks;++chunk)
  
  //now kill the worker tasks
  for(i=1;i<NTasks;++i)
    MPI_Send(&NworkHalos,1,MPI_LONG,i,TAG_KILLWORKER,MPI_COMM_WORLD);
    
  //clean up
  free(halos);
  fclose(fp);
  fclose(fph);
  fclose(fpb);
  free(readQueue);
  if(NRnnDataTot > 0)
    free(rnnData);
}

static void worker_code_rnn(BBox *boxes, int Nb)
{
  MPI_Status status;
  MPI_Request request;
  long NworkHalosCurr = 0;
  long NworkHalos;
  long NRnnSend,i,k;
  int fnum;
  float *fptrs[3];
  Halo *halos = NULL,*tmpHalos;
  float *px,*py,*pz,*rnn,*tmp;
  int NpFile,NpTot,Nrnn;
  long NhBuffMax = 0;
  RnnData *rnnData = NULL,*tmpRnnData;
  fptrs[0] = NULL;
  fptrs[1] = NULL;
  fptrs[2] = NULL;
  BBox *tmp_bboxes;

  //alloc boxes to contain modified halo boxes
  tmp_bboxes = (BBox*)malloc(sizeof(BBox)*allData.Ndiv*allData.Ndiv*allData.Ndiv);
  assert(tmp_bboxes != NULL);

  while(1)
    {
      if(NworkHalosCurr == 0)
	MPI_Send(&NRnnSend,1,MPI_LONG,0,TAG_READ_REQ,MPI_COMM_WORLD);
              
      MPI_Recv(&NworkHalos,1,MPI_LONG,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
              
      if(status.MPI_TAG == TAG_KILLWORKER)
	break;
      else if(status.MPI_TAG == TAG_WORKWORKER)
	{

#ifdef DEBUG_QWCOMM               
	  fprintf(stderr,"%d: going to recv %ld halos.\n",ThisTask,NworkHalos);
#endif

	  NworkHalosCurr = NworkHalos;
	  NRnnSend = NworkHalos;
	  
	  //make buffers for halos, halo positions and Rnn data
	  if(NworkHalosCurr > NhBuffMax)
	    {
	      NhBuffMax = NworkHalosCurr;
	      for(i=0;i<3;++i)
		{
		  tmp = (float*)realloc(fptrs[i],sizeof(float)*NhBuffMax);
		  assert(tmp != NULL);
		  fptrs[i] = tmp;
		}
	      
	      tmpRnnData = (RnnData*)realloc(rnnData,sizeof(RnnData)*NhBuffMax);
	      assert(tmpRnnData != NULL);
	      rnnData = tmpRnnData;
	      
	      tmpHalos = (Halo*)realloc(halos,sizeof(Halo)*NhBuffMax);
	      assert(tmpHalos != NULL);
	      halos = tmpHalos;
	    }
	  
	  //recv halos
	  MPI_Recv(halos,sizeof(Halo)*NworkHalosCurr,MPI_BYTE,0,TAG_HALO_DATA,MPI_COMM_WORLD,&status);
	  
	  //fill halo pos & rnn buffers
	  for(i=0;i<NworkHalosCurr;++i)
	    {
	      fptrs[0][i] = halos[i].px;
	      fptrs[1][i] = halos[i].py;
	      fptrs[2][i] = halos[i].pz;
	      rnnData[i].index = halos[i].sortIndex;
	    }
	  
	  //save old boxes
          fnum = halos[0].bboxFileIndex;
	  k = 0;
          for(i=0;i<Nb;++i)
            if(boxes[i].fnum == fnum)
              {
                tmp_bboxes[k] = boxes[i];
                ++k;
              }
          
          //resize boxes to contain all of the halos
          for(i=0;i<NworkHalosCurr;++i)
            {
              k = halos[i].bboxIndex;
              
              if(boxes[k].bmin[0] >= halos[i].px)
                boxes[k].bmin[0] = halos[i].px;
              if(boxes[k].bmax[0] <= halos[i].px)
                boxes[k].bmax[0] = halos[i].px;
              
              if(boxes[k].bmin[1] >= halos[i].py)
                boxes[k].bmin[1] = halos[i].py;
              if(boxes[k].bmax[1] <= halos[i].py)
                boxes[k].bmax[1] = halos[i].py;
              
              if(boxes[k].bmin[2] >= halos[i].pz)
                boxes[k].bmin[2] = halos[i].pz;
              if(boxes[k].bmax[2] <= halos[i].pz)
                boxes[k].bmax[2] = halos[i].pz;
            }
	  
#ifdef DEBUG_QWCOMM               
	  fprintf(stderr,"%d: read %ld halos for file %ld, going to read parts for halos.\n",ThisTask,NworkHalosCurr,halos[0].bboxFileIndex);
#endif
	  //read parts
	  readparts(fnum,&px,&py,&pz,&NpFile,&NpTot,boxes,Nb);
          
	  //reset boxes to old versions
          k = 0;
          for(i=0;i<Nb;++i)
            if(boxes[i].fnum == fnum)
              {
                boxes[i] = tmp_bboxes[k];
                ++k;
              }
	  
	  //tell root done reading
	  MPI_Isend(&NRnnSend,1,MPI_LONG,0,TAG_READ_DONE,MPI_COMM_WORLD,&request);
                  
#ifdef DEBUG_QWCOMM               
	  fprintf(stderr,"%d: read parts for halos, doing rnn halos for file %d.\n",ThisTask,fnum);
#endif
	  
	  //do rnn calc
	  rnncalc_halos(fptrs[0],fptrs[1],fptrs[2],NworkHalosCurr,px,py,pz,NpTot,&rnn,&Nrnn);
          assert(Nrnn == NworkHalosCurr);
	  free(px);
	  free(py);
	  free(pz);
          
	  //fill RnnData to send back
	  for(i=0;i<NworkHalosCurr;++i)
	    rnnData[i].rnn = rnn[i];
	  free(rnn);
	  
	  //make sure root knows
	  MPI_Wait(&request,&status);
          
#ifdef DEBUG_QWCOMM               
	  fprintf(stderr,"%d: did rnn for halos, send rnn for file %d\n",ThisTask,fnum);
#endif
	  
	  //tell root done with Rnn
	  MPI_Send(&NRnnSend,1,MPI_LONG,0,TAG_RNN_DONE,MPI_COMM_WORLD); 
	  
	  //send root the data
	  MPI_Send(rnnData,sizeof(RnnData)*NworkHalosCurr,MPI_BYTE,0,TAG_RNN_DATA,MPI_COMM_WORLD);
	  
	  NworkHalosCurr = 0;
	}
      else
	{
	  fprintf(stderr,"worker %d recv weird tag! tag,src = %d|%d\n",ThisTask,status.MPI_TAG,status.MPI_SOURCE);
	  MPI_Abort(MPI_COMM_WORLD,999);
	}
      
    } //end of while loop
  
  //free rnn data and halo pos buffers
  free(tmp_bboxes);
  if(NhBuffMax > 0)
    {
      free(halos);
      free(rnnData);
      for(i=0;i<3;++i)
	free(fptrs[i]);
    }
}

static void do_rnn_halos_serial(BBox *boxes, int Nb)
{
  FILE *fph,*fp;
  float *px,*py,*pz,*rnn,*tmp,*fptrs[3];
  long NhBuffMax = 0;
  long NhTot;
  long Nh,NhAlloc,NhRead;
  char fname[MAX_FILENAME],hline[MAX_FILENAME];
  long baseChunkSize,chunk,chunkStart,Nchunks;
  Halo *halos;
  long NumHaloWorkUnits,startWorkHalos,NworkHalos,workHaloBBoxFileIndex;
  long haloWorkUnit,i,k;
  int NpFile,NpTot,Nrnn;
  fptrs[0] = NULL;
  fptrs[1] = NULL;
  fptrs[2] = NULL;
  BBox *tmp_bboxes;
  
  //alloc boxes to contain modified halo boxes
  tmp_bboxes = (BBox*)malloc(sizeof(BBox)*allData.Ndiv*allData.Ndiv*allData.Ndiv);
  assert(tmp_bboxes != NULL);
  
  //open halo file and get # of halos
  fph = fopen_retry(allData.HaloFile,"r");
  assert(fph != NULL);
  NhTot = fnumlines(fph) - 1;
  fgets(hline,MAX_FILENAME,fph); //read the comment line
  
  //get halo chunk size and mem
  NhAlloc = (allData.HaloChunkSizeMB*1024l*1024l)/sizeof(Halo);
  assert(NhAlloc > 0);
  
  halos = (Halo*)malloc(sizeof(Halo)*NhAlloc);
  assert(halos != NULL);
  
  //get # of halo chunks
  baseChunkSize = NhAlloc;
  Nchunks = NhTot/NhAlloc;
  if(NhTot != baseChunkSize*Nchunks)
    ++Nchunks;
  
  //open Rnn file
  sprintf(hline,"%s",allData.HaloFile);
  extract_filename(hline);
  sprintf(fname,"%s/rnn_%s",allData.OutputPath,hline);
  fp = fopen_retry(fname,"w");
  fprintf(fp,"#ID Rnn [for N = %d]\n",allData.NRnn);
  
  for(chunk=0;chunk<Nchunks;++chunk)
    {
      //get parms of chunk
      Nh = baseChunkSize;
      chunkStart = baseChunkSize*chunk;
      if(chunkStart + Nh > NhTot)
	Nh = NhTot - 1 - chunkStart + 1;
      
      //read halos
      NhRead = read_halo_chunk(halos,Nh,fph,NULL);
      assert(NhRead == Nh);
      for(i=0;i<Nh;++i)
	halos[i].haloFileIndex = chunkStart + i;
      
      //sort them on disk into work units
      compute_bboxes_for_halos(halos,Nh,boxes,Nb);
      NumHaloWorkUnits = sort_halo_chunk_into_bboxes(halos,Nh,boxes,Nb);
            
      //do rnn for each work unit
      startWorkHalos = 0;
      for(haloWorkUnit=0;haloWorkUnit<NumHaloWorkUnits;++haloWorkUnit)
	{
	  //get all halos with the same work index
	  workHaloBBoxFileIndex = halos[startWorkHalos].bboxFileIndex;
	  NworkHalos = 0;
	  for(i=startWorkHalos;i<Nh;++i)
	    {
	      if(halos[i].bboxFileIndex == workHaloBBoxFileIndex)
		++NworkHalos;
	      else
		break;
	    }
	  
	  //now make extra halo part buffers
	  if(NworkHalos > NhBuffMax)
	    {
	      NhBuffMax = NworkHalos;
	      for(i=0;i<3;++i)
		{
		  tmp = (float*)realloc(fptrs[i],sizeof(float)*NhBuffMax);
		  assert(tmp != NULL);
		  fptrs[i] = tmp;
		}
	    }
	  
	  //fill the buffer
	  for(i=0;i<NworkHalos;++i)
	    {
	      fptrs[0][i] = halos[startWorkHalos+i].px;
	      fptrs[1][i] = halos[startWorkHalos+i].py;
	      fptrs[2][i] = halos[startWorkHalos+i].pz;
	    }
	  
	  //save old boxes
	  k = 0;
	  for(i=0;i<Nb;++i)
	    if(boxes[i].fnum == workHaloBBoxFileIndex)
	      {
		tmp_bboxes[k] = boxes[i];
		++k;
	      }
	  
	  //resize boxes to contain all of the halos
	  for(i=0;i<NworkHalos;++i)
	    {
	      k = halos[startWorkHalos+i].bboxIndex;
	      
	      if(boxes[k].bmin[0] >= halos[startWorkHalos+i].px)
		boxes[k].bmin[0] = halos[startWorkHalos+i].px;
	      if(boxes[k].bmax[0] <= halos[startWorkHalos+i].px)
		boxes[k].bmax[0] = halos[startWorkHalos+i].px;
	      
	      if(boxes[k].bmin[1] >= halos[startWorkHalos+i].py)
		boxes[k].bmin[1] = halos[startWorkHalos+i].py;
	      if(boxes[k].bmax[1] <= halos[startWorkHalos+i].py)
		boxes[k].bmax[1] = halos[startWorkHalos+i].py;
	      
	      if(boxes[k].bmin[2] >= halos[startWorkHalos+i].pz)
		boxes[k].bmin[2] = halos[startWorkHalos+i].pz;
	      if(boxes[k].bmax[2] <= halos[startWorkHalos+i].pz)
		boxes[k].bmax[2] = halos[startWorkHalos+i].pz;
	    }
	  
	  //read the parts needed
	  readparts(workHaloBBoxFileIndex,&px,&py,&pz,&NpFile,&NpTot,boxes,Nb);
	  
	  //reset boxes to old versions
	  k = 0;
	  for(i=0;i<Nb;++i)
	    if(boxes[i].fnum == workHaloBBoxFileIndex)
	      {
		boxes[i] = tmp_bboxes[k];
		++k;
	      }
	  
	  //finally do Rnn
	  rnncalc_halos(fptrs[0],fptrs[1],fptrs[2],NworkHalos,px,py,pz,NpTot,&rnn,&Nrnn);
	  assert(Nrnn == NworkHalos);
	  free(px);
	  free(py);
	  free(pz);
	  
	  //now store Rnn back in main halo vector
	  for(i=0;i<NworkHalos;++i)
	    halos[startWorkHalos+i].rnn = rnn[i];
	  free(rnn);
	  
	  //clean up and set parms for next chunk
	  startWorkHalos += NworkHalos;
	}
      
      //resort into file index order
      for(i=0;i<Nh;++i)
	halos[i].sortIndex = halos[i].haloFileIndex;
      sort_halos_sortindex(halos,Nh);
      
      //write to disk
      for(i=0;i<Nh;++i)
	fprintf(fp,"%ld %f\n",halos[i].id,halos[i].rnn);
    
    }//for(chunk=0;chunk<Nchunks;++chunk)
  
  //clean up
  free(halos);
  fclose(fp);
  fclose(fph);
  if(NhBuffMax > 0)
    {
      free(fptrs[0]);
      free(fptrs[1]);
      free(fptrs[2]);
    }
  free(tmp_bboxes);
}

static float get_rnn_val_knn(float centPos[3], int periodic, double BoxLength, float *nbrsRad2, int *nbrsInd, long NRnnFind, float *px, float *py, float *pz, kdTreeData *kdTree)
{
  int retNRnn,k,j;
  float max2Rnn[2];
  float maxRad;
  
  //get Rnn
  retNRnn = get_knnbrs_kdtree(centPos,periodic,allData.BoxLength,nbrsRad2,nbrsInd,NRnnFind,px,py,pz,kdTree);

  //count # of nbrs not including the particle
  //keep track of 2 max distances
  k = 0;
  for(j=0;j<retNRnn;++j)
    {
      if(nbrsRad2[j] != 0.0)
	{
	  if(k < 2)
	    max2Rnn[k] = nbrsRad2[j];
	  else
	    {
	      if(nbrsRad2[j] > max2Rnn[0] || nbrsRad2[j] > max2Rnn[1])
		{
		  if(max2Rnn[1] > max2Rnn[0])
		    max2Rnn[0] = nbrsRad2[j];
		  else
		    max2Rnn[1] = nbrsRad2[j];
		}
	    }
	  
	  ++k;
	}
    }
  
  //if particle is in own nnbrs list then keep max
  if(k != retNRnn)
    {
      if(max2Rnn[0] > max2Rnn[1])
	maxRad = max2Rnn[0];
      else
	maxRad = max2Rnn[1];
    }
  else //keep second max
    {
      if(max2Rnn[0] > max2Rnn[1])
	maxRad = max2Rnn[1];
      else
	maxRad = max2Rnn[0];
    }
  
  return sqrt(maxRad);
}

static int compFloat(const void *a, const void *b)
{
  if((*((const float*) a)) > (*((const float*) b)))
    return 1;
  else if((*((const float*) a)) < (*((const float*) b)))
    return -1;
  else
    return 0;
}

static float get_rnn_val_nn_qsort(float centPos[3], int periodic, double BoxLength,
				  float **nbrsRad2, int **nbrsInd, int *maxNNbrs, float *Rmax, long NRnnFind,
				  float *px, float *py, float *pz, kdTreeData *kdTree)
{
  int retNNbrs;
  int j;
  float fac = 1.1;

  //get all neighbors in some distance
  retNNbrs = get_nnbrs_kdtree(centPos,*Rmax,periodic,allData.BoxLength,nbrsRad2,nbrsInd,maxNNbrs,px,py,pz,kdTree);

  while(retNNbrs < NRnnFind && *Rmax < allData.DomainBuffSize/fac)
    {
      *Rmax = (*Rmax)*fac;
      retNNbrs = get_nnbrs_kdtree(centPos,*Rmax,periodic,allData.BoxLength,nbrsRad2,nbrsInd,maxNNbrs,px,py,pz,kdTree);
    }

  //make sure did not exceed buffers
  if(*Rmax > allData.DomainBuffSize)
    {
      fprintf(stderr,"%d: Rmax value greater than domain buffer size! (Rmax = %lf, buff = %lf)\n",ThisTask,*Rmax,allData.DomainBuffSize);
      MPI_Abort(MPI_COMM_WORLD,999);
    }

  //make sure particle is not its own neighbor
  for(j=0;j<retNNbrs;++j)
    if((*nbrsRad2)[j] == 0.0)
      {
        (*nbrsInd)[j] = (*nbrsInd)[retNNbrs-1];
        (*nbrsRad2)[j] = (*nbrsRad2)[retNNbrs-1];
        --retNNbrs;
        break;
      }

  //sort by distance
  qsort(*nbrsRad2,(size_t) retNNbrs,sizeof(float),compFloat);

  return sqrt((*nbrsRad2)[NRnnFind-2]);
}

void rnncalc_halos(float *phx, float *phy, float *phz, int Nh, float *px, float *py, float *pz, int Np, float **rnn, int *Nrnn)
{
  int i;
  int periodic;
  float centPos[3];
  float *nbrsRad2;
  int *nbrsInd;
  int NumToTest,NumTested;
  kdTreeData *kdTree;
  double time[2];
  long *pinds;
  long NRnnFind = allData.NRnn + 1;
  int maxNNbrs;
  float Rmax,maxRad,Rstart;
  double tnn,tknn;

#ifdef USEOPENMP
#define TFUNC omp_get_wtime()
#else
#define TFUNC MPI_Wtime()
#endif

  //alloc part indexes
  pinds = (long*)malloc(sizeof(long)*Np);
  assert(pinds != NULL);
  for(i=0;i<Np;++i)
    pinds[i] = -1;
  
  //presort parts along PH curve
  reorder_partshalos_phcurve(px,py,pz,pinds,Np);
  
  //build the tree
  time[0] = -MPI_Wtime();
  kdTree = buildkdTree(px,py,pz,Np,NULL);
  time[0] += MPI_Wtime();

#ifdef DEBUG
#if DEBUG_LEVEL > 0
  fprintf(stderr,"%d: built kdTree in %f seconds.\n",ThisTask,time[0]);
#endif
#endif
  
  (*rnn) = (float*)malloc(sizeof(float)*Nh);
  assert((*rnn) != NULL);
  *Nrnn = Nh;
  
  for(i=0;i<Nh;++i)
    (*rnn)[i] = 0.0;
  
  //setup
  if(allData.BoxLength > 0)
    periodic = 1;
  else
    periodic = 0;
  time[1] = -MPI_Wtime();
  
#ifdef USEOPENMP    
#pragma omp parallel default(none) \
  shared(ThisTask,phx,phy,phz,kdTree,px,py,pz,allData,periodic,rnn,Nh,NRnnFind)	\
  private(nbrsRad2,nbrsInd,NumToTest,i,NumTested,centPos,maxRad,tnn,tknn,Rmax,maxNNbrs,Rstart) 
#endif
  {
    nbrsRad2 = (float*)malloc(sizeof(float)*NRnnFind);
    assert(nbrsRad2 != NULL);
    nbrsInd = (int*)malloc(sizeof(int)*NRnnFind);
    assert(nbrsInd != NULL);
    maxNNbrs = NRnnFind;
    Rmax = 2.0;
    
#ifdef USEOPENMP
    NumToTest = Nh/omp_get_num_threads();
#else
    NumToTest = Nh;
#endif
    NumTested = 0;

#ifdef USEOPENMP
#pragma omp for schedule(guided)
#endif
    for(i=0;i<Nh;++i)
      {

#ifdef DEBUG
#if DEBUG_LEVEL > 2
	if((NumTested%(NumToTest/25) == 0 && ThisTask == 1 
#ifdef USEOPENMP
	    && omp_get_thread_num() == 0
#endif
	    && NumTested > 0) && 1)
          fprintf(stderr,"%d of %d (%.2lf percent) done in rnn calc.\n",
		  NumTested,NumToTest,((double) NumTested)/((double) NumToTest)*100.0);
#endif
#endif
	++NumTested;
	
	centPos[0] = phx[i];
	centPos[1] = phy[i];
	centPos[2] = phz[i];
	
	//if starting, time both options and pick the fastest                                                                                                                                                                                                                                        
	if(NumTested == 1)
	  {
	    tknn = -TFUNC;
	    maxRad = get_rnn_val_knn(centPos,periodic,allData.BoxLength,nbrsRad2,nbrsInd,NRnnFind,px,py,pz,kdTree);
	    tknn += TFUNC;
	    
	    tnn = -TFUNC;
	    maxRad = get_rnn_val_nn_qsort(centPos,periodic,allData.BoxLength,&nbrsRad2,&nbrsInd,&maxNNbrs,&Rmax,NRnnFind,px,py,pz,kdTree);
	    tnn += TFUNC;
	  }
	
	//always use the fastest
	if(tknn < tnn && NRnnFind < 700)
	  maxRad = get_rnn_val_knn(centPos,periodic,allData.BoxLength,nbrsRad2,nbrsInd,NRnnFind,px,py,pz,kdTree);
	else
	  {
	    Rstart = Rmax;
	    maxRad = get_rnn_val_nn_qsort(centPos,periodic,allData.BoxLength,&nbrsRad2,&nbrsInd,&maxNNbrs,&Rmax,NRnnFind,px,py,pz,kdTree);
	    Rmax = (Rstart*(NumTested-1.0) + Rmax)/NumTested;
	  }
	
#ifdef USEOPENMP
#pragma omp atomic
	(*rnn)[i] += maxRad;
#else
	(*rnn)[i] += maxRad;
#endif
      }
    
    free(nbrsRad2);
    free(nbrsInd);
  } // end of parallel region
  
  //error check rnn values
  for(i=0;i<Nh;++i)
    if((*rnn)[i] > allData.DomainBuffSize)
      {
        fprintf(stderr,"%d: Rnn value greater than domain buffer size! (rnn = %lf, buff = %lf)\n",ThisTask,(*rnn)[i],allData.DomainBuffSize);
        MPI_Abort(MPI_COMM_WORLD,999);
      }
  
  time[1] += MPI_Wtime();
#ifdef DEBUG
#if DEBUG_LEVEL > 0
  fprintf(stderr,"%d: did Rnn calc in %f seconds.\n",ThisTask,time[1]);
#endif
#endif
  
#undef TFUNC
  
  //free mem
  destroykdTree(kdTree);
  free(pinds);
}

static float get_mindist_bbox_point(BBox bx, float px, float py, float pz, float L)
{
  int i,j,start = 1;
  float mind,mindf,d,dl,du;
  float pos[3];
  
  pos[0] = px;
  pos[1] = py;
  pos[2] = pz;
  
  //loop through each face
  //compute min dist from face
  //then take min overall
  for(j=0;j<3;++j)
    {
      /////////////////
      //lower face
      /////////////////
      d = bx.bmin[j]-pos[j];
      mindf = d*d;
      
      for(i=0;i<3;++i)
	if(i != j)
	  {
	    if(pos[i] >= bx.bmin[i] && pos[i] <= bx.bmax[i])
	      d = 0.0;
	    else
	      {
		dl = fabs(bx.bmin[i]-pos[i]);
		du = fabs(bx.bmax[i]-pos[i]);
		if(dl < du)
		  mindf += dl*dl;
		else
		  mindf += du*du;
	      }
	  }
      
      //check (and init)
      if(start)
	{
	  mind = mindf;
	  start = 0;
	}
      else if(mindf < mind)
	mind = mindf;
      
      /////////////////
      //upper face
      /////////////////
      d = bx.bmax[j]-pos[j];
      mindf = d*d;
      
      for(i=0;i<3;++i)
	if(i != j)
	  {
	    if(pos[i] >= bx.bmin[i] && pos[i] <= bx.bmax[i])
	      d = 0.0;
	    else
	      {
		dl = fabs(bx.bmin[i]-pos[i]);
		du = fabs(bx.bmax[i]-pos[i]);
		if(dl < du)
		  mindf += dl*dl;
		else
		  mindf += du*du;
	      }
	  }
      
      //check
      if(mindf < mind)
	mind = mindf;
    }
  
  return mind;
}

static void compute_bboxes_for_halos(Halo *halos, long Nh, BBox *boxes, int Nb)
{
  long i,j;
  int fvmax,fmind,bvmax,bmind;
  float vmax,mind,d;
  
  //get box for each halo
  for(i=0;i<Nh;++i)
    {
      fvmax = -1;
      fmind = -1;
      bvmax = -1;
      bmind = -1;
      halos[i].bboxFileIndex = -1;
      halos[i].bboxIndex = -1;
      
      for(j=0;j<Nb;++j)
	{
	  if(test_point_bboxes_general(boxes[j],halos[i].px,halos[i].py,halos[i].pz,allData.BoxLength,0.0))
	    {
	      //if halo is in a box, make sure it is the one with the most volume
	      if(fvmax == -1)
		{
		  vmax = vol_bbox(boxes[j]);
		  fvmax = boxes[j].fnum;
		  bvmax = j;
		}
	      else if(vol_bbox(boxes[j]) > vmax)
		{
		  vmax = vol_bbox(boxes[j]);
		  fvmax = boxes[j].fnum;
		  bvmax = j;
		}
	    }
	  else if(fvmax == -1)
	    {
	      //if halo is not in a box, then we will put it into the closest box, so keep track of that
	      d = get_mindist_bbox_point(boxes[j],halos[i].px,halos[i].py,halos[i].pz,allData.BoxLength);
	      
	      if(fmind == -1)
		{
		  mind = d;
		  fmind = boxes[j].fnum;
		  bmind = j;
		}
	      else if(d < mind)
		{
		  mind = d;
		  fmind = boxes[j].fnum;
		  bmind = j;
		}
	    }
	}
      
      //what do we do if halo is not in any of the bboxes in the file?
      //place it into the closest box!
      if(fvmax == -1)
	{
	  assert(fmind != -1);
	  halos[i].bboxFileIndex = fmind;
	  halos[i].bboxIndex = bmind;
	}
      else
	{
	  halos[i].bboxFileIndex = fvmax;
	  halos[i].bboxIndex = bvmax;
	}
      
      assert(halos[i].bboxFileIndex >= 0);
      assert(halos[i].bboxIndex >= 0);
    }
}

static long sort_halo_chunk_into_bboxes(Halo *halos, long Nh, BBox *boxes, int Nb)
{
  long i;
  
  //sort halos by bboxes
  fprintf(stderr,"sorting halos by BBoxes.\n");
  
  for(i=0;i<Nh;++i)
    halos[i].sortIndex = halos[i].bboxFileIndex;
  sort_halos_sortindex(halos,Nh);
  
  //sort halos in each bbox by ph curve
  //also get total # of uniqe boxes
  fprintf(stderr,"sorting halos by PH ind for each BBox.\n");
  
  long Nhbx = 0;
  long start,N;
  int sortChunk;
  
  start = 0;
  N = 1;
  sortChunk = 0;
  for(i=1;i<Nh;++i)
    {
      if(halos[i].bboxFileIndex == halos[i-1].bboxFileIndex)
	++N;
      else
	sortChunk = 1;
	
      if(sortChunk)
	{
	  //sort halos by PH inds
	  reorder_halos_phcurve(halos+start,N);
	  
	  start += N;
	  N = 1;
	  
	  ++Nhbx;
	  sortChunk = 0;
	}
    }
  
  //always sort the chunk
  reorder_halos_phcurve(halos+start,N);
  start += N;
  ++Nhbx;
    
  //make sure have covered all of the halos
  assert(start == Nh);
  
  return Nhbx;
}

static void reorder_halos_phcurve(Halo *halos, long Nh)
{
  int i;
  float domainLengths[3],domainStart[3],dmax[3];
  int N,x,y,z;
  double dx,dy,dz;
  
  //get start and length of domain in each direction
  domainStart[0] = halos[0].px;
  domainStart[1] = halos[0].py;
  domainStart[2] = halos[0].pz;
  dmax[0] = domainStart[0];
  dmax[1] = domainStart[1];
  dmax[2] = domainStart[2];
  for(i=0;i<Nh;++i)
    {
      if(halos[i].px > dmax[0])
        dmax[0] = halos[i].px;
      
      if(halos[i].py > dmax[1])
        dmax[1] = halos[i].py;
      
      if(halos[i].pz > dmax[2])
        dmax[2] = halos[i].pz;
      
      if(halos[i].px < domainStart[0])
        domainStart[0] = halos[i].px;
      
      if(halos[i].py < domainStart[1])
        domainStart[1] = halos[i].py;
      
      if(halos[i].pz < domainStart[2])
        domainStart[2] = halos[i].pz;
    }
  domainLengths[0] = dmax[0] - domainStart[0];
  domainLengths[1] = dmax[1] - domainStart[1];
  domainLengths[2] = dmax[2] - domainStart[2];
    
  N = 1;
  N = N << BITS_PER_DIMENSION;
  
  dx = domainLengths[0]/((double) (N));
  dy = domainLengths[1]/((double) (N));
  dz = domainLengths[2]/((double) (N));
  
  for(i=0;i<Nh;++i)
    {
      x = (int) ((halos[i].px-domainStart[0])/dx);
      if(x >= N)
        x = N-1;
      if(x < 0)
        x = 0;
      
      y = (int) ((halos[i].py-domainStart[1])/dy);
      if(y >= N)
        y = N-1;
      if(y < 0)
        y = 0;
      
      z = (int) ((halos[i].pz-domainStart[2])/dz);
      if(z >= N)
        z = N-1;
      if(z < 0)
        z = 0;
      
      halos[i].sortIndex = peano_hilbert_key(x,y,z,BITS_PER_DIMENSION);
    }
  
  sort_halos_sortindex(halos,Nh);
}

static void sort_halos_sortindex(Halo *halos, long Nh)
{
  qsort(halos,(size_t) Nh,sizeof(Halo),compHaloSortIndex);
}

static int compHaloSortIndex(const void *a, const void *b)
{
  if(((const Halo*)a)->sortIndex > ((const Halo*)b)->sortIndex)
    return 1;
  else if(((const Halo*)a)->sortIndex < ((const Halo*)b)->sortIndex)
    return -1;
  else
    return 0;
}

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

//#define DEBUG_QWCOMM

void do_rnn_parts(BBox *boxes, int Nb)
{
  int i;
  float *px,*py,*pz,*rnn;
  int NpTot,NpFile,Nrnn,Nf,Np;
  int startF,endF,dF;
  int Nmod;
  
  Nf = allData.Nf;
  
  if(NTasks == 1)
    {
      dF = Nf/NTasks;
      startF = dF*ThisTask;
      if(ThisTask == NTasks-1)
        endF = Nf-1;
      else
        endF = startF + dF - 1;
      
      Nmod = (endF-startF+1.0)/10.0;
      if(Nmod == 0)
        Nmod = 1;
      for(i=startF;i<=endF;++i)
        {
          if(ThisTask == 0 && (i-startF+1)%(Nmod) == 0)
            fprintf(stderr,"Rnn: task 0 on file %d of %d (%.2f percent).\n",
                    i-startF+1,endF-startF+1,
                    ((double) (i-startF+1))/((double) (endF-startF+1))*100.0);

	  //Check to see if file contains particles
	  Np = read_num_parts_file(allData.fnames[i]);
	  if(Np == 0)
	    {
	      write_rnn(i,rnn,0);
	      continue;
	    }
	    
	  readparts(i,&px,&py,&pz,&NpFile,&NpTot,boxes,Nb);

          //do rnn calc
          rnncalc_parts(px,py,pz,NpFile,NpTot,&rnn,&Nrnn);
          free(px);
          free(py);
          free(pz);

          //write to disc
          write_rnn(i,rnn,Nrnn);
          free(rnn);
        }
    }
  else //do queen-worker comp - should load balance better - uses queues to control I/O load
    {
      int NumTasksDoingIO,NumReadQueue,NumWriteQueue,NumTasksWorking;
      int *readQueue,*writeQueue;
      MPI_Status status;
      int fnumCurr,workerFnum;
      MPI_Request request;
      int NumDone;
      
      if(ThisTask == 0) //queen code
        {
          //alloc queues
          readQueue = (int*)malloc(sizeof(int)*NTasks);
          assert(readQueue != NULL);
          writeQueue = (int*)malloc(sizeof(int)*NTasks);
          assert(writeQueue != NULL);
          
          //init for loop logic
          NumDone = 0;
          fnumCurr = 0;
          NumTasksDoingIO = 0;
          NumWriteQueue = 0;
          NumTasksWorking = 0;
          NumReadQueue = 0;
                  
          while(fnumCurr < Nf || NumTasksWorking > 0)
            {
              //see if anyone needs to write
              if(NumWriteQueue > 0 && NumTasksDoingIO < allData.NumTasksIOInParallel)
                {
                  MPI_Send(&fnumCurr,1,MPI_INT,writeQueue[NumWriteQueue-1],TAG_WRITE_OK,MPI_COMM_WORLD);
                  --NumWriteQueue;
                  ++NumTasksDoingIO;
                }
              
              //see if anyone needs to read
              if(NumReadQueue > 0 && NumTasksDoingIO < allData.NumTasksIOInParallel && fnumCurr < Nf)
                {
                  MPI_Send(&fnumCurr,1,MPI_INT,readQueue[NumReadQueue-1],TAG_WORKWORKER,MPI_COMM_WORLD);
                  --NumReadQueue;
                  ++NumTasksDoingIO;
                  ++NumTasksWorking;
                  ++fnumCurr;
                }
              
              //now recv messages from workers
              MPI_Recv(&workerFnum,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
              
              if(status.MPI_TAG == TAG_READ_REQ)
                {
                  readQueue[NumReadQueue] = status.MPI_SOURCE;
                  ++NumReadQueue;
                }
              else if(status.MPI_TAG == TAG_WRITE_REQ)
                {
                  writeQueue[NumWriteQueue] = status.MPI_SOURCE;
                  ++NumWriteQueue;
                }
              else if(status.MPI_TAG == TAG_READ_DONE)
                {
                  --NumTasksDoingIO;
                }
              else if(status.MPI_TAG == TAG_WRITE_DONE)
                {
                  --NumTasksDoingIO;
                  --NumTasksWorking;
                  ++NumDone;
                  
                  fprintf(stderr,"Rnn: finished %d of %d (%.2f percent) of files.\n",
                          NumDone,Nf,((double) NumDone)/((double) Nf)*100.0);
                }
              else
                {
                  fprintf(stderr,"recv weird tag at root! tag,src = %d|%d\n",status.MPI_TAG,status.MPI_SOURCE);
                  MPI_Abort(MPI_COMM_WORLD,123);
                }
            }
          
          assert(NumDone == Nf);
          assert(fnumCurr == Nf);
	  
	  //make sure to get last of read requests
	  while(NumReadQueue < NTasks-1)
	    {
	      MPI_Recv(&workerFnum,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	      
	      if(status.MPI_TAG == TAG_READ_REQ)
                {
                  readQueue[NumReadQueue] = status.MPI_SOURCE;
                  ++NumReadQueue;
                }
	      else
                {
                  fprintf(stderr,"recv weird tag at root while cleaning up read requests! tag,src = %d|%d\n",status.MPI_TAG,status.MPI_SOURCE);
                  MPI_Abort(MPI_COMM_WORLD,123);
                }
	    }
	  
          //now kill the worker tasks
          for(i=1;i<NTasks;++i)
            MPI_Send(&fnumCurr,1,MPI_INT,i,TAG_KILLWORKER,MPI_COMM_WORLD);
	  
	  free(readQueue);
	  free(writeQueue);
        }
      else //worker code
        {
          fnumCurr = -1;
          
          while(1)
            {
              if(fnumCurr == -1)
                MPI_Send(&workerFnum,1,MPI_INT,0,TAG_READ_REQ,MPI_COMM_WORLD);
              
              MPI_Recv(&workerFnum,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
              
              if(status.MPI_TAG == TAG_KILLWORKER)
                break;
              else if(status.MPI_TAG == TAG_WORKWORKER)
                {
                  fnumCurr = workerFnum;

#ifdef DEBUG_QWCOMM               
                  fprintf(stderr,"%d: read %d\n",ThisTask,fnumCurr);
#endif
		  //Check to see if file contains particles
		  Np = read_num_parts_file(allData.fnames[fnumCurr]);
		  if(Np == 0)
		    {
		      //request permission to write
		      MPI_Send(&workerFnum,1,MPI_INT,0,TAG_WRITE_REQ,MPI_COMM_WORLD);
		      MPI_Recv(&workerFnum,1,MPI_INT,0,TAG_WRITE_OK,MPI_COMM_WORLD,&status);
		      
		      //Write file with just header
		      write_rnn(fnumCurr,rnn,0);
		      MPI_Send(&workerFnum,1,MPI_INT,0,TAG_WRITE_DONE,MPI_COMM_WORLD); 

		      //Continue to next file
		      fnumCurr=-1;
		      continue;
		    }
		  
                  //read parts
                  readparts(fnumCurr,&px,&py,&pz,&NpFile,&NpTot,boxes,Nb);
                  
                  //tell root done reading
                  MPI_Isend(&workerFnum,1,MPI_INT,0,TAG_READ_DONE,MPI_COMM_WORLD,&request);
                  
#ifdef DEBUG_QWCOMM               
                  fprintf(stderr,"%d: rnn %d\n",ThisTask,fnumCurr);
#endif
                  
                  //do rnn calc
                  rnncalc_parts(px,py,pz,NpFile,NpTot,&rnn,&Nrnn);
                  free(px);
                  free(py);
                  free(pz);
                  
                  //make sure root knows
                  MPI_Wait(&request,&status);
                  
                  //request permission to write
                  MPI_Send(&workerFnum,1,MPI_INT,0,TAG_WRITE_REQ,MPI_COMM_WORLD);
                  MPI_Recv(&workerFnum,1,MPI_INT,0,TAG_WRITE_OK,MPI_COMM_WORLD,&status);
                  
#ifdef DEBUG_QWCOMM               
                  fprintf(stderr,"%d: write %d\n",ThisTask,fnumCurr);
#endif
                  
                  //write to disc
                  write_rnn(fnumCurr,rnn,Nrnn);
                  free(rnn);
                  
                  //tell root done writing
                  MPI_Send(&workerFnum,1,MPI_INT,0,TAG_WRITE_DONE,MPI_COMM_WORLD); 
                  
                  fnumCurr = -1;
                }
              else
                {
                  fprintf(stderr,"worker %d recv weird tag! tag,src = %d|%d\n",ThisTask,status.MPI_TAG,status.MPI_SOURCE);
                  MPI_Abort(MPI_COMM_WORLD,999);
                }
            } //end of while loop
        } //end of worker block
    } //end of queen-worker block
}

static float get_rnn_val_knn(int i, float centPos[3], int periodic, double BoxLength, float *nbrsRad2, int *nbrsInd, long NRnnFind, float *px, float *py, float *pz, kdTreeData *kdTree)
{
  int retNRnn,k,j;
  float max2Rnn[2];
  float maxRad;
  
  //do the knn call
  retNRnn = get_knnbrs_kdtree(centPos,periodic,allData.BoxLength,nbrsRad2,nbrsInd,NRnnFind,px,py,pz,kdTree);

  //count # of nbrs not including the particle
  //keep track of 2 max distances
  k = 0;
  for(j=0;j<retNRnn;++j)
    {
      if(nbrsInd[j] != i)
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

static float get_rnn_val_nn_qsort(int i, float centPos[3], int periodic, double BoxLength, 
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
    if((*nbrsInd)[j] == i)
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

void rnncalc_parts(float *px, float *py, float *pz, int NpFile, int NpTot, float **rnn, int *Nrnn)
{
  long i,j;
  int periodic;
  float centPos[3];
  float *nbrsRad2;
  int *nbrsInd;
  int NumToTest,NumTested;
  kdTreeData *kdTree;
  double time[2];
  long *pinds;
  long *rnninds;
  long NRnnFind = allData.NRnn + 1;
  int maxNNbrs;
  float Rmax,maxRad,Rstart;
  double tnn,tknn;
  double tstart;
  int Nmod;
  
#ifdef USEOPENMP
#define TFUNC omp_get_wtime()
#else
#define TFUNC MPI_Wtime()
#endif

  //alloc part indexes
  pinds = (long*)malloc(sizeof(long)*NpTot);
  assert(pinds != NULL);
  for(i=0;i<NpTot;++i)
    pinds[i] = -1;
  for(i=0;i<NpFile;++i)
    pinds[i] = i;
  
  //presort parts along PH curve
  reorder_partshalos_phcurve(px,py,pz,pinds,NpTot);
  
  //get inds of just parts needed for Rnn
  rnninds = (long*)malloc(sizeof(long)*NpFile);
  assert(rnninds != NULL);
  j = 0;
  for(i=0;i<NpTot;++i)
    if(pinds[i] >= 0)
      {
	rnninds[j] = i;
	++j;
      }
    
  //build the tree
  time[0] = -MPI_Wtime();
  kdTree = buildkdTree(px,py,pz,NpTot,NULL);
  time[0] += MPI_Wtime();

#ifdef DEBUG
#if DEBUG_LEVEL > 0
  fprintf(stderr,"%d: built kdTree in %f seconds.\n",ThisTask,time[0]);
#endif
#endif
  
  //setup vec to record rnn
  (*rnn) = (float*)malloc(sizeof(float)*NpFile);
  assert((*rnn) != NULL);
  *Nrnn = NpFile;
    
  //setup
  if(allData.BoxLength > 0)
    periodic = 1;
  else
    periodic = 0;
  time[1] = -MPI_Wtime();
  
#ifdef USEOPENMP    
#pragma omp parallel default(none) \
  shared(ThisTask,NRnnFind,NpFile,NpTot,px,py,pz,kdTree,allData,periodic,rnn,stderr,rnninds,pinds) \
  private(i,j,nbrsRad2,nbrsInd,NumToTest,NumTested,centPos,maxRad,maxNNbrs,Rmax,tknn,tnn,Rstart,tstart,Nmod)
#endif
  {
    nbrsRad2 = (float*)malloc(sizeof(float)*NRnnFind);
    assert(nbrsRad2 != NULL);
    nbrsInd = (int*)malloc(sizeof(int)*NRnnFind);
    assert(nbrsInd != NULL);
    Rmax = 2.0;
    maxNNbrs = NRnnFind;
    tstart = TFUNC;
    
#ifdef USEOPENMP
    NumToTest = NpFile/omp_get_num_threads();
#else
    NumToTest = NpFile;
#endif
    NumTested = 0;
    Nmod = NumToTest/1000;
    if(Nmod <= 0)
      Nmod = 1;

#ifdef USEOPENMP
#pragma omp for schedule(guided)
#endif
    for(j=0;j<NpFile;++j)
      {
	i = rnninds[j];
	
#ifdef DEBUG
#if DEBUG_LEVEL > 2
	//if((NumTested%(NumToTest/25) == 0 && ThisTask == 1 && NumTested > 0
	if(
	   ((NumTested%Nmod == 0 && NumTested > 0) || NumTested == 1)
#ifdef USEOPENMP
	   && omp_get_thread_num() == 0
#endif
	   )
	  //) && 1)
	  fprintf(stderr,"%d: %d of %d (%.2lf percent) done in rnn calc with tknn = %lg, tnn = %lg, Rmax = %lf. %lf hours to finish.\n",
		  ThisTask,NumTested,NumToTest,((double) NumTested)/((double) NumToTest)*100.0,tknn,tnn,Rmax,
		  (TFUNC-tstart)/((double) NumTested)*((double) (NumToTest-NumTested))/60.0/60.0);
#endif
#endif
	
	++NumTested;
	
	centPos[0] = px[i];
	centPos[1] = py[i];
	centPos[2] = pz[i];
	
	//if starting, time both options and pick the fastest
	if(NumTested == 1)
	  {
	    tknn = -TFUNC;
	    maxRad = get_rnn_val_knn(i,centPos,periodic,allData.BoxLength,nbrsRad2,nbrsInd,NRnnFind,px,py,pz,kdTree);
	    tknn += TFUNC;
	    
	    tnn = -TFUNC;
	    maxRad = get_rnn_val_nn_qsort(i,centPos,periodic,allData.BoxLength,&nbrsRad2,&nbrsInd,&maxNNbrs,&Rmax,NRnnFind,px,py,pz,kdTree);
	    tnn += TFUNC;
	  }
	
	//always use the fastest
	if(tknn < tnn && NRnnFind < 700)
	  maxRad = get_rnn_val_knn(i,centPos,periodic,allData.BoxLength,nbrsRad2,nbrsInd,NRnnFind,px,py,pz,kdTree);
	else
	  {
	    Rstart = Rmax;
	    maxRad = get_rnn_val_nn_qsort(i,centPos,periodic,allData.BoxLength,&nbrsRad2,&nbrsInd,&maxNNbrs,&Rmax,NRnnFind,px,py,pz,kdTree);
	    Rmax = (Rstart*(NumTested-1.0) + Rmax)/NumTested;
	  }
	
	(*rnn)[pinds[i]] = maxRad;
      }//for(j=0;j<NpFile;++j)
    
    free(nbrsRad2);
    free(nbrsInd);
  } /* end of parallel region */
  
  //error check rnn values
  for(i=0;i<NpFile;++i)
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
  free(rnninds);
}

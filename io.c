#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort_int.h>

#include "allheader.h"

FILE *fopen_retry(const char *filename, const char *mode)
{
  int try,Ntry = 10;
  FILE *fp;
  
  //try Ntry times, if opens, return fp
  for(try=0;try<Ntry;++try)
    {
      fp = fopen(filename,mode);
      
      if(fp != NULL)
	return fp;
    }

  //if we get to here, return NULL
  return NULL;
}

static void get_bbox_intersect_index(int fnum, BBox *bindex, int Nbindex, int **filesToRead, int *NumFilesToRead, BBox **fbx, int *Nfbx);

float read_period_length_file(char fname[])
{
  float (*read_pl)(char *);
  
  if(strcmp(allData.SimulationType,"LGADGET") == 0)
    read_pl = &get_period_length_LGADGET;
  else if(strcmp(allData.SimulationType,"LGADGETLC") == 0)
    read_pl = &get_period_length_LGADGETLC;
  else if(strcmp(allData.SimulationType,"LGADGETBCCLC") == 0)
    read_pl = &get_period_length_LGADGETBCCLC;
  else
    {
      if(ThisTask == 0)
        fprintf(stderr,"read period length - could not find I/O code for simulation type '%s'!\n",allData.SimulationType);
      MPI_Abort(MPI_COMM_WORLD,666);
    }
  
  return read_pl(fname);
}

int read_num_parts_file(char fname[])
{
  int (*read_np)(char *);
  
  if(strcmp(allData.SimulationType,"LGADGET") == 0)
    read_np = &get_num_parts_LGADGET;
  else if(strcmp(allData.SimulationType,"LGADGETLC") == 0)
    read_np = &get_num_parts_LGADGETLC;
  else if(strcmp(allData.SimulationType,"LGADGETBCCLC") == 0)
    read_np = &get_num_parts_LGADGETBCCLC;
  else
    {
      if(ThisTask == 0)
        fprintf(stderr,"read number of particles - could not find I/O code for simulation type '%s'!\n",allData.SimulationType);
      MPI_Abort(MPI_COMM_WORLD,666);
    }
  
  return read_np(fname);
}


void readparts_file(char fname[], float **px, float **py, float **pz, int *Np)
{
  void (*read_file)(char *, float**, float**, float**, int*);
  
  if(strcmp(allData.SimulationType,"LGADGET") == 0)
    read_file = &read_LGADGET;
  else if(strcmp(allData.SimulationType,"LGADGETLC") == 0)
    read_file = &read_LGADGETLC;
  else if(strcmp(allData.SimulationType,"LGADGETBCCLC") == 0)
    read_file = &read_LGADGETBCCLC;
  else
    {
      if(ThisTask == 0)
        fprintf(stderr,"read parts - could not find I/O code for simulation type '%s'!\n",allData.SimulationType);
      MPI_Abort(MPI_COMM_WORLD,666);
    }
  
  read_file(fname,px,py,pz,Np);
  
  //if have a periodic box, make sure parts are in primary image
  long i;
  if(allData.BoxLength > 0.0)
    {
      for(i=0;i<(*Np);++i)
	{
	  (*px)[i] = wrap_pos((*px)[i],allData.BoxLength);
	  (*py)[i] = wrap_pos((*py)[i],allData.BoxLength);
	  (*pz)[i] = wrap_pos((*pz)[i],allData.BoxLength);
	}
    }
}

//writes rnn to a file
void write_rnn(int fnum, float *rnn, int Nrnn)
{
  FILE *fp;
  char fname[MAX_FILENAME];
  char fout[MAX_FILENAME];
  int dummy;

  //get snap file name
  sprintf(fname,"%s",allData.fnames[fnum]);
  extract_filename(fname);
  
  //write rnn to disk
  sprintf(fout,"%s/rnn_%s",allData.OutputPath,fname);
  fp = fopen_retry(fout,"w");
  assert(fp != NULL);

  dummy = 2*sizeof(int);
  fwrite(&dummy,sizeof(int),(size_t) 1,fp);
  fwrite(&Nrnn,sizeof(int),(size_t) 1,fp);
  fwrite(&(allData.NRnn),sizeof(int),(size_t) 1,fp);
  fwrite(&dummy,sizeof(int),(size_t) 1,fp);

  dummy = sizeof(float)*Nrnn;
  fwrite(&dummy,sizeof(int),(size_t) 1,fp);
  fwrite(rnn,(size_t) Nrnn,sizeof(float),fp);
  fwrite(&dummy,sizeof(int),(size_t) 1,fp);
  
  fclose(fp);
}

//reads all parts
void readparts(int fnum, float **px, float **py, float **pz, int *NpFile, int *NpTot, BBox *bindex, int Nbindex)
{
  int *filesToRead;
  int NumFilesToRead;
  BBox *fbx;
  int Nfbx;
  
#ifdef DEBUG
#if DEBUG_LEVEL > 0
  fprintf(stderr,"\n%d: on file %d\n",ThisTask,fnum);
#endif
#endif
  
  //first get indexes of files to read
  get_bbox_intersect_index(fnum,bindex,Nbindex,&filesToRead,&NumFilesToRead,&fbx,&Nfbx);
  assert(filesToRead[0] == fnum);
  
  float *ppx,*ppy,*ppz;
  int Npp,NpSoFar = 0,NpMax = -1;
  float *tmp;
  int file,i,j,NppKeep;
    
  for(file=0;file<NumFilesToRead;++file)
    {
#ifdef DEBUG
#if DEBUG_LEVEL > 0
      fprintf(stderr,"%d: on file %d of %d for fnum %d - %d '%s'\n",ThisTask,file+1,NumFilesToRead,fnum,filesToRead[file],allData.fnames[filesToRead[file]]);
#endif
#endif

      //read parts and chuck extras
      readparts_file(allData.fnames[filesToRead[file]],&ppx,&ppy,&ppz,&Npp);
      if(filesToRead[file] != fnum)
	{
	  NppKeep = 0;
	  for(i=0;i<Npp;++i)
	    {
	      for(j=0;j<Nfbx;++j)
		{
		  if(test_point_bboxes_general(fbx[j],ppx[i],ppy[i],ppz[i],allData.BoxLength,0.0))
		    {
		      ppx[NppKeep] = ppx[i];
		      ppy[NppKeep] = ppy[i];
		      ppz[NppKeep] = ppz[i];
		      ++NppKeep;
		      break;
		    }
		}
	    }
	}
      else
	{
	  NppKeep = Npp;
	  *NpFile = Npp;
	}
      
      //add to main mem vec
      if(file == 0)  //init main mem
	{
	  NpMax = 4*Npp;
	  NpSoFar = 0;
	  
	  *px = (float*)malloc(sizeof(float)*NpMax);
	  assert((*px) != NULL);
	  *py = (float*)malloc(sizeof(float)*NpMax);
	  assert((*py) != NULL);
	  *pz = (float*)malloc(sizeof(float)*NpMax);
	  assert((*pz) != NULL);
	}
      else if(NppKeep + NpSoFar >= NpMax)
	{
	  NpMax += Npp;
	  
	  tmp = (float*)realloc(*px,sizeof(float)*NpMax);
	  assert(tmp != NULL);
	  *px = tmp;
	  
	  tmp = (float*)realloc(*py,sizeof(float)*NpMax);
	  assert(tmp != NULL);
	  *py = tmp;
	  
	  tmp = (float*)realloc(*pz,sizeof(float)*NpMax);
	  assert(tmp != NULL);
	  *pz = tmp;
	}
      
      for(i=0;i<NppKeep;++i)
	{
	  (*px)[NpSoFar+i] = ppx[i];
	  (*py)[NpSoFar+i] = ppy[i];
	  (*pz)[NpSoFar+i] = ppz[i];
	}
      NpSoFar += NppKeep;
      
      //free mem!
      free(ppx);
      free(ppy);
      free(ppz);
    }
  
  //realloc parts
  *NpTot = NpSoFar;
  
  tmp = (float*)realloc(*px,sizeof(float)*NpSoFar);
  assert(tmp != NULL);
  *px = tmp;
  
  tmp = (float*)realloc(*py,sizeof(float)*NpSoFar);
  assert(tmp != NULL);
  *py = tmp;
  
  tmp = (float*)realloc(*pz,sizeof(float)*NpSoFar);
  assert(tmp != NULL);
  *pz = tmp;
  
  //print emmory usage
  fprintf(stderr,"%d: used %lf MB of memory for paritcles.\n",ThisTask,NpSoFar*3.0*4.0/1024.0/1024.0);
  
  //free mem
  free(fbx);
  free(filesToRead);
}

//computes which files to read for a given fnum
static void get_bbox_intersect_index(int fnum, BBox *bindex, int Nbindex, int **filesToRead, int *NumFilesToRead, BBox **fbx, int *Nfbx)
{
  int i,j,k;
  int *ftr,Nftr,NftrA,*tmp,NfbxAlloc,Nlook;
  long Nlistpix = 0;
  BBox *tmpBox,box;
  BBox *keepBoxes;
  int NkeepBoxes,NkeepBoxesAlloc;
  float qrad,  brad;
  double theta, phi;
  long *listpix = NULL;

  Nftr = 1;
  NftrA = 100;
  ftr = (int*)malloc(sizeof(int)*NftrA);
  assert(ftr != NULL);
  ftr[0] = fnum;
  
  //get decomposed boxes that cover this file
  NfbxAlloc = allData.Ndiv*allData.Ndiv*allData.Ndiv;
  *Nfbx = 0;
  (*fbx) = (BBox*)malloc(sizeof(BBox)*NfbxAlloc);
  assert((*fbx) != NULL);
  for(i=0;i<Nbindex;++i)
    if(bindex[i].fnum == fnum)
      {
	//copy box
	box = bindex[i];
	
	//get range
	box.bmin[0] -= allData.DomainBuffSize;
	box.bmin[1] -= allData.DomainBuffSize;
	box.bmin[2] -= allData.DomainBuffSize;
	
	box.bmax[0] += allData.DomainBuffSize;
	box.bmax[1] += allData.DomainBuffSize;
	box.bmax[2] += allData.DomainBuffSize;
	
	if(allData.BoxLength > 0.0)
	  {
	    //decompose if needed
	    keepBoxes = NULL;
	    NkeepBoxesAlloc = 0;
	    decompose_bbox_periodic(box,&keepBoxes,&NkeepBoxes,&NkeepBoxesAlloc,allData.BoxLength);
	  }
	else
	  {
	    keepBoxes = &box;
	    NkeepBoxes = 1;
	  }
	
	//add to total and realloc if needed
	if(NkeepBoxes + (*Nfbx) > NfbxAlloc)
	  {
	    NfbxAlloc += (*Nfbx)*5;
	    tmpBox = (BBox*)realloc(*fbx,sizeof(BBox)*NfbxAlloc);
	    assert(tmpBox != NULL);
	    *fbx = tmpBox;
	  }
	
	for(j=0;j<NkeepBoxes;++j)
	  (*fbx)[(*Nfbx)+j] = keepBoxes[j];
	(*Nfbx) += NkeepBoxes;
	
	if(allData.BoxLength > 0.0)
	  free(keepBoxes);
      }
  
  //realloc to save mem
  NfbxAlloc = (*Nfbx);
  tmpBox = (BBox*)realloc(*fbx,sizeof(BBox)*NfbxAlloc);
  assert(tmpBox != NULL);
  *fbx = tmpBox;

  Nlook = (allData.DomainBuffSize/25 < 1 ) ? 1 : allData.DomainBuffSize/25;

  //get boxes that intersect
  for(i=0;i<(*Nfbx);++i)
    {
      brad = sqrt(4*M_PI/(12*pow(1<<(*fbx)[i].forder,2))/2)/(((*fbx)[i].rbin+1)*25.0);
      if(strcmp(allData.SimulationType,"LGADGETBCCLC") == 0)
	{
	  //fprintf(stderr, "\n%d: Healpix cell, forder: %d, %d\n", (*fbx)[i].hpix, (*fbx)[i].forder);
	  nest2ang((*fbx)[i].hpix,&theta,&phi,(*fbx)[i].forder);
	}

      for(j=0;j<Nbindex;++j)
	{
	  if(strcmp(allData.SimulationType,"LGADGETBCCLC") == 0)
	    {
	      if (abs((*fbx)[i].rbin-bindex[j].rbin) <= Nlook){
		qrad = atan((brad + allData.DomainBuffSize)/(bindex[j].rbin*25.0));
		//listpix = (long*)malloc(sizeof(long)*12*(1<<(2*(*fbx)[i].forder)));
		//Nlistpix = 12*(1<<(2*(*fbx)[i].forder));
		//#ifdef DEBUG
		//#if DEBUG_LEVEL > 0
		//		fprintf(stderr,"\n%d: About to query disc, hpix, theta, phi, qrad, forder[i], forder[j], rbin[i], rbin[j]: %d, %f, %f, %f, %d, %d, %d, %d\n", ThisTask, (*fbx)[i].hpix,theta,phi,qrad,(*fbx)[i].forder, bindex[j].forder,(*fbx)[i].rbin,bindex[j].rbin);
		//#endif
		//#endif
		query_disc_inclusive_nest(theta, phi, qrad, &listpix, &Nlistpix, bindex[j].forder);
		if (Nlistpix>0){
		  for (k=0; k<Nlistpix; k++)
		    {

		      //Make sure enough memory allocated for files to read list
		      if(Nftr == NftrA)
			{
			  NftrA += 100;
			  tmp = (int*)realloc(ftr,sizeof(int)*NftrA);
			  assert(tmp != NULL);
			  ftr = tmp;
			}

		      //Check to see if this box is in nearby pixels
		      //if so, add to files to read
		      if (bindex[j].hpix==listpix[k])
			{
			  ftr[Nftr] = bindex[j].fnum;
			  ++Nftr;
			}
		    }
		  free(listpix);
		  listpix=NULL;
		  Nlistpix=0;
		}
	      }   
	    } else {
	    if(intersect_bboxes_general((*fbx)[i],bindex[j],allData.BoxLength,0.0,0.0))
	      {
		if(Nftr == NftrA)
		  {
		    NftrA += 100;
		    tmp = (int*)realloc(ftr,sizeof(int)*NftrA);
		    assert(tmp != NULL);
		    ftr = tmp;
		  }
		
		ftr[Nftr] = bindex[j].fnum;
		++Nftr;
	      }
	  }
	}
    }
    
  //get uniq files
  gsl_sort_int(ftr,(size_t) 1,(size_t) Nftr);
  NftrA = Nftr;
  Nftr = 1;
  for(i=1;i<NftrA;++i)
    if(ftr[i] != ftr[Nftr-1])
      {
	ftr[Nftr] = ftr[i];
	++Nftr;
      }
  
  //put base file at front
  for(i=0;i<Nftr;++i)
    {
      if(ftr[i] == fnum)
	{
	  ftr[i] = ftr[0];
	  ftr[0] = fnum;
	  
	  break;
	}
    }
  
  //realloc and return
  tmp = (int*)realloc(ftr,sizeof(int)*Nftr);
  assert(tmp != NULL);
  ftr = tmp;
  *filesToRead = ftr;
  *NumFilesToRead = Nftr;
}

//gets the number of lines in a file
long fnumlines(FILE *fp)
{
  long i=-1;
  char c[5000];
  while(!feof(fp))
    {
      ++i;
      fgets(c,5000,fp);
    }
  rewind(fp);
  return i;
}

//extracts file name at end of path+file string
void extract_filename(char fname[MAX_FILENAME])
{
  int len = strlen(fname);
  int fstart,flen;
  int p;
  
  p = len-1;
  flen = 0;
  while(fname[p] != '/')
    {
      fstart = p;
      ++flen;
      --p;
    }
  
  for(p=0;p<flen;++p)
    fname[p] = fname[fstart+p];
  fname[flen] = '\0';
}

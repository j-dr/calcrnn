#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort_int.h>

#include "allheader.h"

//LGadget-2 I/O header for BCC Lightcones 
struct io_header_1 {
  unsigned long npart;      /*!< npart[1] gives the number of particles in the present file, other particle types are ignored */
  unsigned int nside;
  unsigned int filenside;
  float rmin;
  float rmax;
  unsigned long npartrad;
  float BoxSize;
  double mass;          /*!< mass[1] gives the particle mass */
  double Omega0;           /*!< matter density */
  double OmegaLambda;      /*!< vacuum energy density */
  double HubbleParam;      /*!< little 'h' */
};

float get_period_length_LGADGETBCCLC(char fname[])
{
  return -1.0;
}

int get_num_parts_LGADGETBCCLC(char fname[])
{
  FILE *fp;
  int Np;
  struct io_header_1 header1;
  int dummy,k;
  float *partChunk;

  fp = fopen_retry(fname,"r");
  if(fp == NULL)
    {
      fprintf(stderr,"%d: could not open file '%s'!\n",ThisTask,fname);
      assert(fp != NULL);
    }

  fread(&header1,sizeof(header1),(size_t) 1,fp);
  Np = header1.npart;
  
  fclose(fp);
  
  return Np;
}

void get_fnside_LGADGETBCCLC(char fname[], long *fnside, long *Np)
{
  FILE *fp;
  struct io_header_1 header1;
  int dummy,k;
  float *partChunk;

  fp = fopen_retry(fname,"r");
  if(fp == NULL)
    {
      fprintf(stderr,"%d: could not open file '%s'!\n",ThisTask,fname);
      assert(fp != NULL);
    }

  fread(&header1,sizeof(header1),(size_t) 1,fp);
  *fnside = header1.filenside;
  *Np = header1.npart;
  
  fclose(fp);
}


void read_LGADGETBCCLC(char fname[], float **px, float **py, float **pz, int *Np)
{
  FILE *fp;
  struct io_header_1 header1;
  int dummy,k,temp,npix,order_=0;
  float *partChunk;
  long *idx;

  fp = fopen_retry(fname,"r");
  if(fp == NULL)
    {
      fprintf(stderr,"%d: could not open file '%s'!\n",ThisTask,fname);
      assert(fp != NULL);
    }

  fread(&header1,sizeof(header1),(size_t) 1,fp);

  *Np = header1.npart;
  temp = header1.nside;
  while(temp>>=1) order_++;
  npix = 12<<(2*order_);

  idx = (long*)malloc(sizeof(long)*(npix));
  assert((idx) != NULL);

  fread(idx,sizeof(long),(size_t) npix,fp);

  *px = (float*)malloc(sizeof(float)*(*Np));
  assert((*px) != NULL);
  *py = (float*)malloc(sizeof(float)*(*Np));
  assert((*py) != NULL);
  *pz = (float*)malloc(sizeof(float)*(*Np));
  assert((*pz) != NULL);

  //read positions                                                                                                                                                                                           
  partChunk = (float*)malloc(sizeof(float)*(*Np)*3);
  assert(partChunk != NULL);

#ifdef DEBUG
#if DEBUG_LEVEL > 1
  fprintf(stderr,"\n%d: Allocated memory for particles \n",ThisTask);
#endif
#endif

  fread(partChunk,(size_t) (*Np),3*sizeof(float),fp);
  for(k=0;k<(*Np);++k)
    {
      (*px)[k] = partChunk[k*3 + 0];
      (*py)[k] = partChunk[k*3 + 1];
      (*pz)[k] = partChunk[k*3 + 2];
    }

#ifdef DEBUG
#if DEBUG_LEVEL > 1
  fprintf(stderr,"\n%d: Read in particles \n",ThisTask);
#endif
#endif

  free(partChunk);
  
  fclose(fp);
}

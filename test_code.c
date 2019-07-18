#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_heapsort.h>

#include "allheader.h"

typedef struct {
  int index;
  float rad2;
} KNNbrsSort;

static int compKNNbrsIndexSortRad(const void *a, const void *b)
{
  if(((const KNNbrsSort*)a)->rad2 > ((const KNNbrsSort*)b)->rad2)
    return 1;
  else if(((const KNNbrsSort*)a)->rad2 < ((const KNNbrsSort*)b)->rad2)
    return -1;
  else
    return 0;
}

typedef union {
  float d;
  peanokey key;
} floatKeyUnion;

typedef struct {
  int index;
  floatKeyUnion u;
} SortPart;

static int compPeanoKey(const void *a, const void *b)
{
  if(((const SortPart*)a)->u.key > ((const SortPart*)b)->u.key)
    return 1;
  else if(((const SortPart*)a)->u.key < ((const SortPart*)b)->u.key)
    return -1;
  else
    return 0;
}

static void reorder_parts_phcurve(long Np, float *px, float *py, float *pz, float domainLengths[3])
{
  int i,j;
  float *pos[3];
  SortPart *sp;
  double sortTime = 0.0;
  int N,x,y,z;
  double dx,dy,dz;
  
#ifdef DEBUG
  if(ThisTask == 0)
    fprintf(stderr,"%d: redordering parts according to peano-hilbert index\n",ThisTask);
#endif

  //setup
  pos[0] = px;
  pos[1] = py;
  pos[2] = pz;
  
  N = 1;
  N = N << BITS_PER_DIMENSION;
  
  dx = domainLengths[0]/((double) (N));
  dy = domainLengths[1]/((double) (N));
  dz = domainLengths[2]/((double) (N));
  
#ifdef DEBUG
  if(ThisTask == 0)
    fprintf(stderr,"%d: starting sort of parts\n",ThisTask);
#endif
  sp = (SortPart*)malloc(sizeof(SortPart)*Np);
  assert(sp != NULL);
  for(i=0;i<Np;++i)
    {
      sp[i].index = i;
      x = (int) (px[i]/dx);
      if(x == N)
        --x;
      y = (int) (py[i]/dy);
      if(y == N)
        --y;
      z = (int) (pz[i]/dz);
      if(z == N)
        --z;
      sp[i].u.key = peano_hilbert_key(x,y,z,BITS_PER_DIMENSION);
    }
  sortTime -= MPI_Wtime();
  qsort(sp,(size_t) Np,sizeof(SortPart),compPeanoKey);
  sortTime += MPI_Wtime();
#ifdef DEBUG
  if(ThisTask == 0)  
    fprintf(stderr,"%d: done with sort of parts, time = %f [s]\n",ThisTask,sortTime);
#endif

  for(i=0;i<3;++i)
    {
#ifdef DEBUG
      if(ThisTask == 0)
        fprintf(stderr,"%d: on pos = %d\n",ThisTask,i);
#endif
      
      for(j=0;j<Np;++j)
        sp[j].u.d = pos[i][sp[j].index];      
      for(j=0;j<Np;++j)
        pos[i][j] = sp[j].u.d;
    }
  
  free(sp);
    
#ifdef DEBUG
  if(ThisTask == 0)
    fprintf(stderr,"%d: done redordering parts according to peano-hilbert index\n",ThisTask);
#endif
}

void test_kdtree_knnbrsfind(void)
{
  kdTreeData *td;
  float domainLengths[3],centPos[3];
  float dx,dy,dz;
  int NumKeepParts,NumNbrs,retNumNbrs,periodic,i;
  gsl_rng *rng;
  float *nbrsRad2;
  int *nbrsInd;
  float *px,*py,*pz;
  KNNbrsSort *ns,*nsf;
  float L,halfL;
  double time = 0.0;
  long Ntest=100;
  int k;
  float dlims[3][2];

  L = 1000.0;
  halfL = (float) (L/2.0);
  
  rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  gsl_rng_set(rng,(unsigned long) (ThisTask+1));
  
  NumKeepParts = 10000000;
  domainLengths[0] = L;
  domainLengths[1] = L;
  domainLengths[2] = L;
  NumNbrs = 65;
  for(i=0;i<3;++i)
    {
      dlims[i][0] = 0.0;
      dlims[i][1] = L;
    }
  
  //set max num nbrs
  nbrsRad2 = (float*)malloc(sizeof(float)*NumNbrs);
  assert(nbrsRad2 != NULL);
  nbrsInd = (int*)malloc(sizeof(int)*NumNbrs);
  assert(nbrsInd != NULL);
  nsf = (KNNbrsSort*)malloc(sizeof(KNNbrsSort)*NumNbrs);
  assert(nsf != NULL);
  
  //make parts
  ns = (KNNbrsSort*)malloc(sizeof(KNNbrsSort)*NumKeepParts);
  assert(ns != NULL);
  px = (float*)malloc(sizeof(float)*NumKeepParts);
  assert(px != NULL);
  py = (float*)malloc(sizeof(float)*NumKeepParts);
  assert(py != NULL);
  pz = (float*)malloc(sizeof(float)*NumKeepParts);
  assert(pz != NULL);
  for(i=0;i<NumKeepParts;++i)
    {
      px[i] = (float) (gsl_rng_uniform(rng)*domainLengths[0]);
      py[i] = (float) (gsl_rng_uniform(rng)*domainLengths[1]);
      pz[i] = (float) (gsl_rng_uniform(rng)*domainLengths[2]);
    }
  reorder_parts_phcurve(NumKeepParts,px,py,pz,domainLengths);
  
  for(periodic=0;periodic<2;++periodic)
    {
      td = buildkdTree(px,py,pz,(int) NumKeepParts,dlims);
      fprintf(stderr,"%d: domain = %f|%f|%f to %f|%f|%f\n",ThisTask,td->treeNodes[0].baseLoc[0],td->treeNodes[0].baseLoc[1],td->treeNodes[0].baseLoc[2],
	      td->treeNodes[0].baseLoc[0]+td->treeNodes[0].sideLengths[0],
	      td->treeNodes[0].baseLoc[1]+td->treeNodes[0].sideLengths[1],
	      td->treeNodes[0].baseLoc[2]+td->treeNodes[0].sideLengths[2]);
      
      for(k=0;k<=4;++k)
	{
	  if(periodic && (k < 0 || k > 4))
	    continue;
	  
	  centPos[0] = (float) (k*halfL/2.0);
	  centPos[1] = (float) (k*halfL/2.0);
	  centPos[2] = (float) (k*halfL/2.0);
	  	  
	  for(i=0;i<NumKeepParts;++i)
	    {
	      dx = (float) (fabs(px[i] - centPos[0]));
	      dy = (float) (fabs(py[i] - centPos[1]));
	      dz = (float) (fabs(pz[i] - centPos[2]));
	      if(periodic)
		{
		  if(dx > halfL)
		    dx = L - dx;
		  
		  if(dy > halfL)
		    dy = L - dy;
		  
		  if(dz > halfL)
		    dz = L - dz;
		}
	      
	      ns[i].rad2 = (float) ((dx*dx + dy*dy + dz*dz));
	      ns[i].index = i;
	    }
	  
	  qsort(ns,(size_t) NumKeepParts,sizeof(KNNbrsSort),compKNNbrsIndexSortRad);
	  
	  time -= MPI_Wtime();
	  for(i=0;i<Ntest;++i)
	    retNumNbrs = get_knnbrs_kdtree(centPos,periodic,L,nbrsRad2,nbrsInd,NumNbrs,px,py,pz,td);
	  time += MPI_Wtime();
	  
	  if(NumNbrs > NumKeepParts)
	    {
	      fprintf(stderr,"%d: retNumNbrs = %d (?= %d), time = %lf [s]\n",ThisTask,retNumNbrs,NumKeepParts,time/100.0);
	      assert(retNumNbrs == NumKeepParts);
	    }
	  else
	    {
	      fprintf(stderr,"%d: retNumNbrs = %d (?= %d), time = %lf [s]\n",ThisTask,retNumNbrs,NumNbrs,time/100.0);
	      assert(retNumNbrs == NumNbrs);
	    }
	  
	  for(i=0;i<retNumNbrs;++i)
	    {
	      nsf[i].rad2 = nbrsRad2[i];
	      nsf[i].index = nbrsInd[i];
	    }
	  qsort(nsf,(size_t) retNumNbrs,sizeof(KNNbrsSort),compKNNbrsIndexSortRad);
	  
	  for(i=0;i<retNumNbrs;++i)
	    {
	      if(ns[i].index != nsf[i].index || ns[i].rad2 != nsf[i].rad2)
		fprintf(stderr,"%d: ns,nsf index = %d|%d, ns,nsf rad2 = %f|%f\n",ThisTask,ns[i].index,nsf[i].index,ns[i].rad2,nsf[i].rad2);
	      assert(ns[i].index == nsf[i].index);
	      assert(ns[i].rad2 == nsf[i].rad2);
	    }
	  
	  fprintf(stderr,"%d: all knnbrs passed test!\n",ThisTask);
	}
      
      destroykdTree(td);
    }
  
  free(nbrsRad2);
  free(nbrsInd);
  free(px);
  free(py);
  free(pz);
  free(ns);
  free(nsf);
  gsl_rng_free(rng);
}

void test_kdtree_nbrsfind(void)
{
  kdTreeData *td;
  float domainLengths[3],centPos[3];
  float searchRadius,dx,dy,dz,rad;
  int NumKeepParts,NumNbrs,periodic,i,j,maxNumNbrs,n,found;
  gsl_rng *rng;
  float *nbrsRad2;
  int *nbrsInd;
  float *px,*py,*pz;
  float L,halfL;
  int k;
  
  L = 1000.0;
  halfL = (float) (L/2.0);
    
  rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  gsl_rng_set(rng,(unsigned long) (ThisTask+1));
  
  NumKeepParts = 100000;
  domainLengths[0] = L;
  domainLengths[1] = L;
  domainLengths[2] = L;
  maxNumNbrs = 10;
    
  //set max num nbrs
  maxNumNbrs = 1;
  nbrsRad2 = (float*)malloc(sizeof(float)*maxNumNbrs);
  assert(nbrsRad2 != NULL);
  nbrsInd = (int*)malloc(sizeof(int)*maxNumNbrs);
  assert(nbrsInd != NULL);
  
  //make parts
  px = (float*)malloc(sizeof(float)*NumKeepParts);
  assert(px != NULL);
  py = (float*)malloc(sizeof(float)*NumKeepParts);
  assert(py != NULL);
  pz = (float*)malloc(sizeof(float)*NumKeepParts);
  assert(pz != NULL);
  for(i=0;i<NumKeepParts;++i)
    {
      px[i] = (float) (gsl_rng_uniform(rng)*domainLengths[0]);
      py[i] = (float) (gsl_rng_uniform(rng)*domainLengths[1]);
      pz[i] = (float) (gsl_rng_uniform(rng)*domainLengths[2]);
    }
  
  for(periodic=0;periodic<2;++periodic)
    {
      td = buildkdTree(px,py,pz,(int) NumKeepParts,NULL);
      
      for(k=-6;k<=6;++k)
        {
	  if(periodic && (k < 0 || k > 4))
	    continue;
	     
          centPos[0] = (float) (k*halfL/2.0);
          centPos[1] = (float) (k*halfL/2.0);
          centPos[2] = (float) (k*halfL/2.0);
	  
	  searchRadius = (float) (halfL/2.0);
	  NumNbrs = get_nnbrs_kdtree(centPos,searchRadius,periodic,L,
				     &nbrsRad2,&nbrsInd,&maxNumNbrs,px,py,pz,td);
	  
	  fprintf(stderr,"%d: NumNbrs = %d\n",ThisTask,NumNbrs);
	  
	  for(i=0;i<NumNbrs;++i)
	    {
	      dx = (float) (fabs(px[nbrsInd[i]] - centPos[0]));
	      dy = (float) (fabs(py[nbrsInd[i]] - centPos[1]));
	      dz = (float) (fabs(pz[nbrsInd[i]] - centPos[2]));
	      if(periodic)
		{
		  if(dx > halfL)
		    dx = L - dx;
		  
		  if(dy > halfL)
		    dy = L - dy;
		  
		  if(dz > halfL)
		    dz = L - dz;
		}
	      
	      rad = (float) (sqrt(dx*dx + dy*dy + dz*dz));
	      if(rad > searchRadius || fabs(rad - sqrt(nbrsRad2[i])) >= 1e-5)
		{
		  fprintf(stderr,"%d: search did not work! rad = %f, searchRadius = %f, rad calc = %f\n",ThisTask,rad,searchRadius,sqrt(nbrsRad2[i]));
		  MPI_Abort(MPI_COMM_WORLD,123);
		}
	    }
	  
	  fprintf(stderr,"%d: all nbrs passed test!\n",ThisTask);
	  
	  searchRadius = searchRadius*searchRadius;
	  j = 0;
	  for(i=0;i<NumKeepParts;++i)
	    {
	      dx = (float) (fabs(px[i] - centPos[0]));
	      dy = (float) (fabs(py[i] - centPos[1]));
	      dz = (float) (fabs(pz[i] - centPos[2]));
	      if(periodic)
		{
		  if(dx > halfL)
		    dx = L - dx;
		  
		  if(dy > halfL)
		    dy = L - dy;
		  
		  if(dz > halfL)
		    dz = L - dz;
		}
	      
	      rad = dx*dx + dy*dy + dz*dz;
	      if(rad <= searchRadius)
		{
		  ++j;
		  found = 0;
		  for(n=0;n<NumNbrs;++n)
		    if(nbrsInd[n] == i)
		      {
			found = 1;
			break;
		      }
		  if(!found)
		    fprintf(stderr,"%d: ind = %d, rad = %le, searchRadius = %le, pos = %f|%f|%f\n",ThisTask,i,rad,searchRadius,
			    px[i],
			    py[i],
			    pz[i]
			    );
		}
	    }
	  
	  if(j != NumNbrs)
	    {
	      fprintf(stderr,"%d: search did not work! NumNbrs = %d, NumNbrs found by meshlink = %d\n",ThisTask,j,NumNbrs);
	      MPI_Abort(MPI_COMM_WORLD,123);
	    }
	  
	  fprintf(stderr,"%d: number of nbrs passed test!\n",ThisTask);
	}
      
      destroykdTree(td);
    }
  
  free(nbrsRad2);
  free(nbrsInd);
  free(px);
  free(py);
  free(pz);
  gsl_rng_free(rng);
}


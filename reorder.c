#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allheader.h"

//from gadget-2 peano.c functions
typedef union {
  float d;
  peanokey key;
  long index;
} floatKeyUnion;

typedef struct {
  long index;
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

static int quadrants[24][2][2][2] = {
  /* rotx=0, roty=0-3 */
  {{{0, 7}, {1, 6}}, {{3, 4}, {2, 5}}},
  {{{7, 4}, {6, 5}}, {{0, 3}, {1, 2}}},
  {{{4, 3}, {5, 2}}, {{7, 0}, {6, 1}}},
  {{{3, 0}, {2, 1}}, {{4, 7}, {5, 6}}},
  /* rotx=1, roty=0-3 */
  {{{1, 0}, {6, 7}}, {{2, 3}, {5, 4}}},
  {{{0, 3}, {7, 4}}, {{1, 2}, {6, 5}}},
  {{{3, 2}, {4, 5}}, {{0, 1}, {7, 6}}},
  {{{2, 1}, {5, 6}}, {{3, 0}, {4, 7}}},
  /* rotx=2, roty=0-3 */
  {{{6, 1}, {7, 0}}, {{5, 2}, {4, 3}}},
  {{{1, 2}, {0, 3}}, {{6, 5}, {7, 4}}},
  {{{2, 5}, {3, 4}}, {{1, 6}, {0, 7}}},
  {{{5, 6}, {4, 7}}, {{2, 1}, {3, 0}}},
  /* rotx=3, roty=0-3 */
  {{{7, 6}, {0, 1}}, {{4, 5}, {3, 2}}},
  {{{6, 5}, {1, 2}}, {{7, 4}, {0, 3}}},
  {{{5, 4}, {2, 3}}, {{6, 7}, {1, 0}}},
  {{{4, 7}, {3, 0}}, {{5, 6}, {2, 1}}},
  /* rotx=4, roty=0-3 */
  {{{6, 7}, {5, 4}}, {{1, 0}, {2, 3}}},
  {{{7, 0}, {4, 3}}, {{6, 1}, {5, 2}}},
  {{{0, 1}, {3, 2}}, {{7, 6}, {4, 5}}},
  {{{1, 6}, {2, 5}}, {{0, 7}, {3, 4}}},
  /* rotx=5, roty=0-3 */
  {{{2, 3}, {1, 0}}, {{5, 4}, {6, 7}}},
  {{{3, 4}, {0, 7}}, {{2, 5}, {1, 6}}},
  {{{4, 5}, {7, 6}}, {{3, 2}, {0, 1}}},
  {{{5, 2}, {6, 1}}, {{4, 3}, {7, 0}}}
};


static int rotxmap_table[24] = { 4, 5, 6, 7, 8, 9, 10, 11,
				 12, 13, 14, 15, 0, 1, 2, 3, 17, 18, 19, 16, 23, 20, 21, 22
};

static int rotymap_table[24] = { 1, 2, 3, 0, 16, 17, 18, 19,
				 11, 8, 9, 10, 22, 23, 20, 21, 14, 15, 12, 13, 4, 5, 6, 7
};

static int rotx_table[8] = { 3, 0, 0, 2, 2, 0, 0, 1 };
static int roty_table[8] = { 0, 1, 1, 2, 2, 3, 3, 0 };

static int sense_table[8] = { -1, -1, -1, +1, +1, -1, -1, -1 };

/* vars needed for inverse peano-hilbert construction
  
  static int flag_quadrants_inverse = 1;
  static char quadrants_inverse_x[24][8];
  static char quadrants_inverse_y[24][8];
  static char quadrants_inverse_z[24][8];
*/

/*! This function computes a Peano-Hilbert key for an integer triplet (x,y,z),
 *  with x,y,z in the range between 0 and 2^bits-1.
 */
peanokey peano_hilbert_key(int x, int y, int z, int bits)
{
  int i, quad, bitx, bity, bitz;
  int mask, rotation, rotx, roty, sense;
  peanokey key;
  
  mask = 1 << (bits - 1);
  key = 0;
  rotation = 0;
  sense = 1;
  
  for(i = 0; i < bits; i++, mask >>= 1)
    {
      bitx = (x & mask) ? 1 : 0;
      bity = (y & mask) ? 1 : 0;
      bitz = (z & mask) ? 1 : 0;

      quad = quadrants[rotation][bitx][bity][bitz];

      key <<= 3;
      key += (sense == 1) ? (quad) : (7 - quad);

      rotx = rotx_table[quad];
      roty = roty_table[quad];
      sense *= sense_table[quad];

      while(rotx > 0)
	{
	  rotation = rotxmap_table[rotation];
	  rotx--;
	}

      while(roty > 0)
	{
	  rotation = rotymap_table[rotation];
	  roty--;
	}
    }

  return key;
}

void reorder_partshalos_phcurve(float *px, float *py, float *pz, long *pind, int Np)
{
  int i,j;
  float domainLengths[3],domainStart[3],dmax[3];
  float *pos[3];
  SortPart *sp;
  double sortTime = 0.0;
  int N,x,y,z;
  double dx,dy,dz;
  double time[1];

#ifdef DEBUG
#if DEBUG_LEVEL > 0
  fprintf(stderr,"%d: reordering parts according to peano-hilbert index.\n",ThisTask);
#endif
#endif
  time[0] = -MPI_Wtime();
  
  //setup
  pos[0] = px;
  pos[1] = py;
  pos[2] = pz;
  
  //get start and length of domain in each direction
  domainStart[0] = pos[0][0];
  domainStart[1] = pos[1][0];
  domainStart[2] = pos[2][0];
  dmax[0] = domainStart[0];
  dmax[1] = domainStart[1];
  dmax[2] = domainStart[2];
  for(i=0;i<Np;++i)
    {
      if(pos[0][i] > dmax[0])
	dmax[0] = pos[0][i];
      
      if(pos[1][i] > dmax[1])
	dmax[1] = pos[1][i];
      
      if(pos[2][i] > dmax[2])
	dmax[2] = pos[2][i];
      
      if(pos[0][i] < domainStart[0])
	domainStart[0] = pos[0][i];
      
      if(pos[1][i] < domainStart[1])
	domainStart[1] = pos[1][i];
      
      if(pos[2][i] < domainStart[2])
	domainStart[2] = pos[2][i];
    }
  domainLengths[0] = dmax[0] - domainStart[0];
  domainLengths[1] = dmax[1] - domainStart[1];
  domainLengths[2] = dmax[2] - domainStart[2];
    
  N = 1;
  N = N << BITS_PER_DIMENSION;
  
  dx = domainLengths[0]/((double) (N));
  dy = domainLengths[1]/((double) (N));
  dz = domainLengths[2]/((double) (N));
  
  for(i=0;i<Np;++i)
    {
      if(!(px[i] >= domainStart[0] && px[i] <= dmax[0]))
	fprintf(stderr,"%d: part %d - px,py,pz = %f|%f|%f\n",ThisTask,i,
		px[i],py[i],pz[i]);
      assert(px[i] >= domainStart[0] && px[i] <= dmax[0]);
      
      if(!(py[i] >= domainStart[1] && py[i] <= dmax[1]))
	fprintf(stderr,"%d: part %d - px,py,pz = %f|%f|%f\n",ThisTask,i,
		px[i],py[i],pz[i]);
      assert(py[i] >= domainStart[1] && py[i] <= dmax[1]);
      
      if(!(pz[i] >= domainStart[2] && pz[i] <= dmax[2]))
	fprintf(stderr,"%d: part %d - px,py,pz = %f|%f|%f\n",ThisTask,i,
		px[i],py[i],pz[i]);
      assert(pz[i] >= domainStart[2] && pz[i] <= dmax[2]);
    }
  
#ifdef DEBUG
#if DEBUG_LEVEL > 1
    fprintf(stderr,"%d: starting sort of parts.\n",ThisTask);
#endif
#endif
    
  sp = (SortPart*)malloc(sizeof(SortPart)*Np);
  assert(sp != NULL);
  for(i=0;i<Np;++i)
    {
      sp[i].index = i;
      x = (int) ((px[i]-domainStart[0])/dx);
      if(x == N)
	--x;
      if(x < 0)
	x = 0;
      
      y = (int) ((py[i]-domainStart[1])/dy);
      if(y == N)
	--y;
      if(y < 0)
	y = 0;
      
      z = (int) ((pz[i]-domainStart[2])/dz);
      if(z == N)
	--z;
      if(z < 0)
	z = 0;
      
      sp[i].u.key = peano_hilbert_key(x,y,z,BITS_PER_DIMENSION);
    }
  sortTime -= MPI_Wtime();
  qsort(sp,(size_t) Np,sizeof(SortPart),compPeanoKey);
  sortTime += MPI_Wtime();

#ifdef DEBUG
#if DEBUG_LEVEL > 1  
    fprintf(stderr,"done with sort of parts, time = %f [s].\n",sortTime);
#endif
#endif
    
  for(i=0;i<3;++i)
    {
      for(j=0;j<Np;++j)
	sp[j].u.d = pos[i][sp[j].index];
      
      for(j=0;j<Np;++j)
	pos[i][j] = sp[j].u.d;
    }
  
  for(j=0;j<Np;++j)
    sp[j].u.index = pind[sp[j].index];
  for(j=0;j<Np;++j)
    pind[j] = sp[j].u.index;
  
  free(sp);
    
  time[0] += MPI_Wtime();

#ifdef DEBUG
#if DEBUG_LEVEL > 0
    fprintf(stderr,"redordered parts in %f seconds.\n",time[0]);
#endif
#endif
}


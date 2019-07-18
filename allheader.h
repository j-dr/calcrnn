#include <stdio.h>
#include <mpi.h>

#include "kdtree.h"

#ifndef _CALCRNN_
#define _CALCRNN_

#ifdef MEMWATCH
#include "memwatch.h"
#endif

/* Some debugging macros
   undef DEBUG for no debugging
   DEBUG_LEVEL = 0 is for basic debugging
   DEBUG_LEVEL = 1 is for messages printed by a single task but not as critical
   DEBUG_LEVEL = 2 or above is used for messages that every task will print 
   
   define DEBUG_IO some output as follows
*/

#ifdef NDEBUG
#undef DEBUG
#define DEBUG_LEVEL -1
#undef DEBUG_IO
#endif

#ifdef DEBUG
#ifndef DEBUG_LEVEL
#define DEBUG_LEVEL 0
#endif
#endif

#define MAX_FILENAME 1024

//#define BIG_G 4.304e-9
//#define RHO_CRIT 2.77519737e11 /* Critial mass density  in h^2 M_sun/Mpc^3 */
//#define CSOL 299792.458 /* velocity of light in km/s */

#define TAG_WORKWORKER  10
#define TAG_KILLWORKER  12
#define TAG_WRITE_REQ   13
#define TAG_WRITE_OK    13
#define TAG_WRITE_DONE  14
#define TAG_READ_REQ    15    
#define TAG_READ_DONE   16
#define TAG_RNN_DONE    17
#define TAG_RNN_DATA    18
#define TAG_HALO_DATA   19

typedef struct {
  char SimulationType[MAX_FILENAME];
  char SnapshotFileList[MAX_FILENAME]; /* file(s) to read particle data from */
  char BBoxOutputFile[MAX_FILENAME]; /* file to output bbox index to */
  char HaloFile[MAX_FILENAME];       /* rockstar type halo list*/ 
  char HaloFileFormat[MAX_FILENAME]; /* format for halo file*/
  int NumTasksIOInParallel; /* the number of tasks which are allowed to read from disk at the same time */
  float DomainBuffSize; /* size of buffer between domains for halo finding */
  float BoxLength;
  char OutputPath[MAX_FILENAME]; //output file path
  int NRnn; //# of nbrs to find
  int Ndiv; //# of divisions 
  char **fnames;
  int Nf;
  long HaloChunkSizeMB;
} AllData;

typedef struct {
  float px;           //position x
  float py;           //position x
  float pz;           //position x
  float rnn;          //Rnn value on output
  long id;            //unique ID from the halo file
  long haloFileIndex; //where in the halo file was this halo found
  long bboxFileIndex; //which file in the bbox index is this halo file "in"
  long bboxIndex;     //exact box which contains halo or is nearest to it
  long sortIndex;     //index to use for sorting halos - could be anything
} Halo;

typedef struct {
  int fnum;          //index of file for this BBox
  int rbin;
  long hpix;
  long forder;
  int ind;           //sub-index within each for for this BBox
  float bmin[3];     //min loc of box axes in all three dimensions
  float bmax[3];     //max loc of box axes in all three dimensions
  int np;            //# of objects in this BBox
} BBox;

typedef struct {
  int fnum;          //index of file for this BBox
  int rbin;
  int hpix;
  int forder;
  int np;            //# of objects in this BBox
} Cell;


//from gadget-2 peano.c functions 
typedef long peanokey;
#define BITS_PER_DIMENSION 15 /*!< Bits per dimension available for Peano-Hilbert order.                                                                       
				Note: If peanokey is defined as type int, the allowed maximum is 10.                                                                                             
				If 64-bit integers are used, the maximum is 21. */
//#define  PEANOCELLS (((peanokey)1)<<(3*BITS_PER_DIMENSION))  /*!< The number of different Peano-Hilbert cells */

/* extern defs of global vars in globalvars.c */
extern AllData allData;                              /* global struct with all vars from config file */
extern int ThisTask;                                 /* this task's rank in MPI_COMM_WORLD */
extern int NTasks;                                   /* number of tasks in MPI_COMM_WORLD */

extern int *noParticles;                             /* list of file indicies that do not contain particles */
extern int npInd;                                    /* number of files that don't contain particles */


//////////////////////////////////////////////////
//files for reading specific simulation types go here

/* in read_LGADGET.c */
float get_period_length_LGADGET(char fname[]);
int get_num_parts_LGADGET(char fname[]);
void read_LGADGET(char fname[], float **px, float **py, float **pz, int *Np);

/* in read_LGADGETLC.c */
float get_period_length_LGADGETLC(char fname[]);
int get_num_parts_LGADGETLC(char fname[]);
void read_LGADGETLC(char fname[], float **px, float **py, float **pz, int *Np);

/* in read_LGADGETBCCLC.c */
float get_period_length_LGADGETBCCLC(char fname[]);
int get_num_parts_LGADGETBCCLC(char fname[]);
void read_LGADGETBCCLC(char fname[], float **px, float **py, float **pz, int *Np);

//////////////////////////////////////////////////

/* in bbox.c */
void compute_bbox_index(BBox **boxes, int *Nboxes, int Ndiv);
BBox *compute_bboxes(int Ndiv, float *px, float *py, float *pz, int Np, int *Nb, float L);
void write_bboxes(char fname[], BBox *boxes, int Nb);
int test_point_bboxes_general(BBox box, float px, float py, float pz, float L, float buff);
int intersect_bboxes_general(BBox box1, BBox box2, float L, float buff1, float buff2);
float vol_bbox(BBox bx);
void decompose_bbox_periodic(BBox baseBox, BBox **keepBoxes, int *NkeepBoxes, int *NkeepBoxesAlloc, float L);

/* in rnn_parts.c */
void do_rnn_parts(BBox *boxes, int Nb);
void rnncalc_parts(float *px, float *py, float *pz, int NpFile, int NpTot, float **rnn, int *Nrnn);

/* in rnn_halos.c */
void do_rnn_halos(BBox *boxes, int Nb);
void rnncalc_halos(float *phx, float *phy, float *phz, int Nh, float *px, float *py, float *pz, int Np, float **rnn, int *Nrnn);

/* in config.c */
void begin_run(char config[]);
void read_config_file(char fname[]);
void end_run(void);

/* in io.c */
void readparts(int fnum, float **px, float **py, float **pz, int *NpFile, int *NpTot, BBox *bindex, int Nbindex);
long fnumlines(FILE *fp);
void readparts_file(char fname[], float **px, float **py, float **pz, int *Np);
void write_rnn(int fnum, float *rnn, int Nrnn);
void extract_filename(char fname[MAX_FILENAME]);
float read_period_length_file(char fname[]);
FILE *fopen_retry(const char *filename, const char *mode);

/* reorder.c */
void reorder_partshalos_phcurve(float *px, float *py, float *pz, long *pind, int Np);
peanokey peano_hilbert_key(int x, int y, int z, int bits);

/* utils.c */
double wrap_pos(double pin, double L);

/* healpix_utils.c */
void higher_nest(long pix, long order1_, long order2_, long *hopix);
long lower_nest(long pix, long order1_, long order2_);
void getneighbors_nest(long pix, long *result, long order_);

#ifdef TEST_CODE
/* in test_code.c */
void test_kdtree_knnbrsfind(void);
void test_kdtree_nbrsfind(void);
#endif

#endif /* _CALCRNN_ */

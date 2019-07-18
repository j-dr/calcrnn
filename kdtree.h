#ifndef _KDTREE_
#define _KDTREE_

#ifdef MEMWATCH
#include "memwatch.h"
#endif

#ifndef KDTREE_REFINEPARTLIM
#define KDTREE_REFINEPARTLIM 100
#endif

typedef struct {
  int over;
  int down;
  int splitDim;
  float splitVal;
  float baseLoc[3];
  float sideLengths[3];
  int startPart;
} kdTreeNode;

typedef struct {
  int NumTreeNodesUsed;
  int NumTreeNodes;
  kdTreeNode *treeNodes;
  int *treeParts;
  int NumParts;
} kdTreeData;

/* in kdtree.c */
int get_nnbrs_kdtree(float centPos[3], float searchRadius, int periodic, float periodLength,
                     float **nbrsRad2, int **nbrsInd, int *maxNumNbrs, float *px, float *py, float *pz,
                     kdTreeData *td);
int getNNbrsRadPeriodic(float centPos[3], float searchRadius, float periodLength, float **nbrsRad2, int **nbrsInd, int *maxNumNbrs,
                float *px, float *py, float *pz, kdTreeData *td);
int intersectSpherekdTreeNodePeriodic(kdTreeNode *node, float pos[3], float rad2, float periodLength);
int getNNbrsRad(float centPos[3], float searchRadius, float **nbrsRad2, int **nbrsInd, int *maxNumNbrs,
                float *px, float *py, float *pz, kdTreeData *td);
int intersectSpherekdTreeNode(kdTreeNode *node, float pos[3], float rad2);
kdTreeData *buildkdTree(float *px, float *py, float *pz, int Np, float (*dlims)[2]);
void destroykdTree(kdTreeData *treeData);

/* in kdtree_knnbrs.c */
int get_knnbrs_kdtree(float centPos[3], int periodic, float periodLength,
		      float *nbrsRad2, int *nbrsInd, int NumNbrs, float *px, float *py, float *pz,
		      kdTreeData *td);
int getKNNbrsPeriodic(float centPos[3], float periodLength, float *nbrsRad2, int *nbrsInd, int NumNbrs,
                      float *px, float *py, float *pz, kdTreeData *td);
int getKNNbrs(float centPos[3], float *nbrsRad2, int *nbrsInd, int NumNbrs,
              float *px, float *py, float *pz, kdTreeData *td);


#endif /* _KDTREE_ */

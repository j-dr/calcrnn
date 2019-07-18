#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "kdtree.h"

//does kdtree nearest neighbors searching - 
//can handle periodic spaces, but all points must be in the same image of the periodic space

static int getSplitNode(int node, kdTreeData *td, int *splitDim, float *splitVal, int NumPointsChildren[2], float *px, float *py, float *pz);

int get_nnbrs_kdtree(float centPos[3], float searchRadius, int periodic, float periodLength,
                     float **nbrsRad2, int **nbrsInd, int *maxNumNbrs, float *px, float *py, float *pz, 
                     kdTreeData *td)
{
  if(periodic)
    return (getNNbrsRadPeriodic(centPos,searchRadius,periodLength,nbrsRad2,nbrsInd,maxNumNbrs,px,py,pz,td));
  else
    return (getNNbrsRad(centPos,searchRadius,nbrsRad2,nbrsInd,maxNumNbrs,px,py,pz,td));
}

int getNNbrsRadPeriodic(float centPos[3], float searchRadius, float periodLength, float **nbrsRad2, int **nbrsInd, int *maxNumNbrs, 
			float *px, float *py, float *pz, kdTreeData *td)
{
  int partInd,nodeInd,NumNbrs,*tmpip;
  float dx,dy,dz,searchRadius2,rad2,*tmpfp;
  float halfPeriodLength;
  
  //init
  halfPeriodLength = (float) (periodLength/2.0);
  searchRadius2 = searchRadius*searchRadius;
  NumNbrs = 0;
  nodeInd = 0;
  
  //do search - walk linked list of nodes always going down first if node needs to be opened, else go over
  while(nodeInd >= 0)
    {
      if(intersectSpherekdTreeNodePeriodic(td->treeNodes+nodeInd,centPos,searchRadius2,periodLength))
	{
	  if(td->treeNodes[nodeInd].down >= 0 && td->treeNodes[nodeInd].startPart == -1)
	    nodeInd = td->treeNodes[nodeInd].down;
	  else if(td->treeNodes[nodeInd].startPart >= 0 && td->treeNodes[nodeInd].down == -1)
	    {
	      partInd = td->treeNodes[nodeInd].startPart;
	      while(partInd >= 0)
		{
		  dx = (float) (fabs(px[partInd] - centPos[0]));
		  if(dx > halfPeriodLength)
		    dx = periodLength - dx;
		  dy = (float) (fabs(py[partInd] - centPos[1]));
		  if(dy > halfPeriodLength)
		    dy = periodLength - dy;
		  dz = (float) (fabs(pz[partInd] - centPos[2]));
		  if(dz > halfPeriodLength)
		    dz = periodLength - dz;
		  
		  rad2 = dx*dx + dy*dy + dz*dz;
		  if(rad2 <= searchRadius2)
		    {
		      (*nbrsRad2)[NumNbrs] = (float) (rad2);
		      (*nbrsInd)[NumNbrs] = partInd;
		      ++NumNbrs;
		      
		      if(NumNbrs >= (*maxNumNbrs))
			{
			  tmpfp = (float*)realloc(*nbrsRad2,sizeof(float)*((*maxNumNbrs) + 10000));
			  assert(tmpfp != NULL);
			  *nbrsRad2 = tmpfp;
                              
			  tmpip = (int*)realloc(*nbrsInd,sizeof(int)*((*maxNumNbrs) + 10000));
			  assert(tmpip != NULL);
			  *nbrsInd = tmpip;
                              
			  (*maxNumNbrs) = (*maxNumNbrs) + 10000;
			}
		    }
                      
		  //move to next ind in linked list
		  partInd = td->treeParts[partInd];
		}
	      
	      //move to next node
	      nodeInd = td->treeNodes[nodeInd].over;
	    }
	  else
	    {
	      fprintf(stderr,"weird node - opened it up but could not go down and could not look at parts\n");
	      assert(0);
	    }
	}
      else
	nodeInd = td->treeNodes[nodeInd].over;
    }

  return NumNbrs;
}

int intersectSpherekdTreeNodePeriodic(kdTreeNode *node, float pos[3], float rad2, float periodLength)
{
  int i;
  float dis=0.0,dmin,dmax;
  float halfPeriodLength = (float) (periodLength/2.0);
  
  for(i=0;i<3;++i)
    {
      if(node->baseLoc[i] < pos[i] && pos[i] < node->baseLoc[i] + node->sideLengths[i])
	continue;
      else
	{
	  dmin = (float) (fabs(node->baseLoc[i] - pos[i]));
	  if(dmin > halfPeriodLength)
	    dmin = periodLength - dmin;
	  
	  dmax = (float) (fabs(node->baseLoc[i] + node->sideLengths[i] - pos[i]));
	  if(dmax > halfPeriodLength)
	    dmax = periodLength - dmax;
	  	  
	  if(dmax < dmin)
	    dis += dmax*dmax;
	  else
	    dis += dmin*dmin;
	}
    }
    
  if(dis <= rad2)
    return 1;
  else
    return 0;
}

int getNNbrsRad(float centPos[3], float searchRadius, float **nbrsRad2, int **nbrsInd, int *maxNumNbrs, 
		float *px, float *py, float *pz, kdTreeData *td)
{
  int partInd,nodeInd,NumNbrs,*tmpip;
  float dx,dy,dz,searchRadius2,rad2,*tmpfp;
  
  //init
  searchRadius2 = searchRadius*searchRadius;
  NumNbrs = 0;
  nodeInd = 0;
 
  //do search - walk linked list of nodes always going down first if node needs to be opened, else go over
  while(nodeInd >= 0)
    {
      if(intersectSpherekdTreeNode(td->treeNodes+nodeInd,centPos,searchRadius2))
	{
	  if(td->treeNodes[nodeInd].down >= 0 && td->treeNodes[nodeInd].startPart == -1)
	    nodeInd = td->treeNodes[nodeInd].down;
	  else if(td->treeNodes[nodeInd].startPart >= 0 && td->treeNodes[nodeInd].down == -1)
	    {
	      partInd = td->treeNodes[nodeInd].startPart;
	      while(partInd >= 0)
		{
		  dx = (float) (px[partInd] - centPos[0]);
		  dy = (float) (py[partInd] - centPos[1]);
		  dz = (float) (pz[partInd] - centPos[2]);
                      
		  rad2 = dx*dx + dy*dy + dz*dz;
		  if(rad2 <= searchRadius2)
		    {
		      (*nbrsRad2)[NumNbrs] = (float) (rad2);
		      (*nbrsInd)[NumNbrs] = partInd;
		      ++NumNbrs;
                          
		      if(NumNbrs >= (*maxNumNbrs))
			{
			  tmpfp = (float*)realloc(*nbrsRad2,sizeof(float)*((*maxNumNbrs) + 10000));
			  assert(tmpfp != NULL);
			  *nbrsRad2 = tmpfp;
                              
			  tmpip = (int*)realloc(*nbrsInd,sizeof(int)*((*maxNumNbrs) + 10000));
			  assert(tmpip != NULL);
			  *nbrsInd = tmpip;
                              
			  (*maxNumNbrs) = (*maxNumNbrs) + 10000;
			}
		    }
                      
		  //move to next ind in linked list
		  partInd = td->treeParts[partInd];
		}
	      
	      //move to next node
	      nodeInd = td->treeNodes[nodeInd].over;
	    }
	  else
	    {
	      fprintf(stderr,"weird node - opened it up but could not go down and could not look at parts\n");
	      assert(0);
	    }
	}
      else
	nodeInd = td->treeNodes[nodeInd].over;
    }

  return NumNbrs;
}

int intersectSpherekdTreeNode(kdTreeNode *node, float pos[3], float rad2)
{
  int i;
  float dis = 0.0;

  for(i=0;i<3;++i)
    {
      if(pos[i] < node->baseLoc[i])
	dis += (pos[i] - node->baseLoc[i])*(pos[i] - node->baseLoc[i]);
      else if(pos[i] > node->baseLoc[i] + node->sideLengths[i])
	dis += (pos[i] - node->baseLoc[i] - node->sideLengths[i])*(pos[i] - node->baseLoc[i] - node->sideLengths[i]);
    }
  
  if(dis <= rad2)
    return 1;
  else
    return 0;
}

#define TIMEKDTREE

#ifdef TIMEKDTREE
#include <sys/time.h>
static double wtime(void)
{
  struct timeval tp;
  gettimeofday(&tp,NULL);
  return ((double) (tp.tv_sec)) + ((double) (tp.tv_usec))/1e6;
}
extern int ThisTask;
#endif

kdTreeData *buildkdTree(float *px, float *py, float *pz, int Np, float (*dlims)[2])
{
  int leftNode,rightNode,i,ind,nextInd,NpRight,NpLeft;
  kdTreeData *td;
  kdTreeNode *tmpkdTreeNode;
  int *tmpInt;
  int RefinePartLim;
  int splitDim;
  int NumPointsChildren[2];
  float splitVal,*pos[3];
  int *nodeStack;
  int NumNodesInStack,NumNodeStack,currNode;
  int splitThisNode;
  float minPos[3],maxPos[3];
#ifdef TIMEKDTREE
  double totTime,getSplitDimTime,splitTime;
#endif
  
#ifdef TIMEKDTREE
  totTime = -wtime();
  getSplitDimTime = 0.0;
  splitTime = 0.0;
    
#ifdef DEBUG
  if(ThisTask == 0)
    fprintf(stderr,"%d: building kd-tree...\n",ThisTask);
#endif
#endif
  
  //alloc main data struct
  td = (kdTreeData*)malloc(sizeof(kdTreeData));
  assert(td != NULL);
  
  //init params
  NumNodeStack = 10000;
  RefinePartLim = KDTREE_REFINEPARTLIM;
  pos[0] = px;
  pos[1] = py;
  pos[2] = pz;
  td->NumTreeNodes = 10000;
  td->NumTreeNodesUsed = 1;
  td->NumParts = Np;
  
  //alloc mem
  td->treeNodes = (kdTreeNode*)malloc(sizeof(kdTreeNode)*td->NumTreeNodes);
  assert(td->treeNodes != NULL);
  td->treeParts = (int*)malloc(sizeof(int)*Np);
  assert(td->treeParts != NULL);
  nodeStack = (int*)malloc(sizeof(int)*NumNodeStack);
  assert(nodeStack != NULL);

  //init base node
  if(dlims == NULL)
    {
      minPos[0] = px[0];
      minPos[1] = py[0];
      minPos[2] = pz[0];
      maxPos[0] = minPos[0];
      maxPos[1] = minPos[1];
      maxPos[2] = minPos[2];
      for(i=0;i<Np;++i)
	{
	  if(px[i] < minPos[0])
	    minPos[0] = px[i];
	  if(px[i] > maxPos[0])
	    maxPos[0] = px[i];
	  
	  if(py[i] < minPos[1])
	    minPos[1] = py[i];
	  if(py[i] > maxPos[1])
	    maxPos[1] = py[i];
	  
	  if(pz[i] < minPos[2])
	    minPos[2] = pz[i];
	  if(pz[i] > maxPos[2])
	    maxPos[2] = pz[i];
	}
    }
  else
    {
      for(i=0;i<3;++i)
	{
	  minPos[i] = dlims[i][0];
	  maxPos[i] = dlims[i][1];
	}
    }
  for(i=0;i<3;++i)
    {
      td->treeNodes[0].baseLoc[i] = minPos[i];
      td->treeNodes[0].sideLengths[i] = maxPos[i] - minPos[i];
    }
  td->treeNodes[0].startPart = 0;
  td->treeNodes[0].over = -1;
  td->treeNodes[0].down = -1;
  td->treeNodes[0].splitDim = -1;

  //insert all particles into base node
  for(i=0;i<Np;++i)
    td->treeParts[i] = i+1;
  td->treeParts[Np-1] = -1;
  
  if(Np <= RefinePartLim)
    NumNodesInStack = 0;
  else
    NumNodesInStack = 1;
  
  //now refine the tree until each left has no more than RefinePartLim parts
  /*
    uses a stack/queue and while loop to avoid recursive function calls
    the idea is to start with the first node, if it needs to be refined, then it is and its children are added to the stack/queue
    we continually loop through the queue unitl there are no more nodes left in the queue
    nodes are added at the end of the queue and the queue is traversed from the end to the front
  */
  nodeStack[0] = 0;
  currNode = 0;
  while(NumNodesInStack > 0)
    {
      //if(ThisTask == 0)
      //fprintf(stderr,"NumNodesInStack = %d\n",NumNodesInStack);
      
      //grab last node in stack
      currNode = nodeStack[NumNodesInStack-1];
      --NumNodesInStack;
      
      //get currNode splitting
#ifdef TIMEKDTREE
      getSplitDimTime -= wtime();
#endif
      splitThisNode = getSplitNode(currNode,td,&splitDim,&splitVal,NumPointsChildren,px,py,pz);
#ifdef TIMEKDTREE
      getSplitDimTime += wtime();
#endif
      if(splitThisNode)
	{
#ifdef TIMEKDTREE
	  splitTime -= wtime();
#endif
	  
	  //make sure we have enough memory
	  if(td->NumTreeNodesUsed + 2 >= td->NumTreeNodes)
	    {
	      td->NumTreeNodes = td->NumTreeNodes + 10000;
	      tmpkdTreeNode = (kdTreeNode*)realloc(td->treeNodes,sizeof(kdTreeNode)*td->NumTreeNodes);
	      assert(tmpkdTreeNode != NULL);
	      td->treeNodes = tmpkdTreeNode;
	    }
	  
	  //put plitting into parent node
	  td->treeNodes[currNode].splitDim = splitDim;
	  td->treeNodes[currNode].splitVal = splitVal;
	  
	  //alloc new ndoes
	  leftNode = td->NumTreeNodesUsed;
	  ++(td->NumTreeNodesUsed);
	  rightNode = td->NumTreeNodesUsed;
	  ++(td->NumTreeNodesUsed);
	  
	  td->treeNodes[leftNode].splitDim = -1;
	  td->treeNodes[leftNode].splitVal = -1.0;
	  td->treeNodes[rightNode].splitDim = -1;
	  td->treeNodes[rightNode].splitVal = -1.0;
	  
	  //if(ThisTask == 0)
	  //fprintf(stderr,"NumTreeNodesUsed = %d\n",td->NumTreeNodesUsed);
	  
	  //link up nodes
	  // parent node goes down to left node
	  // left node goes over to right node
	  // right node goes over to the over node of the parent node
	  td->treeNodes[currNode].down = leftNode;
          td->treeNodes[leftNode].down = -1;
          td->treeNodes[leftNode].over = rightNode;
          td->treeNodes[rightNode].down = -1;
          td->treeNodes[rightNode].over = td->treeNodes[currNode].over;
	  
	  //set up domain of each node
	  for(i=0;i<3;++i)
	    {
	      if(i == splitDim)
		{
		  td->treeNodes[leftNode].baseLoc[i] = td->treeNodes[currNode].baseLoc[i];
		  td->treeNodes[leftNode].sideLengths[i] = splitVal - td->treeNodes[currNode].baseLoc[i];
		  
		  td->treeNodes[rightNode].sideLengths[i] = td->treeNodes[currNode].sideLengths[i] - td->treeNodes[leftNode].sideLengths[i];
		  td->treeNodes[rightNode].baseLoc[i] = td->treeNodes[currNode].baseLoc[i] + td->treeNodes[leftNode].sideLengths[i];
		}
	      else
		{
		  td->treeNodes[leftNode].baseLoc[i] = td->treeNodes[currNode].baseLoc[i];
		  td->treeNodes[leftNode].sideLengths[i] = td->treeNodes[currNode].sideLengths[i];
		  
		  td->treeNodes[rightNode].baseLoc[i] = td->treeNodes[currNode].baseLoc[i];
		  td->treeNodes[rightNode].sideLengths[i] = td->treeNodes[currNode].sideLengths[i];
		}
	    }
      
	  //assign parts
	  NpRight = 0;
	  NpLeft = 0;
	  td->treeNodes[leftNode].startPart = -1;
	  td->treeNodes[rightNode].startPart = -1;
	  ind = td->treeNodes[currNode].startPart;
	  while(ind >= 0)
	    {
	      //save next link in lsit since it will be changed
	      nextInd = td->treeParts[ind];
	      
	      //do the split on splitDim with splitVal
	      if(pos[splitDim][ind] < splitVal) //left
		{
		  ++NpLeft;
		  if(td->treeNodes[leftNode].startPart < 0)
		    {
		      td->treeNodes[leftNode].startPart = ind;
		      td->treeParts[ind] = -1;
		    }
		  else
		    {
		      td->treeParts[ind] = td->treeNodes[leftNode].startPart;
		      td->treeNodes[leftNode].startPart = ind;
		    }
		}
	      else //right
		{
		  ++NpRight;
		  if(td->treeNodes[rightNode].startPart < 0)
		    {
		      td->treeNodes[rightNode].startPart = ind;
		      td->treeParts[ind] = -1;
		    }
		  else
		    {
		      td->treeParts[ind] = td->treeNodes[rightNode].startPart;
		      td->treeNodes[rightNode].startPart = ind;
		    }
		}
	      
	      //move to next ind
	      ind = nextInd;
	    }
	  td->treeNodes[currNode].startPart = -1;
	  //DEBUG: fprintf(stderr,"NpLeft,Npright = %d|%d (=? %d|%d)\n",NpLeft,NpRight,NumPointsChildren[0],NumPointsChildren[1]);
	  assert(NpLeft == NumPointsChildren[0]);
	  assert(NpRight == NumPointsChildren[1]);
	  
	  //add nodes to stack if they need to be refined - do left node first and make sure we have room in the stack
	  i = 0;
	  if(NpLeft > RefinePartLim)
	    ++i;
	  if(NpRight > RefinePartLim)
	    ++i;
	  if(NumNodesInStack + i >= NumNodeStack)
	    {
	      NumNodeStack += 10000;
	      tmpInt = (int*)realloc(nodeStack,sizeof(int)*NumNodeStack);
	      assert(tmpInt != NULL);
	      nodeStack = tmpInt;
	    }
	  if(NpRight > RefinePartLim)
	    {
	      nodeStack[NumNodesInStack] = rightNode;
	      ++NumNodesInStack;
	    }
	  if(NpLeft > RefinePartLim)
	    {
	      nodeStack[NumNodesInStack] = leftNode;
	      ++NumNodesInStack;
	    }
#ifdef TIMEKDTREE  
	  splitTime += wtime();
#endif
	}
    }
  
  //free extra mem
  free(nodeStack);
  if(td->NumTreeNodesUsed < td->NumTreeNodes)
    {
      td->NumTreeNodes = td->NumTreeNodesUsed;
      tmpkdTreeNode = (kdTreeNode*)realloc(td->treeNodes,sizeof(kdTreeNode)*td->NumTreeNodes);
      assert(tmpkdTreeNode != NULL);
      td->treeNodes = tmpkdTreeNode;
    }

#ifdef TIMEKDTREE  
  totTime += wtime();

#ifdef DEBUG  
  if(ThisTask == 0)
    fprintf(stderr,"%d: built kd-tree: size of tree in MB = %f, tot,getSplitDim,split times = %f|%f|%f [s]\n",ThisTask,((double) (sizeof(kdTreeNode)*td->NumTreeNodes))/1024.0/1024.0,
	    totTime,getSplitDimTime,splitTime);
#endif
#endif
  
  return td;
}

static int getSplitNode(int node, kdTreeData *td, int *splitDim, float *splitVal, int NumPointsChildren[2], float *px, float *py, float *pz)
{
  int retval,i,j,tmp,ind,NumLeft,NumRight,sindex[3],maxRangeInd,maxSideLengthDim;
  float minPos[3],maxPos[3],range[3],maxRange;
  float slidingFactor,slidingVal,minmaxVal;
  float *pos[3],midPnt;
    
  //some error checking
  assert(td->treeNodes[node].startPart >= 0);
      
  pos[0] = px;
  pos[1] = py;
  pos[2] = pz;
  
  //sort max sideLengths
  //DEBUG: fprintf(stderr,"sideLengths = %f|%f|%f\n",
  //DEBUG: td->treeNodes[node].sideLengths[0],td->treeNodes[node].sideLengths[1],td->treeNodes[node].sideLengths[2]);
  sindex[0] = 0;
  sindex[1] = 1;
  sindex[2] = 2;
  for(j=1;j<3;++j)
    {
      tmp = sindex[j];
      i = j;
      while(i > 0 && td->treeNodes[node].sideLengths[sindex[i-1]] < td->treeNodes[node].sideLengths[tmp])
	{
	  sindex[i] = sindex[i-1];
	  --i;
	}
      sindex[i] = tmp;
    }
  //DEBUG: fprintf(stderr,"sorted sideLengths = %f|%f|%f\n",
  //DEBUG: td->treeNodes[node].sideLengths[sindex[0]],td->treeNodes[node].sideLengths[sindex[1]],td->treeNodes[node].sideLengths[sindex[2]]);
  //DEBUG: fprintf(stderr,"sindex = %d|%d|%d\n",sindex[0],sindex[1],sindex[2]);
  
  //if all sideLengths are the same or top two are the same, choose the one with the maximum range
  if(td->treeNodes[node].sideLengths[sindex[0]] == td->treeNodes[node].sideLengths[sindex[1]])
    {
      //get ranges
      for(i=0;i<3;++i)
	{
	  ind = td->treeNodes[node].startPart;
	  minPos[i] = pos[i][ind];
	  maxPos[i] = pos[i][ind];
	  while(ind >= 0)
	    {
	      if(pos[i][ind] > maxPos[i])
		maxPos[i] = pos[i][ind];
	      if(pos[i][ind] < minPos[i])
		minPos[i] = pos[i][ind];
	      
	      ind = td->treeParts[ind];
	    }
	  range[i] = maxPos[i] - minPos[i];
	}
      //DEBUG: fprintf(stderr,"ranges = %f|%f|%f\n",range[0],range[1],range[2]);
      
      //if all three agree, take max range
      if(td->treeNodes[node].sideLengths[sindex[0]] == td->treeNodes[node].sideLengths[sindex[2]])
	{
	  maxRange = range[0];
	  maxRangeInd = 0;
	  for(i=1;i<3;++i)
	    if(range[i] > maxRange)
	      {
		maxRange = range[i];
		maxRangeInd = i;
	      }
	}
      else //just top two agree, take maxof top two
	{
	  if(range[sindex[0]] > range[sindex[1]])
	    maxRangeInd = sindex[0];
	  else
	    maxRangeInd = sindex[1];
	}
      
      maxSideLengthDim = maxRangeInd;
    }
  else
    {
      maxSideLengthDim = sindex[0];
    }
  *splitDim = maxSideLengthDim;
    
  //do sliding mid point
  midPnt = (float) (td->treeNodes[node].baseLoc[maxSideLengthDim] + td->treeNodes[node].sideLengths[maxSideLengthDim]/2.0);
  NumLeft = 0;
  NumRight = 0;
  ind = td->treeNodes[node].startPart;
  while(ind >= 0)
    {
      if(pos[maxSideLengthDim][ind] < midPnt)
	++NumLeft;
      else
	++NumRight;
      ind = td->treeParts[ind];
    }
  
  //do sliding part if needed
  if(NumRight == 0 || NumLeft == 0)
    {
      if(NumRight == 0)
	{
	  //get min/max pos
	  ind = td->treeNodes[node].startPart;
	  minmaxVal = pos[maxSideLengthDim][ind];
	  while(ind >= 0)
	    {
	      if(pos[maxSideLengthDim][ind] > minmaxVal)
		minmaxVal = pos[maxSideLengthDim][ind];
	      ind = td->treeParts[ind];
	    }
	  
	  slidingFactor = (float) (-0.001*td->treeNodes[node].sideLengths[maxSideLengthDim]);
	}
      else
	{
	  //get min/max pos
	  ind = td->treeNodes[node].startPart;
	  minmaxVal = pos[maxSideLengthDim][ind];
	  while(ind >= 0)
	    {
	      if(pos[maxSideLengthDim][ind] < minmaxVal)
		minmaxVal = pos[maxSideLengthDim][ind];
	      ind = td->treeParts[ind];
	    }
	  
	  slidingFactor = (float) (0.001*td->treeNodes[node].sideLengths[maxSideLengthDim]);
	}
      
      slidingVal = minmaxVal;
      i = 0;
      while((NumRight == 0 || NumLeft == 0) && i < 20)
	{
	  slidingVal = slidingVal + slidingFactor;
	  NumLeft = 0;
	  NumRight = 0;
	  ind = td->treeNodes[node].startPart;
	  while(ind >= 0)
	    {
	      if(pos[maxSideLengthDim][ind] < slidingVal)
		++NumLeft;
	      else
		++NumRight;
	      ind = td->treeParts[ind];
	    }
	  
	  ++i;
	}
      *splitVal = slidingVal;
      
      if(NumLeft == 0 || NumRight == 0)
	retval = 0;
      else
	retval = 1;
    }
  else
    {
      retval = 1;
      *splitVal = midPnt;
    }
  
  //DEBUG: fprintf(stderr,"retval = %d, splitDim = %d, midPnt = %f, Num L/R = %d|%d, splitVal = %f\n",
  //DEBUG: retval,maxSideLengthDim,midPnt,NumLeft,NumRight,*splitVal);
  
  NumPointsChildren[0] = NumLeft;
  NumPointsChildren[1] = NumRight;
  
  return retval;
}

void destroykdTree(kdTreeData *treeData)
{
  free(treeData->treeNodes);
  free(treeData->treeParts);
  free(treeData);
}

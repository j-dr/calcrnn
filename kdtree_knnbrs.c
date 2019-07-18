#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "kdtree.h"

extern int ThisTask;

int get_knnbrs_kdtree(float centPos[3], int periodic, float periodLength,
                     float *nbrsRad2, int *nbrsInd, int NumNbrs, float *px, float *py, float *pz,
		      kdTreeData *td)
{
  if(periodic)
    return (getKNNbrsPeriodic(centPos,periodLength,nbrsRad2,nbrsInd,NumNbrs,px,py,pz,td));
  else
    return (getKNNbrs(centPos,nbrsRad2,nbrsInd,NumNbrs,px,py,pz,td));
}

/*
  find k nearest neighbors like this
  1) find the leaf tree node which has centPos in it
  2) search up the tree from this point adding neighbors until have NumNeighbors
  2a) While adding neighbors, test if each cell can have closer neighbors by looking at 
      intersection of splitDim splitVal with sphere around centPos of radius equal to distance to furthest neighbor found so far.
  3) when you get to root node, search is done
*/

int getKNNbrsPeriodic(float centPos[3], float periodLength, float *nbrsRad2, int *nbrsInd, int NumNbrs, 
		      float *px, float *py, float *pz, kdTreeData *td)
{
  int partInd,nodeInd,NumNbrsFoundSoFar,i,*tmpInt,maxInd;
  float dx,dy,dz,maxSearchRadius2,rad2;
  float halfPeriodLength;
  int *nodeStack,*searchStack;
  int NumNodesInStack,NumNodeStack,NumSearchStack,NumNodesInSearchStack;
  
  //init
  halfPeriodLength = (float) (periodLength/2.0);
  maxSearchRadius2 = (float) (halfPeriodLength*halfPeriodLength*3.0*1.0001);
      
  NumNodeStack = 10000;
  nodeStack = (int*)malloc(sizeof(int)*NumNodeStack);
  assert(nodeStack != NULL);
  
  NumSearchStack = 10000;
  searchStack = (int*)malloc(sizeof(int)*NumSearchStack);
  assert(nodeStack != NULL);
  
  /*
    walk down the tree
    if node is split, then either go down to left child or go down->over to right child
    keep record of all nodes visited
    will use this to go back up the tree
  */
  nodeInd = 0;
  NumNodesInStack = 0;
  while(1)
    {
      if(td->treeNodes[nodeInd].splitDim >= 0)
	{
	  if(NumNodesInStack >= NumNodeStack)
	    {
	      NumNodeStack += 10000;
	      tmpInt = (int*)realloc(nodeStack,sizeof(int)*NumNodeStack);
	      assert(tmpInt != NULL);
	      nodeStack = tmpInt;
	    }
	  nodeStack[NumNodesInStack] = nodeInd;
	  ++NumNodesInStack;
	  
	  if(centPos[td->treeNodes[nodeInd].splitDim] < td->treeNodes[nodeInd].splitVal)
	    nodeInd = td->treeNodes[nodeInd].down;
	  else
	    nodeInd = td->treeNodes[td->treeNodes[nodeInd].down].over;
	}
      else
	break;
    }
  
  //add last node to searchStack
  NumNodesInSearchStack = 0;
  searchStack[NumNodesInSearchStack] = nodeInd;
  ++NumNodesInSearchStack;
  
  //now go through stack and check each node
  NumNbrsFoundSoFar = 0;
  while(NumNodesInStack > 0 || NumNodesInSearchStack > 0)
    {
      if(NumNodesInSearchStack > 0)
	{
	  nodeInd = searchStack[NumNodesInSearchStack-1];
	  --NumNodesInSearchStack;
	  assert(nodeInd != -1);
	  
	  if(intersectSpherekdTreeNodePeriodic(td->treeNodes+nodeInd,centPos,maxSearchRadius2,periodLength))
	    {
	      if(td->treeNodes[nodeInd].startPart >= 0)
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
		      if(rad2 <= maxSearchRadius2)
			{
			  if(NumNbrsFoundSoFar < NumNbrs) // do not have enough neighbors yet so add to list
			    {
			      nbrsRad2[NumNbrsFoundSoFar] = (float) (rad2);
			      nbrsInd[NumNbrsFoundSoFar] = partInd;
			      ++NumNbrsFoundSoFar;
			    }
			  else //have enough neighbors replace to end of list
			    {
			      nbrsRad2[maxInd] = (float) (rad2);
			      nbrsInd[maxInd] = partInd;
			    }
			  
			  //get maxSearchRadius2
			  if(NumNbrsFoundSoFar == NumNbrs)
			    {
			      maxSearchRadius2 = nbrsRad2[0];
			      maxInd = 0;
			      for(i=1;i<NumNbrs;++i)
				if(nbrsRad2[i] > maxSearchRadius2)
				  {
				    maxInd = i;
				    maxSearchRadius2 = nbrsRad2[i];
				  }
			    }
			}
		      
		      //move to next ind in linked list
		      partInd = td->treeParts[partInd];
		    }
		}
	      else //add left and right children to search stack if needed
		{
		  if(NumNodesInSearchStack + 2 >= NumSearchStack)
		    {
		      NumSearchStack += 10000;
		      tmpInt = (int*)realloc(searchStack,sizeof(int)*NumSearchStack);
		      assert(tmpInt != NULL);
		      searchStack = tmpInt;
		    }
		  
		  nodeInd = td->treeNodes[nodeInd].down;
		  searchStack[NumNodesInSearchStack] = nodeInd;
		  ++NumNodesInSearchStack;
		  
		  nodeInd = td->treeNodes[nodeInd].over;
		  searchStack[NumNodesInSearchStack] = nodeInd;
		  ++NumNodesInSearchStack;
		}
	    }
	}
      else if(NumNodesInStack > 0)
	{
	  //get last node in stack
	  nodeInd = nodeStack[NumNodesInStack-1];
	  --NumNodesInStack;
	  assert(nodeInd != -1);
	  
	  //get node with search point and then check other node
	  if(centPos[td->treeNodes[nodeInd].splitDim] < td->treeNodes[nodeInd].splitVal)
	    nodeInd = td->treeNodes[td->treeNodes[nodeInd].down].over;
	  else
	    nodeInd = td->treeNodes[nodeInd].down;
	  
	  if(intersectSpherekdTreeNodePeriodic(td->treeNodes+nodeInd,centPos,maxSearchRadius2,periodLength))
	    {
	      if(NumNodesInSearchStack >= NumSearchStack)
		{
		  NumSearchStack += 10000;
		  tmpInt = (int*)realloc(searchStack,sizeof(int)*NumSearchStack);
		  assert(tmpInt != NULL);
		  searchStack = tmpInt;
		}
	      searchStack[NumNodesInSearchStack] = nodeInd;
	      ++NumNodesInSearchStack;
	    }
	}
    }
  
  free(searchStack);
  free(nodeStack);
  
  return NumNbrsFoundSoFar;
}


/*
  find k nearest neighbors like this
  1) find the leaf tree node which has centPos in it
  2) search up the tree from this point adding neighbors until have NumNeighbors
  2a) While adding neighbors, test if each cell can have closer neighbors by looking at 
      intersection of splitDim splitVal with sphere around centPos of radius equal to distance to furthest neighbor found so far.
  3) when you get to root node, search is done
*/

int getKNNbrs(float centPos[3], float *nbrsRad2, int *nbrsInd, int NumNbrs, 
	      float *px, float *py, float *pz, kdTreeData *td)
{
  int partInd,nodeInd,NumNbrsFoundSoFar,i,*tmpInt,maxInd;
  float dx,dy,dz,maxSearchRadius2,rad2;
  int *nodeStack,*searchStack;
  int NumNodesInStack,NumNodeStack,NumSearchStack,NumNodesInSearchStack;
  
  //init
  maxSearchRadius2 = (float) (td->treeNodes[0].sideLengths[0]*td->treeNodes[0].sideLengths[0]*1.0001);
  for(i=1;i<3;++i)
    maxSearchRadius2 += (float) (td->treeNodes[0].sideLengths[i]*td->treeNodes[0].sideLengths[i]*1.0001);
    
  NumNodeStack = 10000;
  nodeStack = (int*)malloc(sizeof(int)*NumNodeStack);
  assert(nodeStack != NULL);
  
  NumSearchStack = 10000;
  searchStack = (int*)malloc(sizeof(int)*NumSearchStack);
  assert(nodeStack != NULL);
  
  /*
    walk down the tree
    if node is split, then either go down to left child or go down->over to right child
    keep record of all nodes visited
    will use this to go back up the tree
  */
  nodeInd = 0;
  NumNodesInStack = 0;
  while(1)
    {
      if(td->treeNodes[nodeInd].splitDim >= 0)
	{
	  if(NumNodesInStack >= NumNodeStack)
	    {
	      NumNodeStack += 10000;
	      tmpInt = (int*)realloc(nodeStack,sizeof(int)*NumNodeStack);
	      assert(tmpInt != NULL);
	      nodeStack = tmpInt;
	    }
	  nodeStack[NumNodesInStack] = nodeInd;
	  ++NumNodesInStack;
	  
	  if(centPos[td->treeNodes[nodeInd].splitDim] < td->treeNodes[nodeInd].splitVal)
	    nodeInd = td->treeNodes[nodeInd].down;
	  else
	    nodeInd = td->treeNodes[td->treeNodes[nodeInd].down].over;
	}
      else
	break;
    }
  
  //add last node to searchStack
  NumNodesInSearchStack = 0;
  searchStack[NumNodesInSearchStack] = nodeInd;
  ++NumNodesInSearchStack;
  
  //now go through stack and check each node
  NumNbrsFoundSoFar = 0;
  while(NumNodesInStack > 0 || NumNodesInSearchStack > 0)
    {
      if(NumNodesInSearchStack > 0)
	{
	  nodeInd = searchStack[NumNodesInSearchStack-1];
	  --NumNodesInSearchStack;
	  assert(nodeInd != -1);
	  
	  if(intersectSpherekdTreeNode(td->treeNodes+nodeInd,centPos,maxSearchRadius2))
	    {
	      if(td->treeNodes[nodeInd].startPart >= 0)
		{
		  partInd = td->treeNodes[nodeInd].startPart;
		  while(partInd >= 0)
		    {
		      dx = (float) (fabs(px[partInd] - centPos[0]));
		      dy = (float) (fabs(py[partInd] - centPos[1]));
		      dz = (float) (fabs(pz[partInd] - centPos[2]));
		      
		      rad2 = dx*dx + dy*dy + dz*dz;
		      if(rad2 <= maxSearchRadius2)
			{
			  if(NumNbrsFoundSoFar < NumNbrs) // do not have enough neighbors yet so add to list
			    {
			      nbrsRad2[NumNbrsFoundSoFar] = (float) (rad2);
			      nbrsInd[NumNbrsFoundSoFar] = partInd;
			      ++NumNbrsFoundSoFar;
			    }
			  else //have enough neighbors replace to end of list
			    {
			      nbrsRad2[maxInd] = (float) (rad2);
			      nbrsInd[maxInd] = partInd;
			    }
			  
			  //get maxSearchRadius2
			  if(NumNbrsFoundSoFar == NumNbrs)
			    {
			      maxSearchRadius2 = nbrsRad2[0];
                              maxInd = 0;
                              for(i=1;i<NumNbrs;++i)
                                if(nbrsRad2[i] > maxSearchRadius2)
                                  {
                                    maxInd = i;
                                    maxSearchRadius2 = nbrsRad2[i];
                                  }
			    }
			}
		      
		      //move to next ind in linked list
		      partInd = td->treeParts[partInd];
		    }
		}
	      else //add left and right children to search stack if needed
		{
		  if(NumNodesInSearchStack + 2 >= NumSearchStack)
		    {
		      NumSearchStack += 10000;
		      tmpInt = (int*)realloc(searchStack,sizeof(int)*NumSearchStack);
		      assert(tmpInt != NULL);
		      searchStack = tmpInt;
		    }
		  
		  nodeInd = td->treeNodes[nodeInd].down;
		  searchStack[NumNodesInSearchStack] = nodeInd;
		  ++NumNodesInSearchStack;
		  		  
		  nodeInd = td->treeNodes[nodeInd].over;
		  searchStack[NumNodesInSearchStack] = nodeInd;
		  ++NumNodesInSearchStack;
		}
	    }
	}
      else if(NumNodesInStack > 0)
	{
	  //get last node in stack
	  nodeInd = nodeStack[NumNodesInStack-1];
	  --NumNodesInStack;
	  assert(nodeInd != -1);
	  
	  //get node with search point and then check other node
	  if(centPos[td->treeNodes[nodeInd].splitDim] < td->treeNodes[nodeInd].splitVal)
	    nodeInd = td->treeNodes[td->treeNodes[nodeInd].down].over;
	  else
	    nodeInd = td->treeNodes[nodeInd].down;
	  
	  if(intersectSpherekdTreeNode(td->treeNodes+nodeInd,centPos,maxSearchRadius2))
	    {
	      if(NumNodesInSearchStack >= NumSearchStack)
		{
		  NumSearchStack += 10000;
		  tmpInt = (int*)realloc(searchStack,sizeof(int)*NumSearchStack);
		  assert(tmpInt != NULL);
		  searchStack = tmpInt;
		}
	      searchStack[NumNodesInSearchStack] = nodeInd;
	      ++NumNodesInSearchStack;
	    }
	}
    }
  
  free(searchStack);
  free(nodeStack);
  
  return NumNbrsFoundSoFar;
}

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#ifdef USEOPENMP
#include <omp.h>
#endif

#include "allheader.h"

static int bboxes_intersect(BBox box1, BBox box2);

void compute_bbox_index(BBox **boxes, int *Nboxes, int Ndiv)
{
  int Nf,Nbx,Nbxin,i,j,Np;
  BBox *bx,*bxin;
  float *px,*py,*pz;
  int startF,endF,dF;
  int NTasksIO;
  int NbxTot,*recvcnts,*displs;
  int Nmod;
  int rbin;
  long temp, fnside, forder, pix;
  
  //alloc boxes                                                                                                                                                                  
  Nf = allData.Nf;
  Nbx = Nf*allData.Ndiv*allData.Ndiv*allData.Ndiv;
  bx = (BBox*)malloc(sizeof(BBox)*Nbx);
  assert(bx != NULL);
  Nbx = 0;
  
  //get this task's range
  if(allData.NumTasksIOInParallel < NTasks)
    NTasksIO = allData.NumTasksIOInParallel;
  else
    NTasksIO = NTasks;
  dF = Nf/NTasksIO;
  startF = dF*ThisTask;
  if(ThisTask == NTasksIO-1)
    endF = Nf-1;
  else
    endF = startF + dF - 1;
  
  //get boxes for each snap
  if(ThisTask < NTasksIO)
    {
      Nmod = (endF-startF+1.0)/10.0;
      if(Nmod == 0)
	Nmod = 1;
      
      for(i=startF;i<=endF;++i)
	{
	  if(ThisTask == 0 && (i-startF+1)%(Nmod) == 0)
            fprintf(stderr,"Bbox: task 0 on file %d of %d (%.2f percent).\n",
                    i-startF+1,endF-startF+1,
                    ((double) (i-startF+1))/((double) (endF-startF+1))*100.0);
#ifdef DEBUG
#if DEBUG_LEVEL > 1
	  fprintf(stderr, "%d: Reading file number %d\n", ThisTask, startF);
#endif
#endif
	  //if no particles in file, continue to next file
	  Np = read_num_parts_file(allData.fnames[i]);
	  if (Np==0)
	    {
	      continue;
	    }
	  
	  //read the parts
	  readparts_file(allData.fnames[i],&px,&py,&pz,&Np);

	  if(strcmp(allData.SimulationType,"LGADGETBCCLC") == 0)
	    {
	      char tempfname[MAX_FILENAME];
	      char *tokens[2];

	      strncpy( tempfname, allData.fnames[i], strlen(allData.fnames[i]) );
	      tempfname[strlen(allData.fnames[i])] = '\0';
	      get_fnside_LGADGETBCCLC(allData.fnames[i], &fnside, &Np);

	      //Get pixel value and radial bin
	      char *token = strtok(tempfname, "_");
	      tokens[0] = malloc(strlen(token)+1);

	      while (token != NULL)
		{
		  tokens[1] = malloc(strlen(tokens[0])+1);
		  strcpy( tokens[1], tokens[0] );
		  tokens[0] = malloc(strlen(token)+1);
		  strcpy( tokens[0], token );
		  token = strtok(NULL, "_");
		}
	      pix = atoi(tokens[0]);
	      rbin = atoi(tokens[1]);
	      //Get file order
	      forder = 0;
	      temp = fnside;
	      while( temp >>= 1 ) forder+=1;
	      fprintf(stderr,"%d: Pix = %d, rbin = %d, forder = %d \n",ThisTask,pix,rbin,forder);
	    }    
	  //compute bboxes
	  bxin = compute_bboxes(Ndiv,px,py,pz,Np,&Nbxin,allData.BoxLength);
	  free(px);
	  free(py);
	  free(pz);
	  
	  for(j=0;j<Nbxin;++j)
	    {
	      bxin[j].fnum = i;
	      if (strcmp(allData.SimulationType,"LGADGETBCCLC") == 0)
		{
		  bxin[j].rbin = rbin;
		  bxin[j].hpix = pix;
		  bxin[j].forder = forder;
		}
	    }
	  
	  memcpy(bx+Nbx,bxin,sizeof(BBox)*Nbxin);
	  Nbx += Nbxin;
	  free(bxin);
	}
      
      //realloc to save mem
      if (Nbx==0) {
	free(bx);
	bx = NULL;
      } else {
	bxin = (BBox*)realloc(bx,sizeof(BBox)*Nbx);
	assert(bxin != NULL);
	bx = bxin;
      }
    }
  else
    {
      Nbx = 0;
      free(bx);
      bx = NULL;
    }
  
  //get recv counts from each task and displs
  recvcnts = (int*)malloc(sizeof(int)*NTasks);
  assert(recvcnts != NULL);
  displs = (int*)malloc(sizeof(int)*NTasks);
  assert(displs != NULL);
  
  MPI_Allreduce(&Nbx,&NbxTot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); //get total number of boxes
  
  MPI_Gather(&Nbx,1,MPI_INT,recvcnts,1,MPI_INT,0,MPI_COMM_WORLD); 
  for(i=0;i<NTasks;++i)
    recvcnts[i] *= sizeof(BBox);
  displs[0] = 0;
  for(i=1;i<NTasks;++i)
    displs[i] = displs[i-1] + recvcnts[i-1];
  
  //alloc mem for all boxes
  bxin = (BBox*)malloc(sizeof(BBox)*NbxTot);
  assert(bxin != NULL);
  
  //gather all boxes from tasks to root
  MPI_Gatherv(bx,(int) (sizeof(BBox)*Nbx),MPI_BYTE, 
	      bxin,recvcnts,displs,MPI_BYTE,0,MPI_COMM_WORLD);
  
  //broadcast the total box array back to everyone
  MPI_Bcast(bxin,(int) (sizeof(BBox)*NbxTot),MPI_BYTE,0,MPI_COMM_WORLD);
  
  //clean up and exit
  if(bx != NULL)
    free(bx);
  *boxes = bxin;
  *Nboxes = NbxTot;
}

//define this to fprintf debug this routine
//#define PRINT_TEST_BBOX_DECOMP

void decompose_bbox_periodic(BBox baseBox, BBox **keepBoxes, int *NkeepBoxes, int *NkeepBoxesAlloc, float L)
{
  BBox *testBoxes,*tmpBox;
  int NtestBoxes,NtestBoxesAlloc;
  BBox tbox,bx1,bx2;
  int i,axis;
  int itr = 0;
  
  NtestBoxesAlloc = 10;
  testBoxes = (BBox*)malloc(sizeof(BBox)*NtestBoxesAlloc);
  assert(testBoxes != NULL);
  testBoxes[0] = baseBox;
  NtestBoxes = 1;
  
  *NkeepBoxes = 0;
  if((*keepBoxes) == NULL || (*NkeepBoxesAlloc) == 0)
    {
      *NkeepBoxesAlloc = 10;
      *keepBoxes = (BBox*)malloc(sizeof(BBox)*(*NkeepBoxesAlloc));
      assert((*keepBoxes) != NULL);
    }
  
  while(NtestBoxes > 0)
    {
#ifdef PRINT_TEST_BBOX_DECOMP
      fprintf(stderr,"NtestBoxes = %d, NkeepBoxes = %d\n",NtestBoxes,*NkeepBoxes);
#endif
      
      //get box to test
      tbox = testBoxes[NtestBoxes-1];
      --NtestBoxes;
      
      //wrap the test box as close to the primary image as possible
#ifdef PRINT_TEST_BBOX_DECOMP      
      fprintf(stderr,"before wrap tbox range = %f|%f, %f|%f, %f|%f\n",tbox.bmin[0],tbox.bmax[0],tbox.bmin[1],tbox.bmax[1],tbox.bmin[2],tbox.bmax[2]);
#endif
      
      for(i=0;i<3;++i)
	{
	  while(tbox.bmax[i] <= 0.0)
	    {
	      tbox.bmin[i] += L;
	      tbox.bmax[i] += L;
	    }
	  
	  while(tbox.bmax[i] > L)
	    {
	      tbox.bmin[i] -= L;
	      tbox.bmax[i] -= L;
	    }
	}
#ifdef PRINT_TEST_BBOX_DECOMP
      fprintf(stderr,"after wrap tbox range = %f|%f, %f|%f, %f|%f\n",tbox.bmin[0],tbox.bmax[0],tbox.bmin[1],tbox.bmax[1],tbox.bmin[2],tbox.bmax[2]);
#endif
      
      //if box is completely in primary image, add to keepBoxes
      if(tbox.bmin[0] >= 0.0 && tbox.bmax[0] <= L &&
	 tbox.bmin[1] >= 0.0 && tbox.bmax[1] <= L &&
	 tbox.bmin[2] >= 0.0 && tbox.bmax[2] <= L)
	{
	  //make sure have mem
	  if((*NkeepBoxes) + 1 > (*NkeepBoxesAlloc))
	    {
	      (*NkeepBoxesAlloc) += 10;
	      tmpBox = (BBox*)realloc(*keepBoxes,sizeof(BBox)*(*NkeepBoxesAlloc));
	      assert(tmpBox != NULL);
	      *keepBoxes = tmpBox;
	    }
	  
	  //add the box
	  (*keepBoxes)[(*NkeepBoxes)] = tbox;
	  (*NkeepBoxes) += 1;
	}
      else //subdivide and add kids to test boxes
	{
#ifdef PRINT_TEST_BBOX_DECOMP
	  fprintf(stderr,"subdividing the box!\n");
#endif
	  
	  //find first axis which violates criterion
	  axis = -1;
	  for(i=0;i<3;++i)
	    {
	      if(tbox.bmin[i] < 0.0)
		{
		  axis = i;
		  break;
		}
	    }
	  assert(axis >= 0 && axis <= 2);

#ifdef PRINT_TEST_BBOX_DECOMP	  
	  fprintf(stderr,"subdiv axis = %d, tbox range = %f|%f\n",axis,tbox.bmin[axis],tbox.bmax[axis]);
#endif
	  
	  //split on first axis which violates criterion
	  bx1 = tbox;
	  bx2 = tbox;
	  bx1.bmax[axis] = 0.0; //bx1 goes to the left
	  bx2.bmin[axis] = 0.0; //bx2 goes to the right

#ifdef PRINT_TEST_BBOX_DECOMP
	  fprintf(stderr,"bx1 range = %f|%f, %f|%f, %f|%f\n",bx1.bmin[0],bx1.bmax[0],bx1.bmin[1],bx1.bmax[1],bx1.bmin[2],bx1.bmax[2]);
	  fprintf(stderr,"bx2 range = %f|%f, %f|%f, %f|%f\n",bx2.bmin[0],bx2.bmax[0],bx2.bmin[1],bx2.bmax[1],bx2.bmin[2],bx2.bmax[2]);
#endif
	  
	  //add to test boxes
	  if(NtestBoxes + 2 > NtestBoxesAlloc)
            {
              NtestBoxesAlloc += 10;
              tmpBox = (BBox*)realloc(testBoxes,sizeof(BBox)*NtestBoxesAlloc);
              assert(tmpBox != NULL);
              testBoxes = tmpBox;
            }

          //add the boxes
          testBoxes[NtestBoxes] = bx1;
	  ++NtestBoxes;
	  testBoxes[NtestBoxes] = bx2;
          ++NtestBoxes;
	  
#ifdef PRINT_TEST_BBOX_DECOMP
	  fprintf(stderr,"Ntestboxes = %d\n",NtestBoxes);
#endif
	}
      
#ifdef PRINT_TEST_BBOX_DECOMP
      fprintf(stderr,"\n\n");
#endif
      
      ++itr;
    }
  
#ifdef PRINT_TEST_BBOX_DECOMP
  fprintf(stderr,"NtestBoxes = %d, NkeepBoxes = %d\n",NtestBoxes,*NkeepBoxes);
  exit(1);
#endif
  
  //clean up
  free(testBoxes);
}

int test_point_bboxes_general(BBox box, float px, float py, float pz, float L, float buff)
{
  int stat = 1;
  int i;
  float pos[3];
  
  //add buffers
  for(i=0;i<3;++i)
    {
      box.bmin[i] -= buff;
      box.bmax[i] += buff;
    }
  
  pos[0] = px;
  pos[1] = py;
  pos[2] = pz;
  
  if(L <= 0.0) //non-periodic case
    {
      for(i=0;i<3;++i)
	{
	  if(!(box.bmin[i] <= pos[i] && pos[i] <= box.bmax[i]))
	    {
	      stat = 0;
	      break;
	    }
	}
      
      return stat;
    }
  else //periodic case
    {
      //test if periodic box with buffers is in primary image and properly normalized
      // if yes, then it is just like non-periodic case
      //else decompose box into overlaps with primary image and then test like non-periodic case
      if(!(box.bmin[0] < 0.0 || box.bmax[0] > L ||
	   box.bmin[1] < 0.0 || box.bmax[1] > L ||
	   box.bmin[2] < 0.0 || box.bmax[2] > L))
	{
	  for(i=0;i<3;++i)
	    {
	      if(!(box.bmin[i] <= pos[i] && pos[i] <= box.bmax[i]))
		{
		  stat = 0;
		  break;
		}
	    }

	  return stat;
	}
      else //decompose box into overlaps with primary image and then test like non-periodic case 
	{
	  BBox *keepBoxes = NULL;
	  int NkeepBoxes,NkeepBoxesAlloc = 0;
	  int j;
	  
	  decompose_bbox_periodic(box,&keepBoxes,&NkeepBoxes,&NkeepBoxesAlloc,L);
	  
	  for(j=0;j<NkeepBoxes;++j)
	    {
	      //assume it intersects
	      stat = 1;
	      for(i=0;i<3;++i)
		{
		  if(!(keepBoxes[j].bmin[i] <= pos[i] && pos[i] <= keepBoxes[j].bmax[i]))
		    {
		      stat = 0; //does not actually, so break
		      break;
		    }
		}
	      
	      if(stat)  //if it intersects, then just return
		{
		  free(keepBoxes);
		  return stat;
		}
	    }
	  
	  //if we get here, never found an intersection, so just return 0
	  free(keepBoxes);
	  return 0;
	}
    }
}

int intersect_bboxes_general(BBox box1, BBox box2, float L, float buff1, float buff2)
{
  int i;
  
  for(i=0;i<3;++i)
    {
      box1.bmin[i] -= buff1;
      box1.bmax[i] += buff1;
      
      box2.bmin[i] -= buff2;
      box2.bmax[i] += buff2;
    }
  
  if(L <= 0.0) //non-periodic
    {
      return bboxes_intersect(box1,box2);
    }
  else //periodic
    {
      //if properly normalized, then no need to decompose
      if(!(box1.bmin[0] < 0.0 || box1.bmax[0] > L ||
           box1.bmin[1] < 0.0 || box1.bmax[1] > L ||
           box1.bmin[2] < 0.0 || box1.bmax[2] > L)
	 &&
	 !(box2.bmin[0] < 0.0 || box2.bmax[0] > L ||
           box2.bmin[1] < 0.0 || box2.bmax[1] > L ||
           box2.bmin[2] < 0.0 || box2.bmax[2] > L))
	{
	  return bboxes_intersect(box1,box2);
	}
      else //decompose boxes first into primary image, then intersect
	{
	  BBox *keepBoxes1 = NULL;
	  int NkeepBoxes1,NkeepBoxesAlloc1 = 0;
	  BBox *keepBoxes2 = NULL;
	  int NkeepBoxes2,NkeepBoxesAlloc2 = 0;
	  int j;
	  
	  //decompose into boxes in primary image
	  decompose_bbox_periodic(box1,&keepBoxes1,&NkeepBoxes1,&NkeepBoxesAlloc1,L);
	  decompose_bbox_periodic(box2,&keepBoxes2,&NkeepBoxes2,&NkeepBoxesAlloc2,L);
	  
	  //test all boxes
	  for(i=0;i<NkeepBoxes1;++i)
	    for(j=0;j<NkeepBoxes2;++j)
	      if(bboxes_intersect(keepBoxes1[i],keepBoxes2[j]))
		{
		  free(keepBoxes1);
		  free(keepBoxes2);
		  return 1;
		}
	  
	  //if we get here, never found an intersection
	  free(keepBoxes1);
	  free(keepBoxes2);
	  return 0;
	}
    }
}

static int bboxes_intersect(BBox box1, BBox box2)
{
  int stat = 1;
  int i;

  for(i=0;i<3;++i)
    if(!((box1.bmin[i] <= box2.bmin[i] && box2.bmin[i] <= box1.bmax[i]) ||
	 (box1.bmin[i] <= box2.bmax[i] && box2.bmax[i] <= box1.bmax[i]) ||
	 (box2.bmin[i] <= box1.bmin[i] && box1.bmin[i] <= box2.bmax[i]) ||
	 (box2.bmin[i] <= box1.bmax[i] && box1.bmax[i] <= box2.bmax[i])))
      {
	stat = 0;
	break;
      }

  return stat;
}

void write_bboxes(char fname[], BBox *boxes, int Nb)
{
  FILE *fp;
  int i;
  
  fp = fopen_retry(fname,"w");
  fprintf(fp,"# fileNumber boxIndex xmin ymin zmin xmax ymax zmax numParts\n");
  for(i=0;i<Nb;++i)
    fprintf(fp,"%d %d %e %e %e %e %e %e %d\n",boxes[i].fnum,boxes[i].ind,
	    boxes[i].bmin[0],boxes[i].bmin[1],boxes[i].bmin[2],
            boxes[i].bmax[0],boxes[i].bmax[1],boxes[i].bmax[2],boxes[i].np);
  fclose(fp);
}

BBox *compute_bboxes(int Ndiv, float *px, float *py, float *pz, int Np, int *Nb, float L)
{
  float rmin[3],rmax[3],dr[3];
  BBox *boxes,*tmpbx;
  int i,Nbt,xi,yi,zi,ind;
  
  //prep BBoxes
  Nbt = Ndiv*Ndiv*Ndiv;
  boxes = (BBox*)malloc(sizeof(BBox)*Nbt);
  assert(boxes != NULL);
  for(i=0;i<Nbt;++i)
    boxes[i].ind = 0;
  
  //prep range
  rmin[0] = px[0];
  rmin[1] = py[0];
  rmin[2] = pz[0];
  rmax[0] = px[0];
  rmax[1] = py[0];
  rmax[2] = pz[0];
  
  for(i=0;i<Np;++i)
    {
      //min
      if(px[i] < rmin[0])
	rmin[0] = px[i];
      if(py[i] < rmin[1])
	rmin[1] = py[i];
      if(pz[i] < rmin[2])
	rmin[2] = pz[i];
      
      //max
      if(px[i] > rmax[0])
	rmax[0] = px[i];
      if(py[i] > rmax[1])
	rmax[1] = py[i];
      if(pz[i] > rmax[2])
	rmax[2] = pz[i];
    }
  for(i=0;i<3;++i)
    dr[i] = (rmax[i]-rmin[i])/Ndiv;
  
#ifdef DEBUG
  fprintf(stderr,"min,max,dr = (%f|%f|%f) (%f|%f|%f) (%f|%f|%f)\n",rmin[0],rmax[0],dr[0],rmin[1],rmax[1],dr[1],rmin[2],rmax[2],dr[2]);
#endif
  
  for(i=0;i<Np;++i)
    {
      xi = (px[i] - rmin[0])/dr[0];
      if(xi < 0)
	xi = 0;
      if(xi >= Ndiv)
	xi = Ndiv-1;
      
      yi = (py[i] - rmin[1])/dr[1];
      if(yi < 0)
	yi = 0;
      if(yi >= Ndiv)
	yi = Ndiv-1;
      
      zi = (pz[i] - rmin[2])/dr[2];
      if(zi < 0)
	zi = 0;
      if(zi >= Ndiv)
	zi = Ndiv-1;
      
      ind = (xi*Ndiv + yi)*Ndiv + zi;
      
      if(boxes[ind].ind == 0)
	{
	  boxes[ind].np = 1;
	  
	  boxes[ind].bmin[0] = px[i];
	  boxes[ind].bmin[1] = py[i];
	  boxes[ind].bmin[2] = pz[i];
	  
	  boxes[ind].bmax[0] = px[i];
	  boxes[ind].bmax[1] = py[i];
	  boxes[ind].bmax[2] = pz[i];
	  
	  boxes[ind].ind = 1;
	}
      else
	{
	  ++(boxes[ind].np);
	  
	  if(boxes[ind].bmin[0] > px[i])
	    boxes[ind].bmin[0] = px[i];
	  if(boxes[ind].bmin[1] > py[i])
	    boxes[ind].bmin[1] = py[i];
	  if(boxes[ind].bmin[2] > pz[i])
	    boxes[ind].bmin[2] = pz[i];
	  
	  if(boxes[ind].bmax[0] < px[i])
	    boxes[ind].bmax[0] = px[i];
	  if(boxes[ind].bmax[1] < py[i])
	    boxes[ind].bmax[1] = py[i];
	  if(boxes[ind].bmax[2] < pz[i])
	    boxes[ind].bmax[2] = pz[i];
	}
    }
  
  //get boxes with parts in them
  *Nb = 0;
  for(i=0;i<Nbt;++i)
    if(boxes[i].ind == 1)
      {
	boxes[(*Nb)] = boxes[i];
	boxes[(*Nb)].ind = (*Nb);
	++(*Nb);
      }
  
  tmpbx = (BBox*)realloc(boxes,sizeof(BBox)*(*Nb));
  assert(tmpbx != NULL);
  boxes = tmpbx;
  Nbt = *Nb;
  
  //if total volume of overal bbox is within some fraction of summed volume of all bboxes
  //  then just use a single box
  float totvol,sumvol;
  totvol = 
    (rmax[0]-rmin[0])*
    (rmax[1]-rmin[1])*
    (rmax[2]-rmin[2]);
  
  sumvol = 0.0;
  for(i=0;i<Nbt;++i)
    sumvol += vol_bbox(boxes[i]);
  //fprintf(stdout,"Percent volume difference between subdivided and original bboxes: %g",fabs(totvol/sumvol-1.0));

  if(fabs(totvol/sumvol-1.0) < 0.1) //no more than 20% difference in volume
    {
      for(i=0;i<3;++i)
	{
	  boxes[0].bmin[i] = rmin[i];
	  boxes[0].bmax[i] = rmax[i];
	}
      boxes[0].np = Np;
      boxes[0].ind = 0;
      
      *Nb = 1;
      tmpbx = (BBox*)realloc(boxes,sizeof(BBox)*(*Nb));
      assert(tmpbx != NULL);
      boxes = tmpbx;
      Nbt = *Nb;
    }

  //error check boxes for periodic case
  if(L > 0.0)
    {
      for(i=0;i<Nbt;++i)
	if(!(boxes[i].bmin[0] >= 0.0 && boxes[i].bmax[0] <= L &&
	     boxes[i].bmin[1] >= 0.0 && boxes[i].bmax[1] <= L &&
	     boxes[i].bmin[2] >= 0.0 && boxes[i].bmax[2] <= L))
	  {
	    fprintf(stderr,"!!!!!!!!!!!! All particles must be in primary periodic image from [0,L]! !!!!!!!!!!!!\n");
	    MPI_Abort(MPI_COMM_WORLD,999);
	  }
    }
  
  return boxes;
}

float vol_bbox(BBox bx)
{
  return 
    (bx.bmax[0]-bx.bmin[0])*
    (bx.bmax[1]-bx.bmin[1])*
    (bx.bmax[2]-bx.bmin[2]);
}

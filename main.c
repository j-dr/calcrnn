#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort_long.h>

#ifdef USEOPENMP
#include <omp.h>
#endif

#include "allheader.h"

//#define DEBUG_MSCOMM

int main(int argc, char **argv)
{
  //init MPI and get current tasks and number of tasks
  /*#ifdef USEOPENMP
  int provided;  
  int rc = MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  if(rc != MPI_SUCCESS || provided != MPI_THREAD_FUNNELED)
    {
      fprintf(stderr,"Error starting MPI program. Terminating.\n");
      MPI_Abort(MPI_COMM_WORLD,rc);
    }    
    #else*/
  int rc = MPI_Init(&argc,&argv);
  if(rc != MPI_SUCCESS)  
    {
      fprintf(stderr,"Error starting MPI program. Terminating.\n");
      MPI_Abort(MPI_COMM_WORLD,rc);
    }
  //#endif
  
  MPI_Comm_size(MPI_COMM_WORLD,&NTasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&ThisTask);

  
  //test for proper number of threads
#ifdef USEOPENMP
  omp_set_num_threads(atoi(argv[2]));
  if(ThisTask == 0)
    {
      fprintf(stderr,"%d: # of tasks = %d, # of threads = %d, requested %d threads\n",ThisTask,NTasks,omp_get_max_threads(),atoi(argv[2]));
      fflush(stderr);
#pragma omp parallel
      {
	fprintf(stderr,"%d: thread %d\n",ThisTask,omp_get_thread_num());
	fflush(stderr);
      }
    }
  fflush(stderr);
#endif
  
#ifdef MEMWATCH
  mwDoFlush(1);
#endif
  
  BBox *boxes;
  int Nb;
  double time = 0,tb = 0.0,tr = 0.0,th = 0.0;
  
  ///////////////////////////////
  //start 
  ///////////////////////////////
  begin_run(argv[1]);
  
  //override config file if wanted
#ifndef USEOPENMP
  if(argc == 3)
    {
      allData.NumTasksIOInParallel = atoi(argv[2]);
      
      assert(allData.NumTasksIOInParallel > 0);
      assert(allData.NumTasksIOInParallel <= NTasks);
    }
#else
  if(argc == 4)
    {
      allData.NumTasksIOInParallel = atoi(argv[3]);
      
      assert(allData.NumTasksIOInParallel > 0);
      assert(allData.NumTasksIOInParallel <= NTasks);
    }
#endif
  
#ifdef TEST_CODE
  test_kdtree_knnbrsfind();
  test_kdtree_nbrsfind();
  
  MPI_Finalize();
  return 0;
#endif
  
  ///////////////////////////////
  // Make the BBox index 
  ///////////////////////////////
  if(ThisTask == 0)
    {
      time = -MPI_Wtime();
      fprintf(stderr,"maiking bboxes.\n");
    }
  
  compute_bbox_index(&boxes,&Nb,allData.Ndiv);
  
  if(ThisTask == 0)
    {
      tb = time + MPI_Wtime();
      fprintf(stderr,"making bboxes took %f seconds.\n\nfinding Rnn.\n",tb);
      
      write_bboxes(allData.BBoxOutputFile,boxes,Nb);
    }
  
  ///////////////////////////////
  // Do Rnn for particles
  ///////////////////////////////
  if(strlen(allData.HaloFile) == 0)
    {
      if(ThisTask == 0)
	tr = -MPI_Wtime();
      
      do_rnn_parts(boxes,Nb);
      
      if(ThisTask == 0)
	{
	  tr += MPI_Wtime();
	  fprintf(stderr,"finding Rnn for parts took %f seconds.\n\n",tr);
	}
      
      ///////////////////////////////
      MPI_Barrier(MPI_COMM_WORLD);
      ///////////////////////////////
    }
  
  ///////////////////////////////
  // Do Rnn for halos if needed
  ///////////////////////////////
  if(strlen(allData.HaloFile) > 0)
    {
      if(ThisTask == 0)
	th -= MPI_Wtime();
      
      do_rnn_halos(boxes,Nb);
      
      if(ThisTask == 0)
	{
	  th += MPI_Wtime();
	  fprintf(stderr,"finding Rnn for halos took %f seconds.\n\n",th);
	}
    }
  
  ///////////////////////////////
  //clean-up and print timing info
  ///////////////////////////////
  free(boxes);
  
  if(ThisTask == 0)
    {
      time += MPI_Wtime();
      fprintf(stderr,"total time is %lf seconds.\n",time);
    }
  
  ///////////////////////////////
  // finish
  ///////////////////////////////
  end_run();
  
  MPI_Finalize();
  return 0;
}


#include "allheader.h"

/* global vars are defined here except those in 
   profile.c
*/

/* extern defs of global vars in globalvars.c */
AllData allData;                              /* global struct with all vars from config file */
int ThisTask;                                 /* this task's rank in MPI_COMM_WORLD */
int NTasks;                                   /* number of tasks in MPI_COMM_WORLD */

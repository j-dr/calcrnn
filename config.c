#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>

#include "allheader.h"

#define LINELENGTH 1024

static int strcmp_caseinsens(const char *s1, const char *s2);

void begin_run(char config[])
{
  int Nf,i;
  FILE *fp;
  char dname[1024];
  
  //read config file and fnames
  if(ThisTask == 0)
    {
      //config
      read_config_file(config);
      fflush(stderr);
      
      //make dir
      sprintf(dname,"%s",allData.OutputPath);
      mkdir(dname,02755);
      
      //fnames
      fp = fopen_retry(allData.SnapshotFileList,"r");
      Nf = fnumlines(fp);
      allData.Nf = Nf;
      
      allData.fnames = (char**)malloc(sizeof(char*)*Nf);
      assert(allData.fnames != NULL);
      allData.fnames[0] = (char*)malloc(sizeof(char)*MAX_FILENAME*Nf);
      assert(allData.fnames[0] != NULL);
      for(i=1;i<Nf;++i)
	allData.fnames[i] = allData.fnames[0] + i*MAX_FILENAME;
      
      //read the file names
      for(i=0;i<Nf;++i)
	{
	  fgets(allData.fnames[i],MAX_FILENAME,fp);
	  if(allData.fnames[i][strlen(allData.fnames[i])-1] == '\n')
	    allData.fnames[i][strlen(allData.fnames[i])-1] = '\0';
	  
	  //anything that needs to be read from files once is done here 
	  if(i == 0)
	    {
	      allData.BoxLength = read_period_length_file(allData.fnames[i]);
	    }
	}
      
      fclose(fp);
    }
  
  //broadcast it to the world
  MPI_Bcast(&allData,(int) (sizeof(AllData)),MPI_BYTE,0,MPI_COMM_WORLD);
  Nf = allData.Nf;
  

  
  if(ThisTask != 0)
    {
      //now get names to other tasks
      allData.fnames = (char**)malloc(sizeof(char*)*Nf);
      assert(allData.fnames != NULL);
      allData.fnames[0] = (char*)malloc(sizeof(char)*MAX_FILENAME*Nf);
      assert(allData.fnames[0] != NULL);
      for(i=1;i<Nf;++i)
	allData.fnames[i] = allData.fnames[0] +i*MAX_FILENAME;
    }
  MPI_Bcast(allData.fnames[0],(int) (Nf*MAX_FILENAME),MPI_CHAR,0,MPI_COMM_WORLD);
}

void end_run(void)
{
  //free file names
  free(allData.fnames[0]);
  free(allData.fnames);
}

void read_config_file(char fname[])
{
  FILE *fp;
  char line[LINELENGTH];
  char *p,*tag,*value;
    
  //set some defaults
  allData.OutputPath[0] = '\0';
  allData.HaloFile[0] = '\0';
  allData.HaloFileFormat[0] = '\0';
  allData.BBoxOutputFile[0] = '\0';
  allData.NumTasksIOInParallel = -1;
  allData.DomainBuffSize = -1.0;
  allData.NRnn = -1;
  allData.Ndiv = -1;
  allData.BoxLength = -1.0;
  allData.HaloChunkSizeMB = 0;
  
  //error check the file open
  fp = fopen_retry(fname,"r");
  if(fp == NULL)
    {
      fprintf(stderr,"config file '%s' could not be opened!\n",fname);
      assert(fp != NULL);
    }
  
  while(fgets(line,LINELENGTH,fp) != NULL)
    {
      //convert all tabs and newlines to spaces
#ifdef DEBUG
#if DEBUG_LEVEL > 2
      if(ThisTask == 0)
	fprintf(stderr,"base line: '%s'\n",line);
#endif
#endif
      p = line;
      while(*p != '\0')
	{
	  if(*p == '\n' || *p == '\t')
	    *p = ' ';
	  ++p;
	}
      p = line;
#ifdef DEBUG
#if DEBUG_LEVEL > 2
      if(ThisTask == 0)
	fprintf(stderr,"conv line: '%s'\n",line);
#endif
#endif
      //go past any starting white space
      while(*p == ' ')
	++p;
      
      //if have reached end of string or string line is a comment, then go to next line
      if(*p == '\0' || *p == '#')
	continue;
      
      //now at start of tag - go to end and then null temrinate
      tag = p;
      while(*p != ' ')
	++p;
      *p = '\0';
      ++p;
#ifdef DEBUG
#if DEBUG_LEVEL > 2
      if(ThisTask == 0)
	fprintf(stderr,"tag: '%s'\n",tag);
#endif
#endif

      //skip more white space
      while(*p == ' ')
	++p;
      
      //get value
      value = p;
      while(*p != ' ' && *p != '#')
	++p;
      *p = '\0';
      ++p;
#ifdef DEBUG
#if DEBUG_LEVEL > 2
      if(ThisTask == 0)
	fprintf(stderr,"value: '%s'\n",value);
#endif
#endif

      //put tag value in struct
      if(strcmp_caseinsens(tag,"nrnn") == 0)
	{
	  allData.NRnn = atol(value);
	}
      else if(strcmp_caseinsens(tag,"ndiv") == 0)
	{
	  allData.Ndiv = atol(value);
	}
      else if(strcmp_caseinsens(tag,"simulationtype") == 0)
	{
	  strcpy(allData.SimulationType,value);
	}
      else if(strcmp_caseinsens(tag,"halofile") == 0)
	{
	  strcpy(allData.HaloFile,value);
	}
      else if(strcmp_caseinsens(tag,"halofileformat") == 0)
	{
	  strcpy(allData.HaloFileFormat,value);
	}
      else if(strcmp_caseinsens(tag,"outputpath") == 0)
	{
	  strcpy(allData.OutputPath,value);
	}
      else if(strcmp_caseinsens(tag,"BBoxOutputFile") == 0)
	{
	  strcpy(allData.BBoxOutputFile,value);
	}
      else if(strcmp_caseinsens(tag,"snapshotfilelist") == 0)
	{
	  strcpy(allData.SnapshotFileList,value);
	}
      else if(strcmp_caseinsens(tag,"NumTasksIOInParallel") == 0)
	{
	  allData.NumTasksIOInParallel = atoi(value);
	}
      else if(strcmp_caseinsens(tag,"HaloChunkSizeMB") == 0)
	{
	  allData.HaloChunkSizeMB = atol(value);
	}
      else if(strcmp_caseinsens(tag,"domainbuffsize") == 0)
	{
	  allData.DomainBuffSize = (float) (atof(value));
	}
      else
	{
	  fprintf(stderr,"tag '%s' not found!\n",tag);
	  assert(0);
	}
    }
  
  fclose(fp);
  
  /* error check some stuff */
  assert(strlen(allData.OutputPath) != 0);
  assert(strlen(allData.BBoxOutputFile) != 0);
  assert(allData.NumTasksIOInParallel > 0);
  assert(allData.NumTasksIOInParallel <= NTasks);
  assert(allData.DomainBuffSize > 0.0);
  assert(allData.NRnn > 0);
  assert(allData.Ndiv > 0);
    
  if(strlen(allData.HaloFile) > 0)
    {
      assert(allData.HaloChunkSizeMB > 0.0);
      assert(strlen(allData.HaloFileFormat) > 0);
    }
}

static int strcmp_caseinsens(const char *s1, const char *s2)
{
  int N,i,equal;
  
  if(strlen(s1) != strlen(s2))
    return 1;
  N = strlen(s1);

#ifdef DEBUG
#if DEBUG_LEVEL > 2  
  if(ThisTask == 0)
    fprintf(stderr,"s1 = '%s', s2 = '%s'\n",s1,s2);
#endif
#endif

  equal = 0;
  for(i=0;i<N;++i)
    {
#ifdef DEBUG
#if DEBUG_LEVEL > 2
      if(ThisTask == 0)
	fprintf(stderr,"s1[i] = '%c', s2[i] = '%c'\n",tolower(s1[i]),tolower(s2[i]));
#endif
#endif

      if(tolower(s1[i]) != tolower(s2[i]))
	{
	  equal = 1;
	  break;
	}
    }

#ifdef DEBUG
#if DEBUG_LEVEL > 2  
  if(ThisTask == 0)
    fprintf(stderr,"equal %d (0 is true, 1 is false, weird)\n",equal);
#endif
#endif

  return equal;
}

#undef LINELENGTH

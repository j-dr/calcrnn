#compile time options
#OPTS += -DDEBUG_IO #define for some debuggin I/O - none implemented so far
#OPTS += -DDEBUG -DDEBUG_LEVEL=2 #leave undefined for no debugging - 0,1, and 2 give progressively more output to stderr
#OPTS += -DTEST_CODE #define to run some basic test code
#OPTS += -DMEMWATCH -DMEMWATCH_STDIO #define to test for memory leaks, out of bounds, etc. for memory used in this code
#OPTS += -DDSAMPLEFACTOR=16 #set to a factor by which paticle data will be downsampled
#OPTS += -DUSEOPENMP #define to use openmp 


#select your computer
COMP="cori-haswell"
#COMP="edison"
#COMP="orange"
#COMP="ranger"
#COMP="fulla-gcc"
#COMP="fulla-icc"
#COMP="mandor-icc"
#COMP="mandor-gcc"
#COMP="orion"

#edit/add to match your machine
ifeq ($(COMP),"orion")
CC          =  mpicc
OPTIMIZE    =  -g -O3 #-Werror -Wconversion
OPENMPFLAG  =  -fopenmp
GSLI        =  -I/home/beckermr/include
GSLL        =  -L/home/beckermr/lib -lgsl -lgslcblas
EXTRACFLAGS =  -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
	-Wwrite-strings -Wnested-externs -fshort-enums -fno-common -Dinline= #-ansi -std=c99 -pedantic 
EXTRACLIB   =  
endif

ifeq ($(COMP),"orange")
CC          =  mpicc
OPTIMIZE    =  -g -O3 #-O3 #-Wall -wd981 #-wd1419 -wd810
OPENMPFLAG  =  -fopenmp
GSLI        =  -I/afs/slac.stanford.edu/g/ki/software/gsl/1.15/include
GSLL        =  -L/afs/slac.stanford.edu/g/ki/software/gsl/1.15/lib -lgsl -lgslcblas
EXTRACFLAGS =  -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
	-Wwrite-strings -Wnested-externs -fshort-enums -fno-common -Dinline=
EXTRACLIB   =
endif

ifeq ($(COMP),"edison")
CC          =  cc
OPTIMIZE    =  -g -O3 #-O3 #-Wall -wd981 #-wd1419 -wd810
OPENMPFLAG  =  -fopenmp
GSLI        =  -I${GSL_DIR}/include
GSLL        =  -L${GSL_DIR}/lib -lgsl -lgslcblas
EXTRACFLAGS =  -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
	-Wwrite-strings -Wnested-externs -fshort-enums -fno-common -Dinline=
EXTRACLIB   =
endif

ifeq ($(COMP),"cori-haswell")
CC          =  cc
OPTIMIZE    =  -g -O0 #-O3 #-Wall -wd981 #-wd1419 -wd810
OPENMPFLAG  =  -fopenmp
GSLI        =  -I${GSL_DIR}/include
GSLL        =  -L${GSL_DIR}/lib -lgsl -lgslcblas
EXTRACFLAGS =  -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
	-Wwrite-strings -Wnested-externs -fshort-enums -fno-common -Dinline=
EXTRACLIB   =
endif

ifeq ($(COMP),"ranger")
CC          =  mpicc
OPTIMIZE    =  -g -O3 #-Wall -wd981 #-wd1419 -wd810
OPENMPFLAG  =  -fopenmp
GSLI        =  -I${TACC_GSL_INC}
GSLL        =  -L${TACC_GSL_LIB} -lgsl -lgslcblas
EXTRACFLAGS =  -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
        -Wwrite-strings -Wnested-externs -fshort-enums -fno-common -Dinline=
EXTRACLIB   =
endif

ifeq ($(COMP),"fulla-icc")
CC          =  /usr/local/openmpi-intel/bin/mpicc
OPTIMIZE    =  -O3 -ipo #-static -xW
OPENMPFLAG  = -openmp
GSLI        =  
GSLL        =  -lgsl -lgslcblas
EXTRACFLAGS =  -I/home/beckermr/include -Wall -wd981 
EXTRACLIB   =  -L/home/beckermr/lib
endif

ifeq ($(COMP),"fulla-gcc")
CC          =  /usr/local/openmpi-gcc/bin/mpicc
OPTIMIZE    =  -O3
OPENMPFLAG  = 
GSLI        =  
GSLL        =  -lgsl -lgslcblas
EXTRACFLAGS =  -I/home/beckermr/include -ansi -std=c99 -pedantic -Wall -W -Wmissing-prototypes -Wstrict-prototypes \
	-Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
	-Wwrite-strings -Wnested-externs -fshort-enums -fno-common -Dinline=
EXTRACLIB   =  -L/home/beckermr/lib
endif

ifeq ($(COMP),"mandor-icc")
CC          =  /opt/intel/Compiler/11.1/056/bin/intel64/mpicc 
OPTIMIZE    =  -g -O3
OPENMPFLAG  = 
GSLI        =  -I/home/beckermr/include
GSLL        =  -L/home/beckermr/lib -lgsl -lgslcblas
EXTRACFLAGS =  -Wall -wd981 #-wd1419 -wd810
EXTRACLIB   =  
endif

ifeq ($(COMP),"mandor-gcc")
CC          =  /usr/local/bin/mpicc
OPTIMIZE    =  -g -O3 #-Werror
OPENMPFLAG  =  -fopenmp
GSLI        =  -I/home/beckermr/include
GSLL        =  -L/home/beckermr/lib -lgsl -lgslcblas
EXTRACFLAGS =  -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align \
	-Wwrite-strings -Wnested-externs -fshort-enums -fno-common -Dinline= #-ansi -std=c99 -pedantic 
EXTRACLIB   =  
endif

#set it all up
ifneq (USEOPENMP,$(findstring USEOPENMP,$(OPTS)))
OPENMPFLAG= 
endif

CLINK=$(CC)
CFLAGS=$(OPTIMIZE) $(OPENMPFLAG) $(GSLI) $(FFTWI) $(FITSI) $(EXTRACFLAGS) $(OPTS) 
CLIB=$(EXTRACLIB) $(GSLL) $(FFTWL) $(FITSL) $(GSLL) -lm

ifneq (MEMWATCH,$(findstring MEMWATCH,$(CFLAGS)))
MEMWATCH= 
else
MEMWATCH=memwatch.o
endif

ifneq (TEST_CODE,$(findstring TEST_CODE,$(CFLAGS)))
TESTCODE=
else
TESTCODE=test_code.o
endif

OBJS = $(MEMWATCH) $(TESTCODE) globalvars.o bbox.o config.o io.o rnn_parts.o rnn_halos.o kdtree.o kdtree_knnbrs.o reorder.o \
	read_LGADGET.o read_LGADGETLC.o read_LGADGETBCCLC.o healpix_utils.o utils.o

#config.o globalvars.o partio.o domain.o \
#utils.o kdtree.o readLGADGET.o kdtree_knnbrs.o rnncalc.o

#halopercolation_parallel.o

EXEC = calcrnn
TEST =
all: $(EXEC) 
test: $(TEST) clean

OBJS1=$(OBJS) main.o
$(EXEC): $(OBJS1)
	$(CLINK) $(CFLAGS) -o $@ $(OBJS1) $(CLIB)

$(OBJS1): Makefile allheader.h kdtree.h

.PHONY : clean
clean: 
	rm -f *.o

.PHONY : spotless
spotless: 
	rm -f *.o $(EXEC) $(TEST)

.PHONY : pristine
pristine: 
	rm -f *.o $(EXEC) $(TEST) *~


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>
#include <string.h>
#include <gsl/gsl_math.h>

#include "allheader.h"

double wrap_pos(double pin, double L)
{
  double pout = pin;

  while(pout > L)
    pout -= L;
  while(pout < 0.0)
    pout += L;

  return pout;
}


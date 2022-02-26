/***********************************************************
 *  Copyright Univ. of Texas M.D. Anderson Cancer Center
 *  1992.
 *
 *	Some routines modified from Numerical Recipes in C,
 *	including error report, array or matrix declaration
 *	and releasing.
 ****/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mcml.hpp"

void nrerror(const char error_text[])
{
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}
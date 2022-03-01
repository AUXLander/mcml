#pragma once

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <ctype.h>

#include "io/io.hpp"
#include "mcmlnr.hpp"


#define PI 3.1415926
#define WEIGHT 1E-4		/* Critical weight for roulette. */
#define CHANCE 0.1		/* Chance of roulette survival. */

//#define Boolean char

#define SIGN(x) ((x)>=0 ? 1:-1)

#define _CRT_SECURE_NO_WARNINGS



/***********************************************************
 *	Routine prototypes for dynamic memory allocation and
 *	release of arrays and matrices.
 *	Modified from Numerical Recipes in C.
 ****/

void ShowVersion(const char* version);

/****************** Stuctures *****************************/



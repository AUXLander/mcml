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

/***********************************************************
 *	Routine prototypes for dynamic memory allocation and
 *	release of arrays and matrices.
 *	Modified from Numerical Recipes in C.
 ****/

void ShowVersion(const char* version);

/****************** Stuctures *****************************/

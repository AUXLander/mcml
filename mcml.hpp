#pragma once

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <ctype.h>

#include "mcmlnr.hpp"

#include "io/io.hpp"

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

/****
 *	Structure used to describe a photon packet.
 ****/
struct PhotonStruct
{
  double x, y, z;	/* Cartesian coordinates.[cm] */
  double ux, uy, uz;/* directional cosines of a photon. */
  double w;			/* weight. */
  bool dead;		/* 1 if photon is terminated. */
  short layer;		/* index to layer where the photon */
					/* packet resides. */
  double step_size;	/* current step size. [cm]. */
  double sleft;		/* step size left. dimensionless [-]. */
};

/****
 *	Structure used to describe the geometry and optical
 *	properties of a layer.
 *	z0 and z1 are the z coordinates for the upper boundary
 *	and lower boundary respectively.
 *
 *	cos_crit0 and cos_crit1 are the cosines of the 
 *	critical angle of total internal reflection for the
 *	upper boundary and lower boundary respectively.
 *	They are set to zero if no total internal reflection
 *	exists.
 *	They are used for computation speed.
 ****/
struct LayerStruct 
{
  double z0, z1;	/* z coordinates of a layer. [cm] */
  double n;			/* refractive index of a layer. */
  double mua;	    /* absorption coefficient. [1/cm] */
  double mus;	    /* scattering coefficient. [1/cm] */
  double anisotropy;		    /* anisotropy. */
  
  double cos_crit0,	cos_crit1;	
};

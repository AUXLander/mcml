#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <ctype.h>

#define PI 3.1415926
#define WEIGHT 1E-4		/* Critical weight for roulette. */
#define CHANCE 0.1		/* Chance of roulette survival. */
#define STRLEN 256		/* String length. */

//#define Boolean char

#define SIGN(x) ((x)>=0 ? 1:-1)

#define _CRT_SECURE_NO_WARNINGS

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

/****
 *	Input parameters for each independent run.
 *
 *	z and r are for the cylindrical coordinate system. [cm]
 *	a is for the angle alpha between the photon exiting 
 *	direction and the surface normal. [radian]
 *
 *	The grid line separations in z, r, and alpha
 *	directions are dz, dr, and da respectively.  The numbers 
 *	of grid lines in z, r, and alpha directions are
 *	nz, nr, and na respectively.
 *
 *	The member layerspecs will point to an array of 
 *	structures which store parameters of each layer. 
 *	This array has (number_layers + 2) elements. One
 *	element is for a layer.
 *	The layers 0 and (num_layers + 1) are for top ambient 
 *	medium and the bottom ambient medium respectively.
 ****/
struct InputStruct
{
  char	 out_fname[STRLEN];	/* output file name. */
  char	 out_fformat;		/* output file format. */
							/* 'A' for ASCII, */
							/* 'B' for binary. */
  long	 num_photons; 		/* to be traced. */
  double Wth; 				/* play roulette if photon */
							/* weight < Wth.*/
  
  double dz;				/* z grid separation.[cm] */ 
  double dr;				/* r grid separation.[cm] */
  double da;				/* alpha grid separation. */
							/* [radian] */
  short nz;					/* array range 0..nz-1. */
  short nr;					/* array range 0..nr-1. */
  short na;					/* array range 0..na-1. */
  
  short	num_layers;			/* number of layers. */
  LayerStruct * layerspecs;	/* layer parameters. */	
};

/****
 *	Structures for scoring physical quantities. 
 *	z and r represent z and r coordinates of the 
 *	cylindrical coordinate system. [cm]
 *	a is the angle alpha between the photon exiting 
 *	direction and the normal to the surfaces. [radian]
 *	See comments of the InputStruct.
 *	See manual for the physcial quantities.
 ****/
struct OutStruct 
{
  double    Rsp;	/* specular reflectance. [-] */
  double ** Rd_ra;	/* 2D distribution of diffuse */
					/* reflectance. [1/(cm2 sr)] */
  double *  Rd_r;	/* 1D radial distribution of diffuse */
					/* reflectance. [1/cm2] */
  double *  Rd_a;	/* 1D angular distribution of diffuse */
					/* reflectance. [1/sr] */
  double    Rd;		/* total diffuse reflectance. [-] */
  
  double ** A_rz;	/* 2D probability density in turbid */
					/* media over r & z. [1/cm3] */
  double *  A_z;	/* 1D probability density over z. */
					/* [1/cm] */
  double *  A_l;	/* each layer's absorption */
					/* probability. [-] */
  double    A;		/* total absorption probability. [-] */
  
  double ** Tt_ra;	/* 2D distribution of total */
					/* transmittance. [1/(cm2 sr)] */
  double *  Tt_r;	/* 1D radial distribution of */
					/* transmittance. [1/cm2] */
  double *  Tt_a;	/* 1D angular distribution of */
					/* transmittance. [1/sr] */
  double    Tt;		/* total transmittance. [-] */
};

/***********************************************************
 *	Routine prototypes for dynamic memory allocation and 
 *	release of arrays and matrices.
 *	Modified from Numerical Recipes in C.
 ****/
double  *AllocVector(short, short);
double **AllocMatrix(short, short,short, short);
void 	 FreeVector(double *, short, short);
void 	 FreeMatrix(double **, short, short, short, short);
void	 ShowVersion(const char* version);
void	 nrerror(const char error_text[]);

#pragma once
#include <cmath>
#include "mcmlrandom.h"

// #define PARTIALREFLECTION 0
  /* 1=split photon, 0=statistical reflection. */

constexpr double COSZERO = 1.0 - 1.0E-12;
/* cosine of about 1e-6 rad. */

constexpr double COS90D = 1.0E-6;
/* cosine of about 1.57 - 1e-6 rad. */

double SpinTheta(const double anisotropy);

/***********************************************************
 *	Compute the Fresnel reflectance.
 *
 *	Make sure that the cosine of the incident angle a1
 *	is positive, and the case when the angle is greater
 *	than the critical angle is ruled out.
 *
 * 	Avoid trigonometric function operations as much as
 *	possible, because they are computation-intensive.
 ****/
double RFresnel(double n1,	/* incident refractive index.*/
	double n2,	/* transmit refractive index.*/
	double ca1,	/* cosine of the incident. angle. 0<a1<90 degrees*/
	double* ca2_Ptr)  /* pointer to the cosine of the transmission angle. a2 > 0. */
	;
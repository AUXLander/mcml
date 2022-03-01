#pragma once

#include "io/io.hpp"

#include "math/math.h"

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

	double cos_crit0, cos_crit1;
};

double Rspecular(LayerStruct* Layerspecs_Ptr);

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

	InputStruct& In_Ptr;

	PhotonStruct(InputStruct& cfg)
		: In_Ptr(cfg)
	{;}

	void init(const double Rspecular, LayerStruct* Layerspecs_Ptr);

	void spin(const double anisotropy);

	void Hop();

	void StepSizeInGlass();
	bool HitBoundary();
	void Roulette();
	void RecordR(double	Refl /* reflectance. */, OutStruct& Out_Ptr);
	void RecordT(double Refl, OutStruct& Out_Ptr);
	void Drop(OutStruct& Out_Ptr);

	void CrossUpOrNot(OutStruct& Out_Ptr);
	void CrossDnOrNot(OutStruct& Out_Ptr);

	void CrossOrNot(OutStruct& Out_Ptr);

	void HopInGlass(OutStruct& Out_Ptr);
	void HopDropSpin(OutStruct& Out_Ptr);
};
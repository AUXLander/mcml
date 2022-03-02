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

	inline bool is_glass() const
	{
		return mua == 0.0 && mus == 0.0;
	}
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

	void hop();

	void step_size_in_glass();
	bool hit_boundary();
	void roulette();
	void record_r(double	Refl /* reflectance. */, OutStruct& Out_Ptr);
	void record_t(double Refl, OutStruct& Out_Ptr);
	void drop(OutStruct& Out_Ptr);

	void cross_up_or_not(OutStruct& Out_Ptr);
	void cross_down_or_not(OutStruct& Out_Ptr);

	void cross_or_not(OutStruct& Out_Ptr);

	void hop_in_glass(OutStruct& Out_Ptr);
	void hop_drop_spin(OutStruct& Out_Ptr);
};
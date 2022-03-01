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
#include <cmath>

#include "mcml.hpp"

template<typename T, typename Y>
inline T setsign(T value, bool sign)
{
	// true is positive

	constexpr auto signbit_offs = sizeof(T) * 8U - 1U;
	constexpr auto signbit_mask = ~(static_cast<Y>(1U) << signbit_offs);

	Y &reinterpret = *reinterpret_cast<Y*>(&value);

	reinterpret = (reinterpret & signbit_mask) | (static_cast<Y>(!sign) << signbit_offs);

	return value;
}

double Rspecular(LayerStruct* Layerspecs_Ptr)
{
	double r1, r2;
	/* direct reflections from the 1st and 2nd layers. */
	double temp;

	temp = (Layerspecs_Ptr[0].n - Layerspecs_Ptr[1].n) / (Layerspecs_Ptr[0].n + Layerspecs_Ptr[1].n);
	r1 = temp * temp;

	if ((Layerspecs_Ptr[1].mua == 0.0) && (Layerspecs_Ptr[1].mus == 0.0))
	{
		/* glass layer. */
		temp = (Layerspecs_Ptr[1].n - Layerspecs_Ptr[2].n) / (Layerspecs_Ptr[1].n + Layerspecs_Ptr[2].n);

		r2 = temp * temp;

		r1 = r1 + (1 - r1) * (1 - r1) * r2 / (1 - r1 * r2);
	}

	return (r1);
}

/***********************************************************
 *	Initialize a photon packet.
 ****/
void PhotonStruct::init(const double Rspecular, LayerStruct* Layerspecs_Ptr)
{
	w = 1.0 - Rspecular;
	dead = 0;
	layer = 1;
	step_size = 0;
	sleft = 0;

	x = 0.0;
	y = 0.0;
	z = 0.0;
	
	ux = 0.0;
	uy = 0.0;
	uz = 1.0;

	/* glass layer. */
	if ((Layerspecs_Ptr[1].mua == 0.0) && (Layerspecs_Ptr[1].mus == 0.0))
	{
		layer = 2;
		z = Layerspecs_Ptr[2].z0;
	}
}

void PhotonStruct::spin(const double anisotropy)
{
	double cost, sint;	/* cosine and sine of the */
						/* polar deflection angle theta. */
	double cosp, sinp;	/* cosine and sine of the */
						/* azimuthal angle psi. */

	double ux = this->ux;
	double uy = this->uy;
	double uz = this->uz;
	double psi;

	cost = SpinTheta(anisotropy);
	sint = std::sqrt(1.0 - cost * cost); /* sqrt() is faster than sin(). */

	psi = 2.0 * PI * RandomNum(); /* spin psi 0-2pi. */
	cosp = std::cos(psi);

	sinp = setsign<double, uint64_t>(std::sqrt(1.0 - cosp * cosp), psi < PI);

	//sinp = ((psi < PI) ? +1.0 : -1.0) * std::sqrt(1.0 - cosp * cosp); /* sqrt() is faster than sin(). */

	if (fabs(uz) > COSZERO) /* normal incident. */
	{
		this->ux = sint * cosp;
		this->uy = sint * sinp;
		this->uz = cost * SIGN(uz); /* SIGN() is faster than division. */
	}
	else /* regular incident. */
	{
		double temp = std::sqrt(1.0 - uz * uz);

		this->ux = sint * (ux * uz * cosp - uy * sinp) / temp + ux * cost;
		this->uy = sint * (uy * uz * cosp + ux * sinp) / temp + uy * cost;
		this->uz = -sint * cosp * temp + uz * cost;
	}
}

void PhotonStruct::Hop()
{
	double step_size = this->step_size;

	this->x += step_size * this->ux;
	this->y += step_size * this->uy;
	this->z += step_size * this->uz;
}

void PhotonStruct::StepSizeInGlass()
{
	double dl_b;	/* step size to boundary. */
	short  layer = this->layer;
	double uz = this->uz;

	/* Stepsize to the boundary. */
	if (uz > 0.0)
	{
		dl_b = (In_Ptr.layerspecs[layer].z1 - this->z) / uz;
	}
	else if (uz < 0.0)
	{
		dl_b = (In_Ptr.layerspecs[layer].z0 - this->z) / uz;
	}
	else
	{
		dl_b = 0.0;
	}

	this->step_size = dl_b;
}

bool PhotonStruct::HitBoundary()
{
	double dl_b;  /* length to boundary. */
	short  layer = this->layer;
	double uz = this->uz;
	bool hit;

	/* Distance to the boundary. */
	if (uz > 0.0)
	{
		dl_b = (In_Ptr.layerspecs[layer].z1 - this->z) / uz;	/* dl_b>0. */
	}
	else if (uz < 0.0)
	{
		dl_b = (In_Ptr.layerspecs[layer].z0 - this->z) / uz;	/* dl_b>0. */
	}

	if (uz != 0.0 && this->step_size > dl_b)
	{
		/* not horizontal & crossing. */
		double mut = In_Ptr.layerspecs[layer].mua + In_Ptr.layerspecs[layer].mus;

		this->sleft = (this->step_size - dl_b) * mut;
		this->step_size = dl_b;

		hit = 1;
	}
	else
	{
		hit = 0;
	}

	return hit;
}

void PhotonStruct::Roulette()
{
	if (this->w == 0.0)
	{
		this->dead = 1;
	}
	else if (RandomNum() < CHANCE) /* survived the roulette.*/
	{
		this->w /= CHANCE;
	}
	else
	{
		this->dead = 1;
	}
}

void PhotonStruct::RecordR(double Refl, OutStruct& Out_Ptr)
{
	double x = this->x;
	double y = this->y;
	short  ir, ia;	/* index to r & angle. */

	ir = (short)(sqrt(x * x + y * y) / In_Ptr.dr);
	if (ir > In_Ptr.nr - 1)
	{
		ir = In_Ptr.nr - 1;
	}

	ia = (short)(acos(-this->uz) / In_Ptr.da);
	if (ia > In_Ptr.na - 1)
	{
		ia = In_Ptr.na - 1;
	}

	/* assign photon to the reflection array element. */
	Out_Ptr.Rd_ra[ir][ia] += this->w * (1.0 - Refl);

	this->w *= Refl;
}

void PhotonStruct::Drop(OutStruct& Out_Ptr)
{
	double dwa;		/* absorbed weight.*/
	double x = this->x;
	double y = this->y;
	short  iz, ir;	/* index to z & r. */
	short  layer = this->layer;
	double mua, mus;

	/* compute array indices. */
	iz = (short)(this->z / In_Ptr.dz);
	if (iz > In_Ptr.nz - 1) iz = In_Ptr.nz - 1;

	ir = (short)(sqrt(x * x + y * y) / In_Ptr.dr);
	if (ir > In_Ptr.nr - 1) ir = In_Ptr.nr - 1;

	/* update photon weight. */
	mua = In_Ptr.layerspecs[layer].mua;
	mus = In_Ptr.layerspecs[layer].mus;
	dwa = this->w * mua / (mua + mus);
	this->w -= dwa;

	/* assign dwa to the absorption array element. */
	Out_Ptr.A_rz[ir][iz] += dwa;
}


void PhotonStruct::RecordT(double Refl, OutStruct& Out_Ptr)
{
	double x = this->x;
	double y = this->y;
	short  ir, ia;	/* index to r & angle. */

	ir = (short)(sqrt(x * x + y * y) / In_Ptr.dr);
	if (ir > In_Ptr.nr - 1) ir = In_Ptr.nr - 1;

	ia = (short)(acos(this->uz) / In_Ptr.da);
	if (ia > In_Ptr.na - 1) ia = In_Ptr.na - 1;

	/* assign photon to the transmittance array element. */
	Out_Ptr.Tt_ra[ir][ia] += this->w * (1.0 - Refl);

	this->w *= Refl;
}

void PhotonStruct::CrossUpOrNot(OutStruct& Out_Ptr)
{
	double uz = this->uz; /* z directional cosine. */
	double uz1;					/* cosines of transmission alpha. always positive */
	double r = 0.0;				/* reflectance */
	short  layer = this->layer;
	double ni = In_Ptr.layerspecs[layer].n;
	double nt = In_Ptr.layerspecs[layer - 1].n;

	/* Get r. */
	if (-uz <= In_Ptr.layerspecs[layer].cos_crit0)
	{
		r = 1.0; /* total internal reflection. */
	}
	else
	{
		r = RFresnel(ni, nt, -uz, &uz1);
	}

	if (RandomNum() > r) /* transmitted to layer-1. */
	{
		if (layer == 1)
		{
			this->uz = -uz1;
			this->RecordR(0.0, Out_Ptr);
			this->dead = 1;
		}
		else
		{
			this->layer--;
			this->ux *= ni / nt;
			this->uy *= ni / nt;
			this->uz = -uz1;
		}
	}
	else /* reflected. */
	{
		this->uz = -uz;
	}
}

void PhotonStruct::CrossDnOrNot(OutStruct& Out_Ptr)
{
	double uz = this->uz; /* z directional cosine. */
	double uz1;	/* cosines of transmission alpha. */
	double r = 0.0;	/* reflectance */
	short  layer = this->layer;
	double ni = In_Ptr.layerspecs[layer].n;
	double nt = In_Ptr.layerspecs[layer + 1].n;

	/* Get r. */
	if (uz <= In_Ptr.layerspecs[layer].cos_crit1)
	{
		r = 1.0; /* total internal reflection. */
	}
	else
	{
		r = RFresnel(ni, nt, uz, &uz1);
	}

	if (RandomNum() > r) /* transmitted to layer+1. */
	{
		if (layer == In_Ptr.num_layers)
		{
			this->uz = uz1;
			this->RecordT(0.0, Out_Ptr);
			this->dead = 1;
		}
		else
		{
			this->layer++;
			this->ux *= ni / nt;
			this->uy *= ni / nt;
			this->uz = uz1;
		}
	}
	else /* reflected. */
	{
		this->uz = -uz;
	}
}

void PhotonStruct::CrossOrNot(OutStruct& Out_Ptr)
{
	if (this->uz < 0.0)
	{
		this->CrossUpOrNot(Out_Ptr);
	}
	else
	{
		this->CrossDnOrNot(Out_Ptr);
	}
}

void PhotonStruct::HopInGlass(OutStruct& Out_Ptr)
{
	double dl;     /* step size. 1/cm */

	if (this->uz == 0.0)
	{
		/* horizontal photon in glass is killed. */
		this->dead = 1;
	}
	else
	{
		this->StepSizeInGlass();
		this->Hop(); // Move the photon s away in the current layer of medium. 
		this->CrossOrNot(Out_Ptr);
	}
}

void PhotonStruct::HopDropSpin(OutStruct& Out_Ptr)
{
	short layer = this->layer;

	if ((In_Ptr.layerspecs[layer].mua == 0.0) && (In_Ptr.layerspecs[layer].mus == 0.0))
	{
		/* glass layer. */
		this->HopInGlass(Out_Ptr);
	}
	else
	{
		//HopDropSpinInTissue(In_Ptr, Photon_Ptr, Out_Ptr);

		//StepSizeInTissue(Photon_Ptr, In_Ptr);

		short  layer = this->layer;
		double mua = In_Ptr.layerspecs[layer].mua;
		double mus = In_Ptr.layerspecs[layer].mus;

		if (this->sleft == 0.0) /* make a new step. */
		{
			double rnd;

			do
			{
				rnd = RandomNum();
			}
			while (rnd <= 0.0); /* avoid zero. */

			this->step_size = -log(rnd) / (mua + mus);
		}
		else /* take the leftover. */
		{
			this->step_size = this->sleft / (mua + mus);
			this->sleft = 0.0;
		}

		//
		const auto hit = this->HitBoundary();

		this->Hop();	/* move to boundary plane. */

		if (hit)
		{
			this->CrossOrNot(Out_Ptr);
		}
		else
		{
			this->Drop(Out_Ptr);
			this->spin(In_Ptr.layerspecs[this->layer].anisotropy);
		}
	}

	if (this->w < In_Ptr.Wth && !this->dead)
	{
		this->Roulette();
	}
}
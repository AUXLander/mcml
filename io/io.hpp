#pragma once

#include "../mcmlnr.hpp"
#include <memory>
#include <boost/numeric/ublas/matrix.hpp>

#define STRLEN 256		/* String length. */

struct LayerStruct;

#define CPP_COMPILE

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
	LayerStruct* layerspecs;	/* layer parameters. */
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


 //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\


#ifndef CPP_COMPILE
struct OutStruct
{
	double    Rsp;	/* specular reflectance. [-] */
	double** Rd_ra;	/* 2D distribution of diffuse */
					  /* reflectance. [1/(cm2 sr)] */
	double* Rd_r;	/* 1D radial distribution of diffuse */
					  /* reflectance. [1/cm2] */
	double* Rd_a;	/* 1D angular distribution of diffuse */
					  /* reflectance. [1/sr] */
	double    Rd;		/* total diffuse reflectance. [-] */

	double** A_rz;	/* 2D probability density in turbid */
					  /* media over r & z. [1/cm3] */
	double* A_z;	/* 1D probability density over z. */
					  /* [1/cm] */
	double* A_l;	/* each layer's absorption */
					  /* probability. [-] */
	double    A;		/* total absorption probability. [-] */

	double** Tt_ra;	/* 2D distribution of total */
					  /* transmittance. [1/(cm2 sr)] */
	double* Tt_r;	/* 1D radial distribution of */
					  /* transmittance. [1/cm2] */
	double* Tt_a;	/* 1D angular distribution of */
					  /* transmittance. [1/sr] */
	double    Tt;		/* total transmittance. [-] */
};
#else
using namespace boost::numeric::ublas;

struct OutStruct
{
	using value_type = double;

	using matrix = boost::numeric::ublas::matrix<value_type>;
	using vector = std::vector<value_type>;

	using unique_matrix_ptr = std::unique_ptr<matrix>;
	using unique_vector_ptr = std::unique_ptr<vector>;
	
private:
	struct matrix_adaptor
	{
		const unique_matrix_ptr mtx;

		struct matrix_column_adaptor
		{
			const unique_matrix_ptr& mtx;

			const size_t i;

			matrix_column_adaptor(const unique_matrix_ptr& mtx, const size_t i)
				: mtx(mtx), i(i)
			{;}

			value_type& operator[](const size_t& j)
			{
				return mtx->at_element(i, j);
			}

			const value_type& operator[](const size_t& j) const
			{
				return mtx->at_element(i, j);
			}
		};

		matrix_adaptor(unique_matrix_ptr&& mtx_)
			: mtx(std::move(mtx_))
		{
			for (auto i = 0; i < mtx->size1(); ++i)
			{
				for (auto j = 0; j < mtx->size2(); ++j)
				{
					mtx->at_element(i, j) = 0.0;
				}
			}
		}

		matrix_column_adaptor operator[](const size_t& i)
		{
			return matrix_column_adaptor(mtx, i);
		}

		matrix_column_adaptor operator[](const size_t& i) const
		{
			return matrix_column_adaptor(mtx, i);
		}
	};

public:

	value_type Rsp;

	matrix_adaptor Rd_ra;
	vector Rd_r;
	vector Rd_a;
	value_type Rd;

	matrix_adaptor A_rz;
	vector A_z;
	vector A_l;
	value_type A;

	matrix_adaptor Tt_ra;
	vector Tt_r;
	vector Tt_a;
	value_type Tt;

	OutStruct(InputStruct cfg)
		: Rsp(0.0), Rd(0.0), A(0.0), Tt(0.0),
		  /* Allocate the arrays and the matrices. */
		  Rd_ra(std::make_unique<matrix>(cfg.nr, cfg.na)), Rd_r(cfg.nr), Rd_a(cfg.na),
		  A_rz (std::make_unique<matrix>(cfg.nr, cfg.nz)), A_z(cfg.nz),  A_l(cfg.num_layers + 2),
		  Tt_ra(std::make_unique<matrix>(cfg.nr, cfg.na)), Tt_r(cfg.nr), Tt_a(cfg.na)
	{
		/***********************************************************
		 *	Allocate the arrays in OutStruct for one run, and
		 *	array elements are automatically initialized to zeros.
		 ****/

		short nz = cfg.nz;
		short nr = cfg.nr;
		short na = cfg.na;
		short nl = cfg.num_layers;
		/* remember to use nl+2 because of 2 for ambient. */

		if (nz <= 0 || nr <= 0 || na <= 0 || nl <= 0)
		{
			nrerror("Wrong grid parameters.\n");
		}
	}

	~OutStruct()
	{;}
	
	void FreeData(InputStruct In_Parm)
	{
		free(In_Parm.layerspecs);
	}
};
#endif

//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\
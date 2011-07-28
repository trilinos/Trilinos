#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#ifndef DIVERGENCE
#define DIVERGENCE

#define ONE12TH 0.083333333333333333333333
//
//
//
#define HG_X1 0
#define HG_Y1 1
#define HG_Z1 2
#define HG_X2 3
#define HG_Y2 4
#define HG_Z2 5
#define HG_X3 6
#define HG_Y3 7
#define HG_Z3 8
#define HG_X4 9
#define HG_Y4 10
#define HG_Z4 11
//
//  Indexes into a full 3 by 3 tensor stored as a length 9 vector
//
#define K_F_XX 0
#define K_F_YY 1
#define K_F_ZZ 2
#define K_F_XY 3
#define K_F_YZ 4
#define K_F_ZX 5
#define K_F_YX 6
#define K_F_ZY 7
#define K_F_XZ 8


template< typename Scalar , class DeviceType >
struct divergence;

template<typename Scalar>
struct divergence<Scalar, KOKKOS_MACRO_DEVICE>{

	typedef KOKKOS_MACRO_DEVICE     device_type ;
	typedef typename Kokkos::MDArrayView<Scalar,device_type> array_type ;

	const array_type  coords;    
	const array_type  velocity;
	const array_type  force_new;
	const array_type  vorticity;
	const array_type  rotation;
	const array_type  stress_new;
	const array_type  rot_stress;
	const array_type  rot_stret;
	const array_type  gradop12;
	const array_type  elem_mass;
	const array_type  elem_dilmod;
	const array_type  elem_shrmod;
	const array_type  elem_t_step;
	const array_type  intern_energy;
	const array_type  mid_vol;

	const array_type  hgop;
	const array_type  hg_resist;
	const array_type  hg_energy;  

	const array_type  two_mu;
	const array_type  bulk_mod;

	const Scalar     hg_stiff;
	const Scalar     hg_visc;
	const Scalar     linBulkVisc;
	const Scalar     quadBulkVisc;
	const Scalar     dt;

	const bool    scaleHGRotation;

	divergence(	const array_type & arg_p,
				const array_type & arg_v,
				const array_type & arg_f,
				const array_type & arg_vort,
				const array_type & arg_r,
				const array_type & arg_sn,
				const array_type & arg_rstrss,
				const array_type & arg_rstret,
				const array_type & arg_g,
				const array_type & arg_em,
				const array_type & arg_ed,
				const array_type & arg_es,
				const array_type & arg_et,
				const array_type & arg_ie,
				const array_type & arg_mv,
				const array_type & arg_hg,
				const array_type & arg_hgr,
				const array_type & arg_hge,
				const array_type & arg_tmu,
				const array_type & arg_blk,
				const Scalar arg_hgs,
				const Scalar arg_hgv,
				const Scalar arg_l,
				const Scalar arg_q,
				const Scalar delta_t,
				const bool arg_scale       )
		: coords(     arg_p )
		, velocity(     arg_v )
		, force_new(   arg_f )
		, vorticity(   arg_vort )
		, rotation(     arg_r )
		, stress_new(   arg_sn )
		, rot_stress(   arg_rstrss )
		, rot_stret(   arg_rstret )
		, gradop12(     arg_g )
		, elem_mass(   arg_em )
		, elem_dilmod(   arg_ed )
		, elem_shrmod(   arg_es )
		, elem_t_step(   arg_et )
		, intern_energy( arg_ie )
		, mid_vol(     arg_mv )

		, hgop(       arg_hg )
		, hg_resist(   arg_hgr )
		, hg_energy(   arg_hge )

		, two_mu(     arg_tmu )
		, bulk_mod(     arg_blk )

		, hg_stiff(     arg_hgs )
		, hg_visc(     arg_hgv )
		, linBulkVisc(   arg_l )
		, quadBulkVisc(   arg_q )
		, dt(       delta_t )

  , scaleHGRotation( arg_scale )
  {}


  KOKKOS_MACRO_DEVICE_FUNCTION
  void comp_grad(int ielem) const {

		Scalar x[8], y[8], z[8];

		//  load X coordinate information
		x[0] = coords(ielem, 0, 0);
		x[1] = coords(ielem, 0, 1);
		x[2] = coords(ielem, 0, 2);
		x[3] = coords(ielem, 0, 3);
		x[4] = coords(ielem, 0, 4);
		x[5] = coords(ielem, 0, 5);
		x[6] = coords(ielem, 0, 6);
		x[7] = coords(ielem, 0, 7);

		//   calc X difference vectors
		Scalar R42=(x[3] - x[1]);
		Scalar R52=(x[4] - x[1]);
		Scalar R54=(x[4] - x[3]);

		Scalar R63=(x[5] - x[2]);
		Scalar R83=(x[7] - x[2]);
		Scalar R86=(x[7] - x[5]);

		Scalar R31=(x[2] - x[0]);
		Scalar R61=(x[5] - x[0]);
		Scalar R74=(x[6] - x[3]);

		Scalar R72=(x[6] - x[1]);
		Scalar R75=(x[6] - x[4]);
		Scalar R81=(x[7] - x[0]);

		Scalar t1=(R63 + R54);
		Scalar t2=(R61 + R74);
		Scalar t3=(R72 + R81);

		Scalar t4 =(R86 + R42);
		Scalar t5 =(R83 + R52);
		Scalar t6 =(R75 + R31);  

		//  Load Z information
		z[0] = coords(ielem, 2, 0);
		z[1] = coords(ielem, 2, 1);
		z[2] = coords(ielem, 2, 2);
		z[3] = coords(ielem, 2, 3);
		z[4] = coords(ielem, 2, 4);
		z[5] = coords(ielem, 2, 5);
		z[6] = coords(ielem, 2, 6);
		z[7] = coords(ielem, 2, 7);


		//  Calculate Y gradient from X and Z data
		gradop12(ielem, 1, 0) = (z[1] *  t1) - (z[2] * R42) - (z[3] *  t5)  + (z[4] *  t4) + (z[5] * R52) - (z[7] * R54); 
		gradop12(ielem, 1, 1) = (z[2] *  t2) + (z[3] * R31) - (z[0] *  t1)  - (z[5] *  t6) + (z[6] * R63) - (z[4] * R61); 
		gradop12(ielem, 1, 2) = (z[3] *  t3) + (z[0] * R42) - (z[1] *  t2)  - (z[6] *  t4) + (z[7] * R74) - (z[5] * R72); 
		gradop12(ielem, 1, 3) = (z[0] *  t5) - (z[1] * R31) - (z[2] *  t3)  + (z[7] *  t6) + (z[4] * R81) - (z[6] * R83); 
		gradop12(ielem, 1, 4) = (z[5] *  t3) + (z[6] * R86) - (z[7] *  t2)  - (z[0] *  t4) - (z[3] * R81) + (z[1] * R61); 
		gradop12(ielem, 1, 5) = (z[6] *  t5) - (z[4] *  t3)  - (z[7] * R75) + (z[1] *  t6) - (z[0] * R52) + (z[2] * R72);
		gradop12(ielem, 1, 6) = (z[7] *  t1) - (z[5] *  t5)  - (z[4] * R86) + (z[2] *  t4) - (z[1] * R63) + (z[3] * R83);
		gradop12(ielem, 1, 7) = (z[4] *  t2) - (z[6] *  t1)  + (z[5] * R75) - (z[3] *  t6) - (z[2] * R74) + (z[0] * R54); 


		//   calc Z difference vectors
		R42=(z[3] - z[1]);
		R52=(z[4] - z[1]);
		R54=(z[4] - z[3]);

		R63=(z[5] - z[2]);
		R83=(z[7] - z[2]);
		R86=(z[7] - z[5]);

		 R31=(z[2] - z[0]);
		R61=(z[5] - z[0]);
		R74=(z[6] - z[3]);

		R72=(z[6] - z[1]);
		R75=(z[6] - z[4]);
		R81=(z[7] - z[0]);

		t1=(R63 + R54);
		t2=(R61 + R74);
		t3=(R72 + R81);

		t4 =(R86 + R42);
		t5 =(R83 + R52);
		t6 =(R75 + R31);  

		//  Load Y information
		y[0] = coords(ielem, 1, 0);
		y[1] = coords(ielem, 1, 1);
		y[2] = coords(ielem, 1, 2);
		y[3] = coords(ielem, 1, 3);
		y[4] = coords(ielem, 1, 4);
		y[5] = coords(ielem, 1, 5);
		y[6] = coords(ielem, 1, 6);
		y[7] = coords(ielem, 1, 7);


		//  Calculate X gradient from Y and Z data
		gradop12(ielem, 0, 0) = (y[1] *  t1) - (y[2] * R42) - (y[3] *  t5) + (y[4] *  t4) + (y[5] * R52) - (y[7] * R54); 
		gradop12(ielem, 0, 1) = (y[2] *  t2) + (y[3] * R31) - (y[0] *  t1) - (y[5] *  t6) + (y[6] * R63) - (y[4] * R61); 
		gradop12(ielem, 0, 2) = (y[3] *  t3) + (y[0] * R42) - (y[1] *  t2) - (y[6] *  t4) + (y[7] * R74) - (y[5] * R72); 
		gradop12(ielem, 0, 3) = (y[0] *  t5) - (y[1] * R31) - (y[2] *  t3) + (y[7] *  t6) + (y[4] * R81) - (y[6] * R83); 
		gradop12(ielem, 0, 4) = (y[5] *  t3) + (y[6] * R86) - (y[7] *  t2) - (y[0] *  t4) - (y[3] * R81) + (y[1] * R61); 
		gradop12(ielem, 0, 5) = (y[6] *  t5) - (y[4] *  t3) - (y[7] * R75) + (y[1] *  t6) - (y[0] * R52) + (y[2] * R72);
		gradop12(ielem, 0, 6) = (y[7] *  t1) - (y[5] *  t5) - (y[4] * R86) + (y[2] *  t4) - (y[1] * R63) + (y[3] * R83);
		gradop12(ielem, 0, 7) = (y[4] *  t2) - (y[6] *  t1) + (y[5] * R75) - (y[3] *  t6) - (y[2] * R74) + (y[0] * R54); 


		//   calc Y difference vectors
		R42=(y[3] - y[1]);
		R52=(y[4] - y[1]);
		R54=(y[4] - y[3]);

		R63=(y[5] - y[2]);
		R83=(y[7] - y[2]);
		R86=(y[7] - y[5]);

		 R31=(y[2] - y[0]);
		R61=(y[5] - y[0]);
		R74=(y[6] - y[3]);

		R72=(y[6] - y[1]);
		R75=(y[6] - y[4]);
		R81=(y[7] - y[0]);

		t1=(R63 + R54);
		t2=(R61 + R74);
		t3=(R72 + R81);

		t4 =(R86 + R42);
		t5 =(R83 + R52);
		t6 =(R75 + R31);

		//  Calculate Z gradient from X and Y data

		gradop12(ielem, 2, 0) = (x[1] *  t1) - (x[2] * R42) - (x[3] *  t5)  + (x[4] *  t4) + (x[5] * R52) - (x[7] * R54); 
		gradop12(ielem, 2, 1) = (x[2] *  t2) + (x[3] * R31) - (x[0] *  t1)  - (x[5] *  t6) + (x[6] * R63) - (x[4] * R61); 
		gradop12(ielem, 2, 2) = (x[3] *  t3) + (x[0] * R42) - (x[1] *  t2)  - (x[6] *  t4) + (x[7] * R74) - (x[5] * R72); 
		gradop12(ielem, 2, 3) = (x[0] *  t5) - (x[1] * R31) - (x[2] *  t3)  + (x[7] *  t6) + (x[4] * R81) - (x[6] * R83); 
		gradop12(ielem, 2, 4) = (x[5] *  t3) + (x[6] * R86) - (x[7] *  t2)  - (x[0] *  t4) - (x[3] * R81) + (x[1] * R61); 
		gradop12(ielem, 2, 5) = (x[6] *  t5) - (x[4] *  t3)  - (x[7] * R75) + (x[1] *  t6) - (x[0] * R52) + (x[2] * R72);
		gradop12(ielem, 2, 6) = (x[7] *  t1) - (x[5] *  t5)  - (x[4] * R86) + (x[2] *  t4) - (x[1] * R63) + (x[3] * R83);
		gradop12(ielem, 2, 7) = (x[4] *  t2) - (x[6] *  t1)  + (x[5] * R75) - (x[3] *  t6) - (x[2] * R74) + (x[0] * R54); 

		mid_vol(ielem) = ONE12TH * (gradop12(ielem, 0, 0) * x[0] +
				      gradop12(ielem, 0, 1) * x[1] +
				      gradop12(ielem, 0, 2) * x[2] +
				      gradop12(ielem, 0, 3) * x[3] +
				      gradop12(ielem, 0, 4) * x[4] +
				      gradop12(ielem, 0, 5) * x[5] +
				      gradop12(ielem, 0, 6) * x[6] +
				      gradop12(ielem, 0, 7) * x[7] );

	}

	KOKKOS_MACRO_DEVICE_FUNCTION
    void comp_hgop(  int ielem) const {

	// 	KHP: Alternatively, we could have
	// 	hx0,hx1,hx2,hx3,...,hz0,hz1,hz2,hz3
		Scalar hgconst12th[12];
		Scalar inv_vol = 1.0 / (mid_vol(ielem) * 12.0);

		Scalar q0 = coords(ielem, 0, 0) - coords(ielem, 0, 1);
		Scalar q1 = coords(ielem, 0, 2) - coords(ielem, 0, 3);
		Scalar q2 = coords(ielem, 0, 4) - coords(ielem, 0, 5);
		Scalar q3 = coords(ielem, 0, 6) - coords(ielem, 0, 7);

		hgconst12th[0] = (   (coords(ielem, 0, 0)+coords(ielem, 0, 1)) - (coords(ielem, 0, 2)+coords(ielem, 0, 3)) - 
			  (coords(ielem, 0, 4)+coords(ielem, 0, 5)) + (coords(ielem, 0, 6)+coords(ielem, 0, 7)) ) * inv_vol;
		hgconst12th[1] = (  q0 - q1 - q2 + q3 ) * inv_vol;
		hgconst12th[2] = (  q0 + q1 + q2 + q3 ) * inv_vol;
		hgconst12th[3] = ( -q0 - q1 + q2 + q3 ) * inv_vol;

		q0 = (coords(ielem, 1, 0) - coords(ielem, 1, 1));
		q1 = (coords(ielem, 1, 2) - coords(ielem, 1, 3));
		q2 = (coords(ielem, 1, 4) - coords(ielem, 1, 5));
		q3 = (coords(ielem, 1, 6) - coords(ielem, 1, 7));

		hgconst12th[4] = (   (coords(ielem, 1, 0)+coords(ielem, 1, 1)) - (coords(ielem, 1, 2)+coords(ielem, 1, 3)) - 
			  (coords(ielem, 1, 4)+coords(ielem, 1, 5)) + (coords(ielem, 1, 6)+coords(ielem, 1, 7)) ) * inv_vol;
		hgconst12th[5] = (  q0 - q1 - q2 + q3 ) * inv_vol;
		hgconst12th[6] = (  q0 + q1 + q2 + q3 ) * inv_vol;
		hgconst12th[7] = ( -q0 - q1 + q2 + q3 ) * inv_vol;

		q0 = (coords(ielem, 2, 0) - coords(ielem, 2, 1));
		q1 = (coords(ielem, 2, 2) - coords(ielem, 2, 3));
		q2 = (coords(ielem, 2, 4) - coords(ielem, 2, 5));
		q3 = (coords(ielem, 2, 6) - coords(ielem, 2, 7));

		hgconst12th[8]  = ( (coords(ielem, 2, 0)+coords(ielem, 2, 1)) - (coords(ielem, 2, 2)+coords(ielem, 2, 3)) - 
			  (coords(ielem, 2, 4)+coords(ielem, 2, 5)) + (coords(ielem, 2, 6)+coords(ielem, 2, 7)) ) * inv_vol;
		hgconst12th[9]  = (  q0 - q1 - q2 + q3 ) * inv_vol;
		hgconst12th[10] = (  q0 + q1 + q2 + q3 ) * inv_vol;
		hgconst12th[11] = ( -q0 - q1 + q2 + q3 ) * inv_vol;


		hgop(ielem,  0, 1) =  1.0 - (hgconst12th[0] * gradop12(ielem, 0, 0) + hgconst12th[4] * gradop12(ielem, 1, 0) + hgconst12th[ 8] * gradop12(ielem, 2, 0));
		hgop(ielem,  1, 1) =  1.0 - (hgconst12th[0] * gradop12(ielem, 0, 1) + hgconst12th[4] * gradop12(ielem, 1, 1) + hgconst12th[ 8] * gradop12(ielem, 2, 1));
		hgop(ielem,  2, 1) = -1.0 - (hgconst12th[0] * gradop12(ielem, 0, 2) + hgconst12th[4] * gradop12(ielem, 1, 2) + hgconst12th[ 8] * gradop12(ielem, 2, 2));
		hgop(ielem,  3, 1) = -1.0 - (hgconst12th[0] * gradop12(ielem, 0, 3) + hgconst12th[4] * gradop12(ielem, 1, 3) + hgconst12th[ 8] * gradop12(ielem, 2, 3));
		hgop(ielem,  4, 1) = -1.0 - (hgconst12th[0] * gradop12(ielem, 0, 4) + hgconst12th[4] * gradop12(ielem, 1, 4) + hgconst12th[ 8] * gradop12(ielem, 2, 4));
		hgop(ielem,  5, 1) = -1.0 - (hgconst12th[0] * gradop12(ielem, 0, 5) + hgconst12th[4] * gradop12(ielem, 1, 5) + hgconst12th[ 8] * gradop12(ielem, 2, 5));
		hgop(ielem,  6, 1) =  1.0 - (hgconst12th[0] * gradop12(ielem, 0, 6) + hgconst12th[4] * gradop12(ielem, 1, 6) + hgconst12th[ 8] * gradop12(ielem, 2, 6));
		hgop(ielem,  7, 1) =  1.0 - (hgconst12th[0] * gradop12(ielem, 0, 7) + hgconst12th[4] * gradop12(ielem, 1, 7) + hgconst12th[ 8] * gradop12(ielem, 2, 7));
		hgop(ielem,  8, 1) =  1.0 - (hgconst12th[1] * gradop12(ielem, 0, 0) + hgconst12th[5] * gradop12(ielem, 1, 0) + hgconst12th[ 9] * gradop12(ielem, 2, 0));
		hgop(ielem,  9, 1) = -1.0 - (hgconst12th[1] * gradop12(ielem, 0, 1) + hgconst12th[5] * gradop12(ielem, 1, 1) + hgconst12th[ 9] * gradop12(ielem, 2, 1));
		hgop(ielem, 10, 1) = -1.0 - (hgconst12th[1] * gradop12(ielem, 0, 2) + hgconst12th[5] * gradop12(ielem, 1, 2) + hgconst12th[ 9] * gradop12(ielem, 2, 2));
		hgop(ielem, 11, 1) =  1.0 - (hgconst12th[1] * gradop12(ielem, 0, 3) + hgconst12th[5] * gradop12(ielem, 1, 3) + hgconst12th[ 9] * gradop12(ielem, 2, 3));
		hgop(ielem, 12, 1) = -1.0 - (hgconst12th[1] * gradop12(ielem, 0, 4) + hgconst12th[5] * gradop12(ielem, 1, 4) + hgconst12th[ 9] * gradop12(ielem, 2, 4));
		hgop(ielem, 13, 1) =  1.0 - (hgconst12th[1] * gradop12(ielem, 0, 5) + hgconst12th[5] * gradop12(ielem, 1, 5) + hgconst12th[ 9] * gradop12(ielem, 2, 5));
		hgop(ielem, 14, 1) =  1.0 - (hgconst12th[1] * gradop12(ielem, 0, 6) + hgconst12th[5] * gradop12(ielem, 1, 6) + hgconst12th[ 9] * gradop12(ielem, 2, 6));
		hgop(ielem, 15, 1) = -1.0 - (hgconst12th[1] * gradop12(ielem, 0, 7) + hgconst12th[5] * gradop12(ielem, 1, 7) + hgconst12th[ 9] * gradop12(ielem, 2, 7));
		hgop(ielem, 16, 1) =  1.0 - (hgconst12th[2] * gradop12(ielem, 0, 0) + hgconst12th[6] * gradop12(ielem, 1, 0) + hgconst12th[10] * gradop12(ielem, 2, 0));
		hgop(ielem, 17, 1) = -1.0 - (hgconst12th[2] * gradop12(ielem, 0, 1) + hgconst12th[6] * gradop12(ielem, 1, 1) + hgconst12th[10] * gradop12(ielem, 2, 1));
		hgop(ielem, 18, 1) =  1.0 - (hgconst12th[2] * gradop12(ielem, 0, 2) + hgconst12th[6] * gradop12(ielem, 1, 2) + hgconst12th[10] * gradop12(ielem, 2, 2));
		hgop(ielem, 19, 1) = -1.0 - (hgconst12th[2] * gradop12(ielem, 0, 3) + hgconst12th[6] * gradop12(ielem, 1, 3) + hgconst12th[10] * gradop12(ielem, 2, 3));
		hgop(ielem, 20, 1) =  1.0 - (hgconst12th[2] * gradop12(ielem, 0, 4) + hgconst12th[6] * gradop12(ielem, 1, 4) + hgconst12th[10] * gradop12(ielem, 2, 4));
		hgop(ielem, 21, 1) = -1.0 - (hgconst12th[2] * gradop12(ielem, 0, 5) + hgconst12th[6] * gradop12(ielem, 1, 5) + hgconst12th[10] * gradop12(ielem, 2, 5));
		hgop(ielem, 22, 1) =  1.0 - (hgconst12th[2] * gradop12(ielem, 0, 6) + hgconst12th[6] * gradop12(ielem, 1, 6) + hgconst12th[10] * gradop12(ielem, 2, 6));
		hgop(ielem, 23, 1) = -1.0 - (hgconst12th[2] * gradop12(ielem, 0, 7) + hgconst12th[6] * gradop12(ielem, 1, 7) + hgconst12th[10] * gradop12(ielem, 2, 7));
		hgop(ielem, 24, 1) = -1.0 - (hgconst12th[3] * gradop12(ielem, 0, 0) + hgconst12th[7] * gradop12(ielem, 1, 0) + hgconst12th[11] * gradop12(ielem, 2, 0));
		hgop(ielem, 25, 1) =  1.0 - (hgconst12th[3] * gradop12(ielem, 0, 1) + hgconst12th[7] * gradop12(ielem, 1, 1) + hgconst12th[11] * gradop12(ielem, 2, 1));
		hgop(ielem, 26, 1) = -1.0 - (hgconst12th[3] * gradop12(ielem, 0, 2) + hgconst12th[7] * gradop12(ielem, 1, 2) + hgconst12th[11] * gradop12(ielem, 2, 2));
		hgop(ielem, 27, 1) =  1.0 - (hgconst12th[3] * gradop12(ielem, 0, 3) + hgconst12th[7] * gradop12(ielem, 1, 3) + hgconst12th[11] * gradop12(ielem, 2, 3));
		hgop(ielem, 28, 1) =  1.0 - (hgconst12th[3] * gradop12(ielem, 0, 4) + hgconst12th[7] * gradop12(ielem, 1, 4) + hgconst12th[11] * gradop12(ielem, 2, 4));
		hgop(ielem, 29, 1) = -1.0 - (hgconst12th[3] * gradop12(ielem, 0, 5) + hgconst12th[7] * gradop12(ielem, 1, 5) + hgconst12th[11] * gradop12(ielem, 2, 5));
		hgop(ielem, 30, 1) =  1.0 - (hgconst12th[3] * gradop12(ielem, 0, 6) + hgconst12th[7] * gradop12(ielem, 1, 6) + hgconst12th[11] * gradop12(ielem, 2, 6));
		hgop(ielem, 31, 1) = -1.0 - (hgconst12th[3] * gradop12(ielem, 0, 7) + hgconst12th[7] * gradop12(ielem, 1, 7) + hgconst12th[11] * gradop12(ielem, 2, 7));

  	}

	KOKKOS_MACRO_DEVICE_FUNCTION
	void rotate_tensor_backward(int ielem)const {

    // 	t : temporary variables
    // 	s_n : stress_new in local memory space
    // 	r_n : rotation_new in local memory space
    	Scalar t[9], s_n[6], r_n[9];

		s_n[0] = stress_new(ielem, 0);
		s_n[1] = stress_new(ielem, 1);
		s_n[2] = stress_new(ielem, 2);
		s_n[3] = stress_new(ielem, 3);
		s_n[4] = stress_new(ielem, 4);
		s_n[5] = stress_new(ielem, 5);

		r_n[0] = rotation(ielem, 0, 1);
		r_n[1] = rotation(ielem, 1, 1);
		r_n[2] = rotation(ielem, 2, 1);
		r_n[3] = rotation(ielem, 3, 1);
		r_n[4] = rotation(ielem, 4, 1);
		r_n[5] = rotation(ielem, 5, 1);
		r_n[6] = rotation(ielem, 6, 1);
		r_n[7] = rotation(ielem, 7, 1);
		r_n[8] = rotation(ielem, 8, 1);

		t[0] = s_n[K_S_XX]*r_n[K_F_XX]+ s_n[K_S_XY]*r_n[K_F_XY]+ s_n[K_S_XZ]*r_n[K_F_XZ];
		t[1] = s_n[K_S_YX]*r_n[K_F_XX]+ s_n[K_S_YY]*r_n[K_F_XY]+ s_n[K_S_YZ]*r_n[K_F_XZ];
		t[2] = s_n[K_S_ZX]*r_n[K_F_XX]+ s_n[K_S_ZY]*r_n[K_F_XY]+ s_n[K_S_ZZ]*r_n[K_F_XZ];
		t[3] = s_n[K_S_XX]*r_n[K_F_YX]+ s_n[K_S_XY]*r_n[K_F_YY]+ s_n[K_S_XZ]*r_n[K_F_YZ];
		t[4] = s_n[K_S_YX]*r_n[K_F_YX]+ s_n[K_S_YY]*r_n[K_F_YY]+ s_n[K_S_YZ]*r_n[K_F_YZ];
		t[5] = s_n[K_S_ZX]*r_n[K_F_YX]+ s_n[K_S_ZY]*r_n[K_F_YY]+ s_n[K_S_ZZ]*r_n[K_F_YZ];
		t[6] = s_n[K_S_XX]*r_n[K_F_ZX]+ s_n[K_S_XY]*r_n[K_F_ZY]+ s_n[K_S_XZ]*r_n[K_F_ZZ];
		t[7] = s_n[K_S_YX]*r_n[K_F_ZX]+ s_n[K_S_YY]*r_n[K_F_ZY]+ s_n[K_S_YZ]*r_n[K_F_ZZ];
		t[8] = s_n[K_S_ZX]*r_n[K_F_ZX]+ s_n[K_S_ZY]*r_n[K_F_ZY]+ s_n[K_S_ZZ]*r_n[K_F_ZZ];

		rot_stress(ielem, K_S_XX) = r_n[K_F_XX]*t[0] + r_n[K_F_XY]*t[1] + r_n[K_F_XZ]*t[2];
		rot_stress(ielem, K_S_YY) = r_n[K_F_YX]*t[3] + r_n[K_F_YY]*t[4] + r_n[K_F_YZ]*t[5];
		rot_stress(ielem, K_S_ZZ) = r_n[K_F_ZX]*t[6] + r_n[K_F_ZY]*t[7] + r_n[K_F_ZZ]*t[8];

		rot_stress(ielem, K_S_XY) = r_n[K_F_XX]*t[3] + r_n[K_F_XY]*t[4] + r_n[K_F_XZ]*t[5];
		rot_stress(ielem, K_S_YZ) = r_n[K_F_YX]*t[6] + r_n[K_F_YY]*t[7] + r_n[K_F_YZ]*t[8];
		rot_stress(ielem, K_S_ZX) = r_n[K_F_ZX]*t[0] + r_n[K_F_ZY]*t[1] + r_n[K_F_ZZ]*t[2];

	}

	KOKKOS_MACRO_DEVICE_FUNCTION
	Scalar comp_aspect(int ielem)const {

		return 6.0 * 12.0 *  mid_vol(ielem) /
		    (gradop12(ielem, 0, 0) * gradop12(ielem, 0, 0)+
		     gradop12(ielem, 0, 1) * gradop12(ielem, 0, 1)+
		     gradop12(ielem, 0, 2) * gradop12(ielem, 0, 2)+
		     gradop12(ielem, 0, 3) * gradop12(ielem, 0, 3)+
		     gradop12(ielem, 0, 4) * gradop12(ielem, 0, 4)+
		     gradop12(ielem, 0, 5) * gradop12(ielem, 0, 5)+
		     gradop12(ielem, 0, 6) * gradop12(ielem, 0, 6)+
		     gradop12(ielem, 0, 7) * gradop12(ielem, 0, 7)+
		     gradop12(ielem, 1, 0) * gradop12(ielem, 1, 0)+
		     gradop12(ielem, 1, 1) * gradop12(ielem, 1, 1)+
		     gradop12(ielem, 1, 2) * gradop12(ielem, 1, 2)+
		     gradop12(ielem, 1, 3) * gradop12(ielem, 1, 3)+
		     gradop12(ielem, 1, 4) * gradop12(ielem, 1, 4)+
		     gradop12(ielem, 1, 5) * gradop12(ielem, 1, 5)+
		     gradop12(ielem, 1, 6) * gradop12(ielem, 1, 6)+
		     gradop12(ielem, 1, 7) * gradop12(ielem, 1, 7)+
		     gradop12(ielem, 2, 0) * gradop12(ielem, 2, 0)+
		     gradop12(ielem, 2, 1) * gradop12(ielem, 2, 1)+
		     gradop12(ielem, 2, 2) * gradop12(ielem, 2, 2)+
		     gradop12(ielem, 2, 3) * gradop12(ielem, 2, 3)+
		     gradop12(ielem, 2, 4) * gradop12(ielem, 2, 4)+
		     gradop12(ielem, 2, 5) * gradop12(ielem, 2, 5)+
		     gradop12(ielem, 2, 6) * gradop12(ielem, 2, 6)+
		     gradop12(ielem, 2, 7) * gradop12(ielem, 2, 7));

	}

	KOKKOS_MACRO_DEVICE_FUNCTION
	void comp_force(int ielem, Scalar fac1, Scalar fac2, Scalar * total_stress12th)const {

  
    //  NKC, does Presto Scalar need these spin rate terms?  Pronto appears to have dumped them....
		const Scalar dwxy = dt * vorticity(ielem, 0);
		const Scalar dwyz = dt * vorticity(ielem, 1);
		const Scalar dwzx = dt * vorticity(ielem, 2);

    //  Compute new hourglass resitance by the old rotated hourglass resitance plus a hourglass rate term
		Scalar hg_resist_total[12];
		Scalar hg_temp[8];


		if (!scaleHGRotation){

			for(int i = 0; i < 4; ++i) {

				Scalar hg_rate_0 = 0, hg_rate_1 = 0, hg_rate_2 = 0;
				Scalar hg_resist_old_0 = hg_resist(ielem, (3 * i) + 0, 0);
				Scalar hg_resist_old_1 = hg_resist(ielem, (3 * i) + 1, 0);
				Scalar hg_resist_old_2 = hg_resist(ielem, (3 * i) + 2, 0);

				for(int j = 0; j < 8; j++){

					hg_temp[j] = hgop(ielem, i * 8 + j, 1);

				}

				for(int j = 0; j < 8; j++){

					hg_rate_0 += hg_temp[j] * velocity(ielem, 0, j);
					hg_rate_1 += hg_temp[j] * velocity(ielem, 1, j);
					hg_rate_2 += hg_temp[j] * velocity(ielem, 2, j);

				}

				hg_resist(ielem, (i * 3) + 0, 1)   = hg_resist_old_0 + dwxy*hg_resist_old_1 - dwzx*hg_resist_old_2 + fac1* hg_rate_0 ;
				hg_resist_total[(i * 3) + 0]      = hg_resist(ielem, (i * 3) + 0, 1) + fac2* hg_rate_0 ;

				hg_resist(ielem, (i * 3) + 1, 1)   = hg_resist_old_1 - dwxy*hg_resist_old_0 + dwyz*hg_resist_old_2 + fac1* hg_rate_1 ;
				hg_resist_total[(i * 3) + 1]       = hg_resist(ielem, (i * 3) + 1, 1) + fac2* hg_rate_1 ;

				hg_resist(ielem, (i * 3) + 2, 1)   = hg_resist_old_2 + dwzx*hg_resist_old_0 - dwyz*hg_resist_old_1 + fac1* hg_rate_2 ;
				hg_resist_total[(i * 3) + 2]       = hg_resist(ielem, (i * 3) + 2, 1) + fac2* hg_rate_2 ;

			}

		} 

		else {

			for(int i = 0; i < 4; ++i) {

				Scalar hg_rate_0 = 0, hg_rate_1 = 0, hg_rate_2 = 0;
				Scalar hg_resist_old_0 = hg_resist(ielem, (3 * i) + 0, 0);
				Scalar hg_resist_old_1 = hg_resist(ielem, (3 * i) + 1, 0);
				Scalar hg_resist_old_2 = hg_resist(ielem, (3 * i) + 2, 0);

				for(int j = 0; j < 8; j++){

					hg_temp[j] = hgop(ielem, i * 8 + j, 1);

				}

				for(int j = 0; j < 8; j++){

					hg_rate_0 += hg_temp[j] * velocity(ielem, 0, j);
					hg_rate_1 += hg_temp[j] * velocity(ielem, 1, j);
					hg_rate_2 += hg_temp[j] * velocity(ielem, 2, j);

				}

				const Scalar rot_hg_resist_old_0 = hg_resist_old_0 + dwxy*hg_resist_old_1 - dwzx*hg_resist_old_2;
				const Scalar rot_hg_resist_old_1 = hg_resist_old_1 - dwxy*hg_resist_old_0 + dwyz*hg_resist_old_2;
				const Scalar rot_hg_resist_old_2 = hg_resist_old_2 + dwzx*hg_resist_old_0 - dwyz*hg_resist_old_1;

				Scalar fnorm = 	rot_hg_resist_old_0 *rot_hg_resist_old_0  + 
								rot_hg_resist_old_1 *rot_hg_resist_old_1  + 
								rot_hg_resist_old_2 *rot_hg_resist_old_2 ;

				if (fnorm > 1.e-30){

					fnorm = sqrt ( (hg_resist_old_0*hg_resist_old_0 +
					hg_resist_old_1*hg_resist_old_1 +
					hg_resist_old_2*hg_resist_old_2) / fnorm );
					hg_resist(ielem, i * 3 + 0, 1) = fnorm*rot_hg_resist_old_0 + fac1*hg_rate_0 ;
					hg_resist(ielem, i * 3 + 1, 1) = fnorm*rot_hg_resist_old_1 + fac1*hg_rate_1 ;
					hg_resist(ielem, i * 3 + 2, 1) = fnorm*rot_hg_resist_old_2 + fac1*hg_rate_2 ;

				} 

				else {

					hg_resist(ielem, i * 3 + 0, 1) = rot_hg_resist_old_0 + fac1*hg_rate_0 ;
					hg_resist(ielem, i * 3 + 1, 1) = rot_hg_resist_old_1 + fac1*hg_rate_1 ;
					hg_resist(ielem, i * 3 + 2, 1) = rot_hg_resist_old_2 + fac1*hg_rate_2 ;

				}

				hg_resist_total[(i * 3) + 0] = hg_resist(ielem, (i * 3) + 0, 1) + fac2*hg_rate_0 ;
				hg_resist_total[(i * 3) + 1] = hg_resist(ielem, (i * 3) + 1, 1) + fac2*hg_rate_1 ;
				hg_resist_total[(i * 3) + 2] = hg_resist(ielem, (i * 3) + 2, 1) + fac2*hg_rate_2 ;

			}

		}

		Scalar hg_force_0[8];
		Scalar hg_force_1[8];
		Scalar hg_force_2[8];

		hg_energy(ielem) = 0.0;
		intern_energy(ielem) = 0.0;

		for(int i = 0; i < 8; ++i) {

			hg_force_0[i] =
			   (hg_resist_total[HG_X1] * hgop(ielem, i +  0, 0) +
				hg_resist_total[HG_X2] * hgop(ielem, i +  8, 0) +
				hg_resist_total[HG_X3] * hgop(ielem, i + 16, 0) +
				hg_resist_total[HG_X4] * hgop(ielem, i + 24, 0));

			hg_force_1[i] =
			   (hg_resist_total[HG_Y1] * hgop(ielem, i +  0, 0) +
				hg_resist_total[HG_Y2] * hgop(ielem, i +  8, 0) +
				hg_resist_total[HG_Y3] * hgop(ielem, i + 16, 0) +
				hg_resist_total[HG_Y4] * hgop(ielem, i + 24, 0));

			hg_force_2[i] =
			   (hg_resist_total[HG_Z1] * hgop(ielem, i +  0, 0) +
				hg_resist_total[HG_Z2] * hgop(ielem, i +  8, 0) +
				hg_resist_total[HG_Z3] * hgop(ielem, i + 16, 0) +
				hg_resist_total[HG_Z4] * hgop(ielem, i + 24, 0));

		}

		for(int i = 0; i < 8; ++i) {
			force_new(ielem, 0, i) =
				total_stress12th[K_S_XX] * gradop12(ielem, 0, i) +
				total_stress12th[K_S_XY] * gradop12(ielem, 1, i) +
				total_stress12th[K_S_XZ] * gradop12(ielem, 2, i) + hg_force_0[i] ;

			force_new(ielem, 1, i) =
				total_stress12th[K_S_YX] * gradop12(ielem, 0, i) +
				total_stress12th[K_S_YY] * gradop12(ielem, 1, i) +
				total_stress12th[K_S_YZ] * gradop12(ielem, 2, i) + hg_force_1[i] ;

			force_new(ielem, 2, i) =
				total_stress12th[K_S_ZX] * gradop12(ielem, 0, i) +
				total_stress12th[K_S_ZY] * gradop12(ielem, 1, i) +
				total_stress12th[K_S_ZZ] * gradop12(ielem, 2, i) + hg_force_2[i] ;

			hg_energy(ielem)  +=	hg_force_0[i]   *velocity(ielem, 0, i) + \
									hg_force_1[i]   *velocity(ielem, 1, i) + \
									hg_force_2[i]   *velocity(ielem, 2, i);

			intern_energy(ielem) +=	force_new(ielem, 0, i)*velocity(ielem, 0, i) + \
									force_new(ielem, 1, i)*velocity(ielem, 1, i) + \
									force_new(ielem, 2, i)*velocity(ielem, 2, i);

		}

	}

	KOKKOS_MACRO_DEVICE_FUNCTION
    void operator()( int ielem )const {

		comp_grad(ielem);
		comp_hgop(ielem);

		Scalar fac1_pre = dt * hg_stiff * 0.0625;
		Scalar shr = elem_shrmod(ielem) = two_mu(ielem);
		Scalar dil = elem_dilmod(ielem) = ( bulk_mod(ielem) + 2.0*shr ) * (1.0 / 3.0);


		Scalar aspect = comp_aspect(ielem);

		Scalar inv_aspect = 1.0 / aspect;

		Scalar dtrial = sqrt(elem_mass(ielem) * aspect / dil);
		Scalar traced = (rot_stret(ielem, 0) + rot_stret(ielem, 1) + rot_stret(ielem, 2));

		Scalar eps = traced < 0 ? (linBulkVisc - quadBulkVisc * traced * dtrial) : linBulkVisc ;

		Scalar bulkq = eps * dil * dtrial * traced;

		Scalar cur_time_step = dtrial * ( sqrt( 1.0 + eps * eps) - eps);

		elem_t_step(ielem) = cur_time_step;

		rotate_tensor_backward(ielem);		

		Scalar total_stress12th[6];
		total_stress12th[0] = ONE12TH*(rot_stress(ielem, 0) + bulkq);
		total_stress12th[1] = ONE12TH*(rot_stress(ielem, 1) + bulkq);
		total_stress12th[2] = ONE12TH*(rot_stress(ielem, 2) + bulkq);
		total_stress12th[3] = ONE12TH*(rot_stress(ielem, 3));
		total_stress12th[4] = ONE12TH*(rot_stress(ielem, 4));
		total_stress12th[5] = ONE12TH*(rot_stress(ielem, 5));

		Scalar fac1 = fac1_pre * shr * inv_aspect;
		Scalar fac2 = hg_visc * sqrt(shr * elem_mass(ielem) * inv_aspect);

		comp_force(ielem, fac1, fac2, total_stress12th);

	}

};

#endif

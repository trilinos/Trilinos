#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include <impl/Kokkos_Preprocessing_macros.hpp>


/********************************************************/

#define numNodesPerElem 8 /* don't change */
#define numGaussPointsPerDim 2
#define spatialDim 3

#define ONE8TH 0.125

#define X2 0.577350269
#define X3 0.774596669
#define W1 0.555555556
#define W2 0.888888889


/********************************************************/


template< typename Scalar , class DeviceType >
struct assembleFE;

template<typename Scalar>
struct assembleFE<Scalar, KOKKOS_MACRO_DEVICE>{

	typedef KOKKOS_MACRO_DEVICE     					device_type;
  	typedef Kokkos::MDArrayView<Scalar,device_type> 	double_type;
  	typedef Kokkos::MDArrayView<int,device_type> 		int_type;

	double_type element_stiffness;
	double_type element_vectors;
	double_type element_coords;

	int global_x, global_y, global_z;

  	assembleFE( 	double_type arg_element_stiffness, 
					double_type arg_element_vectors,
					double_type arg_element_coords,  
					int x, 
					int y, 
					int z ){
    	
		element_stiffness = arg_element_stiffness;
		element_vectors = arg_element_vectors;
		element_coords = arg_element_coords;

		global_x = x;
		global_y = y;
		global_z = z;

	}


	KOKKOS_MACRO_DEVICE_FUNCTION
  	void load_coordinates(int ielem, Scalar * coords_x, Scalar * coords_y, Scalar * coords_z) const{

		coords_x[0] = element_coords(ielem, 0, 0);
		coords_x[1] = element_coords(ielem, 0, 1);
		coords_x[2] = element_coords(ielem, 0, 2);
		coords_x[3] = element_coords(ielem, 0, 3);
		coords_x[4] = element_coords(ielem, 0, 4);
		coords_x[5] = element_coords(ielem, 0, 5);
		coords_x[6] = element_coords(ielem, 0, 6);
		coords_x[7] = element_coords(ielem, 0, 7);

		coords_y[0] = element_coords(ielem, 1, 0);
		coords_y[1] = element_coords(ielem, 1, 1);
		coords_y[2] = element_coords(ielem, 1, 2);
		coords_y[3] = element_coords(ielem, 1, 3);
		coords_y[4] = element_coords(ielem, 1, 4);
		coords_y[5] = element_coords(ielem, 1, 5);
		coords_y[6] = element_coords(ielem, 1, 6);
		coords_y[7] = element_coords(ielem, 1, 7);

		coords_z[0] = element_coords(ielem, 2, 0);
		coords_z[1] = element_coords(ielem, 2, 1);
		coords_z[2] = element_coords(ielem, 2, 2);
		coords_z[3] = element_coords(ielem, 2, 3);
		coords_z[4] = element_coords(ielem, 2, 4);
		coords_z[5] = element_coords(ielem, 2, 5);
		coords_z[6] = element_coords(ielem, 2, 6);
		coords_z[7] = element_coords(ielem, 2, 7);

	}

	template<int N>
	KOKKOS_MACRO_DEVICE_FUNCTION
  	void gauss_pts(Scalar * pts, Scalar * wts)const {

		switch(N) {
			case 1:
				pts[0] = 0.0; wts[0] = 2.0;
				break;
			case 2:
				pts[0] = -X2; wts[0] = 1.0;
				pts[1] =  X2; wts[1] = 1.0;
				break;
			case 3:
				pts[0] =  -X3;  wts[0] = W1;
				pts[1] =  0.0;  wts[1] = W2;
				pts[2] =   X3;  wts[2] = W1;
				break;
			default:
				break;
		}

	}

	KOKKOS_MACRO_DEVICE_FUNCTION
  	void gradients(Scalar * x, Scalar * values_per_fn)const {

	//	partial derivatives of shape functions with respect
	//	to master element coordinate system
		const Scalar u = 1.0 - x[0];
		const Scalar v = 1.0 - x[1];
		const Scalar w = 1.0 - x[2];

		const Scalar up1 = 1.0 + x[0];
		const Scalar vp1 = 1.0 + x[1];
		const Scalar wp1 = 1.0 + x[2];

		//fn 0
		values_per_fn[ 0] = -ONE8TH *  v   *  wp1;
		values_per_fn[ 1] = -ONE8TH *  u   *  wp1;
		values_per_fn[ 2] =  ONE8TH *  u   *  v;
		//fn 1
		values_per_fn[ 3] =  ONE8TH *  v   *  wp1;
		values_per_fn[ 4] = -ONE8TH *  up1 *  wp1;
		values_per_fn[ 5] =  ONE8TH *  up1 *  v;
		//fn 2
		values_per_fn[ 6] =  ONE8TH *  v   *  w;
		values_per_fn[ 7] = -ONE8TH *  up1 *  w;
		values_per_fn[ 8] = -ONE8TH *  up1 *  v;
		//fn 3
		values_per_fn[ 9] = -ONE8TH *  v   *  w;
		values_per_fn[10] = -ONE8TH *  u   *  w;
		values_per_fn[11] = -ONE8TH *  u   *  v;
		//fn 4
		values_per_fn[12] = -ONE8TH *  vp1 *  wp1;
		values_per_fn[13] =  ONE8TH *  u   *  wp1;
		values_per_fn[14] =  ONE8TH *  u   *  vp1;
		//fn 5
		values_per_fn[15] =  ONE8TH *  vp1 *  wp1;
		values_per_fn[16] =  ONE8TH *  up1 *  wp1;
		values_per_fn[17] =  ONE8TH *  up1 *  vp1;
		//fn 6
		values_per_fn[18] =  ONE8TH *  vp1 *  w;
		values_per_fn[19] =  ONE8TH *  up1 *  w;
		values_per_fn[20] = -ONE8TH *  up1 *  vp1;
		//fn 7
		values_per_fn[21] = -ONE8TH *  vp1 *  w;
		values_per_fn[22] =  ONE8TH *  u   *  w;
		values_per_fn[23] = -ONE8TH *  u   *  vp1;

	}

	KOKKOS_MACRO_DEVICE_FUNCTION
  	void gradients_and_jacobian(	Scalar * pt, 
									Scalar * x, 
									Scalar * y, 
									Scalar * z, 
									Scalar * grad_vals, 
									Scalar * J) const{

		J[0] = 0.0;
		J[1] = 0.0;
		J[2] = 0.0;
		J[3] = 0.0;
		J[4] = 0.0;
		J[5] = 0.0;
		J[6] = 0.0;
		J[7] = 0.0;
		J[8] = 0.0;

		gradients(pt, grad_vals);

		int i_X_spatialDim = 0;

		for(int i = 0; i < 8; ++i) {
			J[0] += grad_vals[i_X_spatialDim+0]*x[i];
			J[1] += grad_vals[i_X_spatialDim+0]*y[i];
			J[2] += grad_vals[i_X_spatialDim+0]*z[i];

			J[3] += grad_vals[i_X_spatialDim+1]*x[i];
			J[4] += grad_vals[i_X_spatialDim+1]*y[i];
			J[5] += grad_vals[i_X_spatialDim+1]*z[i];

			J[6] += grad_vals[i_X_spatialDim+2]*x[i];
			J[7] += grad_vals[i_X_spatialDim+2]*y[i];
			J[8] += grad_vals[i_X_spatialDim+2]*z[i];

			i_X_spatialDim += spatialDim;
		}

	}

	KOKKOS_MACRO_DEVICE_FUNCTION
  	void inverse_and_determinant3x3(	Scalar * J,
										Scalar * invJ, 
										Scalar & detJ)const {

		Scalar J00 = J[0];
		Scalar J01 = J[1];
		Scalar J02 = J[2];

		Scalar J10 = J[3];
		Scalar J11 = J[4];
		Scalar J12 = J[5];

		Scalar J20 = J[6];
		Scalar J21 = J[7];
		Scalar J22 = J[8];

		Scalar term0 = J22*J11 - J21*J12;
		Scalar term1 = J22*J01 - J21*J02;
		Scalar term2 = J12*J01 - J11*J02;

		detJ = J00*term0 - J10*term1 + J20*term2;

		Scalar inv_detJ = 1.0/detJ;

		invJ[0] =  term0*inv_detJ;
		invJ[1] = -term1*inv_detJ;
		invJ[2] =  term2*inv_detJ;

		invJ[3] = -(J22*J10 - J20*J12)*inv_detJ;
		invJ[4] =  (J22*J00 - J20*J02)*inv_detJ;
		invJ[5] = -(J12*J00 - J10*J02)*inv_detJ;

		invJ[6] =  (J21*J10 - J20*J11)*inv_detJ;
		invJ[7] = -(J21*J00 - J20*J01)*inv_detJ;
		invJ[8] =  (J11*J00 - J10*J01)*inv_detJ;

	}

	KOKKOS_MACRO_DEVICE_FUNCTION
  	void matTransMat3x3_X_3xn(const Scalar * A, int n, const Scalar * B, Scalar * C)const {

		//A is 3x3, B is 3xn. So C is also 3xn.
		//A,B,C are all assumed to be ordered such that columns are contiguous.

		Scalar* Cj = C;
		const Scalar* Bj = B;

		for(int j=0; j<n; ++j) {
			Cj[0] = A[0]*Bj[0] + A[1]*Bj[1] + A[2]*Bj[2];
			Cj[1] = A[3]*Bj[0] + A[4]*Bj[1] + A[5]*Bj[2];
			Cj[2] = A[6]*Bj[0] + A[7]*Bj[1] + A[8]*Bj[2];
			Bj += 3;
			Cj += 3;
		}

	}

	KOKKOS_MACRO_DEVICE_FUNCTION
  	Scalar determinant3x3(const Scalar* J) const{	

		//hardwired "3x3" in function-name allows us to assume that J has length 9:

		Scalar J00 = J[0];
		Scalar J01 = J[1];
		Scalar J02 = J[2];

		Scalar J10 = J[3];
		Scalar J11 = J[4];
		Scalar J12 = J[5];

		Scalar J20 = J[6];
		Scalar J21 = J[7];
		Scalar J22 = J[8];

		Scalar term0 = J22*J11 - J21*J12;
		Scalar term1 = J22*J01 - J21*J02;
		Scalar term2 = J12*J01 - J11*J02;

		Scalar detJ = J00*term0 - J10*term1 + J20*term2;

		return detJ;

	}

	KOKKOS_MACRO_DEVICE_FUNCTION
  	void shape_fns(const Scalar* x, Scalar* values_at_nodes)const {

		//dangerous assumption that values_at_nodes has length numNodesPerElem
		//also, non-general assumption that x has length 3 (hard-coded spatialDim)
		const Scalar u = 1.0 - x[0];
		const Scalar v = 1.0 - x[1];
		const Scalar w = 1.0 - x[2];

		const Scalar up1 = 1.0 + x[0];
		const Scalar vp1 = 1.0 + x[1];
		const Scalar wp1 = 1.0 + x[2];

		values_at_nodes[0] = ONE8TH *   u *   v * wp1;
		values_at_nodes[1] = ONE8TH * up1 *   v * wp1;
		values_at_nodes[2] = ONE8TH * up1 *   v *   w;
		values_at_nodes[3] = ONE8TH *   u *   v *   w;
		values_at_nodes[4] = ONE8TH *   u * vp1 * wp1;
		values_at_nodes[5] = ONE8TH * up1 * vp1 * wp1;
		values_at_nodes[6] = ONE8TH * up1 * vp1 *   w;
		values_at_nodes[7] = ONE8TH *   u * vp1 *   w;
	
	}

	KOKKOS_MACRO_DEVICE_FUNCTION
  	void diffusionMatrix(	int ielem , 
							Scalar * x, 
							Scalar * y, 
							Scalar * z)const {

		Scalar gpts[numGaussPointsPerDim];
		Scalar gwts[numGaussPointsPerDim];

		for(int i = 0; i < 8; i++){
			for(int j = 0; j < 8; j++){
				element_stiffness(ielem, i, j) = 0.0;
			}
		}


		gauss_pts<numGaussPointsPerDim>(gpts, gwts);


		Scalar detJ = 0.0;
		const Scalar k = 2.0;

		Scalar pt[spatialDim];

		Scalar J[spatialDim*spatialDim];
		Scalar invJ[spatialDim*spatialDim];

		Scalar grad_vals[numNodesPerElem*spatialDim];
//		Scalar invJ_grad_vals[numNodesPerElem*spatialDim];

		for(size_t ig=0; ig<numGaussPointsPerDim; ++ig) {
			pt[0] = gpts[ig];
			Scalar wi = gwts[ig];

			for(size_t jg=0; jg<numGaussPointsPerDim; ++jg) {
				pt[1] = gpts[jg];
				Scalar wi_wj = wi*gwts[jg];

				for(size_t kg=0; kg<numGaussPointsPerDim; ++kg) {
					pt[2] = gpts[kg];
					Scalar wi_wj_wk = wi_wj*gwts[kg];

					gradients_and_jacobian(pt, x, y, z, grad_vals, J);

					inverse_and_determinant3x3(J, invJ, detJ);

					Scalar k_detJ_wi_wj_wk = k*detJ*wi_wj_wk;


					Scalar dpsidx[8], dpsidy[8], dpsidz[8];


					for(int i = 0; i < 8; i++){

						dpsidx[i] = grad_vals[i * 3 + 0] * invJ[0] + 
									grad_vals[i * 3 + 1] * invJ[1] + 
									grad_vals[i * 3 + 2] * invJ[2];
						dpsidy[i] = grad_vals[i * 3 + 0] * invJ[3] + 
									grad_vals[i * 3 + 1] * invJ[4] + 
									grad_vals[i * 3 + 2] * invJ[5];
						dpsidz[i] = grad_vals[i * 3 + 0] * invJ[6] + 
									grad_vals[i * 3 + 1] * invJ[7] + 
									grad_vals[i * 3 + 2] * invJ[8];

					}


					for(int m = 0; m < numNodesPerElem; m++) {
						for(int n = 0; n < numNodesPerElem; n++) {

							element_stiffness(ielem, m, n) += k_detJ_wi_wj_wk * \
								((dpsidx[m] * dpsidx[n]) + 
								 (dpsidy[m] * dpsidy[n]) +
								 (dpsidz[m] * dpsidz[n]));						

						}
					}

				}//for kg
			}//for jg
		}//for ig


	}


	KOKKOS_MACRO_DEVICE_FUNCTION
  	void sourceVector( int ielem , Scalar * x, Scalar * y, Scalar * z)const {

		Scalar gpts[numGaussPointsPerDim];
		Scalar gwts[numGaussPointsPerDim];

		Scalar pt[spatialDim];
		Scalar psi[numNodesPerElem];
		Scalar J[spatialDim*spatialDim];
		Scalar grad_vals[numNodesPerElem*spatialDim];

		element_vectors(ielem, 0) = 0.0;
		element_vectors(ielem, 1) = 0.0;
		element_vectors(ielem, 2) = 0.0;
		element_vectors(ielem, 3) = 0.0;
		element_vectors(ielem, 4) = 0.0;
		element_vectors(ielem, 5) = 0.0;
		element_vectors(ielem, 6) = 0.0;
		element_vectors(ielem, 7) = 0.0;

		gauss_pts<numGaussPointsPerDim>(gpts, gwts);

		Scalar Q = 1.0;

		for(size_t ig=0; ig<numGaussPointsPerDim; ++ig) {
			pt[0] = gpts[ig];
			Scalar wi = gwts[ig];

			for(size_t jg=0; jg<numGaussPointsPerDim; ++jg) {
				pt[1] = gpts[jg];
				Scalar wj = gwts[jg];

				for(size_t kg=0; kg<numGaussPointsPerDim; ++kg) {
					pt[2] = gpts[kg];
					Scalar wk = gwts[kg];

					shape_fns(pt, psi);
					gradients_and_jacobian(pt, x, y, z, grad_vals, J);
					Scalar detJ = determinant3x3(J);

					Scalar term = Q*detJ*wi*wj*wk;

					for(int i=0; i<numNodesPerElem; ++i) {
					  element_vectors(ielem, i) += psi[i]*term;
					}
				}
			}
		}


	}

	
	KOKKOS_MACRO_DEVICE_FUNCTION
  	void operator()( int ielem )const {

		Scalar coords_x[8], coords_y[8], coords_z[8];

		load_coordinates(ielem, coords_x, coords_y, coords_z);

		diffusionMatrix(ielem, coords_x, coords_y, coords_z);
		sourceVector(ielem, coords_x, coords_y, coords_z);

	}


};// assembleFE


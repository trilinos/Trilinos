/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include <impl/KokkosArray_Preprocessing_macros.hpp>

/********************************************************/

#define numNodesPerElem 8 /* don't change */
#define spatialDim 3

/********************************************************/

//----------------------------------------------------------------------------

template< typename Scalar , unsigned Dim , unsigned N >
struct TensorIntegration ;

template<typename Scalar >
struct TensorIntegration<Scalar,1,1> {
  Scalar pts[1] ;
  Scalar wts[1] ;

  TensorIntegration() { pts[0] = 0 ; wts[0] = 2 ; }
};

template<typename Scalar >
struct TensorIntegration<Scalar,1,2>
{
  Scalar pts[2] ;
  Scalar wts[2] ;

  TensorIntegration()
  {
    const Scalar x2 = 0.577350269 ;
    pts[0] = -x2; wts[0] = 1.0;
    pts[1] =  x2; wts[1] = 1.0;
  }
};

template<typename Scalar >
struct TensorIntegration<Scalar,1,3>
{
  Scalar pts[3] ;
  Scalar wts[3] ;

  TensorIntegration()
  {
    const Scalar x3 = 0.774596669 ;
    const Scalar w1 = 0.555555556 ;
    const Scalar w2 = 0.888888889 ;
    pts[0] =  -x3 ;  wts[0] = w1 ;
    pts[1] =    0 ;  wts[1] = w2 ;
    pts[2] =   x3 ;  wts[2] = w1 ;
  }
};

template< typename Scalar , unsigned Order >
struct TensorIntegration<Scalar,3,Order>
{
  static const unsigned N = Order * Order * Order ;

  Scalar pts[N][3] ;
  Scalar wts[N];

  TensorIntegration()
  {
    TensorIntegration<Scalar,1,Order> oneD ;

    unsigned n = 0 ;
    for ( unsigned k = 0 ; k < Order ; ++k ) {
    for ( unsigned j = 0 ; j < Order ; ++j ) {
    for ( unsigned i = 0 ; i < Order ; ++i , ++n ) {
      pts[n][0] = oneD.pts[i] ;
      pts[n][1] = oneD.pts[j] ;
      pts[n][2] = oneD.pts[k] ;
      wts[n] = oneD.wts[i] * oneD.wts[j] * oneD.wts[k] ;
    }}}
  }
};

//----------------------------------------------------------------------------

template< typename Scalar >
struct ShapeFunctionEvaluation {
  static const unsigned PointCount = 9 ;
  static const unsigned FunctionCount = 8 ;
  static const unsigned SpatialDimension = 3 ;
  static const unsigned IntegrationOrder = 2 ;

  Scalar value   [ PointCount ][ FunctionCount ] ;
  Scalar gradient[ PointCount ][ FunctionCount * SpatialDimension ];
  Scalar weight  [ PointCount ];

  ShapeFunctionEvaluation()
  {
    const TensorIntegration< Scalar , SpatialDimension , IntegrationOrder > 
      integration ;

    const Scalar ONE8TH = 0.125 ;

    for ( unsigned i = 0 ; i < PointCount ; ++i ) {

      const Scalar u = 1.0 - integration.pts[i][0];
      const Scalar v = 1.0 - integration.pts[i][1];
      const Scalar w = 1.0 - integration.pts[i][2];

      const Scalar up1 = 1.0 + integration.pts[i][0];
      const Scalar vp1 = 1.0 + integration.pts[i][1];
      const Scalar wp1 = 1.0 + integration.pts[i][2];

      weight[i] = integration.wts[i] ;

      // Vaues:
      value[i][0] = ONE8TH *   u *   v *  w ;
      value[i][1] = ONE8TH * up1 *   v *  w ;
      value[i][2] = ONE8TH * up1 * vp1 *  w ;
      value[i][3] = ONE8TH *   u * vp1 *  w ;

      value[i][4] = ONE8TH *   u *   v *  wp1 ;
      value[i][5] = ONE8TH * up1 *   v *  wp1 ;
      value[i][6] = ONE8TH * up1 * vp1 *  wp1 ;
      value[i][7] = ONE8TH *   u * vp1 *  wp1 ;

      //fn 0 = u * v * w
      gradient[i][ 0] = ONE8TH * -1  *  v  *  w  ;
      gradient[i][ 1] = ONE8TH *  u  * -1  *  w  ;
      gradient[i][ 2] = ONE8TH *  u  *  v  * -1  ;

      //fn 1 = up1 * v * w
      gradient[i][ 3] = ONE8TH *  1  *  v  *  w  ;
      gradient[i][ 4] = ONE8TH * up1 * -1  *  w  ;
      gradient[i][ 5] = ONE8TH * up1 *  v  * -1  ;

      //fn 2 = up1 * vp1 * w
      gradient[i][ 6] = ONE8TH *  1  * vp1 *  w ;
      gradient[i][ 7] = ONE8TH * up1 *  1  *  w ;
      gradient[i][ 8] = ONE8TH * up1 * vp1 * -1 ;

      //fn 3 = u * vp1 * w
      gradient[i][ 9] = ONE8TH * -1 * vp1 *  w ;
      gradient[i][10] = ONE8TH *  u *  1  *  w ;
      gradient[i][11] = ONE8TH *  u * vp1 * -1 ;

      //fn 4 = u * v * wp1
      gradient[i][12] = ONE8TH * -1  *  v  * wp1 ;
      gradient[i][13] = ONE8TH *  u  * -1  * wp1 ;
      gradient[i][14] = ONE8TH *  u  *  v  *  1  ;

      //fn 5 = up1 * v * wp1
      gradient[i][15] = ONE8TH *  1  *  v  * wp1 ;
      gradient[i][16] = ONE8TH * up1 * -1  * wp1 ;
      gradient[i][17] = ONE8TH * up1 *  v  *  1  ;

      //fn 6 = up1 * vp1 * wp1
      gradient[i][18] = ONE8TH *  1  * vp1 * wp1 ;
      gradient[i][19] = ONE8TH * up1 *  1  * wp1 ;
      gradient[i][20] = ONE8TH * up1 * vp1 *  1 ;

      //fn 7 = u * vp1 * wp1
      gradient[i][21] = ONE8TH * -1 * vp1 * wp1 ;
      gradient[i][22] = ONE8TH *  u *  1  * wp1 ;
      gradient[i][23] = ONE8TH *  u * vp1 *  1 ;
    }
  }
};

//----------------------------------------------------------------------------

template< typename Scalar , typename ScalarCoord , class DeviceType >
struct ElementComp ;

template<typename Scalar , typename ScalarCoord >
struct ElementComp <Scalar, ScalarCoord, KOKKOS_MACRO_DEVICE> {

  typedef KOKKOS_MACRO_DEVICE                        device_type;
  typedef Scalar                                     scalar_type ;
  typedef device_type::size_type                     index_type ;
  typedef KokkosArray::MDArray<index_type,device_type>    index_array ;
  typedef KokkosArray::MDArray<Scalar,      device_type>  scalar_array ;
  typedef KokkosArray::MDArray<ScalarCoord, device_type>  coord_array ;
  typedef ShapeFunctionEvaluation< Scalar > shape_function_data ;

  static const unsigned SpatialDimension = shape_function_data::SpatialDimension ;
  static const unsigned FunctionCount    = shape_function_data::FunctionCount ;

  shape_function_data  shape_eval ;
  index_array          elem_node_ids ;
  coord_array          node_coords ;
  scalar_array         element_stiffness;
  scalar_array         element_vectors;
  Scalar               coeff_K ;
  Scalar               coeff_Q ;

  ElementComp ( const index_array  & arg_elem_node_ids ,
                const coord_array  & arg_node_coords ,
                const scalar_array & arg_element_stiffness , 
                const scalar_array & arg_element_vectors ,
                const Scalar       & arg_coeff_K ,
                const Scalar       & arg_coeff_Q )
  : shape_eval()
  , elem_node_ids( arg_elem_node_ids )
  , node_coords(   arg_node_coords )
  , element_stiffness( arg_element_stiffness )
  , element_vectors( arg_element_vectors )
  , coeff_K( arg_coeff_K )
  , coeff_Q( arg_coeff_Q )
  {}

  static const unsigned FLOPS_jacobian = FunctionCount * SpatialDimension * SpatialDimension * 2 ;

  KOKKOS_MACRO_DEVICE_FUNCTION
  void jacobian( const ScalarCoord * x, 
                 const ScalarCoord * y, 
                 const ScalarCoord * z, 
                 const Scalar * grad_vals, 
                 Scalar * J) const
  {
    int i_grad = 0 ;

    for(int i = 0; i < 8; ++i , i_grad += SpatialDimension ) {
      const Scalar g0 = grad_vals[ i_grad ];
      const Scalar g1 = grad_vals[ i_grad + 1 ];
      const Scalar g2 = grad_vals[ i_grad + 2 ];
      const Scalar x0 = x[i] ;
      const Scalar x1 = y[i] ;
      const Scalar x2 = z[i] ;

      J[0] += g0 * x0 ;
      J[1] += g0 * x1 ;
      J[2] += g0 * x2 ;

      J[3] += g1 * x0 ;
      J[4] += g1 * x1 ;
      J[5] += g1 * x2 ;

      J[6] += g2 * x0 ;
      J[7] += g2 * x1 ;
      J[8] += g2 * x2 ;
    }
  }

  static const unsigned FLOPS_inverse_and_det = 46 ;

  KOKKOS_MACRO_DEVICE_FUNCTION
  Scalar inverse_and_determinant3x3( Scalar * const J ) const
  {
    const Scalar J00 = J[0];
    const Scalar J01 = J[1];
    const Scalar J02 = J[2];

    const Scalar J10 = J[3];
    const Scalar J11 = J[4];
    const Scalar J12 = J[5];

    const Scalar J20 = J[6];
    const Scalar J21 = J[7];
    const Scalar J22 = J[8];

    const Scalar term0 = J22*J11 - J21*J12;
    const Scalar term1 = J22*J01 - J21*J02;
    const Scalar term2 = J12*J01 - J11*J02;

    const Scalar detJ = J00*term0 - J10*term1 + J20*term2;
    const Scalar inv_detJ = 1.0/detJ;

    J[0] =  term0*inv_detJ;
    J[1] = -term1*inv_detJ;
    J[2] =  term2*inv_detJ;

    J[3] = -(J22*J10 - J20*J12)*inv_detJ;
    J[4] =  (J22*J00 - J20*J02)*inv_detJ;
    J[5] = -(J12*J00 - J10*J02)*inv_detJ;

    J[6] =  (J21*J10 - J20*J11)*inv_detJ;
    J[7] = -(J21*J00 - J20*J01)*inv_detJ;
    J[8] =  (J11*J00 - J10*J01)*inv_detJ;

    return detJ ;
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

  static const unsigned FLOPS_contributeDiffusionMatrix = FunctionCount * ( 3 * 5 + FunctionCount * 7 ) ;

  KOKKOS_MACRO_DEVICE_FUNCTION
  void contributeDiffusionMatrix(
    const Scalar weight ,
    const Scalar grad_vals[] ,
    const Scalar invJ[] ,
    Scalar elem_stiff[][8] ) const
  {
    Scalar dpsidx[8], dpsidy[8], dpsidz[8];

    int i_grad = 0 ;
    for( unsigned i = 0; i < FunctionCount ; ++i , i_grad += 3 ) {
      const Scalar g0 = grad_vals[i_grad+0];
      const Scalar g1 = grad_vals[i_grad+1];
      const Scalar g2 = grad_vals[i_grad+2];

      dpsidx[i] = g0 * invJ[0] + g1 * invJ[1] + g2 * invJ[2];
      dpsidy[i] = g0 * invJ[3] + g1 * invJ[4] + g2 * invJ[5];
      dpsidz[i] = g0 * invJ[6] + g1 * invJ[7] + g2 * invJ[8];
    }

    for( unsigned m = 0; m < FunctionCount; m++) {
      for( unsigned n = 0; n < FunctionCount; n++) {

        elem_stiff[m][n] += weight * 
          ((dpsidx[m] * dpsidx[n]) + 
           (dpsidy[m] * dpsidy[n]) +
           (dpsidz[m] * dpsidz[n]));            
      }
    }
  }

  static const unsigned FLOPS_contributeSourceVector = FunctionCount * 2 ;

  KOKKOS_MACRO_DEVICE_FUNCTION
  void contributeSourceVector( const Scalar term ,
                               const Scalar psi[] ,
                               Scalar elem_vec[] ) const
  {
     for( unsigned i=0; i< FunctionCount ; ++i) {
       elem_vec[i] += psi[i] * term ;
     }
  }


  static const unsigned FLOPS_operator =
           shape_function_data::PointCount * ( 3
             + FLOPS_jacobian
             + FLOPS_inverse_and_det
             + FLOPS_contributeDiffusionMatrix
             + FLOPS_contributeSourceVector ) ;

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( int ielem )const {

    Scalar elem_vec[8] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };
    Scalar elem_stiff[8][8] =
      { { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } };

    ScalarCoord x[8], y[8], z[8];

    for ( int i = 0 ; i < 8 ; ++i ) {
      const int node_index = elem_node_ids( ielem , i );
      x[i] = node_coords( node_index , 0 );
      y[i] = node_coords( node_index , 1 );
      z[i] = node_coords( node_index , 2 );
    }

    // This loop could be parallelized; however,
    // it would require additional per-thread temporaries
    // of 'elem_vec' and 'elem_stiff' which would
    // consume more local memory and have to be reduced.

    for ( unsigned i = 0 ; i < shape_function_data::PointCount ; ++i ) {

      Scalar J[spatialDim*spatialDim] = { 0, 0, 0,  0, 0, 0,  0, 0, 0 };

      jacobian( x, y, z, shape_eval.gradient[i] , J );

      // Overwrite J with its inverse to save scratch memory space.
      const Scalar detJ_w   = shape_eval.weight[i] * inverse_and_determinant3x3(J);
      const Scalar k_detJ_w = coeff_K * detJ_w ;
      const Scalar Q_detJ_w = coeff_Q * detJ_w ;

      contributeDiffusionMatrix( k_detJ_w , shape_eval.gradient[i] , J , elem_stiff );

      contributeSourceVector( Q_detJ_w , shape_eval.value[i] , elem_vec );
    }

    for(int i=0; i<numNodesPerElem; ++i) {
      element_vectors(ielem, i) = elem_vec[i] ;
    }

    for(int i = 0; i < 8; i++){
      for(int j = 0; j < 8; j++){
        element_stiffness(ielem, i, j) = elem_stiff[i][j] ;
      }
    }
  }

};// ElementComp

#if 0
  // Possible experiment

  // Split operator to reduce temporaries
  // but increase global memory write & read.

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( int ielem )const {

    ScalarCoord x[8], y[8], z[8];

    for ( int i = 0 ; i < 8 ; ++i ) {
      const int node_index = elem_node_ids( ielem , i );
      x[i] = node_coords( node_index , 0 );
      y[i] = node_coords( node_index , 1 );
      z[i] = node_coords( node_index , 2 );
    }

    // This loop could be parallelized; however,
    // it would require additional per-thread temporaries
    // of 'elem_vec' and 'elem_stiff' which would
    // consume more local memory and have to be reduced.

    for ( int i = 0 ; i < shape_function_data::PointCount ; ++i ) {

      Scalar J[spatialDim*spatialDim] = { 0, 0, 0,  0, 0, 0,  0, 0, 0 };

      jacobian( x, y, z, shape_eval.gradient[i] , J );

      // Overwrite J with its inverse to save scratch memory space.
      jacobian_det(ielem,i) = inverse_and_determinant3x3(J);

      for ( int j = 0 ; j < 9 ; ++j ) {
        jacobian_inv(ielem,i,j) = J[j] ;
      }
    }
  }

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( int ielem )const {

    Scalar elem_vec[8] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };
    Scalar elem_stiff[8][8] =
      { { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } ,
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 } };

    // This loop could be parallelized; however,
    // it would require additional per-thread temporaries
    // of 'elem_vec' and 'elem_stiff' which would
    // consume more local memory and have to be reduced.

    for ( int i = 0 ; i < shape_function_data::PointCount ; ++i ) {

      Scalar J[spatialDim*spatialDim] ;

      for ( int j = 0 ; j < 9 ; ++j ) {
        J[j] = jacobian_inv(ielem,i,j);
      }

      // Overwrite J with its inverse to save scratch memory space.
      const Scalar detJ_w   = shape_eval.weight[i] * jacobian_det(ielem,i);
      const Scalar k_detJ_w = coeff_K * detJ_w ;
      const Scalar Q_detJ_w = coeff_Q * detJ_w ;

      contributeDiffusionMatrix( k_detJ_w , shape_eval.gradient[i] , J , elem_stiff );

      contributeSourceVector( Q_detJ_w , shape_eval.value[i] , elem_vec );
    }

    for(int i=0; i<numNodesPerElem; ++i) {
      element_vectors(ielem, i) = elem_vec[i] ;
    }

    for(int i = 0; i < 8; i++){
      for(int j = 0; j < 8; j++){
        element_stiffness(ielem, i, j) = elem_stiff[i][j] ;
      }
    }
  }

#endif


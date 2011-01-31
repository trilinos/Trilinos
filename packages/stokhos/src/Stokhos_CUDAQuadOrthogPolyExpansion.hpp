// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_CUDAQUADORTHOGPOLYEXPANSION_HPP
#define STOKHOS_CUDAQUADORTHOGPOLYEXPANSION_HPP

#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_STOKHOS_THRUST

#include "Stokhos_QuadOrthogPolyExpansion.hpp"
#include "Stokhos_CUDAStorage.hpp"
#include "Stokhos_Quadrature.hpp"

#include <thrust/device_vector.h>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

namespace Stokhos {

  /*! \brief Orthogonal polynomial expansions based on numerical quadrature
   * specialized to CUDA GPUs
   */
  template <>
  class QuadOrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> >: 
    public OrthogPolyExpansion<int, float, Stokhos::CUDAStorage<int, float> > {
  public:

    typedef Stokhos::CUDAStorage<int, float> node_type;
    typedef node_type::pointer pointer;
    typedef node_type::const_pointer const_pointer;
    typedef node_type::reference reference;
    typedef node_type::const_reference const_reference;

    //! Constructor
    QuadOrthogPolyExpansion(
    const Teuchos::RCP<const OrthogPolyBasis<int, float> >& basis,
    const Teuchos::RCP<const Stokhos::Sparse3Tensor<int, float> >& Cijk,
    const Teuchos::RCP<const Quadrature<int, float> >& quad,
    bool use_quad_for_times = false);

    //! Destructor
    virtual ~QuadOrthogPolyExpansion();

    //! Get expansion size
    int size() const { return sz; }

    //! Get basis
    Teuchos::RCP< const OrthogPolyBasis<int, float> > 
    getBasis() const {return basis; }

    //! Get triple product
    virtual Teuchos::RCP<const Sparse3Tensor<int, float> >
    getTripleProduct() const { return Cijk; }
 
    // Operations
    void unaryMinus(
      OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);

    void plusEqual(
      OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const float& x);
    void minusEqual(
      OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const float& x);
    void timesEqual(
      OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const float& x);
    void divideEqual(
      OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const float& x);

    void plusEqual(
      OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& x);
    void minusEqual(
      OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& x);
    void timesEqual(
      OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& x);
    void divideEqual(
      OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& x);

    void plus(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
              const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
              const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b);
    void plus(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
              const float& a, 
              const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b);
    void plus(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
              const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
              const float& b);
    void minus(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
               const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a,
               const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b);
    void minus(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
               const float& a, 
               const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b);
    void minus(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
               const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
               const float& b);
    void times(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
               const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
               const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b);
    void times(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
               const float& a, 
               const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b);
    void times(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
               const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
               const float& b);
    void divide(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
                const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
                const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b);
    void divide(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
                const float& a, 
                const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b);
    void divide(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
                const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
                const float& b);

    void exp(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void log(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void log10(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
               const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void sqrt(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
              const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void pow(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b);
    void pow(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
             const float& a, 
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b);
    void pow(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
             const float& b);
    void cos(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void sin(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void tan(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void cosh(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
              const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void sinh(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
              const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void tanh(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
              const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void acos(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
              const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void asin(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
              const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void atan(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
              const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void atan2(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
               const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a,
               const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b);
    void atan2(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
               const float& a, 
               const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b);
    void atan2(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
               const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
               const float& b);
    void acosh(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
               const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void asinh(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
               const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void atanh(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
               const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void abs(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void fabs(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
              const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);
    void max(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a,
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b);
    void max(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
             const float& a, 
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b);
    void max(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
             const float& b);
    void min(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a,
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b);
    void min(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
             const float& a, 
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b);
    void min(OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
             const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
             const float& b);

    template <typename FuncT>
    void nary_op(const FuncT& func,
		 OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c,
		 const OrthogPolyApprox<int, float, CUDAStorage<int,float> >** a);

    template <typename ExprT1, typename ExprT2>
    float compute_times_coeff(int k, const ExprT1& a, const ExprT2& b) const;

    template <typename ExprT1, typename ExprT2>
    float fast_compute_times_coeff(int k, const ExprT1& a, const ExprT2& b) const;

  private:

    // Prohibit copying
    QuadOrthogPolyExpansion(const QuadOrthogPolyExpansion&);

    // Prohibit Assignment
    QuadOrthogPolyExpansion& operator=(const QuadOrthogPolyExpansion& b);

  protected:

     //! Basis
    Teuchos::RCP<const OrthogPolyBasis<int, float> > basis;

    //! Short-hand for Cijk
    typedef Stokhos::Sparse3Tensor<int,float> Cijk_type;

    //! Triple-product tensor
    Teuchos::RCP<const Stokhos::Sparse3Tensor<int, float> > Cijk;

    //! Quadrature routine
    Teuchos::RCP<const Quadrature<int, float> > quad;

    //! Use quadrature for times functions
    bool use_quad_for_times;

    //! Expansions size
    int sz;

    //! Array of Quad points
    const Teuchos::Array< Teuchos::Array<float> >& quad_points;

    //! Array of Quad weights
    const Teuchos::Array<float>& quad_weights;

    //! Values of basis at Quad points
    const Teuchos::Array< Teuchos::Array<float> >& quad_values;

    //! Norms of basis vectors
    const Teuchos::Array<float>& norms;

    //! Number of Quad points
    int nqp;

    //! Temporary array for values of first argument at quad points
    thrust::device_vector<float> avals;

    //! Temporary array for values of second argument at quad points
    thrust::device_vector<float> bvals;

    //! Temporary array for values of n-ary arguments at quad points
    typedef Teuchos::Array< thrust::device_vector<float> > navals_type;
    Teuchos::Array< navals_type > navals;

    //! Temporary array for values of operation at quad points
    thrust::device_vector<float> fvals;

    //! Reshaped quad values into 1D array
    thrust::host_vector<float> host_qv;

    //! Quad values scaled by norms
    thrust::host_vector<float> host_sqv;

    //! Reshaped quad values into 1D array
    thrust::device_vector<float> qv;

    //! Quad values scaled by norms
    thrust::device_vector<float> sqv;

  public:

    //! Nonlinear unary function
    template <typename FuncT>
    void unary_op(
      const FuncT& func,
      OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a);

    //! Nonlinear binary function
    template <typename FuncT>
    void binary_op(
      const FuncT& func,
      OrthogPolyApprox<int, float, CUDAStorage<int,float> >& c, 
      const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& a, 
      const OrthogPolyApprox<int, float, CUDAStorage<int,float> >& b);
    
  }; // class QuadOrthogPolyExpansion

} // namespace Stokhos

#include "Stokhos_CUDAQuadOrthogPolyExpansionImp.hpp"

#endif // HAVE_STOKHOS_THRUST

#endif // STOKHOS_CUDAQUADORTHOGPOLYEXPANSION_HPP

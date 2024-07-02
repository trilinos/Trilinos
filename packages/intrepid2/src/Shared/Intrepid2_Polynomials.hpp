// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_Polynomials.hpp
    \brief  Free functions, callable from device code, that implement various polynomials useful in basis definitions.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_Polynomials_h
#define Intrepid2_Polynomials_h

#include "Intrepid2_Polylib.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

namespace Intrepid2
{
  namespace Polynomials
  {
    /*
     These polynomials are supplemental to those defined in the Polylib class; there is some overlap.
     We actually take advantage of the overlap in our verification tests, using the Polylib functions
     to verify the ones defined here.  Our interface here is a little simpler, and the functions are a little less
     general than those in Polylib.
     
     We define some polynomial functions that are useful in a variety of contexts.
     In particular, the (integrated) Legendre and Jacobi polynomials are useful in defining
     hierarchical bases.  See in particular:
     
     Federico Fuentes, Brendan Keith, Leszek Demkowicz, Sriram Nagaraj.
     "Orientation embedded high order shape functions for the exact sequence elements of all shapes."
     Computers & Mathematics with Applications, Volume 70, Issue 4, 2015, Pages 353-458, ISSN 0898-1221.
     https://doi.org/10.1016/j.camwa.2015.04.027.
     
     In this implementation, we take care to make minimal assumptions on both the containers
     and the scalar type.  The containers need to support a one-argument operator() for assignment and/or
     lookup (as appropriate).  The scalar type needs to support a cast from Intrepid2::ordinal_type, as well
     as standard arithmetic operations.  In particular, both 1-rank Kokkos::View and Kokkos::DynRankView are
     supported, as are C++ floating point types and Sacado scalar types.
     */
    
    /** \brief  Evaluate Legendre polynomials up to order n at a specified point.
        \param [out] outputValues - the view into which to place the output values (must have at least n+1 entries)
        \param [in] n - the maximum polynomial order of Legendre polynomials to compute
        \param [in] x - point at which to evaluate the polynomials
     These are defined for x in [-1,1].  See equation (2.7) in Fuentes et al.
     */
    template<typename OutputValueViewType, typename ScalarType>
    KOKKOS_INLINE_FUNCTION void legendreValues(OutputValueViewType outputValues, Intrepid2::ordinal_type n, ScalarType x)
    {
      if (n >= 0) outputValues(0) = 1.0;
      if (n >= 1) outputValues(1) = x;
      for (int i=2; i<=n; i++)
      {
        const ScalarType i_scalar = ScalarType(i);
        outputValues(i) = (2. - 1. / i_scalar) * x * outputValues(i-1) - (1. - 1. / i_scalar) * outputValues(i-2);
      }
    }
    
    /** \brief  Evaluate first derivatives of Legendre polynomials up to order n at a specified point, based on provided container with the values of the Legendre polynomials.
        \param [out] outputValues - the view into which to place the output values (must have at least n+1 entries)
        \param [in] legendreValues - the view containing previously evaluated Legendre polynomial values at the point (must have at least n entries)
        \param [in] n - the maximum polynomial order of Legendre polynomials to compute a derivative for
        \param [in] x - point at which to evaluate the derivatives.
        These are defined for x in [-1,1].
     */
    template<typename OutputValueViewType, typename ScalarType>
    KOKKOS_INLINE_FUNCTION void legendreDerivativeValues(OutputValueViewType outputValues, const OutputValueViewType legendreValues, Intrepid2::ordinal_type n, ScalarType x)
    {
      if (n >= 0) outputValues(0) = 0.0;
      if (n >= 1) outputValues(1) = 1.0;
      for (int i=2; i<=n; i++)
      {
        const ScalarType i_scalar = ScalarType(i);
        outputValues(i) = outputValues(i-2) + (2. * i_scalar - 1.) * legendreValues(i-1);
      }
    }
    
    // derivative values can be computed using the Legendre values
    // n: number of Legendre polynomials' derivative values to compute.  outputValues must have at least n+1 entries
    // x: value in [-1, 1]
    // dn: number of derivatives to take.  Must be >= 1.
    /** \brief  Evaluate the dnth derivative of Legendre polynomials up to order n at a specified point.
        \param [out] outputValues - the view into which to place the output values (must have at least n+1 entries)
        \param [in] n - the maximum polynomial order of Legendre polynomials to compute the derivative for
        \param [in] x - point at which to evaluate the derivatives (should be in [-1,1]).
        \param [in] dn - the order of differentiation.
        These are defined for x in [-1,1].
     */
    template<typename OutputValueViewType, typename PointScalarType>
    KOKKOS_INLINE_FUNCTION void legendreDerivativeValues(OutputValueViewType outputValues, Intrepid2::ordinal_type n, PointScalarType x, Intrepid2::ordinal_type dn)
    {
      const OutputValueViewType nullOutputScalarView;
      const double alpha = 0;
      const double beta = 0;
      const int numPoints = 1;
      
      using Layout = typename NaturalLayoutForType<PointScalarType>::layout;
      
      using UnmanagedPointScalarView = Kokkos::View<PointScalarType*, Layout, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
      UnmanagedPointScalarView pointView = UnmanagedPointScalarView(&x,numPoints);
      
      for (int i=0; i<=n; i++)
      {
        auto jacobiValue = Kokkos::subview(outputValues,Kokkos::pair<Intrepid2::ordinal_type,Intrepid2::ordinal_type>(i,i+1));
        jacobiValue(0) = 0.0;
        Intrepid2::Polylib::Serial::JacobiPolynomial(numPoints, pointView, jacobiValue, nullOutputScalarView, i-dn, alpha+dn, beta+dn);
        
        double scaleFactor = 1.0;
        for (int j=1; j<=dn; j++)
        {
          scaleFactor *= 0.5 * (j+alpha+beta+i);
        }
        
        outputValues(i) = jacobiValue(0) * scaleFactor;
      }
    }
    
    /** \brief  Evaluate shifted Legendre polynomials up to order n at a specified point in [0,1].
        \param [out] outputValues - the view into which to place the output values (must have at least n+1 entries)
        \param [in] n - the maximum polynomial order of Legendre polynomials to compute
        \param [in] x - point at which to evaluate the polynomials
     These are defined for x in [0,1].  See equation (2.12) in Fuentes et al.
     
     If these are defined by P_i, and the (unshifted) Legendre polynomials are given by ~P_i, we have:
     P_i(x) = ~P_i(2x-1)
     */
    template<typename OutputValueViewType, typename ScalarType>
    KOKKOS_INLINE_FUNCTION void shiftedLegendreValues(OutputValueViewType outputValues, Intrepid2::ordinal_type n, ScalarType x)
    {
      legendreValues(outputValues, n, 2.*x-1.);
    }
    
    /** \brief  Evaluate shifted, scaled Legendre polynomials up to order n at a specified point in [0,1].
        \param [out] outputValues - the view into which to place the output values (must have at least n+1 entries)
        \param [in] n - the maximum polynomial order of Legendre polynomials to compute
        \param [in] x - point at which to evaluate the polynomials
        \param [in] t - scaling parameter
     These are defined for x in [0,1].  See equation (2.14) in Fuentes et al.
     
     Shifted, scaled Legendre polynomials are given by
        P_i(x;t) = P_i(x/t) * t^i = ~P_i(2x-t;t).
     */
    template<typename OutputValueViewType, typename ScalarType, typename ScalarTypeForScaling>
    KOKKOS_INLINE_FUNCTION void shiftedScaledLegendreValues(OutputValueViewType outputValues, Intrepid2::ordinal_type n, ScalarType x, ScalarTypeForScaling t)
    {
      using OutputScalar = typename OutputValueViewType::value_type;
      OutputScalar two_x_minus_t = 2. * x - t;
      OutputScalar t_squared = t * t;
      if (n >= 0) outputValues(0) = 1.0;
      if (n >= 1) outputValues(1) = two_x_minus_t;
      for (int i=2; i<=n; i++)
      {
        const ScalarType one_over_i = 1.0 / ScalarType(i);
        outputValues(i) = one_over_i * ( (2. *i - 1. ) * two_x_minus_t * outputValues(i-1) - (i - 1.) * t_squared * outputValues(i-2));
      }
    }
    
    /** \brief  Integrated Legendre polynomials L_i for i>=1, defined for x in [0,1].
        \param [out] outputValues - the view into which to place the output values (must have at least n+1 entries)
        \param [in] n - the maximum polynomial order of integrated Legendre polynomials to compute
        \param [in] x - point at which to evaluate the polynomials
        \param [in] t - scaling parameter; may be of type double or may match the type of x
     These are defined for x in [0,1].  See equation (2.18) in Fuentes et al.  (We additionally define L_0 = 1.)
     
     Shifted, scaled Legendre polynomials are given by
        P_i(x;t) = P_i(x/t) * t^i = ~P_i(2x-t;t).
     
      The formula in Fuentes et al. is defined in terms of P_i and P_{i-2}.  We offer two versions of this computation, one which can
      reuse an existing P_i computation (in the form of a shiftedScaledLegendreValues input container), and one which reuses space in outputValues.
     */
    template<typename OutputValueViewType, typename ScalarType, typename ScalarTypeForScaling>
    KOKKOS_INLINE_FUNCTION void shiftedScaledIntegratedLegendreValues(OutputValueViewType outputValues, Intrepid2::ordinal_type n, ScalarType x, ScalarTypeForScaling t)
    {
      // reduced memory version: compute P_i in outputValues
      shiftedScaledLegendreValues(outputValues,n,x,t);
      // keep a copy of the last two P_i values around; update these before overwriting in outputValues
      ScalarType P_i_minus_two, P_i_minus_one;
      if (n >= 0) P_i_minus_two = outputValues(0);
      if (n >= 1) P_i_minus_one = outputValues(1);
      
      if (n >= 0) outputValues(0) = 1.0;
      if (n >= 1) outputValues(1) = x;
      for (int i=2; i<=n; i++)
      {
        const ScalarType & P_i = outputValues(i); // define as P_i just for clarity of the code below
        const ScalarType i_scalar = ScalarType(i);
        ScalarType L_i = (P_i - t * t * P_i_minus_two) /( 2. * (2. * i_scalar - 1.));
        
        // get the next values of P_{i-1} and P_{i-2} before overwriting the P_i value
        P_i_minus_two = P_i_minus_one;
        P_i_minus_one = P_i;
        
        // overwrite P_i value
        outputValues(i) = L_i;
      }
    }
    
    /** \brief  Integrated Legendre polynomials L_i for i>=1, defined for x in [0,1].
        \param [out] outputValues - the view into which to place the output values (must have at least n+1 entries)
        \param [in] shiftedScaledLegendreValues - view containing the shifted, scaled (unintegrated) Legendre values (must have at least n+1 entries)
        \param [in] n - the maximum polynomial order of integrated Legendre polynomials to compute
        \param [in] x - point at which to evaluate the polynomials
        \param [in] t - scaling parameter
     These are defined for x in [0,1].  See equation (2.18) in Fuentes et al.  (We additionally define L_0 = 1.)
     
     Shifted, scaled Legendre polynomials are given by
        P_i(x;t) = P_i(x/t) * t^i = ~P_i(2x-t;t).
     
      The formula in Fuentes et al. is defined in terms of P_i and P_{i-2}.  We offer two versions of this computation, one which can
      reuse an existing P_i computation (in the form of a shiftedScaledLegendreValues input container), and one which reuses space in outputValues.
     */
    template<typename OutputValueViewType, typename ScalarType, typename ScalarTypeForScaling>
    KOKKOS_INLINE_FUNCTION void shiftedScaledIntegratedLegendreValues(OutputValueViewType outputValues, const OutputValueViewType shiftedScaledLegendreValues,
                                                                      Intrepid2::ordinal_type n, ScalarType x, ScalarTypeForScaling t)
    {
      // reduced flops version: rely on previously computed P_i
      if (n >= 0) outputValues(0) = 1.0;
      if (n >= 1) outputValues(1) = x;
      for (int i=2; i<=n; i++)
      {
        const ScalarType & P_i           = shiftedScaledLegendreValues(i); // define as P_i just for clarity of the code below
        const ScalarType & P_i_minus_two = shiftedScaledLegendreValues(i-2);
        const ScalarType i_scalar = ScalarType(i);
        outputValues(i) = (P_i - t * t * P_i_minus_two) /( 2. * (2. * i_scalar - 1.));
      }
    }
    
    // the below implementation is commented out for now to guard a certain confusion.
    // the integratedLegendreValues() implementation below is implemented in such a way as to agree, modulo a coordinate transformation, with the
    // shiftedScaledIntegratedLegendreValues() above.  Since the latter really is an integral of shiftedScaledLegendreValues(), the former isn't
    // actually an integral of the (unshifted, unscaled) Legendre polynomials.
//    template<typename OutputValueViewType, typename ScalarType>
//    KOKKOS_INLINE_FUNCTION void integratedLegendreValues(OutputValueViewType outputValues, Intrepid2::ordinal_type n, ScalarType x)
//    {
//      // to reduce memory requirements, compute P_i in outputValues
//      legendreValues(outputValues,n,x);
//      // keep a copy of the last two P_i values around; update these before overwriting in outputValues
//      ScalarType P_i_minus_two, P_i_minus_one;
//      if (n >= 0) P_i_minus_two = outputValues(0);
//      if (n >= 1) P_i_minus_one = outputValues(1);
//
//      if (n >= 0) outputValues(0) = 1.0;
//      if (n >= 1) outputValues(1) = (x + 1.0) / 2.0;
//      for (int i=2; i<=n; i++)
//      {
//        const ScalarType & P_i = outputValues(i); // define as P_i just for clarity of the code below
//        const ScalarType i_scalar = ScalarType(i);
//        ScalarType L_i = (P_i - P_i_minus_two) /( 2. * (2. * i_scalar - 1.));
//
//        // get the next values of P_{i-1} and P_{i-2} before overwriting the P_i value
//        P_i_minus_two = P_i_minus_one;
//        P_i_minus_one = P_i;
//
//        // overwrite P_i value
//        outputValues(i) = L_i;
//      }
//    }
    
    /** \brief x derivative of shifted, scaled integrated Legendre polynomials L_i for i>=1, defined for x in [0,1].
        \param [out] outputValues - the view into which to place the output values (must have at least n+1 entries)
        \param [in] n - the maximum polynomial order of integrated Legendre polynomials to compute
        \param [in] x - point at which to evaluate the polynomials
        \param [in] t - scaling parameter
     These are defined for x in [0,1].  The x derivative of integrated Legendre is just Legendre; the only distinction is in the index -- outputValues indices are shifted by 1 relative to shiftedScaledLegendreValues, above.
     */
    template<typename OutputValueViewType, typename ScalarType, typename ScalarTypeForScaling>
    KOKKOS_INLINE_FUNCTION void shiftedScaledIntegratedLegendreValues_dx(OutputValueViewType outputValues, Intrepid2::ordinal_type n, ScalarType x, ScalarTypeForScaling t)
    {
      if (n >= 0) outputValues(0) = 0.0;
      if (n >= 1) outputValues(1) = 1.0;
      if (n >= 2) outputValues(2) = 2. * x - t;
      for (int i=2; i<=n-1; i++)
      {
        const ScalarType one_over_i = 1.0 / ScalarType(i);
        outputValues(i+1) = one_over_i * (2. * i - 1.) * (2. * x - t) * outputValues(i) -  one_over_i * (i - 1.0) * t * t * outputValues(i-1);
      }
    }
    
    /** \brief t derivative of shifted, scaled integrated Legendre polynomials L_i for i>=1, defined for x in [0,1].
        \param [out] outputValues - the view into which to place the output values (must have at least n+1 entries)
        \param [in] n - the maximum polynomial order of integrated Legendre polynomials to compute
        \param [in] x - point at which to evaluate the polynomials
        \param [in] t - scaling parameter
     These are defined for x in [0,1].  See Fuentes et al. (2.20)-(2.22).
     
     This implementation uses less memory than the one below, but depending on the application may introduce some extra computation, in the form of a call to shiftedScaledLegendreValues().
     */
    template<typename OutputValueViewType, typename ScalarType, typename ScalarTypeForScaling>
    KOKKOS_INLINE_FUNCTION void shiftedScaledIntegratedLegendreValues_dt(OutputValueViewType outputValues, Intrepid2::ordinal_type n, ScalarType x, ScalarTypeForScaling t)
    {
      // memory-conserving version -- place the Legendre values in the final output container
      shiftedScaledLegendreValues(outputValues, n, x, t);
      
      ScalarType P_i_minus_2 = outputValues(0);
      ScalarType P_i_minus_1 = outputValues(1);
      
      if (n >= 0) outputValues(0) = 0.0;
      if (n >= 1) outputValues(1) = 0.0;
      
      for (int i=2; i<=n; i++)
      {
        const ScalarType L_i_dt = -0.5 * (P_i_minus_1 + t * P_i_minus_2);
        
        P_i_minus_2 = P_i_minus_1;
        P_i_minus_1 = outputValues(i);
        
        outputValues(i) = L_i_dt;
      }
    }
    
    /** \brief t derivative of shifted, scaled integrated Legendre polynomials L_i for i>=1, defined for x in [0,1].
        \param [out] outputValues - the view into which to place the output values (must have at least n+1 entries)
        \param [in] shiftedScaledLegendreValues - previously computed shifted, scaled Legendre values; must have at least n entries
        \param [in] n - the maximum polynomial order of integrated Legendre polynomials to compute
        \param [in] x - point at which to evaluate the polynomials
        \param [in] t - scaling parameter
     These are defined for x in [0,1].  See Fuentes et al. (2.20)-(2.22).
     
     This implementation uses more memory than the one above, but depending on the application may save some computation, in that it can reuse previously computed shiftedScaledLegendreValues.
     */
    template<typename OutputValueViewType, typename ScalarType, typename ScalarTypeForScaling>
    KOKKOS_INLINE_FUNCTION void shiftedScaledIntegratedLegendreValues_dt(OutputValueViewType outputValues, const OutputValueViewType shiftedScaledLegendreValues,
                                                                         Intrepid2::ordinal_type n, ScalarType x, ScalarTypeForScaling t)
    {
      // reduced flops version: rely on previously computed P_i
      if (n >= 0) outputValues(0) = 0.0;
      if (n >= 1) outputValues(1) = 0.0;
      for (int i=2; i<=n; i++)
      {
        const ScalarType & P_i_minus_1 = shiftedScaledLegendreValues(i-1); // define as P_i just for clarity of the code below
        const ScalarType & P_i_minus_2 = shiftedScaledLegendreValues(i-2);
        outputValues(i) = -0.5 * (P_i_minus_1 + t * P_i_minus_2);
      }
    }
    
    /** \brief Shifted, scaled Jacobi values, defined for x in [0,1].
        \param [out] outputValues - the view into which to place the output values (must have at least n+1 entries)
        \param [in] alpha - Jacobi alpha parameter (beta is taken to be 0)
        \param [in] n - the maximum polynomial order of Jacobi polynomials to compute
        \param [in] x - point at which to evaluate the polynomials
        \param [in] t - scaling parameter
     These are defined for x in [0,1].  See Fuentes et al. (2.24).
     
     In general, Jacobi polynomials have two parameters, alpha and beta.
     
     Following Fuentes et al., we only consider beta=0, and use a domain for x of [0,1] (compared to the classical [-1,1]).
     
     When alpha = 0, Jacobi coincides with Legendre.
    */
    template<typename OutputValueViewType, typename ScalarType, typename ScalarTypeForScaling>
    KOKKOS_INLINE_FUNCTION void shiftedScaledJacobiValues(OutputValueViewType outputValues, double alpha, Intrepid2::ordinal_type n, ScalarType x, ScalarTypeForScaling t)
    {
      ScalarType two_x_minus_t = 2. * x - t;
      ScalarTypeForScaling alpha_squared_t = alpha * alpha * t;
      
      if (n >= 0) outputValues(0) = 1.0;
      if (n >= 1) outputValues(1) = two_x_minus_t + alpha * x;
      
      for (int i=2; i<=n; i++)
      {
        const ScalarType & P_i_minus_one = outputValues(i-1); // define as P_i just for clarity of the code below
        const ScalarType & P_i_minus_two = outputValues(i-2);
        
        double a_i = (2. * i) * (i + alpha) * (2. * i + alpha - 2.);
        double b_i = 2. * i + alpha - 1.;
        double c_i = (2. * i + alpha) * (2. * i + alpha - 2.);
        double d_i = 2. * (i + alpha - 1.) * (i - 1.) * (2. * i + alpha);
        
        outputValues(i) = (b_i / a_i) * (c_i * two_x_minus_t + alpha_squared_t) * P_i_minus_one - (d_i / a_i) * t * t * P_i_minus_two;
      }
    }
    
    /** \brief Integrated Jacobi values, defined for x in [0,1].
     \param [out] outputValues - the view into which to place the output values (must have at least n+1 entries)
     \param [in] jacobiValues - previously computed Jacobi values; must have at least n+1 entries
     \param [in] alpha - Jacobi alpha parameter (beta is taken to be 0)
     \param [in] n - the maximum polynomial order of Jacobi polynomials to compute
     \param [in] x - point at which to evaluate the polynomials
     \param [in] t - scaling parameter
     These are defined for x in [0,1].  See Fuentes et al. (2.27).
     
     In general, Jacobi polynomials have two parameters, alpha and beta.
     
     Following Fuentes et al., we only consider beta=0, and use a domain for x of [0,1] (compared to the classical [-1,1]).
     
     When alpha = 0, Jacobi coincides with Legendre.
     
     Compared with the shiftedScaledIntegratedJacobiValues() below, this version uses more memory, but may require fewer floating point computations by reusing the values in jacobiValues.
     */
    template<typename OutputValueViewType, typename ScalarType, typename ScalarTypeForScaling>
    KOKKOS_INLINE_FUNCTION void shiftedScaledIntegratedJacobiValues(OutputValueViewType outputValues, const OutputValueViewType jacobiValues,
                                                                    double alpha, Intrepid2::ordinal_type n, ScalarType x, ScalarTypeForScaling t)
    {
      // reduced flops version: rely on previously computed P_i
      if (n >= 0) outputValues(0) = 1.0;
      if (n >= 1) outputValues(1) = x;
      
      ScalarType t_squared = t * t;
      for (int i=2; i<=n; i++)
      {
        const ScalarType & P_i         = jacobiValues(i);   // define as P_i just for clarity of the code below
        const ScalarType & P_i_minus_1 = jacobiValues(i-1);
        const ScalarType & P_i_minus_2 = jacobiValues(i-2);
        
        double a_i = (i + alpha) / ((2. * i + alpha - 1.) * (2. * i + alpha     ));
        double b_i = alpha       / ((2. * i + alpha - 2.) * (2. * i + alpha     ));
        double c_i = (i - 1.)    / ((2. * i + alpha - 2.) * (2. * i + alpha - 1.));
        
        outputValues(i) = a_i * P_i + b_i * t * P_i_minus_1 - c_i * t_squared * P_i_minus_2;
      }
    }
    
    /** \brief Integrated Jacobi values, defined for x in [0,1].
     \param [out] outputValues - the view into which to place the output values (must have at least n+1 entries)
     \param [in] jacobiValues - previously computed Jacobi values; must have at least n+1 entries
     \param [in] alpha - Jacobi alpha parameter (beta is taken to be 0)
     \param [in] n - the maximum polynomial order of Jacobi polynomials to compute
     \param [in] x - point at which to evaluate the polynomials
     \param [in] t - scaling parameter
     These are defined for x in [0,1].  See Fuentes et al. (2.27).
     
     In general, Jacobi polynomials have two parameters, alpha and beta.
     
     Following Fuentes et al., we only consider beta=0, and use a domain for x of [0,1] (compared to the classical [-1,1]).
     
     When alpha = 0, Jacobi coincides with Legendre.
     
     Compared with the shiftedScaledIntegratedJacobiValues() above, this version uses less memory, but may require more floating point computations.
     */
    template<typename OutputValueViewType, typename ScalarType, typename ScalarTypeForScaling>
    KOKKOS_INLINE_FUNCTION void shiftedScaledIntegratedJacobiValues(OutputValueViewType outputValues,
                                                                    double alpha, Intrepid2::ordinal_type n, ScalarType x, ScalarTypeForScaling t)
    {
      // memory-conserving version -- place the Jacobi values in the final output container
      shiftedScaledJacobiValues(outputValues, alpha, n, x, t);
      
      ScalarType P_i_minus_2 = outputValues(0);
      ScalarType P_i_minus_1 = outputValues(1);
      
      if (n >= 0) outputValues(0) = 1.0;
      if (n >= 1) outputValues(1) = x;
      
      ScalarType t_squared = t * t;
      for (int i=2; i<=n; i++)
      {
        const ScalarType & P_i = outputValues(i);
        
        double a_i = (i + alpha) / ((2. * i + alpha - 1.) * (2. * i + alpha     ));
        double b_i = alpha       / ((2. * i + alpha - 2.) * (2. * i + alpha     ));
        double c_i = (i - 1.)    / ((2. * i + alpha - 2.) * (2. * i + alpha - 1.));
        
        ScalarType L_i = a_i * P_i + b_i * t * P_i_minus_1 - c_i * t_squared * P_i_minus_2;
        
        P_i_minus_2 = P_i_minus_1;
        P_i_minus_1 = P_i;
        
        outputValues(i) = L_i;
      }
    }
    
    // x derivative of integrated Jacobi is just Jacobi
    // only distinction is in the index -- outputValues indices are shifted by 1 relative to jacobiValues, above
    /** \brief x derivative of integrated Jacobi polynomials L_i for i>=1, defined for x in [0,1].
        \param [out] outputValues - the view into which to place the output values (must have at least n+1 entries)
        \param [in] n - the maximum polynomial order of integrated Legendre polynomials to compute
        \param [in] x - point at which to evaluate the polynomials
        \param [in] t - scaling parameter
     These are defined for x in [0,1].  The x derivative of integrated Jacobi is just Jacobi; the only distinction is in the index -- outputValues indices are shifted by 1 relative to shiftedScaledJacobiValues, above.
     */
    template<typename OutputValueViewType, typename ScalarType, typename ScalarTypeForScaling>
    KOKKOS_INLINE_FUNCTION void shiftedScaledIntegratedJacobiValues_dx(OutputValueViewType outputValues,
                                                                       double alpha, Intrepid2::ordinal_type n, ScalarType x, ScalarTypeForScaling t)
    {
      // rather than repeating the somewhat involved implementation of jacobiValues here,
      // call with (n-1), and then move values accordingly
      shiftedScaledJacobiValues(outputValues, alpha, n-1, x, t);
      
      // forward implementation
      ScalarType nextValue = 0.0;
      ScalarType nextNextValue = 0.0;
      for (int i=0; i<=n-1; i++)
      {
        nextNextValue = outputValues(i);
        outputValues(i) = nextValue;
        nextValue = nextNextValue;
      }
      outputValues(n-1) = nextValue;
    }
    
     /** \brief t derivative of shifted, scaled integrated Jacobi polynomials L_i for i>=1, defined for x in [0,1].
        \param [out] outputValues - the view into which to place the output values (must have at least n+1 entries)
        \param [in] jacobiValues - previously computed shifted, scaled Jacobi values; must have at least n entries
        \param [in] n - the maximum polynomial order of integrated Legendre polynomials to compute
        \param [in] x - point at which to evaluate the polynomials
        \param [in] t - scaling parameter
     These are defined for x in [0,1].  See Fuentes et al. (2.30).
     
     This implementation uses more memory than the one above, but depending on the application may save some computation, in that it can reuse previously computed jacobiValues.
     */
    template<typename OutputValueViewType, typename ScalarType, typename ScalarTypeForScaling>
    KOKKOS_INLINE_FUNCTION void shiftedScaledIntegratedJacobiValues_dt(OutputValueViewType outputValues, const OutputValueViewType jacobiValues,
                                                                       double alpha, Intrepid2::ordinal_type n, ScalarType x, ScalarTypeForScaling t)
    {
      // reduced flops version: rely on previously computed P_i
      if (n >= 0) outputValues(0) = 0.0;
      if (n >= 1) outputValues(1) = 0.0;
      for (int i=2; i<=n; i++)
      {
        const ScalarType & P_i_minus_1 = jacobiValues(i-1); // define as P_i just for clarity of the code below
        const ScalarType & P_i_minus_2 = jacobiValues(i-2);
        outputValues(i) = - (i-1.) / (2. * i - 2. + alpha) * (P_i_minus_1 + t * P_i_minus_2);
      }
    }
    
     /** \brief t derivative of shifted, scaled integrated Jacobi polynomials L_i for i>=1, defined for x in [0,1].
        \param [out] outputValues - the view into which to place the output values (must have at least n+1 entries)
        \param [in] jacobiValues - previously computed shifted, scaled Jacobi values; must have at least n entries
        \param [in] n - the maximum polynomial order of integrated Legendre polynomials to compute
        \param [in] x - point at which to evaluate the polynomials
        \param [in] t - scaling parameter
     These are defined for x in [0,1].  See Fuentes et al. (2.30).
     
     This implementation requires less memory than the one above, but depending on the application may require some extra computation.
     */
    template<typename OutputValueViewType, typename ScalarType, typename ScalarTypeForScaling>
    KOKKOS_INLINE_FUNCTION void shiftedScaledIntegratedJacobiValues_dt(OutputValueViewType outputValues,
                                                                       double alpha, Intrepid2::ordinal_type n, ScalarType x, ScalarTypeForScaling t)
    {
      // memory-conserving version -- place the Jacobi values in the final output container
      shiftedScaledJacobiValues(outputValues, alpha, n, x, t);
      
      ScalarType P_i_minus_2 = outputValues(0);
      ScalarType P_i_minus_1 = outputValues(1);
      
      if (n >= 0) outputValues(0) = 0.0;
      if (n >= 1) outputValues(1) = 0.0;
      
      for (int i=2; i<=n; i++)
      {
        const ScalarType L_i_dt =  - (i-1.) / (2. * i - 2. + alpha) * (P_i_minus_1 + t * P_i_minus_2);
        
        P_i_minus_2 = P_i_minus_1;
        P_i_minus_1 = outputValues(i);
        
        outputValues(i) = L_i_dt;
      }
    }
  } // namespace Polynomials
} // namespace Intrepid2

#endif /* Intrepid2_Polynomials_h */

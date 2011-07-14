// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov) or
//                    Robert Kirby (robert.c.kirby@ttu.edu)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_HGRAD_TET_Cn_FEM.hpp
    \brief  Header file for the Intrepid::HGRAD_TET_Cn_FEM_ORTH class.
    \author Created by Robert C. Kirby
*/

#ifndef INTREPID_HGRAD_TET_Cn_FEM_ORTH_HPP
#define INTREPID_HGRAD_TET_Cn_FEM_ORTH_HPP

#include "Intrepid_Basis.hpp"
#include "Sacado.hpp"

namespace Intrepid {
  
/** \class  Intrepid::Basis_HGRAD_TET_Cn_FEM_ORTH
    \brief  Implementation of the default H(grad)-compatible orthogonal basis of
            arbitrary degree on tetrahedron.
  
    \remarks

    \li   All degrees of freedom are considered to be internal (ie not assembled)
  */
  
template<class Scalar, class ArrayScalar> 
class Basis_HGRAD_TET_Cn_FEM_ORTH: public Basis<Scalar, ArrayScalar> {
private:
  /** \brief  Initializes <var>tagToOrdinal_</var> and
      <var>ordinalToTag_</var> lookup arrays.   */
  void initializeTags();
  
public:
  
  /** \brief  Constructor.
      \param  degree            [in] - the degree of polynomials contained in the basis.
   */
  Basis_HGRAD_TET_Cn_FEM_ORTH( int degree );  
  
  
  /** \brief  Evaluation of a FEM basis on a <strong>reference Tetrahedron</strong> cell. 
  
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Tetrahedron</strong> cell. For rank and dimensions of
              I/O array arguments see Section \ref basis_md_array_sec .

      \param  outputValues      [out] - variable rank array with the basis values
      \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points
      \param  operatorType      [in]  - the operator acting on the basis functions    
  */
  void getValues(ArrayScalar &          outputValues,
                 const ArrayScalar &    inputPoints,
                 const EOperator        operatorType) const;
  
  
  /**  \brief  FVD basis evaluation: invocation of this method throws an exception.
   */
  void getValues(ArrayScalar &          outputValues,
                 const ArrayScalar &    inputPoints,
                 const ArrayScalar &    cellVertices,
                 const EOperator        operatorType = OPERATOR_VALUE) const;

};

  /** \class TabulatorTet
      \brief This is an internal class with a static member function for
      tabulating derivatives of orthogonal expansion functions.  

      It is a separate class to allow recursive templates partially specialized on
      derivative order without throwing the HGRAD_TET_Cn_FEM_ORTH class into
      an infinite compiler loop.

      This class is intended only to be used internally by the HGRAD_TET_Cn_FEM_ORTH basis
      to implement all the derivative orders in the Basis interface, hiding recursion
      and calls to Sacado. 
  */

template<typename Scalar,typename ArrayScalar, unsigned derivOrder>
class TabulatorTet
{
public:
  /** \brief basic tabulate mathod evaluates the derivOrder^th derivatives
      of the basis functions at inputPoints into outputValues.
      
      \param      [out] outputValues - rank 2 (if derivOrder == 0) or rank 3
                                      array holding the result
      \param      [in]  deg - the degree up to which to tabulate the bases
      \param      [in] inputPoints - a rank 2 array containing the
                                    points at which to evaluate the basis functions.
  */
  static void tabulate( ArrayScalar & outputValues ,
                        const int deg ,
                        const ArrayScalar &inputPoints );
};
  

  /** \class TabulatorTet<Scalar,ArrayScalar,0>
      \brief This is specialized on 0th derivatives to make the
      tabulate function run through recurrence relations.
  */
template<typename Scalar,typename ArrayScalar>
class TabulatorTet<Scalar,ArrayScalar,0>
{
public:
  /** \brief basic tabulate mathod evaluates the basis functions at inputPoints into outputValues. 
      
      \param      [out] outputValues - rank 2 array (F,P) holding
                                      the basis functions at points.
      \param      [in]  deg - the degree up to which to tabulate the bases
      \param      [in]  inputPoints - a rank 2 array containing the
                                      points at which to evaluate the basis functions.
  */
  static void tabulate( ArrayScalar & outputValues ,
                        const int deg ,
                        const ArrayScalar &inputPoints );

  /** \brief function for indexing from orthogonal expansion indices into linear space
      p+q+r = the degree of the polynomial.
      \param p [in] - the first index
      \param q [in] - the second index
      \param r [in] - the third index */
  static int idx(int p, int q,int r)
  {
    return (p+q+r)*(p+q+r+1)*(p+q+r+2)/6+(q+r)*(q+r+1)/2+r;
  }

  /** \brief function for computing the Jacobi recurrence coefficients so that
      \param alpha [in] - the first Jacobi weight
      \param beta  [in] - the second Jacobi weight
      \param n     [n]  - the polynomial degree
      \param an    [out] - the a weight for recurrence
      \param bn    [out] - the b weight for recurrence
      \param cn    [out] - the c weight for recurrence

      The recurrence is
      \f[
      P^{\alpha,\beta}_{n+1} = \left( a_n + b_n x\right) P^{\alpha,\beta}_n - c_n P^{\alpha,\beta}_{n-1}
      \f],
      where
      \f[
      P^{\alpha,\beta}_0 = 1
      \f]
  */
  static void jrc( const Scalar &alpha , const Scalar &beta , 
		   const int &n ,
		   Scalar &an , Scalar &bn, Scalar &cn )
  
  {
    an = (2.0 * n + 1.0 + alpha + beta) * ( 2.0 * n + 2.0 + alpha + beta ) 
      / ( 2.0 * ( n + 1 ) * ( n + 1 + alpha + beta ) );
    bn = (alpha*alpha-beta*beta)*(2.0*n+1.0+alpha+beta) 
      / ( 2.0*(n+1.0)*(2.0*n+alpha+beta)*(n+1.0+alpha+beta) );
    cn = (n+alpha)*(n+beta)*(2.0*n+2.0+alpha+beta) 
      / ( (n+1.0)*(n+1.0+alpha+beta)*(2.0*n+alpha+beta) );
  
  return;
  }
  

    
};

  /** \class TabulatorTet<Scalar,ArrayScalar,1>
      \brief This is specialized on 1st derivatives 
      since it recursively calls the 0th derivative class
      with Sacado AD types, and so the outputValues it passes
      to that function needs to have a rank 2 rather than rank 3
  */
template<typename Scalar,typename ArrayScalar>
class TabulatorTet<Scalar,ArrayScalar,1>
{
public:
  /** \brief basic tabulate mathod evaluates the first derivatives
      of the basis functions at inputPoints into outputValues.
      
      \param      [out] outputValues - rank 3 array (F,P,D) holding
                                       the derivatives of basis functions at points.
      \param      [in]  deg - the degree up to which to tabulate the bases
      \param      [in]  inputPoints - a rank 3 array containing the
                                      points at which to evaluate the basis functions.
  */
  static void tabulate( ArrayScalar & outputValues ,
                        const int deg ,
                        const ArrayScalar &inputPoints );
    
};



}// namespace Intrepid

#include "Intrepid_HGRAD_TET_Cn_FEM_ORTHDef.hpp"

#endif


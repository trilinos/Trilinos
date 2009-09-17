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
//                    Robert C. Kirby (robert.c.kirby@ttu.edu)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_OrthogonalBases.hpp
\brief  Header file for orthogonal bases on various cell types
\author Created by R. Kirby
*/

#ifndef INTREPID_ORTHOGONALBASES_HPP
#define INTREPID_ORTHGONALBASES_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Types.hpp"
#include "Intrepid_Utils.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_TestForException.hpp"

namespace Intrepid {
  
  /** \class Intrepid::OrthgonalBases
  \brief Basic implementation of general orthogonal polynomials on a
  range of shapes, including the triangle, and tetrahedron.

  Each basis is templated over Scalar type to allow multiple
  precisions and automatic differentiation.

  All methods of this class are static

  The recurrence relations are formulated in such a way that automatic
  differentiation on collapsed-coordinate bases works at all points in
  the domain.

  Also provided are routines for obtaining nodal and modal
  derivative matrices for each basis.
  */

  class OrthogonalBases {
  public:
    OrthogonalBases() {;}
    ~OrthogonalBases() {;}

    /** \brief Calculates triangular orthogonal expansions
        (e.g. Dubiner basis) at a range of input points 
        
        \param np       [in]    - number of input points
        \param z        [in]    - 2d array of points z(pt,2)
        \param n        [in]    - the maximum polynomial degree tabulated
        \param poly_val [out]   - 2d array poly_val((n+1)(n+2)/2,np)

      \li The ScalarArray types must support (i,j) indexing 
      and a dimension(i) operation.
    */

    template<class Scalar, class ScalarArray1, class ScalarArray2>
    static void tabulateTriangle( const ScalarArray1& z ,
                                  const int n ,
                                  ScalarArray2 & poly_val );

    /** \brief Calculates triangular orthogonal expansions
        (e.g. Dubiner basis) at a range of input points 
        
        \param np       [in]    - number of input points
        \param z        [in]    - 2d array of points z(pt,3)
        \param n        [in]    - the maximum polynomial degree tabulated
        \param poly_val [out]   - 2d array poly_val((n+1)(n+2)(n+3)/6,np)

      \li The ScalarArray types must support (i,j) indexing 
      and a dimension(i) operation.
    */

    template<class Scalar, class ScalarArray1, class ScalarArray2>
    static void tabulateTetrahedron( const ScalarArray1& z ,
                                    const int n ,
                                    ScalarArray2 & poly_val );
    
  private:
    /** \brief computes Jacobi recurrence coefficients of
        order n with weights a,b so that
        P^{alpha,beta}_{n+1}(x) 
        = (an x + bn) P^{alpha,beta}_n(x) - cn P^{alpha,beta}_{n-1}(x)
    */

    template<class Scalar>
    static void jrc( const Scalar &alpha , const Scalar &beta , const int &n ,
                    Scalar &an , Scalar &bn, Scalar &cn );

    /** \brief Given indices p,q, computes the linear index of
        the Dubiner polynomial D^{p,q} */
    static inline int idxtri(int p, int q)
    {
      return (p+q)*(p+q+1)/2+q;
    }

    /** \brief Given indices p,q,r, computes the linear index of the
        tetrahedral polynomial D^{p,q,r} */
    static inline int idxtet(int p, int q, int r)
    {
      return (p+q+r)*(p+q+r+1)*(p+q+r+2)/6+(q+r)*(q+r+1)/2+r;
    }


  }; // class OrthogonalBases
} // namespace Intrepid

#include "Intrepid_OrthogonalBasesDef.hpp"

#endif

    





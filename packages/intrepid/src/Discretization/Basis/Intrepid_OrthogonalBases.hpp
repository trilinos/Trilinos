// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
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

    





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

/** \file   Intrepid_HGRAD_LINE_Hermite_FEM.hpp
    \brief  Header file for the Intrepid::HGRAD_LINE_Hermite_FEM class.
    \author Created by G. von Winckel
*/

#ifndef INTREPID_HGRAD_LINE_HERMITE_FEM_HPP
#define INTREPID_HGRAD_LINE_HERMITE_FEM_HPP

#include "Intrepid_Basis.hpp"
#include "Intrepid_PointTools.hpp"
#include "Intrepid_Polylib.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include <ostream>
#include <cstdint>

namespace Intrepid {

/** \class  Intrepid::Basis_HGRAD_LINE_Hermite_FEM
    \brief  Implements Hermite interpolant basis of degree n on the reference Line cell. 
            The basis has cardinality 2n and spans a COMPLETE linear polynomial space. 

            The Hermite interpolant approximation space is \f$\mathcal{H}_n\subset H_0^2(\Omega)\f$
            where \f$\Omega=[-1,1]\f$

            \f[\mathcal{H}_n = \{h_{2j}(x),h_{2j+1}\}_{j=0}^{n-1} \f]

            The basis functions determined by the interpolation points \f$\{x_i\}_{i=0}^{n-1}\f$
            which satisfy the ordering property

            \f[ -1\leq x_0 \leq \cdots \leq x_i \leq x_{i+1} \leq \cdots \leq x_n \leq 1 \f]

            while the basis functions satisfy

            \f[ h_{2j}(x_i) = \delta_{ij},\quad h'_{2j}(x_i) = 0 \f]
            \f[ h_{2j+1}(x_i) = 0,\quad h'_{2j+1}(x_i) = \delta_{ij} \f]

            Basis functions and their derivatives are evaluated from the Legendre polynomials
            and their derivatives evaluated on the interpolation points and the input points in
            the getValues() method. Legendre polynomials and their derivatives are computed
            using the relations

            \f[ P_{j+1}(x) = x P_j(x) + \frac{x^2-1}{j+1}P_j^{(1)}(x),\quad P_0(x) = 1 \f]
            \f[ P_{j+1}^{(m)(x)} = (j+m)P_j^{(m-1)}(x) + x P_j^{(m)}(x),\quad P_j^{(m)} = 0,\;\forall j<m\f]


*/
template<class Scalar, class ArrayScalar>   
class Basis_HGRAD_LINE_Hermite_FEM: public Basis<Scalar, ArrayScalar> {

  template<typename T> using RCP = Teuchos::RCP<T>;
  using SerialDenseMatrix = Teuchos::SerialDenseMatrix<int,Scalar>;

  using TAGS = std::vector<int>;

private:

  /** \brief Holds the points defining the Hermite basis */
  FieldContainer<Scalar> latticePts_;  

  /** \brief Contains the values of the Legendre polynomials and their derivatives */
  mutable SerialDenseMatrix V_;

  mutable Teuchos::SerialDenseSolver<int,Scalar> solver_;  

  mutable TAGS tags_;

  mutable bool isFactored_;

  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
  */
  void initializeTags();

  /** \brief Evaluates \f$P_j(x_i)\f$ and \f$P_j'(x_i)\f$ at a particular point \f$x_i\f$
  */
  void recurrence( ArrayScalar &P, ArrayScalar &Px, const Scalar x ) const; 

  /** \brief Evaluates \f$P_j^{(m)}(x_i) \f$ and \f$P_j^{(m+1)}(x_i)  \f$ at a 
    *        particular point \f$x_i\f$ 
    */
  void legendre_d( ArrayScalar &Pm, ArrayScalar &Pm1, const int m, const Scalar pt ) const;

  /** \brief Form the Legendre/Derivative Vandermonde matrix at the given lattice points and
             have the linear solver factor with equilibration */                        
  void setupVandermonde( bool factor=true );
 
   
public:

  /** \brief  Default Constructor assumes the two interpolation points are the cell vertices. 
   *          Cubic Hermite Interpolation
  */
  Basis_HGRAD_LINE_Hermite_FEM();  

   /** \brief Constructor.
      \param  int pointType:    [in] array containing interpolation points */
  Basis_HGRAD_LINE_Hermite_FEM( const ArrayScalar &pts );
  
  /** \brief Constructor.
      \param  int numPoints:    [in] number of interpolation points 
      \param  int pointType:    [in] type of points, either POINTTYPE_EQUISPACED or POINTTYPE_SPECTRAL */
  Basis_HGRAD_LINE_Hermite_FEM( int numPoints , const EPointType &pointType ); 


  /** \brief  Evaluation of a FEM basis on a <strong>reference Line</strong> cell. 
  
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Line</strong> cell. For rank and dimensions of
              I/O array arguments see Section \ref basis_md_array_sec .

      \param  outputValues      [out] - variable rank array with the basis values
      \param  inputPoints       [in]  - rank-n array (P,D) with the evaluation points
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

  /** \brief implements the dofcoords interface */
  virtual void getDofCoords( ArrayScalar & DofCoords ) const;


  void printTags( std::ostream &os ); 

}; // class Basis_HGRAD_LINE_Hermite_FEM

}// namespace Intrepid

#include "Intrepid_HGRAD_LINE_Hermite_FEMDef.hpp"

#endif // INTREPID_HGRAD_LINE_HERMITE_FEM_HPP

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif


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
//                    Denis Ridzal (dridzal@sandia.gov)
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_HDIV_QUAD_In_FEM.hpp
    \brief  Header file for the Intrepid::HDIV_QUAD_In_FEM class.
    \author Created by R. Kirby and P. Bochev and D. Ridzal and K. Petrson.
 */

#ifndef INTREPID_HDIV_QUAD_In_FEM_HPP
#define INTREPID_HDIV_QUAD_In_FEM_HPP
#include "Intrepid_Basis.hpp"
#include "Intrepid_ProductTopology.hpp"
#include "Intrepid_HGRAD_LINE_Cn_FEM.hpp"

namespace Intrepid {
  
/** \class  Intrepid::Basis_HDIV_QUAD_In_FEM
    \brief  Implementation of the default H(div)-compatible FEM basis of degree 1 on Quadrilateral cell 
  
            Implements Raviart-Thomas basis of degree n on the reference Quadrilateral cell. The basis has
            cardinality 2(n+1)n and spans a INCOMPLETE polynomial space. 
  \endverbatim
  
 */
  
template<class Scalar, class ArrayScalar> 
class Basis_HDIV_QUAD_In_FEM : 
    public Basis<Scalar, ArrayScalar> ,
    public DofCoordsInterface<ArrayScalar> {
private:
  
  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  void initializeTags();
  
  Basis_HGRAD_LINE_Cn_FEM<Scalar,ArrayScalar> closedBasis_;
  Basis_HGRAD_LINE_Cn_FEM<Scalar,ArrayScalar> openBasis_;

  FieldContainer<double> closedPts_;
  FieldContainer<double> openPts_;


public:
  /** \brief Destructor
   */
  virtual ~Basis_HDIV_QUAD_In_FEM() {;}

  /** \brief  Constructor.
      \param order     [in] - order of polynomial space
      \param ptsClosed [in] - pts that include the endpoints, used in the direction of normal continuity
      \param ptsOpen   [in] - used in "off" direction
    */
  Basis_HDIV_QUAD_In_FEM( int order , 
			  const ArrayScalar &ptsClosed ,
			  const ArrayScalar &ptsOpen );

  /** \brief  Streamlined constructor that allows user to request equispaced points or
      Gauss-Lobatto cross with Gauss-Legendre points in each vector component
    */
  Basis_HDIV_QUAD_In_FEM( int order , const EPointType &pointType );

  
    
  /** \brief  Evaluation of a FEM basis on a <strong>reference Quadrilateral</strong> cell. 
    
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Quadrilateral</strong> cell. For rank and dimensions of
              I/O array arguments see Section \ref basis_md_array_sec.
  
      \param  outputValues      [out] - rank-3 or 4 array with the computed basis values
      \param  inputPoints       [in]  - rank-2 array with dimensions (P,D) containing reference points  
      \param  operatorType      [in]  - operator applied to basis functions    
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

  /** \brief  Returns spatial locations (coordinates) of degrees of freedom on a
              <strong>reference cell</strong>; defined for interpolatory bases.

      \param  DofCoords      [out] - array with the coordinates of degrees of freedom,
                                     dimensioned (F,D)
   */
   virtual void getDofCoords(ArrayScalar & DofCoords) const;
  
};
}// namespace Intrepid

#include "Intrepid_HDIV_QUAD_In_FEMDef.hpp"

#endif

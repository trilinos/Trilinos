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

/** \file   Intrepid_HGRAD_QUAD_Cn_FEM.hpp
    \brief  Header file for the Intrepid::HGRAD_QUAD_Cn_FEM class.
    \author Created by R. Kirby.
 */

#ifndef INTREPID_HGRAD_QUAD_Cn_FEM_HPP
#define INTREPID_HGRAD_QUAD_Cn_FEM_HPP
#include "Intrepid_Basis.hpp"
#include "Intrepid_ProductTopology.hpp"
#include "Intrepid_HGRAD_LINE_Cn_FEM.hpp"
#include "Teuchos_Array.hpp"
#include "Intrepid_TensorBasis.hpp"
#include "Intrepid_Utils.hpp"

namespace Intrepid {
  
/** \class  Intrepid::Basis_HGRAD_QUAD_C1_FEM
    \brief  Implementation of the default H(grad)-compatible FEM basis of degree 1 on Quadrilateral cell
            Implements Lagrangian basis of degree n on the reference Quadrilateral cell using
	    a tensor product of points
*/
template<class Scalar, class ArrayScalar> 
class Basis_HGRAD_QUAD_Cn_FEM : public TensorBasis<Scalar,ArrayScalar>
{
private:
  
  /** \brief Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  void initializeTags();
  
public:

  /** \brief Constructor.
  */
  Basis_HGRAD_QUAD_Cn_FEM( const int orderx , const int ordery,
			   const ArrayScalar &pts_x ,
			   const ArrayScalar &pts_y );

  /** \brief Constructor giving a PointType switch for equispaced or spectral (i.e. Gauss-Lobatto) points */
  Basis_HGRAD_QUAD_Cn_FEM( const int order ,
			   const EPointType &pointType );
  
    
  /** \brief  FEM basis evaluation on a <strong>reference Quadrilateral</strong> cell. 
    
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Quadrilateral</strong> cell. For rank and dimensions of 
              I/O array arguments see Section \ref basis_md_array_sec .
  
      \param  outputValues      [out] - rank-2 or 3 array with the computed basis values
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
};
}// namespace Intrepid

#include "Intrepid_HGRAD_QUAD_Cn_FEMDef.hpp"

#endif

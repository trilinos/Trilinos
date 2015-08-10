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

/** \file   Intrepid_HGRAD_HEX_Cn_FEM.hpp
    \brief  Header file for the Intrepid::HGRAD_HEX_Cn_FEM class.
    \author Created by P. Bochev and D. Ridzal.
 */

#ifndef INTREPID_HGRAD_HEX_Cn_FEM_HPP
#define INTREPID_HGRAD_HEX_Cn_FEM_HPP
#include "Intrepid_Basis.hpp"
#include "Intrepid_ProductTopology.hpp"
#include "Intrepid_HGRAD_LINE_Cn_FEM.hpp"
#include "Teuchos_Array.hpp"
#include "Intrepid_TensorBasis.hpp"

namespace Intrepid {
  
/** \class  Intrepid::Basis_HGRAD_HEX_Cn_FEM
    \brief  Implementation of the default H(grad)-compatible FEM basis of degree 2 on Hexahedron cell 
  
            Implements Lagrangian basis of degree n on the reference Hexahedron cell. The basis has
            cardinality (n+1)^3 and spans a COMPLETE polynomial space. Basis functions are dual 
            to a unisolvent set of degrees-of-freedom (DoF) defined lexicographically on an
	    array of input points.
  
  \endverbatim
  
 */
  
template<class Scalar, class ArrayScalar> 
class Basis_HGRAD_HEX_Cn_FEM : 
    public TensorBasis<Scalar, ArrayScalar> ,
    public DofCoordsInterface<ArrayScalar>
{
private:
  FieldContainer<double> ptsx_;
  FieldContainer<double> ptsy_;
  FieldContainer<double> ptsz_;
  
  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  void initializeTags();

  
public:
  /** \brief Destructor.
   */
  virtual ~Basis_HGRAD_HEX_Cn_FEM( ) {;}

  /** \brief  Constructor.
    */
  Basis_HGRAD_HEX_Cn_FEM( const int orderx , const int ordery, const int orderz ,
			  const ArrayScalar &pts_x ,
			  const ArrayScalar &pts_y ,
			  const ArrayScalar &pts_z );

  /** \brief streamlined constructor giving two choices: equispaced
      lattice points or spectral (i.e. Lobatto) */
  Basis_HGRAD_HEX_Cn_FEM( const int order , const EPointType & pointType );
  
    
  /** \brief  Evaluation of a FEM basis on a <strong>reference Hexahedron</strong> cell. 
    
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference Hexahedron</strong> cell. For rank and dimensions of
              I/O array arguments see Section \ref basis_md_array_sec.
  
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

  /** \brief implement the DofCoordsInterface interface */
  virtual void getDofCoords( ArrayScalar & DofCoords) const;
};
}// namespace Intrepid

#include "Intrepid_HGRAD_HEX_Cn_FEMDef.hpp"

#endif

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
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_HGRAD_HEX_C2_FEMDef.hpp
    \brief  Definition file for bi-linear FEM basis functions for H(grad) functions on Hexahedron cells.
    \author Created by P. Bochev and D. Ridzal.
 */

namespace Intrepid {

  
template<class Scalar, class ArrayScalar>
Basis_HGRAD_HEX_C2_FEM<Scalar,ArrayScalar>::Basis_HGRAD_HEX_C2_FEM()
  {
    this -> basisCardinality_  = 27;
    this -> basisDegree_       = 2;    
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
    this -> basisType_         = BASIS_FEM_DEFAULT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;
  }
  
  
  
template<class Scalar, class ArrayScalar>
void Basis_HGRAD_HEX_C2_FEM<Scalar, ArrayScalar>::initializeTags() {
  
  // Basis-dependent intializations
  int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
  int posScDim = 0;        // poisition in the tag, counting from 0, of the subcell dim 
  int posScOrd = 1;        // poisition in the tag, counting from 0, of the subcell ordinal
  int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

  // An array with local DoF tags assigned to basis functions, in the order of their local enumeration 
  int tags[]  = { 0, 0, 0, 1,     // Nodes 0 to 7 follow vertex order of the topology
                  0, 1, 0, 1,
                  0, 2, 0, 1,
                  0, 3, 0, 1,
                  0, 4, 0, 1,
                  0, 5, 0, 1,
                  0, 6, 0, 1,
                  0, 7, 0, 1,
                  1, 0, 0, 1,      // Node 8  -> edge 0
                  1, 1, 0, 1,      // Node 9  -> edge 1
                  1, 2, 0, 1,      // Node 10 -> edge 2
                  1, 3, 0, 1,      // Node 11 -> edge 3
                  1, 8, 0, 1,      // Node 12 -> edge 8
                  1, 9, 0, 1,      // Node 13 -> edge 9
                  1,10, 0, 1,      // Node 14 -> edge 10
                  1,11, 0, 1,      // Node 15 -> edge 11
                  1, 4, 0, 1,      // Node 16 -> edge 4
                  1, 5, 0, 1,      // Node 17 -> edge 5
                  1, 6, 0, 1,      // Node 18 -> edge 6
                  1, 7, 0, 1,      // Node 19 -> edge 7
                  3, 0, 0, 1,      // Node 20 -> Hexahedron
                  2, 4, 0, 1,      // Node 21 -> face 4
                  2, 5, 0, 1,      // Node 22 -> face 5
                  2, 3, 0, 1,      // Node 23 -> face 3
                  2, 1, 0, 1,      // Node 24 -> face 1
                  2, 0, 0, 1,      // Node 25 -> face 0
                  2, 2, 0, 1,      // Node 26 -> face 2
  };
  
  // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
  Intrepid::setOrdinalTagData(this -> tagToOrdinal_,
                              this -> ordinalToTag_,
                              tags,
                              this -> basisCardinality_,
                              tagSize,
                              posScDim,
                              posScOrd,
                              posDfOrd);
}



template<class Scalar, class ArrayScalar>
void Basis_HGRAD_HEX_C2_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
                                                            const ArrayScalar &  inputPoints,
                                                            const EOperator      operatorType) const {
  
  // Verify arguments
#ifdef HAVE_INTREPID_DEBUG
  Intrepid::getValues_HGRAD_Args<Scalar, ArrayScalar>(outputValues,
                                                      inputPoints,
                                                      operatorType,
                                                      this -> getBaseCellTopology(),
                                                      this -> getCardinality() );
#endif
  
  // Number of evaluation points = dim 0 of inputPoints
  int dim0 = inputPoints.dimension(0);  
  
  // Temporaries: (x,y,z) coordinates of the evaluation point
  Scalar x = 0.0;                                    
  Scalar y = 0.0;                                    
  Scalar z = 0.0;                                    
  
  switch (operatorType) {
    
    case OPERATOR_VALUE:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0, 0);
        y = inputPoints(i0, 1);
        z = inputPoints(i0, 2);
        
        // outputValues is a rank-2 array with dimensions (basisCardinality_, dim0)
        outputValues( 0, i0) = 0.125*(-1. + x)*x*(-1. + y)*y*(-1. + z)*z;
        outputValues( 1, i0) = 0.125*x*(1.+ x)*(-1. + y)*y*(-1. + z)*z;
        outputValues( 2, i0) = 0.125*x*(1.+ x)*y*(1.+ y)*(-1. + z)*z;
        outputValues( 3, i0) = 0.125*(-1. + x)*x*y*(1.+ y)*(-1. + z)*z;
        outputValues( 4, i0) = 0.125*(-1. + x)*x*(-1. + y)*y*z*(1.+ z);
        outputValues( 5, i0) = 0.125*x*(1.+ x)*(-1. + y)*y*z*(1.+ z);
        outputValues( 6, i0) = 0.125*x*(1.+ x)*y*(1.+ y)*z*(1.+ z);
        outputValues( 7, i0) = 0.125*(-1. + x)*x*y*(1.+ y)*z*(1.+ z);
        outputValues( 8, i0) = 0.25*(1. - x)*(1. + x)*(-1. + y)*y*(-1. + z)*z;
        outputValues( 9, i0) = 0.25*x*(1.+ x)*(1. - y)*(1. + y)*(-1. + z)*z;
        outputValues(10, i0) = 0.25*(1. - x)*(1. + x)*y*(1.+ y)*(-1. + z)*z;
        outputValues(11, i0) = 0.25*(-1. + x)*x*(1. - y)*(1. + y)*(-1. + z)*z;
        outputValues(12, i0) = 0.25*(-1. + x)*x*(-1. + y)*y*(1. - z)*(1. + z);
        outputValues(13, i0) = 0.25*x*(1.+ x)*(-1. + y)*y*(1. - z)*(1. + z);
        outputValues(14, i0) = 0.25*x*(1.+ x)*y*(1.+ y)*(1. - z)*(1. + z);
        outputValues(15, i0) = 0.25*(-1. + x)*x*y*(1.+ y)*(1. - z)*(1. + z);
        outputValues(16, i0) = 0.25*(1. - x)*(1. + x)*(-1. + y)*y*z*(1.+ z);
        outputValues(17, i0) = 0.25*x*(1.+ x)*(1. - y)*(1. + y)*z*(1.+ z);
        outputValues(18, i0) = 0.25*(1. - x)*(1. + x)*y*(1.+ y)*z*(1.+ z);
        outputValues(19, i0) = 0.25*(-1. + x)*x*(1. - y)*(1. + y)*z*(1.+ z);
        outputValues(20, i0) = (1. - x)*(1. + x)*(1. - y)*(1. + y)*(1. - z)*(1. + z);
        outputValues(21, i0) = 0.5*(1. - x)*(1. + x)*(1. - y)*(1. + y)*(-1. + z)*z;
        outputValues(22, i0) = 0.5*(1. - x)*(1. + x)*(1. - y)*(1. + y)*z*(1.+ z);
        outputValues(23, i0) = 0.5*(-1. + x)*x*(1. - y)*(1. + y)*(1. - z)*(1. + z);
        outputValues(24, i0) = 0.5*x*(1.+ x)*(1. - y)*(1. + y)*(1. - z)*(1. + z);
        outputValues(25, i0) = 0.5*(1. - x)*(1. + x)*(-1. + y)*y*(1. - z)*(1. + z);
        outputValues(26, i0) = 0.5*(1. - x)*(1. + x)*y*(1.+ y)*(1. - z)*(1. + z);
      }
      break;
      
    case OPERATOR_GRAD:
    case OPERATOR_D1:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        z = inputPoints(i0,2);

        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        outputValues(0, i0, 0) = (-0.125 + 0.25*x)*(-1. + y)*y*(-1. + z)*z;
        outputValues(0, i0, 1) = (-1. + x)*x*(-0.125 + 0.25*y)*(-1. + z)*z;
        outputValues(0, i0, 2) = (-1. + x)*x*(-1. + y)*y*(-0.125 + 0.25*z);
          
        outputValues(1, i0, 0) = (0.125 + 0.25*x)*(-1. + y)*y*(-1. + z)*z;
        outputValues(1, i0, 1) = x*(1. + x)*(-0.125 + 0.25*y)*(-1. + z)*z;
        outputValues(1, i0, 2) = x*(1. + x)*(-1. + y)*y*(-0.125 + 0.25*z);
          
        outputValues(2, i0, 0) = (0.125 + 0.25*x)*y*(1. + y)*(-1. + z)*z;
        outputValues(2, i0, 1) = x*(1. + x)*(0.125 + 0.25*y)*(-1. + z)*z;
        outputValues(2, i0, 2) = x*(1. + x)*y*(1. + y)*(-0.125 + 0.25*z);
          
        outputValues(3, i0, 0) = (-0.125 + 0.25*x)*y*(1. + y)*(-1. + z)*z;
        outputValues(3, i0, 1) = (-1. + x)*x*(0.125 + 0.25*y)*(-1. + z)*z;
        outputValues(3, i0, 2) = (-1. + x)*x*y*(1. + y)*(-0.125 + 0.25*z);
          
        outputValues(4, i0, 0) = (-0.125 + 0.25*x)*(-1. + y)*y*z*(1. + z);
        outputValues(4, i0, 1) = (-1. + x)*x*(-0.125 + 0.25*y)*z*(1. + z);
        outputValues(4, i0, 2) = (-1. + x)*x*(-1. + y)*y*(0.125 + 0.25*z);
          
        outputValues(5, i0, 0) = (0.125 + 0.25*x)*(-1. + y)*y*z*(1. + z);
        outputValues(5, i0, 1) = x*(1. + x)*(-0.125 + 0.25*y)*z*(1. + z);
        outputValues(5, i0, 2) = x*(1. + x)*(-1. + y)*y*(0.125 + 0.25*z);
          
        outputValues(6, i0, 0) = (0.125 + 0.25*x)*y*(1. + y)*z*(1. + z);
        outputValues(6, i0, 1) = x*(1. + x)*(0.125 + 0.25*y)*z*(1. + z);
        outputValues(6, i0, 2) = x*(1. + x)*y*(1. + y)*(0.125 + 0.25*z);
          
        outputValues(7, i0, 0) = (-0.125 + 0.25*x)*y*(1. + y)*z*(1. + z);
        outputValues(7, i0, 1) = (-1. + x)*x*(0.125 + 0.25*y)*z*(1. + z);
        outputValues(7, i0, 2) = (-1. + x)*x*y*(1. + y)*(0.125 + 0.25*z);
          
        outputValues(8, i0, 0) = -0.5*x*(-1. + y)*y*(-1. + z)*z;
        outputValues(8, i0, 1) = (1. - x)*(1. + x)*(-0.25 + 0.5*y)*(-1. + z)*z;
        outputValues(8, i0, 2) = (1. - x)*(1. + x)*(-1. + y)*y*(-0.25 + 0.5*z);
        
        outputValues(9, i0, 0) = (0.25 + 0.5*x)*(1. - y)*(1. + y)*(-1. + z)*z;
        outputValues(9, i0, 1) = x*(1. + x)*(-0.5*y)*(-1. + z)*z;
        outputValues(9, i0, 2) = x*(1. + x)*(1. - y)*(1. + y)*(-0.25 + 0.5*z);
          
        outputValues(10,i0, 0) = -0.5*x*y*(1. + y)*(-1. + z)*z;
        outputValues(10,i0, 1) = (1. - x)*(1. + x)*(0.25 + 0.5*y)*(-1. + z)*z;
        outputValues(10,i0, 2) = (1. - x)*(1. + x)*y*(1. + y)*(-0.25 + 0.5*z);
          
        outputValues(11,i0, 0) = (-0.25 + 0.5*x)*(1. - y)*(1. + y)*(-1. + z)*z;
        outputValues(11,i0, 1) = (-1. + x)*x*(-0.5*y)*(-1. + z)*z;
        outputValues(11,i0, 2) = (-1. + x)*x*(1. - y)*(1. + y)*(-0.25 + 0.5*z);
          
        outputValues(12,i0, 0) = (-0.25 + 0.5*x)*(-1. + y)*y*(1. - z)*(1. + z);
        outputValues(12,i0, 1) = (-1. + x)*x*(-0.25 + 0.5*y)*(1. - z)*(1. + z);
        outputValues(12,i0, 2) = (-1. + x)*x*(-1. + y)*y*(-0.5*z);
          
        outputValues(13,i0, 0) = (0.25 + 0.5*x)*(-1. + y)*y*(1. - z)*(1. + z);
        outputValues(13,i0, 1) = x*(1. + x)*(-0.25 + 0.5*y)*(1. - z)*(1. + z);
        outputValues(13,i0, 2) = x*(1. + x)*(-1. + y)*y*(-0.5*z);
          
        outputValues(14,i0, 0) = (0.25 + 0.5*x)*y*(1. + y)*(1. - z)*(1. + z);
        outputValues(14,i0, 1) = x*(1. + x)*(0.25 + 0.5*y)*(1. - z)*(1. + z);
        outputValues(14,i0, 2) = x*(1. + x)*y*(1. + y)*(-0.5*z);
          
        outputValues(15,i0, 0) = (-0.25 + 0.5*x)*y*(1. + y)*(1. - z)*(1. + z);
        outputValues(15,i0, 1) = (-1. + x)*x*(0.25 + 0.5*y)*(1. - z)*(1. + z);
        outputValues(15,i0, 2) = (-1. + x)*x*y*(1. + y)*(-0.5*z);
          
        outputValues(16,i0, 0) = -0.5*x*(-1. + y)*y*z*(1. + z);
        outputValues(16,i0, 1) = (1. - x)*(1. + x)*(-0.25 + 0.5*y)*z*(1. + z);
        outputValues(16,i0, 2) = (1. - x)*(1. + x)*(-1. + y)*y*(0.25 + 0.5*z);
          
        outputValues(17,i0, 0) = (0.25 + 0.5*x)*(1. - y)*(1. + y)*z*(1. + z);
        outputValues(17,i0, 1) = x*(1. + x)*(-0.5*y)*z*(1. + z);
        outputValues(17,i0, 2) = x*(1. + x)*(1. - y)*(1. + y)*(0.25 + 0.5*z);
          
        outputValues(18,i0, 0) = -0.5*x*y*(1. + y)*z*(1. + z);
        outputValues(18,i0, 1) = (1. - x)*(1. + x)*(0.25 + 0.5*y)*z*(1. + z);
        outputValues(18,i0, 2) = (1. - x)*(1. + x)*y*(1. + y)*(0.25 + 0.5*z);
          
        outputValues(19,i0, 0) = (-0.25 + 0.5*x)*(1. - y)*(1. + y)*z*(1. + z);
        outputValues(19,i0, 1) = (-1. + x)*x*(-0.5*y)*z*(1. + z);
        outputValues(19,i0, 2) = (-1. + x)*x*(1. - y)*(1. + y)*(0.25 + 0.5*z);
          
        outputValues(20,i0, 0) = -2.*x*(1. - y)*(1. + y)*(1. - z)*(1. + z);
        outputValues(20,i0, 1) = (1. - x)*(1. + x)*(-2.*y)*(1. - z)*(1. + z);
        outputValues(20,i0, 2) = (1. - x)*(1. + x)*(1. - y)*(1. + y)*(-2.*z);
          
        outputValues(21,i0, 0) = -x*(1. - y)*(1. + y)*(-1. + z)*z;
        outputValues(21,i0, 1) = (1. - x)*(1. + x)*(-y)*(-1. + z)*z;
        outputValues(21,i0, 2) = (1. - x)*(1. + x)*(1. - y)*(1. + y)*(-0.5 + z);
          
        outputValues(22,i0, 0) = -x*(1. - y)*(1. + y)*z*(1. + z);
        outputValues(22,i0, 1) = (1. - x)*(1. + x)*(-y)*z*(1. + z);
        outputValues(22,i0, 2) = (1. - x)*(1. + x)*(1. - y)*(1. + y)*(0.5 + z);
          
        outputValues(23,i0, 0) = (-0.5 + x)*(1. - y)*(1. + y)*(1. - z)*(1. + z);
        outputValues(23,i0, 1) = (-1. + x)*x*(-y)*(1. - z)*(1. + z);
        outputValues(23,i0, 2) = (-1. + x)*x*(1. - y)*(1. + y)*(-z);
          
        outputValues(24,i0, 0) = (0.5 + x)*(1. - y)*(1. + y)*(1. - z)*(1. + z);
        outputValues(24,i0, 1) = x*(1. + x)*(-y)*(1. - z)*(1. + z);
        outputValues(24,i0, 2) = x*(1. + x)*(1. - y)*(1. + y)*(-z);
          
        outputValues(25,i0, 0) = -x*(-1. + y)*y*(1. - z)*(1. + z);
        outputValues(25,i0, 1) = (1. - x)*(1. + x)*(-0.5 + y)*(1. - z)*(1. + z);
        outputValues(25,i0, 2) = (1. - x)*(1. + x)*(-1. + y)*y*(-z);
          
        outputValues(26,i0, 0) = -x*y*(1. + y)*(1. - z)*(1. + z);
        outputValues(26,i0, 1) = (1. - x)*(1. + x)*(0.5 + y)*(1. - z)*(1. + z);
        outputValues(26,i0, 2) = (1. - x)*(1. + x)*y*(1. + y)*(-z);
      }
      break;
      
    case OPERATOR_CURL:
      TEST_FOR_EXCEPTION( (operatorType == OPERATOR_CURL), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_HEX_C2_FEM): CURL is invalid operator for rank-0 (scalar) functions in 3D");
      break;
      
    case OPERATOR_DIV:
      TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_HEX_C2_FEM): DIV is invalid operator for rank-0 (scalar) functions in 3D");
      break;
      
    case OPERATOR_D2:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        z = inputPoints(i0,2);
        
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, D2Cardinality = 6) 
        outputValues(0, i0, 0) = 0.25*(-1. + y)*y*(-1. + z)*z; 
        outputValues(0, i0, 1) = (-0.125 + y*(0.25 - 0.25*z) + x*(0.25 + y*(-0.5 + 0.5*z) - 0.25*z) + 0.125*z)*z; 
        outputValues(0, i0, 2) = y*(-0.125 + x*(0.25 + y*(-0.25 + 0.5*z) - 0.5*z) + y*(0.125 - 0.25*z) + 0.25*z); 
        outputValues(0, i0, 3) = 0.25*(-1. + x)*x*(-1. + z)*z; 
        outputValues(0, i0, 4) = x*(-0.125 + y*(0.25 - 0.5*z) + x*(0.125 + y*(-0.25 + 0.5*z) - 0.25*z) + 0.25*z); 
        outputValues(0, i0, 5) = 0.25*(-1. + x)*x*(-1. + y)*y; 
        
        outputValues(1, i0, 0) = 0.25*(-1. + y)*y*(-1. + z)*z; 
        outputValues(1, i0, 1) = (0.125 + x*(0.25 + y*(-0.5 + 0.5*z) - 0.25*z) + y*(-0.25 + 0.25*z) - 0.125*z)*z; 
        outputValues(1, i0, 2) = y*(0.125 + x*(0.25 + y*(-0.25 + 0.5*z) - 0.5*z) + y*(-0.125 + 0.25*z) - 0.25*z); 
        outputValues(1, i0, 3) = 0.25*x*(1 + x)*(-1. + z)*z; 
        outputValues(1, i0, 4) = x*(1. + x)*(0.125 + y*(-0.25 + 0.5*z) - 0.25*z); 
        outputValues(1, i0, 5) = 0.25*x*(1 + x)*(-1. + y)*y; 
        
        outputValues(2, i0, 0) = 0.25*y*(1 + y)*(-1. + z)*z; 
        outputValues(2, i0, 1) = (0.125 + x*(0.25 + 0.5*y) + 0.25*y)*(-1. + z)*z; 
        outputValues(2, i0, 2) = y*(1. + y)*(-0.125 + x*(-0.25 + 0.5*z) + 0.25*z); 
        outputValues(2, i0, 3) = 0.25*x*(1 + x)*(-1. + z)*z; 
        outputValues(2, i0, 4) = x*(1. + x)*(-0.125 + y*(-0.25 + 0.5*z) + 0.25*z); 
        outputValues(2, i0, 5) = 0.25*x*(1 + x)*y*(1 + y); 
        
        outputValues(3, i0, 0) = 0.25*y*(1 + y)*(-1. + z)*z; 
        outputValues(3, i0, 1) = (0.125 + y*(0.25 - 0.25*z) + x*(-0.25 + y*(-0.5 + 0.5*z) + 0.25*z) - 0.125*z)*z; 
        outputValues(3, i0, 2) = y*(1. + y)*(0.125 + x*(-0.25 + 0.5*z) - 0.25*z); 
        outputValues(3, i0, 3) = 0.25*(-1. + x)*x*(-1. + z)*z; 
        outputValues(3, i0, 4) = x*(0.125 + y*(0.25 - 0.5*z) + x*(-0.125 + y*(-0.25 + 0.5*z) + 0.25*z) - 0.25*z); 
        outputValues(3, i0, 5) = 0.25*(-1. + x)*x*y*(1 + y); 
        
        outputValues(4, i0, 0) = 0.25*(-1. + y)*y*z*(1 + z); 
        outputValues(4, i0, 1) = (0.125 + x*(-0.25 + 0.5*y) - 0.25*y)*z*(1. + z); 
        outputValues(4, i0, 2) = y*(0.125 + x*(-0.25 + y*(0.25 + 0.5*z) - 0.5*z) + y*(-0.125 - 0.25*z) + 0.25*z); 
        outputValues(4, i0, 3) = 0.25*(-1. + x)*x*z*(1 + z); 
        outputValues(4, i0, 4) = x*(0.125 + y*(-0.25 - 0.5*z) + x*(-0.125 + y*(0.25 + 0.5*z) - 0.25*z) + 0.25*z); 
        outputValues(4, i0, 5) = 0.25*(-1. + x)*x*(-1. + y)*y; 
        
        outputValues(5, i0, 0) = 0.25*(-1. + y)*y*z*(1 + z); 
        outputValues(5, i0, 1) = (-0.125 + x*(-0.25 + 0.5*y) + 0.25*y)*z*(1. + z); 
        outputValues(5, i0, 2) = (-1. + y)*y*(0.125 + x*(0.25 + 0.5*z) + 0.25*z); 
        outputValues(5, i0, 3) = 0.25*x*(1 + x)*z*(1 + z); 
        outputValues(5, i0, 4) = x*(1. + x)*(-0.125 + y*(0.25 + 0.5*z) - 0.25*z); 
        outputValues(5, i0, 5) = 0.25*x*(1 + x)*(-1. + y)*y; 
        
        outputValues(6, i0, 0) = 0.25*y*(1 + y)*z*(1 + z); 
        outputValues(6, i0, 1) = (0.125 + x*(0.25 + 0.5*y) + 0.25*y)*z*(1. + z); 
        outputValues(6, i0, 2) = y*(1. + y)*(0.125 + x*(0.25 + 0.5*z) + 0.25*z); 
        outputValues(6, i0, 3) = 0.25*x*(1 + x)*z*(1 + z); 
        outputValues(6, i0, 4) = x*(1. + x)*(0.125 + y*(0.25 + 0.5*z) + 0.25*z); 
        outputValues(6, i0, 5) = 0.25*x*(1 + x)*y*(1 + y); 
        
        outputValues(7, i0, 0) = 0.25*y*(1 + y)*z*(1 + z); 
        outputValues(7, i0, 1) = (-0.125 + x*(0.25 + 0.5*y) - 0.25*y)*z*(1. + z); 
        outputValues(7, i0, 2) = y*(1. + y)*(-0.125 + x*(0.25 + 0.5*z) - 0.25*z); 
        outputValues(7, i0, 3) = 0.25*(-1. + x)*x*z*(1 + z); 
        outputValues(7, i0, 4) = (-1. + x)*x*(0.125 + y*(0.25 + 0.5*z) + 0.25*z); 
        outputValues(7, i0, 5) = 0.25*(-1. + x)*x*y*(1 + y); 
        
        outputValues(8, i0, 0) = -0.5*(-1. + y)*y*(-1. + z)*z; 
        outputValues(8, i0, 1) = (0.  + x*(-0.5 + y))*z + (x*(0.5 - y) )*(z*z); 
        outputValues(8, i0, 2) = (y*y)*(x*(0.5 - z) ) + y*(x*(-0.5 + z)); 
        outputValues(8, i0, 3) = 0.5*(1. - x)*(1. + x)*(-1. + z)*z; 
        outputValues(8, i0, 4) = 0.25 + (x*x)*(-0.25 + y*(0.5 - z) + 0.5*z) - 0.5*z + y*(-0.5 + z); 
        outputValues(8, i0, 5) = 0.5*(1. - x)*(1. + x)*(-1. + y)*y; 
        
        outputValues(9, i0, 0) = 0.5*(1. - y)*(1. + y)*(-1. + z)*z; 
        outputValues(9, i0, 1) = (0.5*y + x*(y))*z + (x*(-y) - 0.5*y)*(z*z); 
        outputValues(9, i0, 2) = -0.25 + (y*y)*(0.25 - 0.5*z) + 0.5*z + x*(-0.5 + (y*y)*(0.5 - z) + z); 
        outputValues(9, i0, 3) = -0.5*x*(1 + x)*(-1. + z)*z; 
        outputValues(9, i0, 4) = x*(y*(0.5 - z) ) + (x*x)*(y*(0.5 - z) ); 
        outputValues(9, i0, 5) = 0.5*x*(1 + x)*(1. - y)*(1. + y); 
        
        outputValues(10,i0, 0) = -0.5*y*(1 + y)*(-1. + z)*z; 
        outputValues(10,i0, 1) = (0.  + x*(0.5 + y))*z + (x*(-0.5 - y) )*(z*z); 
        outputValues(10,i0, 2) = y*(x*(0.5 - z) ) + (y*y)*(x*(0.5 - z) ); 
        outputValues(10,i0, 3) = 0.5*(1. - x)*(1. + x)*(-1. + z)*z; 
        outputValues(10,i0, 4) = -0.25 + (x*x)*(0.25 + y*(0.5 - z) - 0.5*z) + 0.5*z + y*(-0.5 + z); 
        outputValues(10,i0, 5) = 0.5*(1. - x)*(1. + x)*y*(1 + y); 
        
        outputValues(11,i0, 0) = 0.5*(1. - y)*(1. + y)*(-1. + z)*z; 
        outputValues(11,i0, 1) = (-0.5*y + x*(y))*z + (x*(-y) + 0.5*y)*(z*z); 
        outputValues(11,i0, 2) = 0.25 + (y*y)*(-0.25 + 0.5*z) - 0.5*z + x*(-0.5 + (y*y)*(0.5 - z) + z); 
        outputValues(11,i0, 3) = -0.5*(-1. + x)*x*(-1. + z)*z; 
        outputValues(11,i0, 4) = (x*x)*(y*(0.5 - z) ) + x*(y*(-0.5 + z)); 
        outputValues(11,i0, 5) = 0.5*(-1. + x)*x*(1. - y)*(1. + y); 
        
        outputValues(12,i0, 0) = 0.5*(-1. + y)*y*(1. - z)*(1. + z); 
        outputValues(12,i0, 1) = 0.25  - 0.25*(z*z) + y*(-0.5  + 0.5*(z*z)) + x*(-0.5  + 0.5*(z*z) + y*(1.  - (z*z))); 
        outputValues(12,i0, 2) = (y*y)*(x*(-z) + 0.5*z) + y*(-0.5*z + x*(z)); 
        outputValues(12,i0, 3) = 0.5*(-1. + x)*x*(1. - z)*(1. + z); 
        outputValues(12,i0, 4) = (x*x)*(y*(-z) + 0.5*z) + x*(-0.5*z + y*(z)); 
        outputValues(12,i0, 5) = -0.5*(-1. + x)*x*(-1. + y)*y; 
        
        outputValues(13,i0, 0) = 0.5*(-1. + y)*y*(1. - z)*(1. + z); 
        outputValues(13,i0, 1) = -0.25  + 0.25*(z*z) + y*(0.5  - 0.5*(z*z)) + x*(-0.5  + 0.5*(z*z) + y*(1.  - (z*z))); 
        outputValues(13,i0, 2) = (y*y)*(x*(-z) - 0.5*z) + y*(0.5*z + x*(z)); 
        outputValues(13,i0, 3) = 0.5*x*(1 + x)*(1. - z)*(1. + z); 
        outputValues(13,i0, 4) = x*(y*(-z) + 0.5*z) + (x*x)*(y*(-z) + 0.5*z); 
        outputValues(13,i0, 5) = -0.5*x*(1 + x)*(-1. + y)*y; 

        outputValues(14,i0, 0) = 0.5*y*(1 + y)*(1. - z)*(1. + z); 
        outputValues(14,i0, 1) = 0.25  - 0.25*(z*z) + y*(0.5  - 0.5*(z*z)) + x*(0.5  - 0.5*(z*z) + y*(1.  - (z*z))); 
        outputValues(14,i0, 2) = y*(x*(-z) - 0.5*z) + (y*y)*(x*(-z) - 0.5*z); 
        outputValues(14,i0, 3) = 0.5*x*(1 + x)*(1. - z)*(1. + z); 
        outputValues(14,i0, 4) = x*(y*(-z) - 0.5*z) + (x*x)*(y*(-z) - 0.5*z); 
        outputValues(14,i0, 5) = -0.5*x*(1 + x)*y*(1 + y); 
        
        outputValues(15,i0, 0) = 0.5*y*(1 + y)*(1. - z)*(1. + z); 
        outputValues(15,i0, 1) = -0.25  + 0.25*(z*z) + y*(-0.5  + 0.5*(z*z)) + x*(0.5  - 0.5*(z*z) + y*(1.  - (z*z))); 
        outputValues(15,i0, 2) = y*(x*(-z) + 0.5*z) + (y*y)*(x*(-z) + 0.5*z); 
        outputValues(15,i0, 3) = 0.5*(-1. + x)*x*(1. - z)*(1. + z); 
        outputValues(15,i0, 4) = (x*x)*(y*(-z) - 0.5*z) + x*(0.5*z + y*(z)); 
        outputValues(15,i0, 5) = -0.5*(-1. + x)*x*y*(1 + y); 
        
        outputValues(16,i0, 0) = -0.5*(-1. + y)*y*z*(1 + z); 
        outputValues(16,i0, 1) = (x*(0.5 - y) )*z + (x*(0.5 - y) )*(z*z); 
        outputValues(16,i0, 2) = (y*y)*(x*(-0.5 - z) ) + y*(x*(0.5 + z)); 
        outputValues(16,i0, 3) = 0.5*(1. - x)*(1. + x)*z*(1 + z); 
        outputValues(16,i0, 4) = -0.25 + (x*x)*(0.25 + y*(-0.5 - z) + 0.5*z) - 0.5*z + y*(0.5 + z); 
        outputValues(16,i0, 5) = 0.5*(1. - x)*(1. + x)*(-1. + y)*y; 
        
        outputValues(17,i0, 0) = 0.5*(1. - y)*(1. + y)*z*(1 + z); 
        outputValues(17,i0, 1) = (x*(-y) - 0.5*y)*z + (x*(-y) - 0.5*y)*(z*z); 
        outputValues(17,i0, 2) = 0.25 + (y*y)*(-0.25 - 0.5*z) + 0.5*z + x*(0.5 + (y*y)*(-0.5 - z) + z); 
        outputValues(17,i0, 3) = -0.5*x*(1 + x)*z*(1 + z); 
        outputValues(17,i0, 4) = x*(y*(-0.5 - z) ) + (x*x)*(y*(-0.5 - z) ); 
        outputValues(17,i0, 5) = 0.5*x*(1 + x)*(1. - y)*(1. + y); 
        
        outputValues(18,i0, 0) = -0.5*y*(1 + y)*z*(1 + z); 
        outputValues(18,i0, 1) = (x*(-0.5 - y) )*z + (x*(-0.5 - y) )*(z*z); 
        outputValues(18,i0, 2) = y*(x*(-0.5 - z) ) + (y*y)*(x*(-0.5 - z) ); 
        outputValues(18,i0, 3) = 0.5*(1. - x)*(1. + x)*z*(1 + z); 
        outputValues(18,i0, 4) = 0.25 + (x*x)*(-0.25 + y*(-0.5 - z) - 0.5*z) + 0.5*z + y*(0.5 + z); 
        outputValues(18,i0, 5) = 0.5*(1. - x)*(1. + x)*y*(1 + y); 
        
        outputValues(19,i0, 0) = 0.5*(1. - y)*(1. + y)*z*(1 + z); 
        outputValues(19,i0, 1) = (x*(-y) + 0.5*y)*z + (x*(-y) + 0.5*y)*(z*z); 
        outputValues(19,i0, 2) = -0.25 + (y*y)*(0.25 + 0.5*z) - 0.5*z + x*(0.5 + (y*y)*(-0.5 - z) + z); 
        outputValues(19,i0, 3) = -0.5*(-1. + x)*x*z*(1 + z); 
        outputValues(19,i0, 4) = (x*x)*(y*(-0.5 - z) ) + x*(y*(0.5 + z)); 
        outputValues(19,i0, 5) = 0.5*(-1. + x)*x*(1. - y)*(1. + y); 
        
        outputValues(20,i0, 0) = -2.*(1. - y)*(1. + y)*(1. - z)*(1. + z); 
        outputValues(20,i0, 1) = -4.*x*y*(-1. + z*z);         
        outputValues(20,i0, 2) = x*((y*y)*(-4.*z) + 4.*z); 
        outputValues(20,i0, 3) = -2.*(1. - x)*(1. + x)*(1. - z)*(1. + z); 
        outputValues(20,i0, 4) = (x*x)*(y*(-4.*z) ) + y*(4.*z); 
        outputValues(20,i0, 5) = -2.*(1. - x)*(1. + x)*(1. - y)*(1. + y); 
        
        outputValues(21,i0, 0) = -(1. - y)*(1. + y)*(-1. + z)*z; 
        outputValues(21,i0, 1) = (x*(-2.*y) )*z + (0.  + x*(2.*y))*(z*z); 
        outputValues(21,i0, 2) =  x*(1. - 2.*z + (y*y)*(-1. + 2.*z)); 
        outputValues(21,i0, 3) = -(1. - x)*(1. + x)*(-1. + z)*z; 
        outputValues(21,i0, 4) = y*(1. - 2.*z)  + (x*x)*(y*(-1. + 2.*z)); 
        outputValues(21,i0, 5) = (1. - x)*(1. + x)*(1. - y)*(1. + y); 
        
        outputValues(22,i0, 0) = -(1. - y)*(1. + y)*z*(1 + z); 
        outputValues(22,i0, 1) = (0.  + x*(2.*y))*z + (0.  + x*(2.*y))*(z*z); 
        outputValues(22,i0, 2) = x*(-1. - 2.*z + (y*y)*(1. + 2.*z)); 
        outputValues(22,i0, 3) = -(1. - x)*(1. + x)*z*(1 + z); 
        outputValues(22,i0, 4) = y*(-1. - 2.*z) + (x*x)*(y*(1. + 2.*z)); 
        outputValues(22,i0, 5) = (1. - x)*(1. + x)*(1. - y)*(1. + y); 
        
        outputValues(23,i0, 0) = (1. - y)*(1. + y)*(1. - z)*(1. + z); 
        outputValues(23,i0, 1) = (-1. + 2.*x)*y*(-1. + z*z);
        outputValues(23,i0, 2) = (-1. + 2.*x)*(-1. + y*y)*z; 
        outputValues(23,i0, 3) =-(-1. + x)*x*(1. - z)*(1. + z); 
        outputValues(23,i0, 4) =  2.*(-1. + x)*x*y*z; 
        outputValues(23,i0, 5) =-(-1. + x)*x*(1. - y)*(1. + y); 
        
        outputValues(24,i0, 0) = (1. - y)*(1. + y)*(1. - z)*(1. + z); 
        outputValues(24,i0, 1) = (1. + 2.*x)*y*(-1. + z*z); 
        outputValues(24,i0, 2) = (1. + 2.*x)*(-1. + y*y)*z; 
        outputValues(24,i0, 3) = x*(1. + x)*(-1. + z)*(1. + z); 
        outputValues(24,i0, 4) = 2.*x*(1. + x)*y*z; 
        outputValues(24,i0, 5) = x*(1. + x)*(-1. + y)*(1. + y); 
        
        outputValues(25,i0, 0) = -(-1. + y)*y*(1. - z)*(1. + z); 
        outputValues(25,i0, 1) = x*(-1. + 2.*y)*(-1. + z*z); 
        outputValues(25,i0, 2) = 2.*x*(-1. + y)*y*z; 
        outputValues(25,i0, 3) = (1. - x)*(1. + x)*(1. - z)*(1. + z); 
        outputValues(25,i0, 4) = (-1. + x*x)*(-1. + 2.*y)*z; 
        outputValues(25,i0, 5) =-(1. - x)*(1. + x)*(-1. + y)*y; 

        outputValues(26,i0, 0) =  y*(1. + y)*(-1. + z)*(1. + z);
        outputValues(26,i0, 1) =  x*(1. + 2.*y)*(-1. + z*z);
        outputValues(26,i0, 2) =  2.*x*y*(1. + y)*z;
        outputValues(26,i0, 3) =  (-1. + x)*(1. + x)*(-1. + z)*(1. + z);
        outputValues(26,i0, 4) =  (-1. + x*x)*(1. + 2.*y)*z;
        outputValues(26,i0, 5) =  (-1. + x)*(1. + x)*y*(1. + y);
      }
      break;
      
    case OPERATOR_D3:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        z = inputPoints(i0,2);
                
        outputValues(0,i0, 0) = 0.;
        outputValues(0,i0, 1) = ((-1.+ 2.*y)*(-1.+ z)*z)/4.;
        outputValues(0,i0, 2) = ((-1.+ y)*y*(-1.+ 2.*z))/4.;
        outputValues(0,i0, 3) = ((-1.+ 2.*x)*(-1.+ z)*z)/4.;
        outputValues(0,i0, 4) = ((-1.+ 2.*x)*(-1.+ 2.*y)*(-1.+ 2.*z))/8.;
        outputValues(0,i0, 5) = ((-1.+ 2.*x)*(-1.+ y)*y)/4.;
        outputValues(0,i0, 6) = 0.;    
        outputValues(0,i0, 7) = ((-1.+ x)*x*(-1.+ 2.*z))/4.;
        outputValues(0,i0, 8) = ((-1.+ x)*x*(-1.+ 2.*y))/4.;
        outputValues(0,i0, 9) = 0.;
                
        outputValues(1, i0, 0) = 0.;
        outputValues(1, i0, 1) = ((-1.+ 2.*y)*(-1.+ z)*z)/4.;
        outputValues(1, i0, 2) = ((-1.+ y)*y*(-1.+ 2.*z))/4.;
        outputValues(1, i0, 3) = ((1.+ 2.*x)*(-1.+ z)*z)/4.;
        outputValues(1, i0, 4) = ((1.+ 2.*x)*(-1.+ 2.*y)*(-1.+ 2.*z))/8.;
        outputValues(1, i0, 5) = ((1.+ 2.*x)*(-1.+ y)*y)/4.;
        outputValues(1, i0, 6) = 0.;
        outputValues(1, i0, 7) = (x*(1.+ x)*(-1.+ 2.*z))/4.;
        outputValues(1, i0, 8) = (x*(1.+ x)*(-1.+ 2.*y))/4.;
        outputValues(1, i0, 9) = 0.;
        
        outputValues(2, i0, 0) = 0.;
        outputValues(2, i0, 1) = ((1.+ 2.*y)*(-1.+ z)*z)/4.;
        outputValues(2, i0, 2) = (y*(1.+ y)*(-1.+ 2.*z))/4.;
        outputValues(2, i0, 3) = ((1.+ 2.*x)*(-1.+ z)*z)/4.;
        outputValues(2, i0, 4) = ((1.+ 2.*x)*(1.+ 2.*y)*(-1.+ 2.*z))/8.;
        outputValues(2, i0, 5) = ((1.+ 2.*x)*y*(1.+ y))/4.;
        outputValues(2, i0, 6) = 0.;
        outputValues(2, i0, 7) = (x*(1.+ x)*(-1.+ 2.*z))/4.;
        outputValues(2, i0, 8) = (x*(1.+ x)*(1.+ 2.*y))/4.;
        outputValues(2, i0, 9) = 0.;
        
        outputValues(3, i0, 0) = 0.;
        outputValues(3, i0, 1) = ((1.+ 2.*y)*(-1.+ z)*z)/4.;
        outputValues(3, i0, 2) = (y*(1.+ y)*(-1.+ 2.*z))/4.;
        outputValues(3, i0, 3) = ((-1.+ 2.*x)*(-1.+ z)*z)/4.;
        outputValues(3, i0, 4) = ((-1.+ 2.*x)*(1.+ 2.*y)*(-1.+ 2.*z))/8.;
        outputValues(3, i0, 5) = ((-1.+ 2.*x)*y*(1.+ y))/4.;
        outputValues(3, i0, 6) = 0.;
        outputValues(3, i0, 7) = ((-1.+ x)*x*(-1.+ 2.*z))/4.;
        outputValues(3, i0, 8) = ((-1.+ x)*x*(1.+ 2.*y))/4.;
        outputValues(3, i0, 9) = 0.;
        
        outputValues(4, i0, 0) = 0.;
        outputValues(4, i0, 1) = ((-1.+ 2.*y)*z*(1.+ z))/4.;
        outputValues(4, i0, 2) = ((-1.+ y)*y*(1.+ 2.*z))/4.;
        outputValues(4, i0, 3) = ((-1.+ 2.*x)*z*(1.+ z))/4.;
        outputValues(4, i0, 4) = ((-1.+ 2.*x)*(-1.+ 2.*y)*(1.+ 2.*z))/8.;
        outputValues(4, i0, 5) = ((-1.+ 2.*x)*(-1.+ y)*y)/4.;
        outputValues(4, i0, 6) = 0.;
        outputValues(4, i0, 7) = ((-1.+ x)*x*(1.+ 2.*z))/4.;
        outputValues(4, i0, 8) = ((-1.+ x)*x*(-1.+ 2.*y))/4.;
        outputValues(4, i0, 9) = 0.;
        
        outputValues(5, i0, 0) = 0.;
        outputValues(5, i0, 1) = ((-1.+ 2.*y)*z*(1.+ z))/4.;
        outputValues(5, i0, 2) = ((-1.+ y)*y*(1.+ 2.*z))/4.;
        outputValues(5, i0, 3) = ((1.+ 2.*x)*z*(1.+ z))/4.;
        outputValues(5, i0, 4) = ((1.+ 2.*x)*(-1.+ 2.*y)*(1.+ 2.*z))/8.;
        outputValues(5, i0, 5) = ((1.+ 2.*x)*(-1.+ y)*y)/4.;
        outputValues(5, i0, 6) = 0.;
        outputValues(5, i0, 7) = (x*(1.+ x)*(1.+ 2.*z))/4.;
        outputValues(5, i0, 8) = (x*(1.+ x)*(-1.+ 2.*y))/4.;
        outputValues(5, i0, 9) = 0.;
        
        outputValues(6, i0, 0) = 0.;
        outputValues(6, i0, 1) = ((1.+ 2.*y)*z*(1.+ z))/4.;
        outputValues(6, i0, 2) = (y*(1.+ y)*(1.+ 2.*z))/4.;
        outputValues(6, i0, 3) = ((1.+ 2.*x)*z*(1.+ z))/4.;
        outputValues(6, i0, 4) = ((1.+ 2.*x)*(1.+ 2.*y)*(1.+ 2.*z))/8.;
        outputValues(6, i0, 5) = ((1.+ 2.*x)*y*(1.+ y))/4.;
        outputValues(6, i0, 6) = 0.;
        outputValues(6, i0, 7) = (x*(1.+ x)*(1.+ 2.*z))/4.;
        outputValues(6, i0, 8) = (x*(1.+ x)*(1.+ 2.*y))/4.;
        outputValues(6, i0, 9) = 0.;
        
        outputValues(7, i0, 0) = 0.;
        outputValues(7, i0, 1) = ((1.+ 2.*y)*z*(1.+ z))/4.;
        outputValues(7, i0, 2) = (y*(1.+ y)*(1.+ 2.*z))/4.;
        outputValues(7, i0, 3) = ((-1.+ 2.*x)*z*(1.+ z))/4.;
        outputValues(7, i0, 4) = ((-1.+ 2.*x)*(1.+ 2.*y)*(1.+ 2.*z))/8.;
        outputValues(7, i0, 5) = ((-1.+ 2.*x)*y*(1.+ y))/4.;
        outputValues(7, i0, 6) = 0.;
        outputValues(7, i0, 7) = ((-1.+ x)*x*(1.+ 2.*z))/4.;
        outputValues(7, i0, 8) = ((-1.+ x)*x*(1.+ 2.*y))/4.;
        outputValues(7, i0, 9) = 0.; 
        
        outputValues(8, i0, 0) = 0.;
        outputValues(8, i0, 1) = -((-1.+ 2.*y)*(-1.+ z)*z)/2.;
        outputValues(8, i0, 2) = -((-1.+ y)*y*(-1.+ 2.*z))/2.;
        outputValues(8, i0, 3) = -(x*(-1.+ z)*z);
        outputValues(8, i0, 4) = -(x*(-1.+ 2.*y)*(-1.+ 2.*z))/2.;
        outputValues(8, i0, 5) = -(x*(-1.+ y)*y);
        outputValues(8, i0, 6) = 0.;
        outputValues(8, i0, 7) = -((-1.+ (x*x))*(-1.+ 2.*z))/2.;    
        outputValues(8, i0, 8) = -((-1.+ (x*x))*(-1.+ 2.*y))/2.;
        outputValues(8, i0, 9) = 0.;
        
        outputValues(9, i0, 0) = 0.;
        outputValues(9, i0, 1) = -(y*(-1.+ z)*z);
        outputValues(9, i0, 2) = -((-1.+ (y*y))*(-1.+ 2.*z))/2.;
        outputValues(9, i0, 3) = -((1.+ 2.*x)*(-1.+ z)*z)/2.;
        outputValues(9, i0, 4) = -((1.+ 2.*x)*y*(-1.+ 2.*z))/2.;
        outputValues(9, i0, 5) = -((1.+ 2.*x)*(-1.+ (y*y)))/2.;
        outputValues(9, i0, 6) = 0.;
        outputValues(9, i0, 7) = -(x*(1.+ x)*(-1.+ 2.*z))/2.;
        outputValues(9, i0, 8) = -(x*(1.+ x)*y);
        outputValues(9, i0, 9) = 0.; 
        
        outputValues(10,i0, 0) = 0.;
        outputValues(10,i0, 1) = -((1.+ 2.*y)*(-1.+ z)*z)/2.;
        outputValues(10,i0, 2) = -(y*(1.+ y)*(-1.+ 2.*z))/2.;
        outputValues(10,i0, 3) = -(x*(-1.+ z)*z);
        outputValues(10,i0, 4) = -(x*(1.+ 2.*y)*(-1.+ 2.*z))/2.;
        outputValues(10,i0, 5) = -(x*y*(1.+ y));
        outputValues(10,i0, 6) = 0.;
        outputValues(10,i0, 7) =  -((-1.+ (x*x))*(-1.+ 2.*z))/2.;    
        outputValues(10,i0, 8) = -((-1.+ (x*x))*(1.+ 2.*y))/2.;
        outputValues(10,i0, 9) = 0.;
        
        outputValues(11,i0, 0) = 0.;
        outputValues(11,i0, 1) = -(y*(-1.+ z)*z);
        outputValues(11,i0, 2) = -((-1.+ (y*y))*(-1.+ 2.*z))/2.;
        outputValues(11,i0, 3) = -((-1.+ 2.*x)*(-1.+ z)*z)/2.;
        outputValues(11,i0, 4) = -((-1.+ 2.*x)*y*(-1.+ 2.*z))/2.;
        outputValues(11,i0, 5) = -((-1.+ 2.*x)*(-1.+ (y*y)))/2.;
        outputValues(11,i0, 6) = 0.;
        outputValues(11,i0, 7) = -((-1.+ x)*x*(-1.+ 2.*z))/2.;
        outputValues(11,i0, 8) = -((-1.+ x)*x*y);
        outputValues(11,i0, 9) = 0.;   
        
        outputValues(12,i0, 0) = 0.;
        outputValues(12,i0, 1) = -((-1.+ 2.*y)*(-1.+ (z*z)))/2.;
        outputValues(12,i0, 2) = -((-1.+ y)*y*z);
        outputValues(12,i0, 3) = -((-1.+ 2.*x)*(-1.+ (z*z)))/2.;
        outputValues(12,i0, 4) = -((-1.+ 2.*x)*(-1.+ 2.*y)*z)/2.;
        outputValues(12,i0, 5) = -((-1.+ 2.*x)*(-1.+ y)*y)/2.;
        outputValues(12,i0, 6) = 0.;
        outputValues(12,i0, 7) = -((-1.+ x)*x*z);    
        outputValues(12,i0, 8) = -((-1.+ x)*x*(-1.+ 2.*y))/2.;
        outputValues(12,i0, 9) = 0.;
        
        outputValues(13,i0, 0) = 0.;
        outputValues(13,i0, 1) = -((-1.+ 2.*y)*(-1.+ (z*z)))/2.;
        outputValues(13,i0, 2) = -((-1.+ y)*y*z);
        outputValues(13,i0, 3) = -((1.+ 2.*x)*(-1.+ (z*z)))/2.;
        outputValues(13,i0, 4) = -((1.+ 2.*x)*(-1.+ 2.*y)*z)/2.;
        outputValues(13,i0, 5) = -((1.+ 2.*x)*(-1.+ y)*y)/2.;
        outputValues(13,i0, 6) = 0.;
        outputValues(13,i0, 7) = -(x*(1.+ x)*z);
        outputValues(13,i0, 8) = -(x*(1.+ x)*(-1.+ 2.*y))/2.;
        outputValues(13,i0, 9) = 0.; 
        
        outputValues(14,i0, 0) = 0.;
        outputValues(14,i0, 1) = -((1.+ 2.*y)*(-1.+ (z*z)))/2.;
        outputValues(14,i0, 2) = -(y*(1.+ y)*z);
        outputValues(14,i0, 3) = -((1.+ 2.*x)*(-1.+ (z*z)))/2.;
        outputValues(14,i0, 4) = -((1.+ 2.*x)*(1.+ 2.*y)*z)/2.;
        outputValues(14,i0, 5) = -((1.+ 2.*x)*y*(1.+ y))/2.;
        outputValues(14,i0, 6) = 0.;
        outputValues(14,i0, 7) = -(x*(1.+ x)*z);    
        outputValues(14,i0, 8) = -(x*(1.+ x)*(1.+ 2.*y))/2.;
        outputValues(14,i0, 9) = 0.;
        
        outputValues(15,i0, 0) = 0.;
        outputValues(15,i0, 1) = -((1.+ 2.*y)*(-1.+ (z*z)))/2.;
        outputValues(15,i0, 2) = -(y*(1.+ y)*z);
        outputValues(15,i0, 3) = -((-1.+ 2.*x)*(-1.+ (z*z)))/2.;
        outputValues(15,i0, 4) = -((-1.+ 2.*x)*(1.+ 2.*y)*z)/2.;
        outputValues(15,i0, 5) = -((-1.+ 2.*x)*y*(1.+ y))/2.;
        outputValues(15,i0, 6) = 0.;
        outputValues(15,i0, 7) = -((-1.+ x)*x*z);
        outputValues(15,i0, 8) = -((-1.+ x)*x*(1.+ 2.*y))/2.;
        outputValues(15,i0, 9) = 0.; 
        
        outputValues(16,i0, 0) = 0.;
        outputValues(16,i0, 1) = -((-1.+ 2.*y)*z*(1.+ z))/2.;
        outputValues(16,i0, 2) = -((-1.+ y)*y*(1.+ 2.*z))/2.;
        outputValues(16,i0, 3) = -(x*z*(1.+ z));
        outputValues(16,i0, 4) = -(x*(-1.+ 2.*y)*(1.+ 2.*z))/2.;
        outputValues(16,i0, 5) = -(x*(-1.+ y)*y);
        outputValues(16,i0, 6) = 0.;
        outputValues(16,i0, 7) = -((-1.+ (x*x))*(1.+ 2.*z))/2.;    
        outputValues(16,i0, 8) = -((-1.+ (x*x))*(-1.+ 2.*y))/2.;
        outputValues(16,i0, 9) = 0.;
        
        outputValues(17,i0, 0) = 0.;
        outputValues(17,i0, 1) = -(y*z*(1.+ z));
        outputValues(17,i0, 2) = -((-1.+ (y*y))*(1.+ 2.*z))/2.;
        outputValues(17,i0, 3) = -((1.+ 2.*x)*z*(1.+ z))/2.;
        outputValues(17,i0, 4) = -((1.+ 2.*x)*y*(1.+ 2.*z))/2.;
        outputValues(17,i0, 5) = -((1.+ 2.*x)*(-1.+ (y*y)))/2.;
        outputValues(17,i0, 6) = 0.;
        outputValues(17,i0, 7) = -(x*(1.+ x)*(1.+ 2.*z))/2.;
        outputValues(17,i0, 8) = -(x*(1.+ x)*y);
        outputValues(17,i0, 9) = 0.;
        
        outputValues(18,i0, 0) = 0.;
        outputValues(18,i0, 1) = -((1.+ 2.*y)*z*(1.+ z))/2.;
        outputValues(18,i0, 2) = -(y*(1.+ y)*(1.+ 2.*z))/2.;
        outputValues(18,i0, 3) = -(x*z*(1.+ z));
        outputValues(18,i0, 4) = -(x*(1.+ 2.*y)*(1.+ 2.*z))/2.;
        outputValues(18,i0, 5) = -(x*y*(1.+ y));
        outputValues(18,i0, 6) = 0.;
        outputValues(18,i0, 7) = -((-1.+ (x*x))*(1.+ 2.*z))/2.;    
        outputValues(18,i0, 8) = -((-1.+ (x*x))*(1.+ 2.*y))/2.;
        outputValues(18,i0, 9) = 0.;
        
        outputValues(19,i0, 0) = 0.;
        outputValues(19,i0, 1) = -(y*z*(1.+ z));
        outputValues(19,i0, 2) = -((-1.+ (y*y))*(1.+ 2.*z))/2.;
        outputValues(19,i0, 3) = -((-1.+ 2.*x)*z*(1.+ z))/2.;
        outputValues(19,i0, 4) = -((-1.+ 2.*x)*y*(1.+ 2.*z))/2.;
        outputValues(19,i0, 5) = -((-1.+ 2.*x)*(-1.+ (y*y)))/2.;
        outputValues(19,i0, 6) = 0.;
        outputValues(19,i0, 7) = -((-1.+ x)*x*(1.+ 2.*z))/2.;
        outputValues(19,i0, 8) = -((-1.+ x)*x*y);
        outputValues(19,i0, 9) = 0.;
        
        outputValues(20,i0, 0) = 0.;
        outputValues(20,i0, 1) = -4*y*(-1.+ (z*z));
        outputValues(20,i0, 2) = -4*(-1.+ (y*y))*z;
        outputValues(20,i0, 3) = -4*x*(-1.+ (z*z));
        outputValues(20,i0, 4) = -8*x*y*z;
        outputValues(20,i0, 5) = -4*x*(-1.+ (y*y));
        outputValues(20,i0, 6) = 0.;
        outputValues(20,i0, 7) = -4*(-1.+ (x*x))*z;
        outputValues(20,i0, 8) = -4*(-1.+ (x*x))*y;
        outputValues(20,i0, 9) = 0.;    
        
        outputValues(21,i0, 0) = 0.;
        outputValues(21,i0, 1) = 2.*y*(-1.+ z)*z;
        outputValues(21,i0, 2) = (-1.+ (y*y))*(-1.+ 2.*z);
        outputValues(21,i0, 3) = 2.*x*(-1.+ z)*z;
        outputValues(21,i0, 4) = 2.*x*y*(-1.+ 2.*z);
        outputValues(21,i0, 5) = 2.*x*(-1.+ (y*y));
        outputValues(21,i0, 6) = 0.;
        outputValues(21,i0, 7) = (-1.+ (x*x))*(-1.+ 2.*z);
        outputValues(21,i0, 8) = 2.*(-1.+ (x*x))*y;
        outputValues(21,i0, 9) = 0.;    
        
        outputValues(22,i0, 0) = 0.;
        outputValues(22,i0, 1) = 2.*y*z*(1.+ z);
        outputValues(22,i0, 2) = (-1.+ (y*y))*(1.+ 2.*z);
        outputValues(22,i0, 3) = 2.*x*z*(1.+ z);
        outputValues(22,i0, 4) = 2.*x*y*(1.+ 2.*z);
        outputValues(22,i0, 5) = 2.*x*(-1.+ (y*y));
        outputValues(22,i0, 6) = 0.;
        outputValues(22,i0, 7) = (-1.+ (x*x))*(1.+ 2.*z);
        outputValues(22,i0, 8) = 2.*(-1.+ (x*x))*y;
        outputValues(22,i0, 9) = 0.;    
        
        outputValues(23,i0, 0) = 0.;
        outputValues(23,i0, 1) = 2.*y*(-1.+ (z*z));
        outputValues(23,i0, 2) = 2.*(-1.+ (y*y))*z;
        outputValues(23,i0, 3) = (-1.+ 2.*x)*(-1.+ (z*z));
        outputValues(23,i0, 4) = 2.*(-1.+ 2.*x)*y*z;
        outputValues(23,i0, 5) = (-1.+ 2.*x)*(-1.+ (y*y));
        outputValues(23,i0, 6) = 0.;
        outputValues(23,i0, 7) = 2.*(-1.+ x)*x*z;
        outputValues(23,i0, 8) = 2.*(-1.+ x)*x*y;
        outputValues(23,i0, 9) = 0.;
        
        outputValues(24,i0, 0) = 0.;
        outputValues(24,i0, 1) = 2.*y*(-1.+ (z*z));
        outputValues(24,i0, 2) = 2.*(-1.+ (y*y))*z;
        outputValues(24,i0, 3) = (1.+ 2.*x)*(-1.+ (z*z));
        outputValues(24,i0, 4) = 2.*(1.+ 2.*x)*y*z;
        outputValues(24,i0, 5) = (1.+ 2.*x)*(-1.+ (y*y));
        outputValues(24,i0, 6) = 0.;
        outputValues(24,i0, 7) = 2.*x*(1.+ x)*z;
        outputValues(24,i0, 8) = 2.*x*(1.+ x)*y;
        outputValues(24,i0, 9) = 0.;   
        
        outputValues(25,i0, 0) = 0.;
        outputValues(25,i0, 1) = (-1.+ 2.*y)*(-1.+ (z*z));
        outputValues(25,i0, 2) = 2.*(-1.+ y)*y*z;
        outputValues(25,i0, 3) = 2.*x*(-1.+ (z*z));
        outputValues(25,i0, 4) = 2.*x*(-1.+ 2.*y)*z;
        outputValues(25,i0, 5) = 2.*x*(-1.+ y)*y;
        outputValues(25,i0, 6) = 0.;
        outputValues(25,i0, 7) = 2.*(-1.+ (x*x))*z;
        outputValues(25,i0, 8) = (-1.+ (x*x))*(-1.+ 2.*y);
        outputValues(25,i0, 9) = 0.;  
        
        outputValues(26,i0, 0) = 0.;
        outputValues(26,i0, 1) = (1.+ 2.*y)*(-1.+ (z*z));    
        outputValues(26,i0, 2) = 2.*y*(1.+ y)*z;
        outputValues(26,i0, 3) = 2.*x*(-1.+ (z*z));
        outputValues(26,i0, 4) = 2.*x*(1.+ 2.*y)*z;
        outputValues(26,i0, 5) = 2.*x*y*(1.+ y);
        outputValues(26,i0, 6) = 0.;
        outputValues(26,i0, 7) = 2.*(-1.+ (x*x))*z;
        outputValues(26,i0, 8) = (-1.+ (x*x))*(1.+ 2.*y);
        outputValues(26,i0, 9) = 0.;
      }
      break;
      
    case OPERATOR_D4:
      {
        // Non-zero entries have Dk (derivative cardinality) indices {3,4,5,7,8,12}, all other entries are 0.
        // Intitialize array by zero and then fill only non-zero entries.
        int DkCardinality = Intrepid::getDkCardinality(operatorType, this -> basisCellTopology_.getDimension() );
        for(int dofOrd = 0; dofOrd < this -> basisCardinality_; dofOrd++) {
          for (int i0 = 0; i0 < dim0; i0++) {
            for(int dkOrd = 0; dkOrd < DkCardinality; dkOrd++){
              outputValues(dofOrd, i0, dkOrd) = 0.0;
            }
          }
        }
        
        for (int i0 = 0; i0 < dim0; i0++) {
          x = inputPoints(i0,0);
          y = inputPoints(i0,1);
          z = inputPoints(i0,2);
          
          outputValues(0, i0, 3) = ((-1.+ z)*z)/2.;
          outputValues(0, i0, 4) = ((-1.+ 2.*y)*(-1.+ 2.*z))/4.; 
          outputValues(0, i0, 5) = ((-1.+ y)*y)/2.;
          outputValues(0, i0, 7) = ((-1.+ 2.*x)*(-1.+ 2.*z))/4.;
          outputValues(0, i0, 8) = ((-1.+ 2.*x)*(-1.+ 2.*y))/4.;
          outputValues(0, i0, 12)= ((-1.+ x)*x)/2.;
          
          outputValues(1, i0, 3) = ((-1.+ z)*z)/2.;
          outputValues(1, i0, 4) = ((-1.+ 2.*y)*(-1.+ 2.*z))/4.;
          outputValues(1, i0, 5) = ((-1.+ y)*y)/2.;
          outputValues(1, i0, 7) = ((1. + 2.*x)*(-1.+ 2.*z))/4.;
          outputValues(1, i0, 8) = ((1. + 2.*x)*(-1.+ 2.*y))/4.;
          outputValues(1, i0, 12)= (x*(1. + x))/2.;
          
          outputValues(2, i0, 3) = ((-1.+ z)*z)/2.;
          outputValues(2, i0, 4) = ((1. + 2.*y)*(-1.+ 2.*z))/4.;
          outputValues(2, i0, 5) = (y*(1. + y))/2.;
          outputValues(2, i0, 7) = ((1. + 2.*x)*(-1.+ 2.*z))/4.;
          outputValues(2, i0, 8) =  ((1. + 2.*x)*(1. + 2.*y))/4.;
          outputValues(2, i0, 12)= (x*(1. + x))/2.;
          
          outputValues(3, i0, 3) = ((-1.+ z)*z)/2.;
          outputValues(3, i0, 4) = ((1. + 2.*y)*(-1.+ 2.*z))/4.;
          outputValues(3, i0, 5) = (y*(1. + y))/2.;
          outputValues(3, i0, 7) = ((-1.+ 2.*x)*(-1.+ 2.*z))/4.;
          outputValues(3, i0, 8) = ((-1.+ 2.*x)*(1. + 2.*y))/4.;
          outputValues(3, i0, 12)= ((-1.+ x)*x)/2.;
          
          outputValues(4, i0, 3) = (z*(1. + z))/2.;
          outputValues(4, i0, 4) = ((-1.+ 2.*y)*(1. + 2.*z))/4.;
          outputValues(4, i0, 5) = ((-1.+ y)*y)/2.;
          outputValues(4, i0, 7) = ((-1.+ 2.*x)*(1. + 2.*z))/4.;
          outputValues(4, i0, 8) = ((-1.+ 2.*x)*(-1.+ 2.*y))/4.;
          outputValues(4, i0, 12)= ((-1.+ x)*x)/2.;
          
          outputValues(5, i0, 3) = (z*(1. + z))/2.;
          outputValues(5, i0, 4) = ((-1.+ 2.*y)*(1. + 2.*z))/4.;
          outputValues(5, i0, 5) = ((-1.+ y)*y)/2.;
          outputValues(5, i0, 7) = ((1. + 2.*x)*(1. + 2.*z))/4.; 
          outputValues(5, i0, 8) = ((1. + 2.*x)*(-1.+ 2.*y))/4.;
          outputValues(5, i0, 12)= (x*(1. + x))/2.;
          
          outputValues(6, i0, 3) = (z*(1. + z))/2.;
          outputValues(6, i0, 4) = ((1. + 2.*y)*(1. + 2.*z))/4.;
          outputValues(6, i0, 5) = (y*(1. + y))/2.;
          outputValues(6, i0, 7) = ((1. + 2.*x)*(1. + 2.*z))/4.;
          outputValues(6, i0, 8) = ((1. + 2.*x)*(1. + 2.*y))/4.;
          outputValues(6, i0, 12)=  (x*(1. + x))/2.; 
          
          outputValues(7, i0, 3) = (z*(1. + z))/2.;
          outputValues(7, i0, 4) = ((1. + 2.*y)*(1. + 2.*z))/4.;
          outputValues(7, i0, 5) = (y*(1. + y))/2.;
          outputValues(7, i0, 7) = ((-1.+ 2.*x)*(1. + 2.*z))/4.;
          outputValues(7, i0, 8) = ((-1.+ 2.*x)*(1. + 2.*y))/4.; 
          outputValues(7, i0, 12)= ((-1.+ x)*x)/2.;
          
          outputValues(8, i0, 3) = -((-1.+ z)*z);
          outputValues(8, i0, 4) = -0.5 + y + z - 2.*y*z;
          outputValues(8, i0, 5) = -((-1.+ y)*y);
          outputValues(8, i0, 7) = x - 2.*x*z;  
          outputValues(8, i0, 8) = x - 2.*x*y;
          outputValues(8, i0, 12)= 1. - x*x;  
          
          outputValues(9, i0, 3) = -((-1.+ z)*z);
          outputValues(9, i0, 4) = y - 2.*y*z;   
          outputValues(9, i0, 5) = 1 - y*y;
          outputValues(9, i0, 7) = 0.5 + x - z - 2.*x*z;
          outputValues(9, i0, 8) = -((1. + 2.*x)*y);
          outputValues(9, i0, 12)= -(x*(1. + x)); 
          
          outputValues(10,i0, 3) = -((-1.+ z)*z);
          outputValues(10,i0, 4) = 0.5 + y - z - 2.*y*z; 
          outputValues(10,i0, 5) = -(y*(1. + y));
          outputValues(10,i0, 7) = x - 2.*x*z; 
          outputValues(10,i0, 8) = -(x*(1. + 2.*y)); 
          outputValues(10,i0, 12)=  1. - x*x;   
          
          outputValues(11,i0, 3) = -((-1.+ z)*z);
          outputValues(11,i0, 4) =  y - 2.*y*z;   
          outputValues(11,i0, 5) =  1. - y*y;   
          outputValues(11,i0, 7) = -0.5 + x + z - 2.*x*z;
          outputValues(11,i0, 8) =  y - 2.*x*y;
          outputValues(11,i0, 12)= -((-1.+ x)*x);
          
          outputValues(12,i0, 3) = 1. - z*z;
          outputValues(12,i0, 4) = z - 2.*y*z;
          outputValues(12,i0, 5) = -((-1.+ y)*y);
          outputValues(12,i0, 7) =  z - 2.*x*z;
          outputValues(12,i0, 8) = -0.5 + x + y - 2.*x*y;
          outputValues(12,i0, 12)= -((-1.+ x)*x); 
          
          outputValues(13,i0, 3) =  1. - z*z;  
          outputValues(13,i0, 4) = z - 2.*y*z;
          outputValues(13,i0, 5) = -((-1.+ y)*y);
          outputValues(13,i0, 7) =  -((1. + 2.*x)*z); 
          outputValues(13,i0, 8) = 0.5 + x - y - 2.*x*y; 
          outputValues(13,i0, 12)= -(x*(1. + x));
          
          outputValues(14,i0, 3) = 1. - z*z;
          outputValues(14,i0, 4) = -((1. + 2.*y)*z);
          outputValues(14,i0, 5) = -(y*(1. + y));
          outputValues(14,i0, 7) = -((1. + 2.*x)*z);
          outputValues(14,i0, 8) = -((1. + 2.*x)*(1. + 2.*y))/2.;
          outputValues(14,i0, 12)= -(x*(1. + x)); 
          
          outputValues(15,i0, 3) =  1. - z*z;
          outputValues(15,i0, 4) = -((1. + 2.*y)*z);
          outputValues(15,i0, 5) = -(y*(1. + y));
          outputValues(15,i0, 7) = z - 2.*x*z; 
          outputValues(15,i0, 8) = 0.5 + y - x*(1. + 2.*y); 
          outputValues(15,i0, 12)= -((-1.+ x)*x);
          
          outputValues(16,i0, 3) = -(z*(1. + z)); 
          outputValues(16,i0, 4) = 0.5 + z - y*(1. + 2.*z);
          outputValues(16,i0, 5) = -((-1.+ y)*y);
          outputValues(16,i0, 7) = -(x*(1. + 2.*z));
          outputValues(16,i0, 8) = x - 2.*x*y;
          outputValues(16,i0, 12)= 1. - x*x;  
          
          outputValues(17,i0, 3) = -(z*(1. + z));
          outputValues(17,i0, 4) = -(y*(1. + 2.*z));
          outputValues(17,i0, 5) = 1. - y*y;
          outputValues(17,i0, 7) = -((1. + 2.*x)*(1. + 2.*z))/2.;
          outputValues(17,i0, 8) = -((1. + 2.*x)*y); 
          outputValues(17,i0, 12)= -(x*(1. + x));
          
          outputValues(18,i0, 3) = -(z*(1. + z));
          outputValues(18,i0, 4) = -((1. + 2.*y)*(1. + 2.*z))/2.;
          outputValues(18,i0, 5) = -(y*(1. + y));
          outputValues(18,i0, 7) =  -(x*(1. + 2.*z)); 
          outputValues(18,i0, 8) =  -(x*(1. + 2.*y)); 
          outputValues(18,i0, 12)= 1. - x*x; 
          
          outputValues(19,i0, 3) = -(z*(1. + z));
          outputValues(19,i0, 4) = -(y*(1. + 2.*z));
          outputValues(19,i0, 5) = 1. - y*y; 
          outputValues(19,i0, 7) = 0.5 + z - x*(1. + 2.*z);
          outputValues(19,i0, 8) = y - 2.*x*y;
          outputValues(19,i0, 12)= -((-1.+ x)*x); 
          
          outputValues(20,i0, 3) = 4. - 4.*z*z;
          outputValues(20,i0, 4) = -8.*y*z;
          outputValues(20,i0, 5) = 4. - 4.*y*y; 
          outputValues(20,i0, 7) = -8.*x*z;
          outputValues(20,i0, 8) = -8.*x*y; 
          outputValues(20,i0, 12)= 4. - 4.*x*x; 
          
          outputValues(21,i0, 3) = 2.*(-1.+ z)*z;
          outputValues(21,i0, 4) = 2.*y*(-1.+ 2.*z);  
          outputValues(21,i0, 5) = 2.*(-1.+ y*y);
          outputValues(21,i0, 7) = 2.*x*(-1.+ 2.*z);  
          outputValues(21,i0, 8) = 4.*x*y;  
          outputValues(21,i0, 12)= 2.*(-1.+ x*x); 
          
          outputValues(22,i0, 3) = 2.*z*(1. + z);
          outputValues(22,i0, 4) = 2.*(y + 2.*y*z);
          outputValues(22,i0, 5) = 2.*(-1.+ y*y);
          outputValues(22,i0, 7) = 2.*(x + 2.*x*z); 
          outputValues(22,i0, 8) = 4.*x*y;
          outputValues(22,i0, 12)= 2.*(-1.+ x*x);
          
          outputValues(23,i0, 3) = 2.*(-1.+ z*z); 
          outputValues(23,i0, 4) = 4.*y*z; 
          outputValues(23,i0, 5) = 2.*(-1.+ y*y);
          outputValues(23,i0, 7) = 2.*(-1.+ 2.*x)*z;
          outputValues(23,i0, 8) = 2.*(-1.+ 2.*x)*y; 
          outputValues(23,i0, 12)= 2.*(-1.+ x)*x; 
          
          outputValues(24,i0, 3) = 2.*(-1.+ z*z);
          outputValues(24,i0, 4) = 4.*y*z; 
          outputValues(24,i0, 5) = 2.*(-1.+ y*y);
          outputValues(24,i0, 7) = 2.*(z + 2.*x*z); 
          outputValues(24,i0, 8) = 2.*(y + 2.*x*y);
          outputValues(24,i0, 12)= 2.*x*(1. + x);
          
          outputValues(25,i0, 3) =  2.*(-1.+ z*z); 
          outputValues(25,i0, 4) = 2.*(-1.+ 2.*y)*z;
          outputValues(25,i0, 5) = 2.*(-1.+ y)*y;
          outputValues(25,i0, 7) = 4.*x*z;  
          outputValues(25,i0, 8) = 2.*x*(-1.+ 2.*y);  
          outputValues(25,i0, 12)= 2.*(-1.+ x*x);
          
          outputValues(26,i0, 3) = 2.*(-1.+ z*z);
          outputValues(26,i0, 4) = 2.*(z + 2.*y*z);
          outputValues(26,i0, 5) = 2.*y*(1. + y);
          outputValues(26,i0, 7) =  4.*x*z;
          outputValues(26,i0, 8) = 2.*(x + 2.*x*y); 
          outputValues(26,i0, 12)= 2.*(-1.+ x*x);   
        }
      }
      break;
      
    case OPERATOR_D5:
    case OPERATOR_D6:
      TEST_FOR_EXCEPTION( true, std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_HEX_C2_FEM): operator not supported");
      break;
      
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
      {
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, DkCardinality)
        int DkCardinality = Intrepid::getDkCardinality(operatorType, 
                                                       this -> basisCellTopology_.getDimension() );
        for(int dofOrd = 0; dofOrd < this -> basisCardinality_; dofOrd++) {
          for (int i0 = 0; i0 < dim0; i0++) {
            for(int dkOrd = 0; dkOrd < DkCardinality; dkOrd++){
              outputValues(dofOrd, i0, dkOrd) = 0.0;
            }
          }
        }
      }
      break;
      
    default:
      TEST_FOR_EXCEPTION( !( Intrepid::isValidOperator(operatorType) ), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_HEX_C2_FEM): Invalid operator type");
  }
}



template<class Scalar, class ArrayScalar>
void Basis_HGRAD_HEX_C2_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                            const ArrayScalar &    inputPoints,
                                                            const ArrayScalar &    cellVertices,
                                                            const EOperator        operatorType) const {
  TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (Basis_HGRAD_HEX_C2_FEM): FEM Basis calling an FVD member function");
                                                            }

}// namespace Intrepid

#ifndef INTREPID_HGRAD_QUAD_C2_FEMDEF_HPP
#define INTREPID_HGRAD_QUAD_C2_FEMDEF_HPP
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

/** \file   Intrepid_HGRAD_QUAD_C2_FEMDef.hpp
    \brief  Definition file for bi-linear FEM basis functions for H(grad) functions on QUAD cells.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {

  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_QUAD_C2_FEM<Scalar, ArrayScalar>::Basis_HGRAD_QUAD_C2_FEM()
  {
    this -> basisCardinality_  = 9;
    this -> basisDegree_       = 2;    
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
    this -> basisType_         = BASIS_FEM_DEFAULT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;
  }
    
  
template<class Scalar, class ArrayScalar>
void Basis_HGRAD_QUAD_C2_FEM<Scalar, ArrayScalar>::initializeTags() {
  
  // Basis-dependent intializations
  int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
  int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
  int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
  int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

  // An array with local DoF tags assigned to basis functions, in the order of their local enumeration 
  int tags[]  = { 0, 0, 0, 1,
                  0, 1, 0, 1,
                  0, 2, 0, 1,
                  0, 3, 0, 1,
                  // edge midpoints
                  1, 0, 0, 1,
                  1, 1, 0, 1,
                  1, 2, 0, 1,
                  1, 3, 0, 1,
                  // quad center
                  2, 0, 0, 1};
  
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
void Basis_HGRAD_QUAD_C2_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
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
  
  // Temporaries: (x,y) coordinates of the evaluation point
  Scalar x = 0.0;                                    
  Scalar y = 0.0;                                    
  
  switch (operatorType) {
    
    case OPERATOR_VALUE:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0, 0);
        y = inputPoints(i0, 1);
        
        // outputValues is a rank-2 array with dimensions (basisCardinality_, dim0)
        outputValues(0, i0) = x*(x - 1.0)*y*(y - 1.0)/4.0;
        outputValues(1, i0) = x*(x + 1.0)*y*(y - 1.0)/4.0;
        outputValues(2, i0) = x*(x + 1.0)*y*(y + 1.0)/4.0;
        outputValues(3, i0) = x*(x - 1.0)*y*(y + 1.0)/4.0;
        // edge midpoints basis functions
        outputValues(4, i0) = (1.0 - x)*(1.0 + x)*y*(y - 1.0)/2.0;
        outputValues(5, i0) = x*(x + 1.0)*(1.0 - y)*(1.0 + y)/2.0;
        outputValues(6, i0) = (1.0 - x)*(1.0 + x)*y*(y + 1.0)/2.0;
        outputValues(7, i0) = x*(x - 1.0)*(1.0 - y)*(1.0 + y)/2.0;
        // quad bubble basis function
        outputValues(8, i0) = (1.0 - x)*(1.0 + x)*(1.0 - y)*(1.0 + y); 
      }
      break;
      
    case OPERATOR_GRAD:
    case OPERATOR_D1:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        outputValues(0, i0, 0) = (-0.25 + 0.5*x)*(-1. + y)*y;
        outputValues(0, i0, 1) = (-1.0 + x)*x*(-0.25 + 0.5*y);
        
        outputValues(1, i0, 0) = (0.25 + 0.5*x)*(-1. + y)*y;
        outputValues(1, i0, 1) = x*(1. + x)*(-0.25 + 0.5*y);
        
        outputValues(2, i0, 0) = (0.25 + 0.5*x)*y*(1. + y);
        outputValues(2, i0, 1) = x*(1. + x)*(0.25 + 0.5*y);
 
        outputValues(3, i0, 0) = (-0.25 + 0.5*x)*y*(1. + y);
        outputValues(3, i0, 1) = (-1. + x)*x*(0.25 + 0.5*y);

        outputValues(4, i0, 0) = x*(1.0 - y)*y;
        outputValues(4, i0, 1) = 0.5*(1.0 - x)*(1.0 + x)*(-1.0 + 2.0*y);
          
        outputValues(5, i0, 0) = 0.5*(1.0 - y)*(1.0 + y)*(1.0 + 2.0*x);
        outputValues(5, i0, 1) =-x*(1.0 + x)*y;
          
        outputValues(6, i0, 0) =-y*(1.0 + y)*x;
        outputValues(6, i0, 1) = 0.5*(1.0 - x)*(1.0 + x)*(1.0 + 2.0*y);
          
        outputValues(7, i0, 0) = 0.5*(1.0 - y)*(1.0+ y)*(-1.0 + 2.0*x);
        outputValues(7, i0, 1) = (1.0 - x)*x*y;
 
        outputValues(8, i0, 0) =-2.0*(1.0 - y)*(1.0 + y)*x;
        outputValues(8, i0, 1) =-2.0*(1.0 - x)*(1.0 + x)*y;          
      }
      break;
      
    case OPERATOR_CURL:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
        
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
        // CURL(u) = (u_y, -u_x), is rotated GRAD
        outputValues(0, i0, 1) =-(-0.25 + 0.5*x)*(-1. + y)*y;
        outputValues(0, i0, 0) = (-1.0 + x)*x*(-0.25 + 0.5*y);
        
        outputValues(1, i0, 1) =-(0.25 + 0.5*x)*(-1. + y)*y;
        outputValues(1, i0, 0) = x*(1. + x)*(-0.25 + 0.5*y);
        
        outputValues(2, i0, 1) =-(0.25 + 0.5*x)*y*(1. + y);
        outputValues(2, i0, 0) = x*(1. + x)*(0.25 + 0.5*y);
        
        outputValues(3, i0, 1) =-(-0.25 + 0.5*x)*y*(1. + y);
        outputValues(3, i0, 0) = (-1. + x)*x*(0.25 + 0.5*y);
        
        outputValues(4, i0, 1) =-x*(1.0 - y)*y;
        outputValues(4, i0, 0) = 0.5*(1.0 - x)*(1.0 + x)*(-1.0 + 2.0*y);
        
        outputValues(5, i0, 1) =-0.5*(1.0 - y)*(1.0 + y)*(1.0 + 2.0*x);
        outputValues(5, i0, 0) =-x*(1.0 + x)*y;
        
        outputValues(6, i0, 1) = y*(1.0 + y)*x;
        outputValues(6, i0, 0) = 0.5*(1.0 - x)*(1.0 + x)*(1.0 + 2.0*y);
        
        outputValues(7, i0, 1) =-0.5*(1.0 - y)*(1.0 + y)*(-1.0 + 2.0*x);
        outputValues(7, i0, 0) = (1.0 - x)*x*y;
        
        outputValues(8, i0, 1) = 2.0*(1.0 - y)*(1.0 + y)*x;
        outputValues(8, i0, 0) =-2.0*(1.0 - x)*(1.0 + x)*y;          
      }
      break;
      
    case OPERATOR_DIV:
      TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_QUAD_C2_FEM): DIV is invalid operator for rank-0 (scalar) functions in 2D");
      break;
      
    case OPERATOR_D2:
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);
       
        // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, D2Cardinality=3) 
        outputValues(0, i0, 0) = 0.5*(-1.0 + y)*y;
        outputValues(0, i0, 1) = 0.25 - 0.5*y + x*(-0.5 + 1.*y);
        outputValues(0, i0, 2) = 0.5*(-1.0 + x)*x;

        outputValues(1, i0, 0) = 0.5*(-1.0 + y)*y;
        outputValues(1, i0, 1) =-0.25 + 0.5*y + x*(-0.5 + 1.*y);
        outputValues(1, i0, 2) = 0.5*x*(1.0 + x);
        
        outputValues(2, i0, 0) = 0.5*y*(1.0 + y);
        outputValues(2, i0, 1) = 0.25 + 0.5*y + x*(0.5 + 1.*y);
        outputValues(2, i0, 2) = 0.5*x*(1.0 + x);
        
        outputValues(3, i0, 0) = 0.5*y*(1.0 + y);
        outputValues(3, i0, 1) =-0.25 - 0.5*y + x*(0.5 + 1.*y);
        outputValues(3, i0, 2) = 0.5*(-1.0 + x)*x;
        
        outputValues(4, i0, 0) = (1.0 - y)*y;
        outputValues(4, i0, 1) = x*(1. - 2.*y);
        outputValues(4, i0, 2) = (1.0 - x)*(1.0 + x);

        outputValues(5, i0, 0) = (1.0 - y)*(1.0 + y);
        outputValues(5, i0, 1) = x*(0. - 2.*y) - 1.*y;
        outputValues(5, i0, 2) =-x*(1.0 + x);

        outputValues(6, i0, 0) =-y*(1.0 + y);
        outputValues(6, i0, 1) = x*(-1. - 2.*y);
        outputValues(6, i0, 2) = (1.0 - x)*(1.0 + x);

        outputValues(7, i0, 0) = (1.0 - y)*(1.0 + y);
        outputValues(7, i0, 1) = x*(0. - 2.*y) + 1.*y;
        outputValues(7, i0, 2) = (1.0 - x)*x;

        outputValues(8, i0, 0) =-2.0 + 2.0*y*y;
        outputValues(8, i0, 1) = 4*x*y;
        outputValues(8, i0, 2) =-2.0 + 2.0*x*x;
        
      }
      break;
      
    case OPERATOR_D3:
      // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, D3Cardinality=4) 
      for (int i0 = 0; i0 < dim0; i0++) {
        x = inputPoints(i0,0);
        y = inputPoints(i0,1);

        outputValues(0, i0, 0) = 0.0;
        outputValues(0, i0, 1) =-0.5 + y;
        outputValues(0, i0, 2) =-0.5 + x;
        outputValues(0, i0, 3) = 0.0;

        outputValues(1, i0, 0) = 0.0;
        outputValues(1, i0, 1) =-0.5 + y;
        outputValues(1, i0, 2) = 0.5 + x;
        outputValues(1, i0, 3) = 0.0;

        outputValues(2, i0, 0) = 0.0;
        outputValues(2, i0, 1) = 0.5 + y;
        outputValues(2, i0, 2) = 0.5 + x;
        outputValues(2, i0, 3) = 0.0;

        outputValues(3, i0, 0) = 0.0;
        outputValues(3, i0, 1) = 0.5 + y;
        outputValues(3, i0, 2) =-0.5 + x;
        outputValues(3, i0, 3) = 0.0;

        outputValues(4, i0, 0) = 0.0;
        outputValues(4, i0, 1) = 1.0 - 2.0*y;
        outputValues(4, i0, 2) =-2.0*x;
        outputValues(4, i0, 3) = 0.0;

        outputValues(5, i0, 0) = 0.0;
        outputValues(5, i0, 1) =-2.0*y;
        outputValues(5, i0, 2) =-1.0 - 2.0*x;
        outputValues(5, i0, 3) = 0.0;

        outputValues(6, i0, 0) = 0.0;
        outputValues(6, i0, 1) =-1.0 - 2.0*y;
        outputValues(6, i0, 2) =-2.0*x;
        outputValues(6, i0, 3) = 0.0;

        outputValues(7, i0, 0) = 0.0;
        outputValues(7, i0, 1) =-2.0*y;
        outputValues(7, i0, 2) = 1.0 - 2.0*x;
        outputValues(7, i0, 3) = 0.0;        
        
        outputValues(8, i0, 0) = 0.0;
        outputValues(8, i0, 1) = 4.0*y;
        outputValues(8, i0, 2) = 4.0*x;
        outputValues(8, i0, 3) = 0.0;                
      }
      break;
    
    case OPERATOR_D4:
      // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, D4Cardinality=5) 
      for (int i0 = 0; i0 < dim0; i0++) {
        
        outputValues(0, i0, 0) = 0.0;
        outputValues(0, i0, 1) = 0.0;
        outputValues(0, i0, 2) = 1.0;
        outputValues(0, i0, 3) = 0.0;                
        outputValues(0, i0, 4) = 0.0;                

        outputValues(1, i0, 0) = 0.0;
        outputValues(1, i0, 1) = 0.0;
        outputValues(1, i0, 2) = 1.0;
        outputValues(1, i0, 3) = 0.0;                
        outputValues(1, i0, 4) = 0.0;                
        
        outputValues(2, i0, 0) = 0.0;
        outputValues(2, i0, 1) = 0.0;
        outputValues(2, i0, 2) = 1.0;
        outputValues(2, i0, 3) = 0.0;                
        outputValues(2, i0, 4) = 0.0;                
        
        outputValues(3, i0, 0) = 0.0;
        outputValues(3, i0, 1) = 0.0;
        outputValues(3, i0, 2) = 1.0;
        outputValues(3, i0, 3) = 0.0;                
        outputValues(3, i0, 4) = 0.0;                
        
        outputValues(4, i0, 0) = 0.0;
        outputValues(4, i0, 1) = 0.0;
        outputValues(4, i0, 2) =-2.0;
        outputValues(4, i0, 3) = 0.0;                
        outputValues(4, i0, 4) = 0.0;                

        outputValues(5, i0, 0) = 0.0;
        outputValues(5, i0, 1) = 0.0;
        outputValues(5, i0, 2) =-2.0;
        outputValues(5, i0, 3) = 0.0;                
        outputValues(5, i0, 4) = 0.0;                
        
        outputValues(6, i0, 0) = 0.0;
        outputValues(6, i0, 1) = 0.0;
        outputValues(6, i0, 2) =-2.0;
        outputValues(6, i0, 3) = 0.0;                
        outputValues(6, i0, 4) = 0.0;                
        
        outputValues(7, i0, 0) = 0.0;
        outputValues(7, i0, 1) = 0.0;
        outputValues(7, i0, 2) =-2.0;
        outputValues(7, i0, 3) = 0.0;                
        outputValues(7, i0, 4) = 0.0;                
        
        outputValues(8, i0, 0) = 0.0;
        outputValues(8, i0, 1) = 0.0;
        outputValues(8, i0, 2) = 4.0;
        outputValues(8, i0, 3) = 0.0;                
        outputValues(8, i0, 4) = 0.0;                
      }
      break;
      
    case OPERATOR_D5:
    case OPERATOR_D6:
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
                          ">>> ERROR (Basis_HGRAD_QUAD_C2_FEM): Invalid operator type");
  }
}


  
template<class Scalar, class ArrayScalar>
void Basis_HGRAD_QUAD_C2_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                             const ArrayScalar &    inputPoints,
                                                             const ArrayScalar &    cellVertices,
                                                             const EOperator        operatorType) const {
  TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (Basis_HGRAD_QUAD_C2_FEM): FEM Basis calling an FVD member function");
}



template<class Scalar, class ArrayScalar>
void Basis_HGRAD_QUAD_C2_FEM<Scalar, ArrayScalar>::getDofCoords(ArrayScalar & DofCoords) const {
#ifdef HAVE_INTREPID_DEBUG
  // Verify rank of output array.
  TEST_FOR_EXCEPTION( !(DofCoords.rank() == 2), std::invalid_argument,
                      ">>> ERROR: (Intrepid::Basis_HGRAD_QUAD_C2_FEM::getDofCoords) rank = 2 required for DofCoords array");
  // Verify 0th dimension of output array.
  TEST_FOR_EXCEPTION( !( DofCoords.dimension(0) == this -> basisCardinality_ ), std::invalid_argument,
                      ">>> ERROR: (Intrepid::Basis_HGRAD_QUAD_C2_FEM::getDofCoords) mismatch in number of DoF and 0th dimension of DofCoords array");
  // Verify 1st dimension of output array.
  TEST_FOR_EXCEPTION( !( DofCoords.dimension(1) == (int)(this -> basisCellTopology_.getDimension()) ), std::invalid_argument,
                      ">>> ERROR: (Intrepid::Basis_HGRAD_QUAD_C2_FEM::getDofCoords) incorrect reference cell (1st) dimension in DofCoords array");
#endif

  DofCoords(0,0) = -1.0;   DofCoords(0,1) = -1.0;
  DofCoords(1,0) =  1.0;   DofCoords(1,1) = -1.0;
  DofCoords(2,0) =  1.0;   DofCoords(2,1) =  1.0;
  DofCoords(3,0) = -1.0;   DofCoords(3,1) =  1.0;

  DofCoords(4,0) =  0.0;   DofCoords(4,1) = -1.0;
  DofCoords(5,0) =  1.0;   DofCoords(5,1) =  0.0;
  DofCoords(6,0) =  0.0;   DofCoords(6,1) =  1.0;
  DofCoords(7,0) = -1.0;   DofCoords(7,1) =  0.0;

  DofCoords(8,0) =  0.0;   DofCoords(8,1) =  0.0;

}

}// namespace Intrepid
#endif

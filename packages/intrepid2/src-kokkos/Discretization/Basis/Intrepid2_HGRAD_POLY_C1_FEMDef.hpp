/*
// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include "Intrepid2_FieldContainer.hpp"
#include "Sacado.hpp"


namespace Intrepid2{

  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_POLY_C1_FEM<Scalar, ArrayScalar>::Basis_HGRAD_POLY_C1_FEM(const shards::CellTopology& cellTopology){
    this -> basisCardinality_  = cellTopology.getNodeCount();
    this -> basisDegree_       = 1;
    this -> basisCellTopology_ = cellTopology;
    this -> basisType_         = BASIS_FEM_DEFAULT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;

    initializeTags();
    this->basisTagsAreSet_ = true;
  }

  
  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_POLY_C1_FEM<Scalar, ArrayScalar>::initializeTags(){
    // Basis-dependent initializations
    int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
    int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
    int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
    int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
  
    // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 
    
    int *tags = new int[tagSize * this->getCardinality()];
    for (int i=0;i<this->getCardinality();i++){
      tags[4*i] = 0;
      tags[4*i+1] = i;
      tags[4*i+2] = 0;
      tags[4*i+3] = 1;
    }
    
    
    // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
    Intrepid2::setOrdinalTagData(this -> tagToOrdinal_,
				this -> ordinalToTag_,
				tags,
				this -> basisCardinality_,
				tagSize,
				posScDim,
				posScOrd,
				posDfOrd);

    delete [] tags;
  }  
  
  // this is the FEM reference element version, this should not be called for polygonal basis
  // since polygonal basis is defined on physical element.
  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_POLY_C1_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar& outputValues,
							       const ArrayScalar& inputPoints,
							       const EOperator operatorType) const{
    TEUCHOS_TEST_FOR_EXCEPTION ( true, std::logic_error,
			 ">>>ERROR (Basis_HGRAD_POLY_C1_FEM): Polygonal basis calling FEM member function");
  }


  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_POLY_C1_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar& outputValues,
							       const ArrayScalar& inputPoints,
							       const ArrayScalar& cellVertices,
							       const EOperator operatorType) const{
    
    
    // implement wachspress basis functions
    switch (operatorType) {
    case OPERATOR_VALUE:
      {
	shapeFunctions<Scalar,ArrayScalar>(outputValues,inputPoints,cellVertices);
      }
      break;
    case OPERATOR_GRAD:
      {
	FieldContainer<Sacado::Fad::SFad<Scalar, 2> > dInput(inputPoints.dimension(0),inputPoints.dimension(1));
	for (int i=0;i<inputPoints.dimension(0);i++){
	  for (int j=0;j<2;j++){
	    dInput(i,j) = Sacado::Fad::SFad<Scalar, 2>( inputPoints(i,j));
	    dInput(i,j).diff(j,2);
	  }
	}
	FieldContainer<Sacado::Fad::SFad<Scalar> > cellVerticesFad(cellVertices.dimension(0),cellVertices.dimension(1));
	for (int i=0;i<cellVertices.dimension(0);i++){
	  for (int j=0;j<cellVertices.dimension(1);j++){
	    cellVerticesFad(i,j) = Sacado::Fad::SFad<Scalar,2>( cellVertices(i,j) );
	  }
	}
	
	FieldContainer<Sacado::Fad::SFad<Scalar,2> > dOutput(outputValues.dimension(0),outputValues.dimension(1));
	
	shapeFunctions<Sacado::Fad::SFad<Scalar,2>, FieldContainer<Sacado::Fad::SFad<Scalar,2> > >(dOutput,dInput,cellVerticesFad);
	
	for (int i=0;i<outputValues.dimension(0);i++){
	  for (int j=0;j<outputValues.dimension(1);j++){
	    for (int k=0;k<outputValues.dimension(2);k++){
	      outputValues(i,j,k) = dOutput(i,j).dx(k);
	    }
	  }
	}
      }
      break;
    case OPERATOR_D1:
    case OPERATOR_D2:
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
      {
	TEUCHOS_TEST_FOR_EXCEPTION ( true, std::invalid_argument, 
			     ">>> ERROR (Basis_HGRAD_POLY_C1_FEM): operator not implemented yet");
      }
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION( !( Intrepid2::isValidOperator(operatorType)), std::invalid_argument,
			  ">>> ERROR (Basis_HGRAD_POLY_C1_FEM): Invalid operator type");
      break;
    }

  }

  
  template<class Scalar, class ArrayScalar>
  template<class Scalar1, class ArrayScalar1>
  void Basis_HGRAD_POLY_C1_FEM<Scalar, ArrayScalar>::shapeFunctions(ArrayScalar1& outputValues,
								    const ArrayScalar1& inputPoints,
								    const ArrayScalar1& cellVertices)const{
    int numPoints = inputPoints.dimension(0);
    FieldContainer<Scalar1> weightFunctions( this->basisCardinality_,numPoints );
    evaluateWeightFunctions<Scalar1, ArrayScalar1>(weightFunctions,inputPoints,cellVertices);
    for (int pointIndex = 0;pointIndex<outputValues.dimension(1);pointIndex++){
      Scalar1 denominator=0;
      
      for (int k=0;k<weightFunctions.dimension(0);k++){
	denominator += weightFunctions(k,pointIndex);
      }
      for (int k=0;k<outputValues.dimension(0);k++){
	outputValues(k,pointIndex) = weightFunctions(k,pointIndex)/denominator;
      }
    }
  }
							    


  template<class Scalar, class ArrayScalar>
  template<class Scalar1, class ArrayScalar1>
  Scalar1 Basis_HGRAD_POLY_C1_FEM<Scalar, ArrayScalar>::computeArea(const ArrayScalar1& p1,
								    const ArrayScalar1& p2,
								    const ArrayScalar1& p3) const{
    Scalar1 area = 0.5*(p1(0)*(p2(1)-p3(1))
			-p2(0)*(p1(1)-p3(1))
			+p3(0)*(p1(1)-p2(1)));
    return area;
  }
  

  template<class Scalar, class ArrayScalar>
  template<class Scalar1, class ArrayScalar1>
  void Basis_HGRAD_POLY_C1_FEM<Scalar, ArrayScalar>::evaluateWeightFunctions(ArrayScalar1& outputValues,
									     const ArrayScalar1& inputValues,
									     const ArrayScalar1& cellVertices)const{
    
    
    int spaceDim = this->basisCellTopology_.getDimension();
    for (int k=0;k<outputValues.dimension(0);k++){
      
      // compute a_k for each weight function
      // need a_k = A(p_i,p_j,p_k) where nodes i and j are adjacent to node k
      int adjIndex1 = -1, adjIndex2 = -1;
      for (int i = 0;i < this->basisCellTopology_.getEdgeCount();i++){
	if ( this->basisCellTopology_.getNodeMap(1,i,0) == k )
	  adjIndex1 = this->basisCellTopology_.getNodeMap(1,i,1);
	else if ( this->basisCellTopology_.getNodeMap(1,i,1) == k )
	  adjIndex2 = this->basisCellTopology_.getNodeMap(1,i,0);
      }
      TEUCHOS_TEST_FOR_EXCEPTION( (adjIndex1 == -1 || adjIndex2 == -1), std::invalid_argument,
			  ">>> ERROR (Intrepid_HGRAD_POLY_C1_FEM): cannot find adjacent nodes when evaluating Wachspress weight function.");
      FieldContainer<Scalar1> p1(spaceDim);
      FieldContainer<Scalar1> p2(spaceDim);
      FieldContainer<Scalar1> p3(spaceDim);
      for (int i=0;i<spaceDim;i++){
	p1(i) = cellVertices(adjIndex1,i);
	p2(i) = cellVertices(k,i);
	p3(i) = cellVertices(adjIndex2,i);
      }
      Scalar1 a_k = computeArea<Scalar1, ArrayScalar1>(p1,p2,p3);
      // now compute prod_{ij!=k} r_ij
      for (int pointIndex = 0;pointIndex < inputValues.dimension(0);pointIndex++){
	Scalar1 product = a_k;
	for (int edgeIndex = 0;edgeIndex < this->basisCellTopology_.getEdgeCount();edgeIndex++){
	  int edgeNode1 = this->basisCellTopology_.getNodeMap(1,edgeIndex,0);
	  int edgeNode2 = this->basisCellTopology_.getNodeMap(1,edgeIndex,1);
	  if ( edgeNode1 != k && edgeNode2 != k ){
	    for (int i=0;i<spaceDim;i++){
	      p1(i) = inputValues(pointIndex,i);
	      p2(i) = cellVertices(edgeNode1,i);
	      p3(i) = cellVertices(edgeNode2,i);
	    }
	    product *= computeArea<Scalar1, ArrayScalar1>(p1,p2,p3);
	  }
	}
	outputValues(k,pointIndex) = product;
      }
    }
  }
} // namespace Intrepid2
	   
    

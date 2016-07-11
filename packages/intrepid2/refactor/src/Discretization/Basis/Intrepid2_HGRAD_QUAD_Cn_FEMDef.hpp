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

/** \file   Intrepid_HGRAD_QUAD_Cn_FEMDef.hpp
    \brief  Definition file for the Intrepid2::HGRAD_QUAD_Cn_FEM class.
    \author Created by R. Kirby.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_QUAD_CN_FEM_DEF_HPP__
#define __INTREPID2_HGRAD_QUAD_CN_FEM_DEF_HPP__

namespace Intrepid2 {

  // -------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------
  template<typename SpT, typename OT, typename PT>
  Basis_HGRAD_QUAD_Cn_FEM<SpT,OT,PT>::
  Basis_HGRAD_QUAD_Cn_FEM( const ordinal_type order,
                           const EPointType   pointType )
    : lineX_(order, pointType),
      lineY_(order, pointType) {

    this->basisCardinality_  = linX_.getCardinality()*linY_.getCardinality();
    this->basisDegree_       = order;
    this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
    this->basisType_         = BASIS_FEM_FIAT;
    this->basisCoordinates_  = COORDINATES_CARTESIAN;

    const auto card = this->basisCardinality_;

    // initialize tags
    {
      // Basis-dependent initializations
      const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
      
      // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 
      ordinal_type tags[(MaxOrder+1)*(MaxOrder+1)][4];
      
      // four vertices
      ordinal_type idx = 0;
      for (auto i=0;i<4;++i,++idx) {
        tags[idx][0] = 0; // vertex dof
        tags[idx][1] = i; // vertex id
        tags[idx][2] = 0; // local dof id
        tags[idx][3] = 1; // total number of dofs in this vertex
      }

      // four edges
      ordinal_type dofs[4] = {};
      dofs[0] = linX_.getCardinality();
      dofs[1] = linY_.getCardinality();
      dofs[2] = linX_.getCardinality();
      dofs[3] = linY_.getCardinality();
      for (auto i=0;i<4;++i) 
        for (auto j=0;j<dofs[i];++j,++idx) {
          tags[idx][0] = 1; // edge dof
          tags[idx][1] = i; // edge id
          tags[idx][2] = j; // local dof id
          tags[idx][3] = dofs[i]; // total number of dofs in this edge
        }
      
      // interior
      const auto intr = card - idx;
      for (auto i=0;i<intr;++i,++idx) {
        tags[idx][0] = 2; // face dof
        tags[idx][1] = 0; // edge id
        tags[idx][2] = i; // local dof id
        tags[idx][3] = intr; // total number of dofs in this interior
      }

      ordinal_type_array_1d_host tagView(&tag[0][0], card*4);

      // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
      // tags are constructed on host
      this->setOrdinalTagData(this->tagToOrdinal_,
                              this->ordinalToTag_,
                              tagView,
                              this->basisCardinality_,
                              tagSize,
                              posScDim,
                              posScOrd,
                              posDfOrd);
    }

    // dofCoords on host and create its mirror view to device
    Kokkos::DynRankView<PT,typename SpT::array_layout,Kokkos::HostSpace>
      dofCoords("dofCoordsHost", this->basisCardinality_,this->basisCellTopology_.getDimension());

    lineX_.impl_.getDofCoords(dofCoordsLinX);
    lineY_.impl_.getDofCoords(dofCoordsLinY);

    // four vertices
    begin = 0; end = 4;
    dofCoords(0,0) = -1.0;   dofCoords(0,1) = -1.0;
    dofCoords(1,0) =  1.0;   dofCoords(1,1) = -1.0;
    dofCoords(2,0) =  1.0;   dofCoords(2,1) =  1.0;
    dofCoords(3,0) = -1.0;   dofCoords(3,1) =  1.0;

    // four edges
    begin = end; end += lineX_.getCardinality() - 2;
    {
      auto dofCoordsEdge0 = Kokkos::subdynrankview(dofCoords, range_type(begin, end), 0);
      auto dofCoordsEdge1 = Kokkos::subdynrankview(dofCoords, range_type(begin, end), 1);
      Kokkos::deep_copy(dofCoordsEdge0, internalDofCoordsX);
      Kokkos::deep_copy(dofCoordsEdge1, -1.0);
    }
    begin = end; end += lineY_.getCardinality() - 2;
    {
      auto dofCoordsEdge0 = Kokkos::subdynrankview(dofCoords, range_type(begin, end), 0);
      auto dofCoordsEdge1 = Kokkos::subdynrankview(dofCoords, range_type(begin, end), 1);
      Kokkos::deep_copy(dofCoordsEdge0, 1.0);
      Kokkos::deep_copy(dofCoordsEdge1, internalDofCoordsY);
    }
    begin = end; end += lineX_.getCardinality() - 2;
    {
      auto dofCoordsEdge0 = Kokkos::subdynrankview(dofCoords, range_type(begin, end), 0);
      auto dofCoordsEdge1 = Kokkos::subdynrankview(dofCoords, range_type(begin, end), 1);
      Kokkos::deep_copy(dofCoordsEdge0, internalDofCoordsX);
      Kokkos::deep_copy(dofCoordsEdge1, 1.0);
    }
    begin = end; end += lineY_.getCardinality() - 2;
    {
      auto dofCoordsEdge0 = Kokkos::subdynrankview(dofCoords, range_type(begin, end), 0);
      auto dofCoordsEdge1 = Kokkos::subdynrankview(dofCoords, range_type(begin, end), 1);
      Kokkos::deep_copy(dofCoordsEdge0, -1.0);
      Kokkos::deep_copy(dofCoordsEdge1, internalDofCoordsY);
    }

    // interior
    begin = end; end += edgeDofX*edgeDofY;
    auto dofCoordsIntr = Kokkos::subdynrankview(dofCoords, range_type(begin, end), Kokkos::ALL());
    
    for (auto j=0;j<edgeDofY;++j)  {
      for (auto i=0;i<edgeDofX;++i,++idx)
        dofCoordsIntr(idx, 0) = internalDofCoordsX(i);
        dofCoordsIntr(idx, 1) = internalDofCoordsY(j);
      }
    
    this->dofCoords_ = Kokkos::create_mirror_view(typename SpT::memory_space(), dofCoords);
    Kokkos::deep_copy(this->dofCoords_, dofCoords);
  }

  template<typename SpT, typename OT, typename PT>
  template<typename outputValueValueType, class ...outputValueProperties,
           typename inputPointValueType,  class ...inputPointProperties>
  void
  Basis_HGRAD_QUAD_C1_FEM<SpT,OT,PT>::Internal::
  getValues( /**/  Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
             const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
             const EOperator operatorType ) const {
#ifdef HAVE_INTREPID2_DEBUG
    getValues_HGRAD_Args<Scalar, ArrayScalar>(outputValues,
					      inputPoints,
					      operatorType,
					      obj_->getBaseCellTopology(),
					      obj_->getCardinality() );
#endif

    ArrayScalar xInputPoints(inputPoints.dimension(0),1);
    ArrayScalar yInputPoints(inputPoints.dimension(0),1);

    switch (operatorType) {
    case OPERATOR_VALUE:
      {
	ArrayScalar xBasisValues(xBasis_.getCardinality(),xInputPoints.dimension(0));
	ArrayScalar yBasisValues(yBasis_.getCardinality(),yInputPoints.dimension(0));



	xBasis_.getValues(xBasisValues,xInputPoints,OPERATOR_VALUE);
	yBasis_.getValues(yBasisValues,yInputPoints,OPERATOR_VALUE);

	int bfcur = 0;
	for (int j=0;j<yBasis_.getCardinality();j++) {
	  for (int i=0;i<xBasis_.getCardinality();i++) {
	    for (int k=0;k<inputPoints.dimension(0);k++) {
	      outputValues(bfcur,k) = xBasisValues(i,k) * yBasisValues(j,k);
	    }
	    bfcur++;
	  }
	}
      }
      break;
    case OPERATOR_GRAD:
    case OPERATOR_D1:
      {
	ArrayScalar xBasisValues(xBasis_.getCardinality(),xInputPoints.dimension(0));
	ArrayScalar yBasisValues(yBasis_.getCardinality(),yInputPoints.dimension(0));
	ArrayScalar xBasisDerivs(xBasis_.getCardinality(),xInputPoints.dimension(0),1);
	ArrayScalar yBasisDerivs(yBasis_.getCardinality(),yInputPoints.dimension(0),1);

	xBasis_.getValues(xBasisValues,xInputPoints,OPERATOR_VALUE);
	yBasis_.getValues(yBasisValues,yInputPoints,OPERATOR_VALUE);
	xBasis_.getValues(xBasisDerivs,xInputPoints,OPERATOR_D1);
	yBasis_.getValues(yBasisDerivs,yInputPoints,OPERATOR_D1);	

	// there are two multiindices: I need the (1,0) and (0,1) derivatives
	int bfcur = 0;

	for (int j=0;j<yBasis_.getCardinality();j++) {
	  for (int i=0;i<xBasis_.getCardinality();i++) {
	    for (int k=0;k<inputPoints.dimension(0);k++) {
	      outputValues(bfcur,k,0) = xBasisDerivs(i,k,0) * yBasisValues(j,k);
	      outputValues(bfcur,k,1) = xBasisValues(i,k) * yBasisDerivs(j,k,0);
	    }
	    bfcur++;
	  }
	}
      }
      break;
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
	ArrayScalar xBasisValues(xBasis_.getCardinality(),xInputPoints.dimension(0));
	ArrayScalar yBasisValues(yBasis_.getCardinality(),yInputPoints.dimension(0));

	Teuchos::Array<int> partialMult;

	for (int d=0;d<getDkCardinality(operatorType,2);d++) {
	  getDkMultiplicities( partialMult , d , operatorType , 2 );
	  if (partialMult[0] == 0) {
	    xBasisValues.resize(xBasis_.getCardinality(),xInputPoints.dimension(0));
	    xBasis_.getValues( xBasisValues , xInputPoints, OPERATOR_VALUE );
	  }
	  else {
	    xBasisValues.resize(xBasis_.getCardinality(),xInputPoints.dimension(0),1);
	    EOperator xop = (EOperator) ( (int) OPERATOR_D1 + partialMult[0] - 1 );
	    xBasis_.getValues( xBasisValues , xInputPoints, xop );
	    xBasisValues.resize(xBasis_.getCardinality(),xInputPoints.dimension(0));
	  }
	  if (partialMult[1] == 0) {
	    yBasisValues.resize(yBasis_.getCardinality(),yInputPoints.dimension(0));
	    yBasis_.getValues( yBasisValues , yInputPoints, OPERATOR_VALUE );
	  }
	  else {
	    yBasisValues.resize(yBasis_.getCardinality(),yInputPoints.dimension(0),1);
	    EOperator yop = (EOperator) ( (int) OPERATOR_D1 + partialMult[1] - 1 );
	    yBasis_.getValues( yBasisValues , yInputPoints, yop );
	    yBasisValues.resize(yBasis_.getCardinality(),yInputPoints.dimension(0));
	  }


	  int bfcur = 0;
	  for (int j=0;j<yBasis_.getCardinality();j++) {
	    for (int i=0;i<xBasis_.getCardinality();i++) {
	      for (int k=0;k<inputPoints.dimension(0);k++) {
		outputValues(bfcur,k,d) = xBasisValues(i,k) * yBasisValues(j,k);
	      }
	      bfcur++;
	    }
	  }
	}
      }
      break;
    case OPERATOR_CURL:
      {
	ArrayScalar xBasisValues(xBasis_.getCardinality(),xInputPoints.dimension(0));
	ArrayScalar yBasisValues(yBasis_.getCardinality(),yInputPoints.dimension(0));
	ArrayScalar xBasisDerivs(xBasis_.getCardinality(),xInputPoints.dimension(0),1);
	ArrayScalar yBasisDerivs(yBasis_.getCardinality(),yInputPoints.dimension(0),1);

	xBasis_.getValues(xBasisValues,xInputPoints,OPERATOR_VALUE);
	yBasis_.getValues(yBasisValues,yInputPoints,OPERATOR_VALUE);
	xBasis_.getValues(xBasisDerivs,xInputPoints,OPERATOR_D1);
	yBasis_.getValues(yBasisDerivs,yInputPoints,OPERATOR_D1);	

	// there are two multiindices: I need the (1,0) and (0,1) derivatives
	int bfcur = 0;

	for (int j=0;j<yBasis_.getCardinality();j++) {
	  for (int i=0;i<xBasis_.getCardinality();i++) {
	    for (int k=0;k<inputPoints.dimension(0);k++) {
	      outputValues(bfcur,k,0) = xBasisValues(i,k) * yBasisDerivs(j,k,0);
	      outputValues(bfcur,k,1) = -xBasisDerivs(i,k,0) * yBasisValues(j,k);
	    }
	    bfcur++;
	  }
	}
      }
      break;      
    default:
        TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                            ">>> ERROR (Basis_HGRAD_QUAD_Cn_FEM): Operator type not implemented");
        break;
    }
  }

  template<class Scalar,class ArrayScalar>
  void Basis_HGRAD_QUAD_Cn_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
							       const ArrayScalar &    inputPoints,
							       const ArrayScalar &    cellVertices,
							       const EOperator        operatorType) const {
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
			">>> ERROR (Basis_HGRAD_QUAD_Cn_FEM): FEM Basis calling an FVD member function");
  }

  template<class Scalar,class ArrayScalar>
  void Basis_HGRAD_QUAD_Cn_FEM<Scalar, ArrayScalar>::getDofCoords( ArrayScalar & dofCoords ) const
  {
    int cur = 0;
    for (int j=0;j<ptsy_.dimension(0);j++)
      {
	for (int i=0;i<ptsx_.dimension(0);i++)
	  {
	    dofCoords(cur,0) = ptsx_(i,0);
	    dofCoords(cur,1) = ptsy_(j,0);
	    cur++;
	  }
      }
  }

  
}// namespace Intrepid2

#endif

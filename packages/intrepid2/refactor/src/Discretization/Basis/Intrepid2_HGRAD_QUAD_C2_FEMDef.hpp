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

/** \file   Intrepid_HGRAD_QUAD_C2_FEMDef.hpp
    \brief  Definition file for bi-linear FEM basis functions for H(grad) functions on QUAD cells.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_QUAD_C2_FEM_DEF_HPP__
#define __INTREPID2_HGRAD_QUAD_C2_FEM_DEF_HPP__

namespace Intrepid2 {

  // -------------------------------------------------------------------------------------

  template<typename SpT, typename OT, typename PT>
  template<EOperator opType>
  template<typename outputValueValueType, class ...outputValueProperties,
           typename inputPointValueType,  class ...inputPointProperties>
  KOKKOS_INLINE_FUNCTION
  void
  Basis_HGRAD_QUAD_C2_FEM<SpT,OT,PT>::Serial<opType>::
  getValues( /**/  Kokkos::DynRankView<outputValueValueType,outputValueProperties...> output,
             const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  input ) {
    switch (opType) {
    case OPERATOR_VALUE : {
      const auto x = input(0);
      const auto y = input(1);

      // output is a rank-2 array with dimensions (basisCardinality_, dim0)
      output(0) = x*(x - 1.0)*y*(y - 1.0)/4.0;
      output(1) = x*(x + 1.0)*y*(y - 1.0)/4.0;
      output(2) = x*(x + 1.0)*y*(y + 1.0)/4.0;
      output(3) = x*(x - 1.0)*y*(y + 1.0)/4.0;
      // edge midpoints basis functions
      output(4) = (1.0 - x)*(1.0 + x)*y*(y - 1.0)/2.0;
      output(5) = x*(x + 1.0)*(1.0 - y)*(1.0 + y)/2.0;
      output(6) = (1.0 - x)*(1.0 + x)*y*(y + 1.0)/2.0;
      output(7) = x*(x - 1.0)*(1.0 - y)*(1.0 + y)/2.0;
      // quad bubble basis function
      output(8) = (1.0 - x)*(1.0 + x)*(1.0 - y)*(1.0 + y); 
      break;
    }
    case OPERATOR_D1 :
    case OPERATOR_GRAD : {
      const auto x = input(0);
      const auto y = input(1);

      // output is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
      output(0, 0) = (-0.25 + 0.5*x)*(-1. + y)*y;
      output(0, 1) = (-1.0 + x)*x*(-0.25 + 0.5*y);
      
      output(1, 0) = (0.25 + 0.5*x)*(-1. + y)*y;
      output(1, 1) = x*(1. + x)*(-0.25 + 0.5*y);
      
      output(2, 0) = (0.25 + 0.5*x)*y*(1. + y);
      output(2, 1) = x*(1. + x)*(0.25 + 0.5*y);
 
      output(3, 0) = (-0.25 + 0.5*x)*y*(1. + y);
      output(3, 1) = (-1. + x)*x*(0.25 + 0.5*y);

      output(4, 0) = x*(1.0 - y)*y;
      output(4, 1) = 0.5*(1.0 - x)*(1.0 + x)*(-1.0 + 2.0*y);
        
      output(5, 0) = 0.5*(1.0 - y)*(1.0 + y)*(1.0 + 2.0*x);
      output(5, 1) =-x*(1.0 + x)*y;
        
      output(6, 0) =-y*(1.0 + y)*x;
      output(6, 1) = 0.5*(1.0 - x)*(1.0 + x)*(1.0 + 2.0*y);
        
      output(7, 0) = 0.5*(1.0 - y)*(1.0+ y)*(-1.0 + 2.0*x);
      output(7, 1) = (1.0 - x)*x*y;
 
      output(8, 0) =-2.0*(1.0 - y)*(1.0 + y)*x;
      output(8, 1) =-2.0*(1.0 - x)*(1.0 + x)*y;          
      break;
    }
    case OPERATOR_CURL : {
      const auto x = input(0);
      const auto y = input(1);

      // output is a rank-3 array with dimensions (basisCardinality_, dim0, spaceDim)
      // CURL(u) = (u_y, -u_x), is rotated GRAD
      output(0, 1) =-(-0.25 + 0.5*x)*(-1. + y)*y;
      output(0, 0) = (-1.0 + x)*x*(-0.25 + 0.5*y);
      
      output(1, 1) =-(0.25 + 0.5*x)*(-1. + y)*y;
      output(1, 0) = x*(1. + x)*(-0.25 + 0.5*y);
      
      output(2, 1) =-(0.25 + 0.5*x)*y*(1. + y);
      output(2, 0) = x*(1. + x)*(0.25 + 0.5*y);
      
      output(3, 1) =-(-0.25 + 0.5*x)*y*(1. + y);
      output(3, 0) = (-1. + x)*x*(0.25 + 0.5*y);
      
      output(4, 1) =-x*(1.0 - y)*y;
      output(4, 0) = 0.5*(1.0 - x)*(1.0 + x)*(-1.0 + 2.0*y);
      
      output(5, 1) =-0.5*(1.0 - y)*(1.0 + y)*(1.0 + 2.0*x);
      output(5, 0) =-x*(1.0 + x)*y;
      
      output(6, 1) = y*(1.0 + y)*x;
      output(6, 0) = 0.5*(1.0 - x)*(1.0 + x)*(1.0 + 2.0*y);
      
      output(7, 1) =-0.5*(1.0 - y)*(1.0 + y)*(-1.0 + 2.0*x);
      output(7, 0) = (1.0 - x)*x*y;
      
      output(8, 1) = 2.0*(1.0 - y)*(1.0 + y)*x;
      output(8, 0) =-2.0*(1.0 - x)*(1.0 + x)*y;          
      break;
    }
    case OPERATOR_D2 : {
      const auto x = input(0);
      const auto y = input(1);
      // output is a rank-3 array with dimensions (basisCardinality_, dim0, D2Cardinality=3) 
      output(0, 0) = 0.5*(-1.0 + y)*y;
      output(0, 1) = 0.25 - 0.5*y + x*(-0.5 + 1.*y);
      output(0, 2) = 0.5*(-1.0 + x)*x;

      output(1, 0) = 0.5*(-1.0 + y)*y;
      output(1, 1) =-0.25 + 0.5*y + x*(-0.5 + 1.*y);
      output(1, 2) = 0.5*x*(1.0 + x);
      
      output(2, 0) = 0.5*y*(1.0 + y);
      output(2, 1) = 0.25 + 0.5*y + x*(0.5 + 1.*y);
      output(2, 2) = 0.5*x*(1.0 + x);
      
      output(3, 0) = 0.5*y*(1.0 + y);
      output(3, 1) =-0.25 - 0.5*y + x*(0.5 + 1.*y);
      output(3, 2) = 0.5*(-1.0 + x)*x;
      
      output(4, 0) = (1.0 - y)*y;
      output(4, 1) = x*(1. - 2.*y);
      output(4, 2) = (1.0 - x)*(1.0 + x);

      output(5, 0) = (1.0 - y)*(1.0 + y);
      output(5, 1) = x*(0. - 2.*y) - 1.*y;
      output(5, 2) =-x*(1.0 + x);

      output(6, 0) =-y*(1.0 + y);
      output(6, 1) = x*(-1. - 2.*y);
      output(6, 2) = (1.0 - x)*(1.0 + x);

      output(7, 0) = (1.0 - y)*(1.0 + y);
      output(7, 1) = x*(0. - 2.*y) + 1.*y;
      output(7, 2) = (1.0 - x)*x;

      output(8, 0) =-2.0 + 2.0*y*y;
      output(8, 1) = 4*x*y;
      output(8, 2) =-2.0 + 2.0*x*x;
      break;
    }
    case OPERATOR_D3 : {
      const auto x = input(0);
      const auto y = input(1);
      output(0, 0) = 0.0;
      output(0, 1) =-0.5 + y;
      output(0, 2) =-0.5 + x;
      output(0, 3) = 0.0;

      output(1, 0) = 0.0;
      output(1, 1) =-0.5 + y;
      output(1, 2) = 0.5 + x;
      output(1, 3) = 0.0;

      output(2, 0) = 0.0;
      output(2, 1) = 0.5 + y;
      output(2, 2) = 0.5 + x;
      output(2, 3) = 0.0;

      output(3, 0) = 0.0;
      output(3, 1) = 0.5 + y;
      output(3, 2) =-0.5 + x;
      output(3, 3) = 0.0;

      output(4, 0) = 0.0;
      output(4, 1) = 1.0 - 2.0*y;
      output(4, 2) =-2.0*x;
      output(4, 3) = 0.0;

      output(5, 0) = 0.0;
      output(5, 1) =-2.0*y;
      output(5, 2) =-1.0 - 2.0*x;
      output(5, 3) = 0.0;

      output(6, 0) = 0.0;
      output(6, 1) =-1.0 - 2.0*y;
      output(6, 2) =-2.0*x;
      output(6, 3) = 0.0;

      output(7, 0) = 0.0;
      output(7, 1) =-2.0*y;
      output(7, 2) = 1.0 - 2.0*x;
      output(7, 3) = 0.0;        
      
      output(8, 0) = 0.0;
      output(8, 1) = 4.0*y;
      output(8, 2) = 4.0*x;
      output(8, 3) = 0.0;                
      break;
    }
    case OPERATOR_D4 : {
      output(0, 0) = 0.0;
      output(0, 1) = 0.0;
      output(0, 2) = 1.0;
      output(0, 3) = 0.0;                
      output(0, 4) = 0.0;                

      output(1, 0) = 0.0;
      output(1, 1) = 0.0;
      output(1, 2) = 1.0;
      output(1, 3) = 0.0;                
      output(1, 4) = 0.0;                
      
      output(2, 0) = 0.0;
      output(2, 1) = 0.0;
      output(2, 2) = 1.0;
      output(2, 3) = 0.0;                
      output(2, 4) = 0.0;                
      
      output(3, 0) = 0.0;
      output(3, 1) = 0.0;
      output(3, 2) = 1.0;
      output(3, 3) = 0.0;                
      output(3, 4) = 0.0;                
      
      output(4, 0) = 0.0;
      output(4, 1) = 0.0;
      output(4, 2) =-2.0;
      output(4, 3) = 0.0;                
      output(4, 4) = 0.0;                

      output(5, 0) = 0.0;
      output(5, 1) = 0.0;
      output(5, 2) =-2.0;
      output(5, 3) = 0.0;                
      output(5, 4) = 0.0;                
      
      output(6, 0) = 0.0;
      output(6, 1) = 0.0;
      output(6, 2) =-2.0;
      output(6, 3) = 0.0;                
      output(6, 4) = 0.0;                
      
      output(7, 0) = 0.0;
      output(7, 1) = 0.0;
      output(7, 2) =-2.0;
      output(7, 3) = 0.0;                
      output(7, 4) = 0.0;                
      
      output(8, 0) = 0.0;
      output(8, 1) = 0.0;
      output(8, 2) = 4.0;
      output(8, 3) = 0.0;                
      output(8, 4) = 0.0;                
      break;
    }
    case OPERATOR_MAX : {
      const auto jend = output.dimension(1);
      const auto iend = output.dimension(0);

      for (size_type j=0;j<jend;++j)
        for (size_type i=0;i<iend;++i)
          output(i, j) = 0.0;
      break;
    }
    default: {
      INTREPID2_TEST_FOR_ABORT( opType != OPERATOR_VALUE &&
                                opType != OPERATOR_GRAD &&
                                opType != OPERATOR_CURL &&
                                opType != OPERATOR_D1 &&
                                opType != OPERATOR_D2 &&
                                opType != OPERATOR_D3 &&
                                opType != OPERATOR_D4 &&
                                opType != OPERATOR_MAX,
                                ">>> ERROR: (Intrepid2::Basis_HGRAD_QUAD_C2_FEM::Serial::getValues) operator is not supported");

    }
    }
  }

  // -------------------------------------------------------------------------------------


  template<typename SpT, typename OT, typename PT>
  Basis_HGRAD_QUAD_C2_FEM<SpT,OT,PT>::
  Basis_HGRAD_QUAD_C2_FEM()
    : impl_(this)
  {
    this -> basisCardinality_  = 9;
    this -> basisDegree_       = 2;    
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
    this -> basisType_         = BASIS_FEM_DEFAULT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;

    {
      // Basis-dependent intializations
      const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
  
      // An array with local DoF tags assigned to basis functions, in the order of their local enumeration 
      ordinal_type tags[36]  = { 0, 0, 0, 1,
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
  
      //host view
      ordinal_type_array_1d_host tagView(&tags[0],36);
    
      // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
      this->setOrdinalTagData(this -> tagToOrdinal_,
                              this -> ordinalToTag_,
                              tagView,
                              this -> basisCardinality_,
                              tagSize,
                              posScDim,
                              posScOrd,
                              posDfOrd);
    }

    // dofCoords on host and create its mirror view to device
    Kokkos::DynRankView<PT,typename SpT::array_layout,Kokkos::HostSpace>
      dofCoords("dofCoordsHost", this->basisCardinality_,this->basisCellTopology_.getDimension());
    
    dofCoords(0,0) = -1.0;   dofCoords(0,1) = -1.0;
    dofCoords(1,0) =  1.0;   dofCoords(1,1) = -1.0;
    dofCoords(2,0) =  1.0;   dofCoords(2,1) =  1.0;
    dofCoords(3,0) = -1.0;   dofCoords(3,1) =  1.0;
  
    dofCoords(4,0) =  0.0;   dofCoords(4,1) = -1.0;
    dofCoords(5,0) =  1.0;   dofCoords(5,1) =  0.0;
    dofCoords(6,0) =  0.0;   dofCoords(6,1) =  1.0;
    dofCoords(7,0) = -1.0;   dofCoords(7,1) =  0.0;
  
    dofCoords(8,0) =  0.0;   dofCoords(8,1) =  0.0;

    this->dofCoords_ = Kokkos::create_mirror_view(typename SpT::memory_space(), dofCoords);
    Kokkos::deep_copy(this->dofCoords_, dofCoords);   
  }



  template<typename SpT, typename OT, typename PT>
  template<typename outputValueValueType, class ...outputValueProperties,
           typename inputPointValueType,  class ...inputPointProperties>
  void 
  Basis_HGRAD_QUAD_C2_FEM<SpT,OT,PT>::Internal::
  getValues( /**/  Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
             const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
             const EOperator operatorType ) const {
  
    // Verify arguments
#ifdef HAVE_INTREPID2_DEBUG
    Intrepid2::getValues_HGRAD_Args(outputValues,
                                    inputPoints,
                                    operatorType,
                                    obj_->getBaseCellTopology(),
                                    obj_->getCardinality() );
#endif
  
    typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
    typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
    typedef typename ExecSpace<typename inputPointViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

    // Number of evaluation points = dim 0 of inputPoints
    const auto loopSize = inputPoints.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
  
    switch (operatorType) {
    
    case OPERATOR_VALUE: {
      typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_VALUE> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
      break;
    } 
    case OPERATOR_GRAD:
    case OPERATOR_D1: {
      typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_GRAD> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
      break;
    } 
    case OPERATOR_CURL: {
      typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_CURL> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
      break;
    } 
    case OPERATOR_DIV: {
      INTREPID2_TEST_FOR_EXCEPTION( (operatorType == OPERATOR_DIV), std::invalid_argument,
                                    ">>> ERROR (Basis_HGRAD_QUAD_C2_FEM): DIV is invalid operator for rank-0 (scalar) functions in 2D");
      break;
    } 
    case OPERATOR_D2: {
      typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D2> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
      break;
    } 
    case OPERATOR_D3: {
      // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, D3Cardinality=4) 
      typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D3> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
      break;
    } 
    case OPERATOR_D4: {
      // outputValues is a rank-3 array with dimensions (basisCardinality_, dim0, D4Cardinality=5) 
      typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_D4> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
      break;
    } 
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10: {
      typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_MAX> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints) );
      break;
    }
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( !( Intrepid2::isValidOperator(operatorType) ), std::invalid_argument,
                                    ">>> ERROR (Basis_HGRAD_QUAD_C2_FEM): Invalid operator type");
    }
    }
  }

  /* 
     template<class Scalar, class ArrayScalar>
     void Basis_HGRAD_QUAD_C2_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
     const ArrayScalar &    inputPoints,
     const ArrayScalar &    cellVertices,
     const EOperator        operatorType) const {
     TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
     ">>> ERROR (Basis_HGRAD_QUAD_C2_FEM): FEM Basis calling an FVD member function");
     }
  */

  template<typename SpT, typename OT, typename PT>
  template<typename dofCoordValueType, class ...dofCoordProperties>
  void
  Basis_HGRAD_QUAD_C2_FEM<SpT,OT,PT>::Internal::
  getDofCoords( Kokkos::DynRankView<dofCoordValueType,dofCoordProperties...> dofCoords ) const {
#ifdef HAVE_INTREPID2_DEBUG
    // Verify rank of output array.
    INTREPID2_TEST_FOR_EXCEPTION( dofCoords.rank() != 2, std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_QUAD_C2_FEM::getDofCoords) rank = 2 required for dofCoords array");
    // Verify 0th dimension of output array.
    INTREPID2_TEST_FOR_EXCEPTION( dofCoords.dimension(0) != obj_->basisCardinality_, std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_QUAD_C2_FEM::getDofCoords) mismatch in number of dof and 0th dimension of dofCoords array");
    // Verify 1st dimension of output array.
    INTREPID2_TEST_FOR_EXCEPTION( dofCoords.dimension(1) != obj_->basisCellTopology_.getDimension(), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_QUAD_C2_FEM::getDofCoords) incorrect reference cell (1st) dimension in dofCoords array");
#endif
    Kokkos::deep_copy(dofCoords, obj_->dofCoords_);
  }

}// namespace Intrepid2
#endif

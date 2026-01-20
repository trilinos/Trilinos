// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HGRAD_TRI_Cn_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree n for H(grad) functions on TRI cells.
    \author Created by R. Kirby and P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
 */

#ifndef __INTREPID2_HGRAD_TRI_CN_FEM_DEF_HPP__
#define __INTREPID2_HGRAD_TRI_CN_FEM_DEF_HPP__

#include "Intrepid2_HGRAD_TRI_Cn_FEM_ORTH.hpp"

namespace Intrepid2 {

// -------------------------------------------------------------------------------------
namespace Impl {

template<EOperator OpType>
template<typename OutputViewType,
typename InputViewType,
typename WorkViewType,
typename VinvViewType>
KOKKOS_INLINE_FUNCTION
void
Basis_HGRAD_TRI_Cn_FEM::Serial<OpType>::
getValues(       OutputViewType output,
    const InputViewType  input,
    WorkViewType         work,
    const VinvViewType   vinv,
    const ordinal_type   order ) {

  constexpr ordinal_type spaceDim = 2;
  const ordinal_type
  card = vinv.extent(0),
  npts = input.extent(0);

  typedef typename Kokkos::DynRankView<typename InputViewType::value_type, typename WorkViewType::memory_space> ViewType;
  auto ptr = work.data();

  switch (OpType) {
  case OPERATOR_VALUE: {
    const ViewType phis = createMatchingUnmanagedView<ViewType>(input, ptr, card, npts);
    ViewType dummyView;

    Impl::Basis_HGRAD_TRI_Cn_FEM_ORTH::
    Serial<OpType>::getValues(phis, input, dummyView, order);

    for (ordinal_type i=0;i<card;++i)
      for (ordinal_type j=0;j<npts;++j) {
        output.access(i,j) = 0.0;
        for (ordinal_type k=0;k<card;++k)
          output.access(i,j) += vinv(k,i)*phis.access(k,j);
      }
    break;
  }
  case OPERATOR_GRAD:
  case OPERATOR_D1: {
    const ViewType phis = createMatchingUnmanagedView<ViewType>(input, ptr, card, npts, spaceDim);
    ptr += card*npts*spaceDim*get_dimension_scalar(input);
    const ViewType workView = createMatchingUnmanagedView<ViewType>(input, ptr, card, npts, spaceDim+1);
    Impl::Basis_HGRAD_TRI_Cn_FEM_ORTH::
    Serial<OpType>::getValues(phis, input, workView, order);

    for (ordinal_type i=0;i<card;++i)
      for (ordinal_type j=0;j<npts;++j)
        for (ordinal_type k=0;k<spaceDim;++k) {
          output.access(i,j,k) = 0.0;
          for (ordinal_type l=0;l<card;++l)
            output.access(i,j,k) += vinv(l,i)*phis.access(l,j,k);
        }
    break;
  }
  case OPERATOR_D2:
  case OPERATOR_D3:
  case OPERATOR_D4:
  case OPERATOR_D5:
  case OPERATOR_D6:
  case OPERATOR_D7:
  case OPERATOR_D8:
  case OPERATOR_D9:
  case OPERATOR_D10: {
    const ordinal_type dkcard = getDkCardinality<OpType,spaceDim>(); //(orDn + 1);
    const ViewType phis = createMatchingUnmanagedView<ViewType>(input, ptr, card, npts, dkcard);
    ViewType dummyView;

    Impl::Basis_HGRAD_TRI_Cn_FEM_ORTH::
    Serial<OpType>::getValues(phis, input, dummyView, order);

    for (ordinal_type i=0;i<card;++i)
      for (ordinal_type j=0;j<npts;++j)
        for (ordinal_type k=0;k<dkcard;++k) {
          output.access(i,j,k) = 0.0;
          for (ordinal_type l=0;l<card;++l)
            output.access(i,j,k) += vinv(l,i)*phis.access(l,j,k);
        }
    break;
  }
  case OPERATOR_CURL: { // only works in 2d. first component is -d/dy, second is d/dx
    const ViewType phis = createMatchingUnmanagedView<ViewType>(input, ptr, card, npts, spaceDim);
    ptr += card*npts*spaceDim*get_dimension_scalar(input);
    const ViewType workView = createMatchingUnmanagedView<ViewType>(input, ptr, card, npts, spaceDim+1);


    Impl::Basis_HGRAD_TRI_Cn_FEM_ORTH::
    Serial<OPERATOR_D1>::getValues(phis, input, workView, order);

    for (ordinal_type i=0;i<card;++i)
      for (ordinal_type j=0;j<npts;++j) {
        output.access(i,j,0) = 0.0;
        for (ordinal_type l=0;l<card;++l)
          output.access(i,j,0) += vinv(l,i)*phis.access(l,j,1);
        output.access(i,j,1) = 0.0;
        for (ordinal_type l=0;l<card;++l)
          output.access(i,j,1) -= vinv(l,i)*phis.access(l,j,0);
      }
    break;
  }
  default: {
    INTREPID2_TEST_FOR_ABORT( true,
        ">>> ERROR (Basis_HGRAD_TRI_Cn_FEM): Operator type not implemented");
  }
  }
}

template<typename DT, ordinal_type numPtsPerEval,
typename outputValueValueType, class ...outputValueProperties,
typename inputPointValueType,  class ...inputPointProperties,
typename vinvValueType,        class ...vinvProperties>
void
Basis_HGRAD_TRI_Cn_FEM::
getValues(
    const typename DT::execution_space& space,
          Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
    const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
    const Kokkos::DynRankView<vinvValueType,       vinvProperties...>        vinv,
    const ordinal_type order,
    const EOperator operatorType) {
  typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
  typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
  typedef          Kokkos::DynRankView<vinvValueType,       vinvProperties...>                vinvViewType;
  typedef typename ExecSpace<typename inputPointViewType::execution_space,typename DT::execution_space>::ExecSpaceType ExecSpaceType;

  // loopSize corresponds to cardinality
  const auto loopSizeTmp1 = (inputPoints.extent(0)/numPtsPerEval);
  const auto loopSizeTmp2 = (inputPoints.extent(0)%numPtsPerEval != 0);
  const auto loopSize = loopSizeTmp1 + loopSizeTmp2;
  Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(space, 0, loopSize);

  const ordinal_type cardinality = outputValues.extent(0);
  const ordinal_type spaceDim = 2;

  typedef typename DeduceDynRankView<inputPointViewType>::type workViewType;

  switch (operatorType) {
  case OPERATOR_VALUE: {
    workViewType work = createMatchingView<workViewType>(inputPoints, "Basis_HGRAD_TRI_Cn_FEM::getValues::work", cardinality, inputPoints.extent(0));
    typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
        OPERATOR_VALUE,numPtsPerEval> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinv, work, order) );
    break;
  }
  case OPERATOR_GRAD:
  case OPERATOR_D1: {
    workViewType work = createMatchingView<workViewType>(inputPoints, "Basis_HGRAD_TRI_Cn_FEM::getValues::work", cardinality*(2*spaceDim+1), inputPoints.extent(0));
    typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
        OPERATOR_D1,numPtsPerEval> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinv, work, order) );
    break;
  }
  case OPERATOR_CURL: {
    workViewType work = createMatchingView<workViewType>(inputPoints, "Basis_HGRAD_TRI_Cn_FEM::getValues::work", cardinality*(2*spaceDim+1), inputPoints.extent(0));
    typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
        OPERATOR_CURL,numPtsPerEval> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinv, work, order) );
    break;
  }
  case OPERATOR_D2: {
    typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
        OPERATOR_D2,numPtsPerEval> FunctorType;
    workViewType work = createMatchingView<workViewType>(inputPoints, "Basis_HGRAD_TRI_Cn_FEM::getValues::work", cardinality*outputValues.extent(2), inputPoints.extent(0));
    Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinv, work, order) );
    break;
  }
  default: {
    INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument,
        ">>> ERROR (Basis_HGRAD_TRI_Cn_FEM): Operator type not implemented" );
  }
  }
}
}

// -------------------------------------------------------------------------------------
template<typename DT, typename OT, typename PT>
Basis_HGRAD_TRI_Cn_FEM<DT,OT,PT>::
Basis_HGRAD_TRI_Cn_FEM( const ordinal_type order,
    const EPointType   pointType ) {
  constexpr ordinal_type spaceDim = 2;

  this->basisCardinality_     = Intrepid2::getPnCardinality<spaceDim>(order); // bigN
  this->basisDegree_          = order; // small n
  this->basisCellTopologyKey_ = shards::Triangle<3>::key;
  this->basisType_            = BASIS_FEM_LAGRANGIAN;
  this->basisCoordinates_     = COORDINATES_CARTESIAN;
  this->functionSpace_        = FUNCTION_SPACE_HGRAD;

  pointType_ = (pointType == POINTTYPE_DEFAULT) ? POINTTYPE_EQUISPACED : pointType;
  const ordinal_type card = this->basisCardinality_;

  // points are computed in the host and will be copied
  Kokkos::DynRankView<scalarType,typename DT::execution_space::array_layout,Kokkos::HostSpace>
  dofCoords("Hgrad::Tri::Cn::dofCoords", card, spaceDim);

  // construct lattice
  const shards::CellTopology cellTopo(shards::getCellTopologyData<shards::Triangle<3>>());
  const ordinal_type offset = 0;
  PointTools::getLattice( dofCoords,
      cellTopo,
      order, offset,
      this->pointType_ );

  this->dofCoords_ = Kokkos::create_mirror_view(typename DT::memory_space(), dofCoords);
  Kokkos::deep_copy(this->dofCoords_, dofCoords);

  // form Vandermonde matrix.  Actually, this is the transpose of the VDM,
  // so we transpose on copy below.
  const ordinal_type lwork = card*card;
  Kokkos::DynRankView<scalarType,Kokkos::LayoutLeft,Kokkos::HostSpace>
  vmat("Hgrad::Tri::Cn::vmat", card, card),
  work("Hgrad::Tri::Cn::work", lwork),
  ipiv("Hgrad::Tri::Cn::ipiv", card);

  Impl::Basis_HGRAD_TRI_Cn_FEM_ORTH::getValues<Kokkos::HostSpace::execution_space,Parameters::MaxNumPtsPerBasisEval>(typename Kokkos::HostSpace::execution_space{},
                                                                                                                     vmat,
                                                                                                                     dofCoords,
                                                                                                                     order,
                                                                                                                     OPERATOR_VALUE);

  ordinal_type info = 0;
  Teuchos::LAPACK<ordinal_type,scalarType> lapack;

  lapack.GETRF(card, card,
      vmat.data(), vmat.stride(1),
      (ordinal_type*)ipiv.data(),
      &info);

  INTREPID2_TEST_FOR_EXCEPTION( info != 0,
      std::runtime_error ,
      ">>> ERROR: (Intrepid2::Basis_HGRAD_TRI_Cn_FEM) lapack.GETRF returns nonzero info." );

  lapack.GETRI(card,
      vmat.data(), vmat.stride(1),
      (ordinal_type*)ipiv.data(),
      work.data(), lwork,
      &info);

  INTREPID2_TEST_FOR_EXCEPTION( info != 0,
      std::runtime_error ,
      ">>> ERROR: (Intrepid2::Basis_HGRAD_TRI_Cn_FEM) lapack.GETRI returns nonzero info." );

  // create host mirror
  Kokkos::DynRankView<scalarType,typename DT::execution_space::array_layout,Kokkos::HostSpace>
  vinv("Hgrad::Line::Cn::vinv", card, card);

  for (ordinal_type i=0;i<card;++i)
    for (ordinal_type j=0;j<card;++j)
      vinv(i,j) = vmat(j,i);

  this->vinv_ = Kokkos::create_mirror_view(typename DT::memory_space(), vinv);
  Kokkos::deep_copy(this->vinv_ , vinv);

  // initialize tags
  {
    // Basis-dependent initializations
    constexpr ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
    const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
    const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
    const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

    // Note: the only reason why equispaced can't support higher order than Parameters::MaxOrder appears to be the fact that the tags below get stored into a fixed-length array.
    // TODO: relax the maximum order requirement by setting up tags in a different container, perhaps directly into an OrdinalTypeArray1DHost (tagView, below).  (As of this writing (1/25/22), looks like other nodal bases do this in a similar way -- those should be fixed at the same time; maybe search for Parameters::MaxOrder.)
    INTREPID2_TEST_FOR_EXCEPTION( order > Parameters::MaxOrder, std::invalid_argument, "polynomial order exceeds the max supported by this class");
    constexpr ordinal_type maxCard = Intrepid2::getPnCardinality<spaceDim, Parameters::MaxOrder>();
    ordinal_type tags[maxCard][tagSize];

    const ordinal_type
    numEdgeDof = Intrepid2::getPnCardinality<1>(order-2),
    numElemDof = (order > 2 ? Intrepid2::getPnCardinality<2>(order-3) : 0);

    scalarType xi0, xi1, xi2;
    const scalarType eps = threshold();

    ordinal_type edgeId[3] = {}, elemId = 0;
    for (ordinal_type i=0;i<card;++i) {

      // compute barycentric coordinates
      const auto x = dofCoords(i,0);
      const auto y = dofCoords(i,1);
      xi0 = 1.0 - x - y;
      xi1= x;
      xi2= y;

      // vertex
      if ((1.0 - xi0) < eps) { // vert 0
        tags[i][0] = 0; // vertex dof
        tags[i][1] = 0; // vertex id
        tags[i][2] = 0; // local dof id
        tags[i][3] = 1; // total vert dof
      }
      else if ((1.0 - xi1) < eps) { // vert 1
        tags[i][0] = 0; // vertex dof
        tags[i][1] = 1; // vertex id
        tags[i][2] = 0; // local dof id
        tags[i][3] = 1; // total vert dof
      }
      else if ((1.0 - xi2) < eps) { // vert 2
        tags[i][0] = 0; // vertex dof
        tags[i][1] = 2; // vertex id
        tags[i][2] = 0; // local dof id
        tags[i][3] = 1; // total vert dof
      }
      else if (xi2 < eps) { // edge 0
        tags[i][0] = 1; // edge dof
        tags[i][1] = 0; // edge id
        tags[i][2] = edgeId[0]++; // local dof id
        tags[i][3] = numEdgeDof; // total vert dof
      }
      else if (xi0 < eps) { // edge 1
        tags[i][0] = 1; // edge dof
        tags[i][1] = 1; // edge id
        tags[i][2] = edgeId[1]++; // local dof id
        tags[i][3] = numEdgeDof; // total vert dof
      }
      else if (xi1 < eps) { // edge 2
        tags[i][0] = 1; // edge dof
        tags[i][1] = 2; // edge id
        tags[i][2] = edgeId[2]++; // local dof id
        tags[i][3] = numEdgeDof; // total vert dof
      }
      else { // elem
        tags[i][0] = 2; // intr dof
        tags[i][1] = 0; // intr id
        tags[i][2] = elemId++; // local dof id
        tags[i][3] = numElemDof; // total vert dof
      }
    }

    OrdinalTypeArray1DHost tagView(&tags[0][0], card*tagSize);

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
}

  template<typename DT, typename OT, typename PT>
  void 
  Basis_HGRAD_TRI_Cn_FEM<DT,OT,PT>::getScratchSpaceSize(       
                                    ordinal_type& perTeamSpaceSize,
                                    ordinal_type& perThreadSpaceSize,
                              const PointViewType inputPoints,
                              const EOperator operatorType) const {
    perTeamSpaceSize = 0;
    perThreadSpaceSize = getWorkSizePerPoint(operatorType)*get_dimension_scalar(inputPoints)*sizeof(typename BasisBase::scalarType);
  }

  template<typename DT, typename OT, typename PT>
  KOKKOS_INLINE_FUNCTION
  void 
  Basis_HGRAD_TRI_Cn_FEM<DT,OT,PT>::getValues(       
          OutputViewType outputValues,
      const PointViewType  inputPoints,
      const EOperator operatorType,
      const typename Kokkos::TeamPolicy<typename DT::execution_space>::member_type& team_member,
      const typename DT::execution_space::scratch_memory_space & scratchStorage, 
      const ordinal_type subcellDim,
      const ordinal_type subcellOrdinal) const {
      
      INTREPID2_TEST_FOR_ABORT( !((subcellDim == -1) && (subcellOrdinal == -1)),
        ">>> ERROR: (Intrepid2::Basis_HGRAD_TRI_Cn_FEM::getValues), The capability of selecting subsets of basis functions has not been implemented yet.");

      const int numPoints = inputPoints.extent(0);
      using ScalarType = typename ScalarTraits<typename PointViewType::value_type>::scalar_type;
      using WorkViewType = Kokkos::DynRankView< ScalarType,typename DT::execution_space::scratch_memory_space,Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
      constexpr ordinal_type spaceDim = 2;
      auto sizePerPoint = (operatorType==OPERATOR_VALUE) ? 
                          this->vinv_.extent(0)*get_dimension_scalar(inputPoints) : 
                          (2*spaceDim+1)*this->vinv_.extent(0)*get_dimension_scalar(inputPoints);
      WorkViewType workView(scratchStorage, sizePerPoint*team_member.team_size());
      using range_type = Kokkos::pair<ordinal_type,ordinal_type>;
      switch(operatorType) {
        case OPERATOR_VALUE:
          Kokkos::parallel_for (Kokkos::TeamThreadRange (team_member, numPoints), [=, &vinv_ = this->vinv_, basisDegree_ = this->basisDegree_] (ordinal_type& pt) {
            auto       output = Kokkos::subview( outputValues, Kokkos::ALL(), range_type  (pt,pt+1), Kokkos::ALL() );
            const auto input  = Kokkos::subview( inputPoints,                 range_type(pt, pt+1), Kokkos::ALL() );
            WorkViewType  work(workView.data() + sizePerPoint*team_member.team_rank(), sizePerPoint);
            Impl::Basis_HGRAD_TRI_Cn_FEM::Serial<OPERATOR_VALUE>::getValues( output, input, work, vinv_, basisDegree_);
          });
          break;
        case OPERATOR_GRAD:
          Kokkos::parallel_for (Kokkos::TeamThreadRange (team_member, numPoints), [=, &vinv_ = this->vinv_, basisDegree_ = this->basisDegree_] (ordinal_type& pt) {
            auto       output = Kokkos::subview( outputValues, Kokkos::ALL(), range_type(pt,pt+1), Kokkos::ALL() );
            const auto input  = Kokkos::subview( inputPoints,                 range_type(pt,pt+1), Kokkos::ALL() );
            WorkViewType  work(workView.data() + sizePerPoint*team_member.team_rank(), sizePerPoint);
            Impl::Basis_HGRAD_TRI_Cn_FEM::Serial<OPERATOR_GRAD>::getValues( output, input, work, vinv_, basisDegree_);
          });
          break;
          case OPERATOR_CURL:
          Kokkos::parallel_for (Kokkos::TeamThreadRange (team_member, numPoints), [=, &vinv_ = this->vinv_, basisDegree_ = this->basisDegree_] (ordinal_type& pt) {
            auto       output = Kokkos::subview( outputValues, Kokkos::ALL(), range_type(pt,pt+1), Kokkos::ALL() );
            const auto input  = Kokkos::subview( inputPoints,                 range_type(pt,pt+1), Kokkos::ALL() );
            WorkViewType  work(workView.data() + sizePerPoint*team_member.team_rank(), sizePerPoint);
            Impl::Basis_HGRAD_TRI_Cn_FEM::Serial<OPERATOR_CURL>::getValues( output, input, work, vinv_, basisDegree_);
          });
          break;
        default: {          
          INTREPID2_TEST_FOR_ABORT( true,
            ">>> ERROR (Basis_HGRAD_TRI_Cn_FEM): getValues not implemented for this operator");
          }
    }
  }

} // namespace Intrepid2

#endif

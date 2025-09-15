// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file   Intrepid2_OrientationToolsDefModifyBasis.hpp
    \brief  Definition file for the Intrepid2::OrientationTools class.
    \author Created by Kyungjoo Kim
*/
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_MODIFY_BASIS_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_MODIFY_BASIS_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

#include "Intrepid2_Orientation.hpp"

namespace Intrepid2 {

  template<typename DT>
  template<typename elemOrtValueType, class ...elemOrtProperties,
           typename elemNodeValueType, class ...elemNodeProperties>
  void
  OrientationTools<DT>::
  getOrientation(      Kokkos::DynRankView<elemOrtValueType,elemOrtProperties...> elemOrts,
                 const Kokkos::DynRankView<elemNodeValueType,elemNodeProperties...> elemNodes,
                 const shards::CellTopology cellTopo,
                 bool isSide) {
    // small meta data modification and it uses shards; let's do this on host
    auto elemOrtsHost = Kokkos::create_mirror_view(elemOrts);
    auto elemNodesHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), elemNodes);
    
    const ordinal_type numCells = elemNodes.extent(0);
    for (auto cell=0;cell<numCells;++cell) {
      const auto nodes = Kokkos::subview(elemNodesHost, cell, Kokkos::ALL());
      elemOrtsHost(cell) = Orientation::getOrientation(cellTopo, nodes, isSide);
    }

    Kokkos::deep_copy(elemOrts, elemOrtsHost);
  }

  template<typename ortViewType,
           typename OutputViewType,
           typename inputViewType,
           typename o2tViewType,
           typename t2oViewType,
           typename dataViewType>
  struct F_modifyBasisByOrientation {
    ortViewType orts;
    OutputViewType output;
    inputViewType input;
    o2tViewType ordinalToTag;
    t2oViewType tagToOrdinal;

    const dataViewType matData;
    const ordinal_type cellDim, numVerts, numEdges, numFaces, numPoints, dimBasis;
    const bool leftMultiply;
    // for simple left-multiplied basis value modification, numPoints is the dimension after the field dimension
    // for matrix value modification (C,F1,F2), numPoints is F2 when left multiplied, and F1 when right multiplied
    const bool transpose; // when true, multiply by the transpose of the matrix

    F_modifyBasisByOrientation(ortViewType orts_,
                               OutputViewType output_,
                               inputViewType input_,
                               o2tViewType ordinalToTag_,
                               t2oViewType tagToOrdinal_,
                               const dataViewType matData_,
                               const ordinal_type cellDim_,
                               const ordinal_type numVerts_,
                               const ordinal_type numEdges_,
                               const ordinal_type numFaces_,
                               const ordinal_type numPoints_,
                               const ordinal_type dimBasis_,
                               const bool leftMultiply_ = true,
                               const bool transpose_ = false)
    : orts(orts_),
      output(output_),
      input(input_),
      ordinalToTag(ordinalToTag_),
      tagToOrdinal(tagToOrdinal_),
      matData(matData_),
      cellDim(cellDim_),
      numVerts(numVerts_),
      numEdges(numEdges_),
      numFaces(numFaces_),
      numPoints(numPoints_),
      dimBasis(dimBasis_),
      leftMultiply(leftMultiply_),
      transpose(transpose_)
    {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const ordinal_type cell) const {
      typedef typename inputViewType::non_const_value_type input_value_type;

      auto out = Kokkos::subview(output, cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
      auto in  = (input.rank() == output.rank()) ?
                 Kokkos::subview(input,  cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL())
               : Kokkos::subview(input,        Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

      // edge transformation
      ordinal_type existEdgeDofs = 0;
      if (numEdges > 0) {
        ordinal_type ortEdges[12];
        orts(cell).getEdgeOrientation(ortEdges, numEdges);

        // apply coeff matrix
        for (ordinal_type edgeId=0;edgeId<numEdges;++edgeId) {
          const ordinal_type ordEdge = (1 < tagToOrdinal.extent(0) ? (static_cast<size_type>(edgeId) < tagToOrdinal.extent(1) ? tagToOrdinal(1, edgeId, 0) : -1) : -1);

          if (ordEdge != -1) {
            existEdgeDofs = 1;
            const ordinal_type ndofEdge = ordinalToTag(ordEdge, 3);
            const auto mat = Kokkos::subview(matData,
                                             edgeId, ortEdges[edgeId],
                                             Kokkos::ALL(), Kokkos::ALL());

            for (ordinal_type j=0;j<numPoints;++j)
              for (ordinal_type i=0;i<ndofEdge;++i) {
                const ordinal_type ii = tagToOrdinal(1, edgeId, i);

                for (ordinal_type k=0;k<dimBasis;++k) {
                  input_value_type temp = 0.0;
                  for (ordinal_type l=0;l<ndofEdge;++l) {
                    const ordinal_type ll = tagToOrdinal(1, edgeId, l);
                    auto & input_ = leftMultiply ? in.access(ll, j, k) : in.access(j, ll, k);
                    auto & mat_il = transpose ? mat(l,i) : mat(i,l);
                    temp += mat_il*input_;
                  }
                  auto & output_ = leftMultiply ? out.access(ii, j, k) : out.access(j, ii, k);
                  output_ = temp;
                }
              }
          }
        }
      }

      // face transformation
      if (numFaces > 0) {
        ordinal_type ortFaces[12];
        orts(cell).getFaceOrientation(ortFaces, numFaces);

        // apply coeff matrix
        for (ordinal_type faceId=0;faceId<numFaces;++faceId) {
          const ordinal_type ordFace = (2 < tagToOrdinal.extent(0) ? (static_cast<size_type>(faceId) < tagToOrdinal.extent(1) ? tagToOrdinal(2, faceId, 0) : -1) : -1);

          if (ordFace != -1) {
            const ordinal_type ndofFace = ordinalToTag(ordFace, 3);
            const auto mat = Kokkos::subview(matData,
                                             numEdges*existEdgeDofs+faceId, ortFaces[faceId],
                                             Kokkos::ALL(), Kokkos::ALL());
            
            for (ordinal_type j=0;j<numPoints;++j)
              for (ordinal_type i=0;i<ndofFace;++i) {
                const ordinal_type ii = tagToOrdinal(2, faceId, i);

                for (ordinal_type k=0;k<dimBasis;++k) {
                  input_value_type temp = 0.0;
                  for (ordinal_type l=0;l<ndofFace;++l) {
                    const ordinal_type ll = tagToOrdinal(2, faceId, l);
                    auto & input_ = leftMultiply ? in.access(ll, j, k) : in.access(j, ll, k);
                    auto & mat_il = transpose ? mat(l,i) : mat(i,l);
                    temp += mat_il*input_;
                  }
                  
                  auto & output_ = leftMultiply ? out.access(ii, j, k) : out.access(j, ii, k);
                  output_ = temp;
                }
              }
          }
        }
      }

      //side orientations
      ordinal_type faceOrt(0), edgeOrt(0);
      if(cellDim == 2) orts(cell).getFaceOrientation(&faceOrt, 1);
      if (faceOrt != 0) {
        const ordinal_type ordFace = (2 < tagToOrdinal.extent(0) ? (static_cast<size_type>(0) < tagToOrdinal.extent(1) ? tagToOrdinal(2, 0, 0) : -1) : -1);

        if (ordFace != -1) {
          const ordinal_type ndofFace = ordinalToTag(ordFace, 3);
          const auto mat = Kokkos::subview(matData,
                                           numEdges*existEdgeDofs, faceOrt,
                                           Kokkos::ALL(), Kokkos::ALL());

          for (ordinal_type j=0;j<numPoints;++j)
            for (ordinal_type i=0;i<ndofFace;++i) {
              const ordinal_type ii = tagToOrdinal(2, 0, i);

              for (ordinal_type k=0;k<dimBasis;++k) {
                input_value_type temp = 0.0;
                for (ordinal_type l=0;l<ndofFace;++l) {
                  const ordinal_type ll = tagToOrdinal(2, 0, l);
                  auto & input_ = leftMultiply ? in.access(ll, j, k) : in.access(j, ll, k);
                  auto & mat_il = transpose ? mat(l,i) : mat(i,l);
                  temp += mat_il*input_;
                }
                auto & output_ = leftMultiply ? out.access(ii, j, k) : out.access(j, ii, k);
                output_ = temp;
              }
            }
        }
      }

      if(cellDim == 1) orts(cell).getEdgeOrientation(&edgeOrt, 1);
      if (edgeOrt != 0) {
        const ordinal_type ordEdge = (1 < tagToOrdinal.extent(0) ? (static_cast<size_type>(0) < tagToOrdinal.extent(1) ? tagToOrdinal(1, 0, 0) : -1) : -1);

        if (ordEdge != -1) {
          const ordinal_type ndofEdge = ordinalToTag(ordEdge, 3);
          const auto mat = Kokkos::subview(matData,
                                           0, edgeOrt,
                                           Kokkos::ALL(), Kokkos::ALL());

          for (ordinal_type j=0;j<numPoints;++j)
            for (ordinal_type i=0;i<ndofEdge;++i) {
              const ordinal_type ii = tagToOrdinal(1, 0, i);

              for (ordinal_type k=0;k<dimBasis;++k) {
                input_value_type temp = 0.0;
                for (ordinal_type l=0;l<ndofEdge;++l) {
                  const ordinal_type ll = tagToOrdinal(1, 0, l);
                  auto & input_ = leftMultiply ? in.access(ll, j, k) : in.access(j, ll, k);
                  auto & mat_il = transpose ? mat(l,i) : mat(i,l);
                  temp += mat_il*input_;
                }
                auto & output_ = leftMultiply ? out.access(ii, j, k) : out.access(j, ii, k);
                output_ = temp;
              }
            }
        }
      }
    }
  };

  template<typename OrtViewType,
             typename OutputViewType,
             typename InputViewType,
             typename OperatorViewType>
    struct F_modifyBasisByOrientationOperator {
      OrtViewType orts;
      OutputViewType output;
      InputViewType input;

      const OperatorViewType edgeOperatorData, faceOperatorData;
      const ordinal_type numEdges, numFaces, numPoints, dimBasis;
      const bool leftMultiply;
      const ordinal_type transpose;
      // for simple left-multiplied basis value modification, numPoints is the dimension after the field dimension
      // for matrix value modification (C,F1,F2), numPoints is F2 when left multiplied, and F1 when right multiplied

      F_modifyBasisByOrientationOperator(OrtViewType orts_,
                                         OutputViewType output_,
                                         InputViewType input_,
                                         const OperatorViewType edgeOperatorData_,
                                         const OperatorViewType faceOperatorData_,
                                         const ordinal_type numEdges_,
                                         const ordinal_type numFaces_,
                                         const ordinal_type numPoints_,
                                         const ordinal_type dimBasis_,
                                         const bool leftMultiply_ = true,
                                         const bool transpose_ = false)
      : orts(orts_),
        output(output_),
        input(input_),
        edgeOperatorData(edgeOperatorData_),
        faceOperatorData(faceOperatorData_),
        numEdges(numEdges_),
        numFaces(numFaces_),
        numPoints(numPoints_),
        dimBasis(dimBasis_),
        leftMultiply(leftMultiply_),
        transpose(static_cast<ordinal_type>(transpose_))
      {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cell) const {
        using input_value_type =  typename InputViewType::non_const_value_type;

        auto out = Kokkos::subview(output, cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto in  = (input.rank() == output.rank()) ?
                   Kokkos::subview(input,  cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL())
                 : Kokkos::subview(input,        Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

        // edge transformation
        if (numEdges > 0)
        {
          ordinal_type ortEdges[12];
          orts(cell).getEdgeOrientation(ortEdges, numEdges);
          
          if ((numEdges != 1) || (ortEdges[0] != 0)) // numEdges == 1 is "side" orientation case; can skip if ort is 0.
          {
            // apply operators on each edge
            for (ordinal_type edgeId=0;edgeId<numEdges;++edgeId)
            {
              //            bool havePrinted = false;
              auto & edgeOp = edgeOperatorData(edgeId, ortEdges[edgeId], transpose);
              auto & rowIndices = edgeOp.rowIndices;
              auto & colIndices = edgeOp.packedColumnIndices;
              auto & weights    = edgeOp.packedWeights;
              
              const ordinal_type numRowIndices = rowIndices.extent_int(0);
              if (edgeOp.isWeightedPermutation)
              {
                for (ordinal_type pointOrdinal=0;pointOrdinal<numPoints;pointOrdinal++)
                {
                  for (ordinal_type rowOrdinal=0; rowOrdinal<numRowIndices; rowOrdinal++)
                  {
                    const ordinal_type & i = rowIndices(rowOrdinal);
                    const ordinal_type & j = colIndices(rowOrdinal);
                    const auto & mat_ij = weights(rowOrdinal);
                    
                    for (int d=0; d<dimBasis; d++)
                    {
                      auto & input_ = leftMultiply ? in.access(j, pointOrdinal, d) : in.access(pointOrdinal, j, d);
                      auto & output_ = leftMultiply ? out.access(i, pointOrdinal, d) : out.access(pointOrdinal, i, d);
                      output_ = mat_ij * input_;
                    }
                  }
                }
              }
              else
              {
                auto & rowOffsets = edgeOp.offsetsForRowOrdinal;
                for (ordinal_type pointOrdinal=0;pointOrdinal<numPoints;pointOrdinal++)
                {
                  for (ordinal_type rowOrdinal=0; rowOrdinal<numRowIndices; rowOrdinal++)
                  {
                    const ordinal_type & i = rowIndices[rowOrdinal];
                    const ordinal_type & thisRowOffset = rowOffsets(rowOrdinal);
                    const ordinal_type & nextRowOffset = rowOffsets(rowOrdinal+1);
                    
                    for (int d=0; d<dimBasis; d++)
                    {
                      input_value_type temp = 0.0;
                      for (ordinal_type rowOffset=thisRowOffset; rowOffset<nextRowOffset; rowOffset++)
                      {
                        const ordinal_type & j = colIndices(rowOffset);
                        const auto & mat_ij = weights(rowOffset);
                        
                        auto & input_ = leftMultiply ? in.access(j, pointOrdinal, d) : in.access(pointOrdinal, j, d);
                        temp += mat_ij * input_;
                      }
                      auto & output_ = leftMultiply ? out.access(i, pointOrdinal, d) : out.access(pointOrdinal, i, d);
                      output_ = temp;
                    }
                  }
                }
              }
            }
          }
        }
        
        // face transformation
        if (numFaces > 0)
        {
          ordinal_type ortFaces[12];
          orts(cell).getFaceOrientation(ortFaces, numFaces);

          if ((numFaces != 1) || (ortFaces[0] != 0)) // numFaces == 1 is "side" orientation case; can skip if ort is 0.
          {
            // apply operators on each face
            for (ordinal_type faceId=0;faceId<numFaces;++faceId)
            {
              auto & faceOp = faceOperatorData(faceId, ortFaces[faceId], transpose);
              auto & rowIndices = faceOp.rowIndices;
              auto & colIndices = faceOp.packedColumnIndices;
              auto & weights    = faceOp.packedWeights;
              
              ordinal_type numRowIndices = rowIndices.extent_int(0);
              if (faceOp.isWeightedPermutation)
              {
                for (ordinal_type pointOrdinal=0;pointOrdinal<numPoints;pointOrdinal++)
                {
                  for (ordinal_type rowOrdinal=0; rowOrdinal<numRowIndices; rowOrdinal++)
                  {
                    const ordinal_type & i = rowIndices[rowOrdinal];
                    const ordinal_type & j = colIndices(rowOrdinal);
                    const auto & mat_ij = weights(rowOrdinal);
                    
                    for (int d=0; d<dimBasis; d++)
                    {
                      auto & input_ = leftMultiply ? in.access(j, pointOrdinal, d) : in.access(pointOrdinal, j, d);
                      auto & output_ = leftMultiply ? out.access(i, pointOrdinal, d) : out.access(pointOrdinal, i, d);
                      output_ = mat_ij * input_;
                    }
                  }
                }
              }
              else
              {
                auto & rowOffsets = faceOp.offsetsForRowOrdinal;
                for (ordinal_type pointOrdinal=0;pointOrdinal<numPoints;pointOrdinal++)
                {
                  for (ordinal_type rowOrdinal=0; rowOrdinal<numRowIndices; rowOrdinal++)
                  {
                    const ordinal_type & i = rowIndices[rowOrdinal];
                    const ordinal_type & thisRowOffset = rowOffsets(rowOrdinal);
                    const ordinal_type & nextRowOffset = rowOffsets(rowOrdinal+1);
                    
                    for (int d=0; d<dimBasis; d++)
                    {
                      input_value_type temp = 0.0;
                      for (ordinal_type rowOffset=thisRowOffset; rowOffset<nextRowOffset; rowOffset++)
                      {
                        const ordinal_type & j = colIndices(rowOffset);
                        const auto & mat_ij = weights(rowOffset);
                        
                        auto & input_ = leftMultiply ? in.access(j, pointOrdinal, d) : in.access(pointOrdinal, j, d);
                        temp += mat_ij * input_;
                      }
                      auto & output_ = leftMultiply ? out.access(i, pointOrdinal, d) : out.access(pointOrdinal, i, d);
                      output_ = temp;
                    }
                  }
                }
              }
            }
          }
        }
      }
    };

  template<typename DT>
  template<typename outputValueType, class ...outputProperties,
           typename inputValueType,  class ...inputProperties,
           typename OrientationViewType,
           typename BasisType>
  void
  OrientationTools<DT>::
  modifyBasisByOrientation(      Kokkos::DynRankView<outputValueType,outputProperties...> output,
                           const Kokkos::DynRankView<inputValueType, inputProperties...>  input,
                           const OrientationViewType orts,
                           const BasisType* basis,
                           const bool transpose) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      if (input.rank() == output.rank())
      {
        for (size_type i=0;i<input.rank();++i)
          INTREPID2_TEST_FOR_EXCEPTION( input.extent(i) != output.extent(i), std::invalid_argument,
                                        ">>> ERROR (OrientationTools::modifyBasisByOrientation): Input and output dimensions do not match.");
      }
      else if (input.rank() == output.rank() - 1)
      {
        for (size_type i=0;i<input.rank();++i)
          INTREPID2_TEST_FOR_EXCEPTION( input.extent(i) != output.extent(i+1), std::invalid_argument,
                                       ">>> ERROR (OrientationTools::modifyBasisByOrientation): Input dimensions must match output dimensions exactly, or else match all but the first dimension (in the case that input does not have a 'cell' dimension).");
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument,
                                     ">>> ERROR (OrientationTools::modifyBasisByOrientation): input and output ranks must either match, or input rank must be one less than that of output.")
      }

      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(output.extent(1)) != basis->getCardinality(), std::invalid_argument,
                                    ">>> ERROR (OrientationTools::modifyBasisByOrientation): Field dimension of input/output does not match to basis cardinality.");
    }
#endif

    const ordinal_type NODE_DIM = 0;
    const ordinal_type EDGE_DIM = 1;
    const ordinal_type FACE_DIM = 2;
    
    const shards::CellTopology cellTopo = basis->getBaseCellTopology();
    const ordinal_type  cellDim = cellTopo.getDimension();
    
    //Initialize output with values from input
    if(input.rank() == output.rank())
      Kokkos::deep_copy(output, input);
    else
      RealSpaceTools<DT>::clone(output, input);
    
    if ((cellDim < 3) || basis->requireOrientation()) {
      const ordinal_type
      numCells  = output.extent(0),
      //numBasis  = output.extent(1),
      numPoints = output.extent(2),
      dimBasis  = output.extent(3); //returns 1 when output.rank() < 4;
      
      const auto op_tuple = createOperators(basis);
      
      ordinal_type numEdges(0), numFaces(0);
      if (basis->requireOrientation()) {
        numEdges = cellTopo.getEdgeCount()*ordinal_type(basis->getDofCount(EDGE_DIM, 0) > 0);
        numFaces = cellTopo.getFaceCount()*ordinal_type(basis->getDofCount(FACE_DIM, 0) > 0);
      }
      else
      {
        // "side" orientations
        if      (cellDim == 1) numEdges = 1;
        else if (cellDim == 2) numFaces = 1;
      }
      
      bool leftMultiply = true;
      
      const Kokkos::RangePolicy<typename DT::execution_space> policy(0, numCells);
      using FunctorType = F_modifyBasisByOrientationOperator
      <decltype(orts),
      decltype(output),decltype(input),
      decltype(std::get<0>(op_tuple))>;
      Kokkos::parallel_for
      (policy,
       FunctorType(orts,
                   output, input,
                   std::get<0>(op_tuple), std::get<1>(op_tuple), numEdges, numFaces,
                   numPoints, dimBasis, leftMultiply, transpose));
    }
  }

  template<typename DT>
  template<typename outputValueType, class ...outputProperties,
           typename inputValueType,  class ...inputProperties,
           typename OrientationViewType,
           typename BasisType>
  void
  OrientationTools<DT>::
  modifyBasisByOrientationTranspose(      Kokkos::DynRankView<outputValueType,outputProperties...> output,
                                    const Kokkos::DynRankView<inputValueType, inputProperties...>  input,
                                    const OrientationViewType orts,
                                    const BasisType* basis ) {
    bool transpose = true; 
    modifyBasisByOrientation(output, input, orts, basis, transpose);
  }

  template<typename DT>
  template<typename outputValueType, class ...outputProperties,
           typename inputValueType,  class ...inputProperties,
           typename OrientationViewType,
           typename BasisType>
  void
  OrientationTools<DT>::
  modifyBasisByOrientationInverse(      Kokkos::DynRankView<outputValueType,outputProperties...> output,
                                    const Kokkos::DynRankView<inputValueType, inputProperties...>  input,
                                    const OrientationViewType orts,
                                    const BasisType* basis,
                                    const bool transpose ) {
  #ifdef HAVE_INTREPID2_DEBUG
    {
      if (input.rank() == output.rank())
      {
        for (size_type i=0;i<input.rank();++i)
          INTREPID2_TEST_FOR_EXCEPTION( input.extent(i) != output.extent(i), std::invalid_argument,
                                        ">>> ERROR (OrientationTools::modifyBasisByOrientation): Input and output dimensions do not match.");
      }
      else if (input.rank() == output.rank() - 1)
      {
        for (size_type i=0;i<input.rank();++i)
          INTREPID2_TEST_FOR_EXCEPTION( input.extent(i) != output.extent(i+1), std::invalid_argument,
                                       ">>> ERROR (OrientationTools::modifyBasisByOrientation): Input dimensions must match output dimensions exactly, or else match all but the first dimension (in the case that input does not have a 'cell' dimension).");
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument,
                                     ">>> ERROR (OrientationTools::modifyBasisByOrientation): input and output ranks must either match, or input rank must be one less than that of output.")
      }

      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(output.extent(1)) != basis->getCardinality(), std::invalid_argument,
                                    ">>> ERROR (OrientationTools::modifyBasisByOrientation): Field dimension of input/output does not match to basis cardinality.");
    }
  #endif
    const ordinal_type EDGE_DIM = 1;
    const ordinal_type FACE_DIM = 2;
    
    const shards::CellTopology cellTopo = basis->getBaseCellTopology();
    const ordinal_type  cellDim = cellTopo.getDimension();

    //Initialize output with values from input
    if(input.rank() == output.rank())
      Kokkos::deep_copy(output, input);
    else
      RealSpaceTools<DT>::clone(output, input);

    if ((cellDim < 3) || basis->requireOrientation()) {
      auto ordinalToTag = Kokkos::create_mirror_view_and_copy(typename DT::memory_space(), basis->getAllDofTags());
      auto tagToOrdinal = Kokkos::create_mirror_view_and_copy(typename DT::memory_space(), basis->getAllDofOrdinal());

      const ordinal_type
        numCells  = output.extent(0),
        //numBasis  = output.extent(1),
        numPoints = output.extent(2),
        dimBasis  = output.extent(3); //returns 1 when output.rank() < 4;

      const CoeffMatrixDataViewType matData = createInvCoeffMatrix(basis);
      auto opData = createInvOperators(basis);
      auto edgeOpData = std::get<0>(opData);
      auto faceOpData = std::get<1>(opData);
      
      ordinal_type numEdges(0), numFaces(0);

      if (basis->requireOrientation()) {
        numEdges = cellTopo.getEdgeCount()*ordinal_type(basis->getDofCount(1, 0) > 0);
        numFaces = cellTopo.getFaceCount()*ordinal_type(basis->getDofCount(2, 0) > 0);
      }
      else
      {
        // "side" orientations
        if      (cellDim == 1) numEdges = 1;
        else if (cellDim == 2) numFaces = 1;
      }
      
      bool leftMultiply = true;
      
      const Kokkos::RangePolicy<typename DT::execution_space> policy(0, numCells);
      typedef F_modifyBasisByOrientationOperator
        <decltype(orts),
         decltype(output),decltype(input),
         decltype(edgeOpData)> FunctorType;
      Kokkos::parallel_for
        (policy,
         FunctorType(orts,
                     output, input,
                     edgeOpData, faceOpData, numEdges, numFaces,
                     numPoints, dimBasis, leftMultiply, transpose));
    }
  }

  template<typename DT>
  template<typename outputValueType, class ...outputProperties,
           typename inputValueType,  class ...inputProperties,
           typename OrientationViewType,
           typename BasisTypeLeft,
           typename BasisTypeRight>
  void
  OrientationTools<DT>::
  modifyMatrixByOrientation(Kokkos::DynRankView<outputValueType,outputProperties...> output,
                            const Kokkos::DynRankView<inputValueType, inputProperties...>  input,
                            const OrientationViewType orts,
                            const BasisTypeLeft* basisLeft,
                            const BasisTypeRight* basisRight)
  {
    const ordinal_type NODE_DIM = 0;
    const ordinal_type EDGE_DIM = 1;
    const ordinal_type FACE_DIM = 2;
    
    const shards::CellTopology cellTopo = basisLeft->getBaseCellTopology();
    const ordinal_type  cellDim = cellTopo.getDimension();
    
    const ordinal_type numCells       = output.extent(0);
    const ordinal_type numFieldsLeft  = basisLeft->getCardinality();
    const ordinal_type numFieldsRight = basisRight->getCardinality();
    
	// apply orientations on left
    decltype(output) outputLeft("temp view - output from left application", numCells, numFieldsLeft, numFieldsRight);
    
    //Initialize outputLeft with values from input
    if(input.rank() == output.rank())
      Kokkos::deep_copy(outputLeft, input);
    else
      RealSpaceTools<DT>::clone(outputLeft, input);
	
#ifdef HAVE_INTREPID2_DEBUG
    {
      if (input.rank() == output.rank())
      {
        for (size_type i=0;i<input.rank();++i)
          INTREPID2_TEST_FOR_EXCEPTION( input.extent(i) != output.extent(i), std::invalid_argument,
                                        ">>> ERROR (OrientationTools::modifyMatrixByOrientation): Input and output dimensions do not match.");
      }
      else if (input.rank() == output.rank() - 1)
      {
        for (size_type i=0;i<input.rank();++i)
          INTREPID2_TEST_FOR_EXCEPTION( input.extent(i) != output.extent(i+1), std::invalid_argument,
                                       ">>> ERROR (OrientationTools::modifyMatrixByOrientation): Input dimensions must match output dimensions exactly, or else match all but the first dimension (in the case that input does not have a 'cell' dimension).");
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument,
                                     ">>> ERROR (OrientationTools::modifyMatrixByOrientation): input and output ranks must either match, or input rank must be one less than that of output.")
      }

      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(output.extent(1)) != numFieldsLeft, std::invalid_argument,
                                    ">>> ERROR (OrientationTools::modifyMatrixByOrientation): First field dimension of input/output does not match left basis cardinality.");
      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(output.extent(2)) != numFieldsRight, std::invalid_argument,
                                    ">>> ERROR (OrientationTools::modifyMatrixByOrientation): Second field dimension of input/output does not match right basis cardinality.");
      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(output.extent(3)) != 1, std::invalid_argument,
                                    ">>> ERROR (OrientationTools::modifyMatrixByOrientation): Third dimension of output must be 1.");
      
    }
#endif
    
    if ((cellDim < 3) || basisLeft->requireOrientation()) {
      const bool leftMultiply = true;
      const bool transpose    = false;
      
      const ordinal_type
        numOtherFields = numFieldsRight,
        dimBasis       = output.extent(3); //returns 1 when output.rank() < 4;

      ordinal_type numEdges(0), numFaces(0);
      if (basisLeft->requireOrientation()) {
        numEdges = cellTopo.getEdgeCount()*ordinal_type(basisLeft->getDofCount(EDGE_DIM, 0) > 0);
        numFaces = cellTopo.getFaceCount()*ordinal_type(basisLeft->getDofCount(FACE_DIM, 0) > 0);
      }
      else
      {
        // "side" orientations
        if      (cellDim == 1) numEdges = 1;
        else if (cellDim == 2) numFaces = 1;
      }
      
      const auto op_tuple = createOperators(basisLeft);
      
      const Kokkos::RangePolicy<typename DT::execution_space> policy(0, numCells);
      using FunctorType = F_modifyBasisByOrientationOperator
      <decltype(orts),
      decltype(outputLeft),decltype(input),
      decltype(std::get<0>(op_tuple))>;
      Kokkos::parallel_for
      (policy,
       FunctorType(orts,
                   outputLeft, input,
                   std::get<0>(op_tuple), std::get<1>(op_tuple), numEdges, numFaces,
                   numOtherFields, dimBasis, leftMultiply, transpose));
    }
    
    // apply orientations on right
    //Initialize output with values from outputLeft
    Kokkos::deep_copy(output, outputLeft);
    if ((cellDim < 3) || basisRight->requireOrientation()) {
      const bool leftMultiply = false;
      const bool transpose = true; // NVR: added this September 2025; I think to date we have only ever run this with symmetric orientation operations, so that it has not mattered to date.
      auto ordinalToTag = Kokkos::create_mirror_view_and_copy(typename DT::memory_space(), basisRight->getAllDofTags());
      auto tagToOrdinal = Kokkos::create_mirror_view_and_copy(typename DT::memory_space(), basisRight->getAllDofOrdinal());

      const ordinal_type
        numOtherFields = numFieldsLeft,
        dimBasis       = output.extent(3); //returns 1 when output.rank() < 4;

      ordinal_type numEdges(0), numFaces(0);
      if (basisRight->requireOrientation()) {
        numEdges = cellTopo.getEdgeCount()*ordinal_type(basisRight->getDofCount(1, 0) > 0);
        numFaces = cellTopo.getFaceCount()*ordinal_type(basisRight->getDofCount(2, 0) > 0);
      }
      else
      {
        // "side" orientations
        if      (cellDim == 1) numEdges = 1;
        else if (cellDim == 2) numFaces = 1;
      }
      
      const auto op_tuple = createOperators(basisRight);
      
      const Kokkos::RangePolicy<typename DT::execution_space> policy(0, numCells);
      using FunctorType = F_modifyBasisByOrientationOperator
      <decltype(orts),
      decltype(output),decltype(input),
      decltype(std::get<0>(op_tuple))>;
      Kokkos::parallel_for
      (policy,
       FunctorType(orts,
                   output, outputLeft,
                   std::get<0>(op_tuple), std::get<1>(op_tuple), numEdges, numFaces,
                   numOtherFields, dimBasis, leftMultiply, transpose));
    }
  }

} // namespace Intrepid2

#endif

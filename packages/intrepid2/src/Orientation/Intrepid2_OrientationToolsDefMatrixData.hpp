// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file   Intrepid2_OrientationToolsDefMatrixData.hpp
    \brief  Definition file for matrix data in the Intrepid2::OrientationTools class.
    \author Created by Kyungjoo Kim
*/
#ifndef __INTREPID2_ORIENTATIONTOOLS_DEF_MATRIX_DATA_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_DEF_MATRIX_DATA_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {

  template<typename DT>
  template<typename BasisHostType>
  typename OrientationTools<DT>::CoeffMatrixDataViewType
  OrientationTools<DT>::createCoeffMatrixInternal(const BasisHostType* basis, const bool inverse) {
    const std::string name(basis->getName());
    CoeffMatrixDataViewType matData;

    const auto cellTopo = basis->getBaseCellTopology();
    const ordinal_type numEdges = cellTopo.getSubcellCount(1);
    const ordinal_type numFaces = cellTopo.getSubcellCount(2);
    ordinal_type matDim = 0, matDim1 = 0, matDim2 = 0, numOrts = 0, numSubCells;
    for(ordinal_type i=0; i<numEdges; ++i) {
      matDim1 = std::max(matDim1, basis->getDofCount(1,i));
      numOrts = std::max(numOrts,2);
    }
    for(ordinal_type i=0; i<numFaces; ++i) {
      matDim2 = std::max(matDim2, basis->getDofCount(2,i));
      numOrts = std::max(numOrts,2*ordinal_type(cellTopo.getSideCount(2,i)));
    }
    matDim = std::max(matDim1,matDim2);
    numSubCells = (matDim1>0)*numEdges + (matDim2>0)*numFaces;


    matData = CoeffMatrixDataViewType("Orientation::CoeffMatrix::"+name,
                                      numSubCells,
                                      numOrts,
                                      matDim,
                                      matDim);

    if(basis->getFunctionSpace() == FUNCTION_SPACE_HGRAD) {
      init_HGRAD(matData, basis, inverse);
    } else if (basis->getFunctionSpace() == FUNCTION_SPACE_HCURL) {
      init_HCURL(matData, basis, inverse);
    } else if (basis->getFunctionSpace() == FUNCTION_SPACE_HDIV) {
      init_HDIV(matData, basis, inverse);
    } else if (basis->getFunctionSpace() == FUNCTION_SPACE_HVOL) {
      init_HVOL(matData, basis, inverse);
    }
    return matData;
  }

  template<typename DT>
  template<typename BasisHostType>
  typename OrientationTools<DT>::OperatorViewType
  OrientationTools<DT>::createEdgeOperatorsInternal(const BasisHostType* basis, CoeffMatrixDataViewType matData)
  {
    const std::string name(basis->getName());
    OperatorViewType operators;

    auto matDataHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), matData);

    ordinal_type matDim = 0, matDim1 = 0, matDim2 = 0, numOrts = 0, numEdgeOrts = 0, numSubCells;
    const auto cellTopo = basis->getBaseCellTopology();
    {
      const ordinal_type numEdges = cellTopo.getSubcellCount(1);
      const ordinal_type numFaces = cellTopo.getSubcellCount(2);
      for(ordinal_type i=0; i<numEdges; ++i) {
        matDim1 = std::max(matDim1, basis->getDofCount(1,i));
        numEdgeOrts = std::max(numOrts,2);
        numOrts     = numEdgeOrts;
      }
      for(ordinal_type i=0; i<numFaces; ++i) {
        matDim2 = std::max(matDim2, basis->getDofCount(2,i));
        numOrts = std::max(numOrts,2*ordinal_type(cellTopo.getSideCount(2,i)));
      }
      matDim = std::max(matDim1,matDim2);
      numSubCells = (matDim1>0)*numEdges + (matDim2>0)*numFaces;
      
      operators = OperatorViewType("Orientation::EdgeOperators::"+name,
                                   numEdges, numOrts);
    }
    
    // NOTE: the OrientationOperators within operatorsHost contain *device* views.
    auto operatorsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), operators);

    if (matDim1 == 0)
    {
      // no non-trivial edge operators (no edge dofs)
      return operators;
    }
    
    {
      auto ordinalToTag = basis->getAllDofTags();
      auto tagToOrdinal = basis->getAllDofOrdinal();
      
      ordinal_type numVerts(0), numEdges(0), numFaces(0);
      
      if (basis->requireOrientation()) {
        numVerts = cellTopo.getVertexCount()*ordinal_type(basis->getDofCount(0, 0) > 0);
        numEdges = cellTopo.getEdgeCount()*ordinal_type(basis->getDofCount(1, 0) > 0);
        numFaces = cellTopo.getFaceCount()*ordinal_type(basis->getDofCount(2, 0) > 0);
      }
      
      ordinal_type rowOffset = 0;
      if (numEdges > 0)
      {
        for (ordinal_type edgeId=0;edgeId<numEdges;++edgeId)
        {
          for (ordinal_type edgeOrt=0; edgeOrt<numEdgeOrts; edgeOrt++)
          {
            std::vector<ordinal_type> nonIdentityDofs;
            std::vector<ordinal_type> rowOffsets; // within the column storage
            std::vector<ordinal_type> colIDs;
            std::vector<double> weights;
            
            const ordinal_type ordEdge = (1 < tagToOrdinal.extent(0) ? (static_cast<size_type>(edgeId) < tagToOrdinal.extent(1) ? tagToOrdinal(1, edgeId, 0) : -1) : -1);
            
            if (ordEdge != -1) {
              const ordinal_type ndofEdge = ordinalToTag(ordEdge, 3);
              const auto mat = Kokkos::subview(matData,
                                               edgeId, edgeOrt,
                                               Kokkos::ALL(), Kokkos::ALL());
              
              for (ordinal_type i=0;i<ndofEdge;++i) {
                const ordinal_type ii = tagToOrdinal(1, edgeId, i);
                
                // first pass for ii:
                // check whether this is different from the identity
                // count number of nonzeros
                int nnz = 0;
                bool deviatesFromIdentity = false;
                for (ordinal_type l=0;l<ndofEdge;++l) {
                  const ordinal_type ll = tagToOrdinal(1, edgeId, l);
                  //                      auto & input_ = in.access(ll);
                  auto & mat_il = mat(i,l);
                  if (mat_il != 0.0)
                  {
                    nnz++;
                    if ((mat_il != 1.0) || (ii != ll))
                    {
                      deviatesFromIdentity = true;
                    }
                  }
                }
                INTREPID2_TEST_FOR_EXCEPTION(nnz != 0, std::invalid_argument, "Each dof should have *some* nonzero weight");
                if (deviatesFromIdentity)
                {
                  // then we store the nonzeros for ii
                  nonIdentityDofs.push_back(ii);
                  rowOffsets.push_back(rowOffset);
                  rowOffset += nnz;
                  
                  for (ordinal_type l=0;l<ndofEdge;++l)
                  {
                    const ordinal_type ll = tagToOrdinal(1, edgeId, l);
                    //                      auto & input_ = in.access(ll);
                    auto & mat_il = mat(i,l);
                    if (mat_il != 0.0)
                    {
                      colIDs.push_back(ll);
                      weights.push_back(mat_il);
                    }
                  }
                }
              }
            }
            // initialize operator
            const int numRows = static_cast<int>(nonIdentityDofs.size());
            Kokkos::View<ordinal_type*,DT> rowIndices("OrientationOperator: rowIndices", numRows);
            Kokkos::View<ordinal_type*,DT> offsetForRowOrdinal("OrientationOperator: rowOffsets", numRows+1); // convenient to be able to check the "next" row offset to get a column count, even when on the last row
            Kokkos::View<ordinal_type*,DT> packedColumnIndices("OrientationOperator: packedColumnIndices", rowOffset);
            Kokkos::View<double*,      DT> packedWeights("OrientationOperator: packedWeights", rowOffset);
            
            auto rowIndicesHost          = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rowIndices);
            auto offsetForRowOrdinalHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), offsetForRowOrdinal);
            auto packedColumnIndicesHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), packedColumnIndices);
            auto packedWeightsHost       = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), packedWeights);
            
            // for convenience of the code below, put the "next" rowOffset at the end:
            rowOffsets.push_back(rowOffset);
            for (int rowOrdinal=0; rowOrdinal<numRows; rowOrdinal++)
            {
              rowIndicesHost(rowOrdinal)          = nonIdentityDofs[rowOrdinal];
              int thisRowOffset                   = rowOffsets[rowOrdinal];
              offsetForRowOrdinalHost(rowOrdinal) = thisRowOffset;
              int nextRowOffset                   = rowOffsets[rowOrdinal+1];
              for (int i=thisRowOffset; i<nextRowOffset; i++)
              {
                packedColumnIndicesHost(i) = colIDs[i];
                packedWeightsHost(i)       = weights[i];
              }
            }
            offsetForRowOrdinalHost(numRows) = rowOffset;
            
            Kokkos::deep_copy(rowIndices, rowIndicesHost);
            Kokkos::deep_copy(offsetForRowOrdinal, offsetForRowOrdinalHost);
            Kokkos::deep_copy(packedColumnIndices, packedColumnIndicesHost);
            Kokkos::deep_copy(packedWeights, packedWeightsHost);
            
            OrientationOperator<DT> orientationOperator(rowIndices, offsetForRowOrdinal, packedColumnIndices, packedWeights);
            operatorsHost(edgeId, edgeOrt) = orientationOperator;
          }
        }
      }
    }
    
    Kokkos::deep_copy(operators, operatorsHost);
    
    return operators;
  }

//template<typename DT>
//  template<typename BasisHostType>
//  typename OrientationTools<DT>::OperatorViewType
//  OrientationTools<DT>::createFaceOperatorsInternal(const BasisHostType* basis, CoeffMatrixDataViewType matData) {
//    const std::string name(basis->getName());
//    OperatorViewType operators;
//
//    const auto cellTopo = basis->getBaseCellTopology();
//    const ordinal_type numEdges = cellTopo.getSubcellCount(1);
//    const ordinal_type numFaces = cellTopo.getSubcellCount(2);
//    ordinal_type matDim = 0, matDim1 = 0, matDim2 = 0, numOrts = 0, numSubCells;
//    for(ordinal_type i=0; i<numEdges; ++i) {
//      matDim1 = std::max(matDim1, basis->getDofCount(1,i));
//      numOrts = std::max(numOrts,2);
//    }
//    for(ordinal_type i=0; i<numFaces; ++i) {
//      matDim2 = std::max(matDim2, basis->getDofCount(2,i));
//      numOrts = std::max(numOrts,2*ordinal_type(cellTopo.getSideCount(2,i)));
//    }
//    matDim = std::max(matDim1,matDim2);
//    numSubCells = (matDim1>0)*numEdges + (matDim2>0)*numFaces;
//
//
//    operators = OperatorViewType("Orientation::Operators::"+name,
//                                 numSubCells,
//                                 numOrts);
//
//    
//    auto ordinalToTag = Kokkos::create_mirror_view_and_copy(typename DT::memory_space(), basis->getAllDofTags());
//    auto tagToOrdinal = Kokkos::create_mirror_view_and_copy(typename DT::memory_space(), basis->getAllDofOrdinal());
//
//    ordinal_type numVerts(0), numEdges(0), numFaces(0);
//
//    if (basis->requireOrientation()) {
//      numVerts = cellTopo.getVertexCount()*ordinal_type(basis->getDofCount(0, 0) > 0);
//      numEdges = cellTopo.getEdgeCount()*ordinal_type(basis->getDofCount(1, 0) > 0);
//      numFaces = cellTopo.getFaceCount()*ordinal_type(basis->getDofCount(2, 0) > 0);
//    }
//
//    bool leftMultiply = true;
//
//    const Kokkos::RangePolicy<typename DT::execution_space> policy(0, numCells);
//    typedef F_modifyBasisByOrientation
//      <decltype(orts),
//       decltype(output),decltype(input),
//       decltype(ordinalToTag),decltype(tagToOrdinal),
//       decltype(matData)> FunctorType;
//    Kokkos::parallel_for
//      (policy,
//       FunctorType(orts,
//                   output, input,
//                   ordinalToTag, tagToOrdinal,
//                   matData,
//                   cellDim, numVerts, numEdges, numFaces,
//                   numPoints, dimBasis, leftMultiply, transpose));
//
//    
//    
//    
//    
//    
//    
//    
//    return operators;
//  }

  //
  // HGRAD elements
  //
  template<typename DT>
  template<typename BasisHostType>
  void
  OrientationTools<DT>::
  init_HGRAD(typename OrientationTools<DT>::CoeffMatrixDataViewType matData,
                                 BasisHostType const *cellBasis,
                                 const bool inverse) {

    const auto cellTopo = cellBasis->getBaseCellTopology();
    const ordinal_type numEdges = cellTopo.getSubcellCount(1);
    const ordinal_type numFaces = cellTopo.getSubcellCount(2);
    Intrepid2::BasisPtr<typename BasisHostType::DeviceType,
        typename BasisHostType::OutputValueType,
        typename BasisHostType::PointValueType>
    basisPtr;
    BasisHostType const *subcellBasis;
    
    { //edges
      subcellBasis = cellBasis; // if (dim==1)
      const ordinal_type numOrt = 2;
      for (ordinal_type edgeId=0;edgeId<numEdges;++edgeId) {
        if(cellBasis->getDofCount(1, edgeId) < 2) continue;
        if(cellTopo.getDimension()!=1) {
          basisPtr = cellBasis->getSubCellRefBasis(1,edgeId);
          subcellBasis = basisPtr.get();
        }

        for (ordinal_type edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
          auto mat = Kokkos::subview(matData,
              edgeId, edgeOrt,
              Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HGRAD
            (mat,
              *subcellBasis, *cellBasis,
              edgeId, edgeOrt, inverse);
        }
      }
    }
    { //faces
      subcellBasis = cellBasis; // if(dim==2)
      for (ordinal_type faceId=0;faceId<numFaces;++faceId) {
        // this works for triangles (numOrt=6) and quadrilaterals (numOrt=8)
        const ordinal_type numOrt = 2*cellTopo.getSideCount(2,faceId);
        if(cellBasis->getDofCount(2, faceId) < 1) continue;
        if(cellTopo.getDimension()!=2) {
          basisPtr = cellBasis->getSubCellRefBasis(2,faceId);
          subcellBasis = basisPtr.get();
        }
        for (ordinal_type faceOrt=0;faceOrt<numOrt;++faceOrt) {
          auto mat = Kokkos::subview(matData,
              numEdges+faceId, faceOrt,
              Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HGRAD
            (mat,
             *subcellBasis, *cellBasis,
             faceId, faceOrt, inverse);
        }
      }
    }
  }

  //
  // HCURL elements
  //
  template<typename DT>
  template<typename BasisHostType>
  void
  OrientationTools<DT>::
  init_HCURL(typename OrientationTools<DT>::CoeffMatrixDataViewType matData,
      BasisHostType const *cellBasis, const bool inverse) {
    const auto cellTopo = cellBasis->getBaseCellTopology();
    const ordinal_type numEdges = cellTopo.getSubcellCount(1);
    const ordinal_type numFaces = cellTopo.getSubcellCount(2);
    Intrepid2::BasisPtr<typename BasisHostType::DeviceType,
        typename BasisHostType::OutputValueType,
        typename BasisHostType::PointValueType>
    basisPtr;
    BasisHostType const* subcellBasis;

    { // edges
      subcellBasis = cellBasis; // if (dim==1)
      const ordinal_type numOrt = 2;
      for (ordinal_type edgeId=0;edgeId<numEdges;++edgeId) {
        if(cellBasis->getDofCount(1, edgeId) < 1) continue;
        if(cellTopo.getDimension()!=1) {
          basisPtr = cellBasis->getSubCellRefBasis(1,edgeId);
          subcellBasis = basisPtr.get();
        }
        for (ordinal_type edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
          auto mat = Kokkos::subview(matData,
                                     edgeId, edgeOrt,
                                     Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HCURL(mat,
              *subcellBasis, *cellBasis,
              edgeId, edgeOrt, inverse);
        }
      }
    }
    { //faces
      subcellBasis = cellBasis; // if (dim==2)
      for (ordinal_type faceId=0;faceId<numFaces;++faceId) {
        // this works for triangles (numOrt=6) and quadratures (numOrt=8)
        const ordinal_type numOrt = 2*cellTopo.getSideCount(2,faceId);
        if(cellBasis->getDofCount(2, faceId) < 1) continue;
        if(cellTopo.getDimension()!=2) {
          basisPtr = cellBasis->getSubCellRefBasis(2,faceId);
          subcellBasis = basisPtr.get();
        }
        for (ordinal_type faceOrt=0;faceOrt<numOrt;++faceOrt) {
          auto mat = Kokkos::subview(matData,
                                     numEdges+faceId, faceOrt,
                                     Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HCURL
            (mat,
             *subcellBasis, *cellBasis,
             faceId, faceOrt, inverse);
        }
      }
    }
  }

  //
  // HDIV elements
  //
  template<typename DT>
  template<typename BasisHostType>
  void
  OrientationTools<DT>::
  init_HDIV(typename OrientationTools<DT>::CoeffMatrixDataViewType matData,
      BasisHostType const *cellBasis, const bool inverse) {
    const auto cellTopo = cellBasis->getBaseCellTopology();
    const ordinal_type numSides = cellTopo.getSideCount();
    const ordinal_type sideDim = cellTopo.getDimension()-1;
    Intrepid2::BasisPtr<typename BasisHostType::DeviceType,
        typename BasisHostType::OutputValueType,
        typename BasisHostType::PointValueType>
    subcellBasisPtr;

    {
      for (ordinal_type sideId=0;sideId<numSides;++sideId) {
        if(cellBasis->getDofCount(sideDim, sideId) < 1) continue;
        const ordinal_type numOrt = (sideDim == 1) ? 2 : 2*cellTopo.getSideCount(sideDim,sideId);
        subcellBasisPtr = cellBasis->getSubCellRefBasis(sideDim,sideId);
        for (ordinal_type faceOrt=0;faceOrt<numOrt;++faceOrt) {
          auto mat = Kokkos::subview(matData, 
                                     sideId, faceOrt,
                                     Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HDIV(mat,
              *subcellBasisPtr, *cellBasis,
              sideId, faceOrt, inverse);
        }
      }
    }
  }

  //
  // HVOL elements (used for 2D and 1D side cells)
  //
  template<typename DT>
  template<typename BasisHostType>
  void
  OrientationTools<DT>::
  init_HVOL(typename OrientationTools<DT>::CoeffMatrixDataViewType matData,
                                 BasisHostType const *cellBasis, const bool inverse) {

    const auto cellTopo = cellBasis->getBaseCellTopology();
    const ordinal_type numEdges = (cellTopo.getDimension()==1);
    const ordinal_type numFaces = (cellTopo.getDimension()==2);

    {
      const ordinal_type numOrt = 2;
      for (ordinal_type edgeId=0;edgeId<numEdges;++edgeId) {
        if(cellBasis->getDofCount(1, edgeId) < 1) continue;
        for (ordinal_type edgeOrt=0;edgeOrt<numOrt;++edgeOrt) {
          auto mat = Kokkos::subview(matData,
              edgeId, edgeOrt,
              Kokkos::ALL(), Kokkos::ALL());
          Impl::OrientationTools::getCoeffMatrix_HVOL
            (mat, *cellBasis, edgeOrt, inverse);
        }
      }
    }
    {
      for (ordinal_type faceId=0;faceId<numFaces;++faceId) {
        // this works for triangles (numOrt=6) and quadratures (numOrt=8)
        const ordinal_type numOrt = 2*cellTopo.getSideCount(2,faceId);
        if(cellBasis->getDofCount(2, faceId) < 1) continue;
        for (ordinal_type faceOrt=0;faceOrt<numOrt;++faceOrt) {
          auto mat = Kokkos::subview(matData,
              numEdges+faceId, faceOrt,
              Kokkos::ALL(), Kokkos::ALL());
            Impl::OrientationTools::getCoeffMatrix_HVOL
              (mat, *cellBasis, faceOrt, inverse);
        }
      }
    }
  }

  template<typename DT>
  template<typename BasisType>
  typename OrientationTools<DT>::CoeffMatrixDataViewType
  OrientationTools<DT>::createCoeffMatrix(const BasisType* basis) {
    Kokkos::push_finalize_hook( [=] {
      ortCoeffData=OrientationTools<DT>::OrtCoeffDataType();
    });

    const KeyType key(basis->getName(), basis->getDegree());
    const auto found = ortCoeffData.find(key);
    
    CoeffMatrixDataViewType matData;
    if (found == ortCoeffData.end()) {
      {
        auto basis_host = basis->getHostBasis();
        matData = createCoeffMatrixInternal(basis_host.getRawPtr());
        ortCoeffData.insert(std::make_pair(key, matData));
      }
    } else {
      matData = found->second;
    }
    
    return matData;
  }

  template<typename DT>
  template<typename BasisType>
  typename OrientationTools<DT>::CoeffMatrixDataViewType
  OrientationTools<DT>::createInvCoeffMatrix(const BasisType* basis) {
    Kokkos::push_finalize_hook( [=] {
      ortInvCoeffData=OrientationTools<DT>::OrtCoeffDataType();
    });

    const KeyType key(basis->getName(), basis->getDegree());
    const auto found = ortInvCoeffData.find(key);
    
    CoeffMatrixDataViewType matData;
    if (found == ortInvCoeffData.end()) {
      {
        auto basis_host = basis->getHostBasis();
        matData = createCoeffMatrixInternal(basis_host.getRawPtr(),true);
        ortInvCoeffData.insert(std::make_pair(key, matData));
      }
    } else {
      matData = found->second;
    }
    
    return matData;
  }
  
  template<typename DT>
  template<typename BasisType>
  std::tuple<typename OrientationTools<DT>::OperatorViewType, typename OrientationTools<DT>::OperatorViewType>
  OrientationTools<DT>::createOperators(const BasisType* basis) {
    Kokkos::push_finalize_hook( [=] {
      edgeOperatorData=std::map<std::pair<std::string,ordinal_type>, typename OrientationTools<DT>::OperatorViewType>();
      faceOperatorData=std::map<std::pair<std::string,ordinal_type>, typename OrientationTools<DT>::OperatorViewType>();
    });

    const std::pair<std::string,ordinal_type> key(basis->getName(), basis->getDegree());
    const auto edgeFound = edgeOperatorData.find(key);
    const auto faceFound = faceOperatorData.find(key);
    
    OperatorViewType edgeOperator, faceOperator;
    if (edgeFound == edgeOperatorData.end()) {
      {
        auto basis_host = basis->getHostBasis();
        auto matData = createCoeffMatrix(basis);
        edgeOperator = createEdgeOperatorsInternal(basis_host.get(), matData);
        edgeOperatorData.insert(std::make_pair(key, edgeOperator));
        
        faceOperator = createFaceOperatorsInternal(basis_host.get(), matData);
        faceOperatorData.insert(std::make_pair(key, faceOperator));
      }
    } else {
      edgeOperator = edgeFound->second;
      INTREPID2_TEST_FOR_EXCEPTION(faceFound == faceOperatorData.end(), std::invalid_argument, "edgeOperator found while faceOperator was not");
      faceOperator = faceFound->second;
    }
      
    return std::make_tuple(edgeOperator,faceOperator);
  }

  template<typename DT>
  template<typename BasisType>
  std::tuple<typename OrientationTools<DT>::OperatorViewType, typename OrientationTools<DT>::OperatorViewType>
  OrientationTools<DT>::createInvOperators(const BasisType* basis) {
    Kokkos::push_finalize_hook( [=] {
      invEdgeOperatorData=std::map<std::pair<std::string,ordinal_type>, typename OrientationTools<DT>::OperatorViewType>();
      invFaceOperatorData=std::map<std::pair<std::string,ordinal_type>, typename OrientationTools<DT>::OperatorViewType>();
    });

    const std::pair<std::string,ordinal_type> key(basis->getName(), basis->getDegree());
    const auto edgeFound = invEdgeOperatorData.find(key);
    const auto faceFound = invFaceOperatorData.find(key);
    
    OperatorViewType edgeOperator, faceOperator;
    if (edgeFound == edgeOperatorData.end()) {
      {
        auto basis_host = basis->getHostBasis();
        auto matData = createInvCoeffMatrix(basis);
        edgeOperator = createEdgeOperatorsInternal(basis_host, matData);
        invEdgeOperatorData.insert(std::make_pair(key, edgeOperator));
        
        faceOperator = createFaceOperatorsInternal(basis_host, matData);
        invFaceOperatorData.insert(std::make_pair(key, faceOperator));
      }
    } else {
      edgeOperator = edgeFound->second;
      INTREPID2_TEST_FOR_EXCEPTION(faceFound == faceOperatorData.end(), std::invalid_argument, "edgeOperator found while faceOperator was not");
      faceOperator = faceFound->second;
    }
      
    return std::make_tuple(edgeOperator,faceOperator);
  }

  template<typename DT>
  void OrientationTools<DT>::clearCoeffMatrix() {
    ortCoeffData.clear();
    ortInvCoeffData.clear();
  }
}

#endif

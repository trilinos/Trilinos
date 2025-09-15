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
    const int EDGE_DIM = 1;
    const int FACE_DIM = 2;
    
    const std::string name(basis->getName());
    OperatorViewType operators;

    auto matDataHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), matData);

    ordinal_type matDim1 = 0, matDim2 = 0, numOrts = 0, numEdgeOrts = 0;
    const auto cellTopo = basis->getBaseCellTopology();
    {
      const ordinal_type numEdges = cellTopo.getSubcellCount(EDGE_DIM);
      const ordinal_type numFaces = cellTopo.getSubcellCount(FACE_DIM);
      for(ordinal_type i=0; i<numEdges; ++i) {
        matDim1 = std::max(matDim1, basis->getDofCount(EDGE_DIM,i));
        numEdgeOrts = std::max(numOrts,2);
        numOrts     = numEdgeOrts;
      }
      for(ordinal_type i=0; i<numFaces; ++i) {
        matDim2 = std::max(matDim2, basis->getDofCount(2,i));
        numOrts = std::max(numOrts,2*ordinal_type(cellTopo.getSideCount(FACE_DIM,i)));
      }
      
      operators = OperatorViewType("Orientation::EdgeOperators::"+name,
                                   numEdges, numOrts, 2);
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
      
      const ordinal_type numEdges = cellTopo.getSubcellCount(1);
      const ordinal_type numFaces = cellTopo.getSubcellCount(2);
      
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
            
            const ordinal_type ordEdge = (1 < tagToOrdinal.extent(0) ? (static_cast<size_type>(edgeId) < tagToOrdinal.extent(1) ? tagToOrdinal(EDGE_DIM, edgeId, 0) : -1) : -1);
            
            ordinal_type rowOffset = 0;
            if (ordEdge != -1) {
              const ordinal_type ndofEdge = ordinalToTag(ordEdge, 3);
              const auto mat = Kokkos::subview(matDataHost,
                                               edgeId, edgeOrt,
                                               Kokkos::ALL(), Kokkos::ALL());
              
              for (ordinal_type i=0;i<ndofEdge;++i) {
                const ordinal_type ii = tagToOrdinal(EDGE_DIM, edgeId, i);
                
                // first pass for ii:
                // check whether this is different from the identity
                // count number of nonzeros in this row
                int nnz = 0;
                bool deviatesFromIdentity = false;
                for (ordinal_type l=0;l<ndofEdge;++l) {
                  const ordinal_type ll = tagToOrdinal(EDGE_DIM, edgeId, l);
                  auto & mat_il = mat(i,l);
                  if (mat_il != 0.0)
                  {
                    nnz++;
                    if ((mat_il != 1.0) || (ii != ll))
                    {
                      deviatesFromIdentity = true;
                    }
                  }
                  else if (ii == ll)
                  {
                    // zero entry on the diagonal is also a deviation from the identity
                    deviatesFromIdentity = true;
                  }
                } // column
                INTREPID2_TEST_FOR_EXCEPTION(nnz == 0, std::invalid_argument, "Each dof should have *some* nonzero weight");
                if (deviatesFromIdentity)
                {
                  // then we store the nonzeros for ii
                  nonIdentityDofs.push_back(ii);
                  rowOffsets.push_back(rowOffset);
                  rowOffset += nnz;
                  
                  for (ordinal_type l=0;l<ndofEdge;++l)
                  {
                    const ordinal_type ll = tagToOrdinal(EDGE_DIM, edgeId, l);
                    auto & mat_il = mat(i,l);
                    if (mat_il != 0.0)
                    {
                      colIDs.push_back(ll);
                      weights.push_back(mat_il);
                    }
                  }
                }
              } // row
              rowOffsets.push_back(rowOffset);
            }
            
            std::vector<bool> transposeVector {false, true};
            for (const bool transpose : transposeVector)
            {
              ordinal_type transposeInt = transpose ? 1 : 0;
              OrientationOperator<DT> orientationOperator = constructOrientationOperatorInternal(nonIdentityDofs, rowOffsets, colIDs, weights, transpose);
              operatorsHost(edgeId, edgeOrt, transposeInt) = orientationOperator;
            }
          }
        }
      }
    }
    
    Kokkos::deep_copy(operators, operatorsHost);
    
    return operators;
  }

  template<typename DT>
  template<typename BasisHostType>
  typename OrientationTools<DT>::OperatorViewType
  OrientationTools<DT>::createFaceOperatorsInternal(const BasisHostType* basis, CoeffMatrixDataViewType matData)
  {
    const int EDGE_DIM = 1;
    const int FACE_DIM = 2;
    
    const std::string name(basis->getName());
    OperatorViewType operators;

    auto matDataHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), matData);

    ordinal_type matDim1 = 0, matDim2 = 0, numOrts = 0, numEdgeOrts = 0, maxFaceOrts = 0;
    const auto cellTopo = basis->getBaseCellTopology();
    {
      const ordinal_type numEdges = cellTopo.getSubcellCount(EDGE_DIM);
      const ordinal_type numFaces = cellTopo.getSubcellCount(FACE_DIM);
      for(ordinal_type i=0; i<numEdges; ++i) {
        matDim1 = std::max(matDim1, basis->getDofCount(EDGE_DIM,i));
        numEdgeOrts = std::max(numOrts,2);
        numOrts     = numEdgeOrts;
      }
      for(ordinal_type i=0; i<numFaces; ++i) {
        matDim2 = std::max(matDim2, basis->getDofCount(FACE_DIM,i));
        const ordinal_type faceEdgeCount = static_cast<ordinal_type>(cellTopo.getSideCount(FACE_DIM,i));
        maxFaceOrts = std::max(maxFaceOrts,2*faceEdgeCount);
        // 2*(#face edges): a formula that happens to work for triangles and quads: 6 triangle orientations, 8 quad orientations.
        numOrts = std::max(numOrts,maxFaceOrts);
        
      }
      
      operators = OperatorViewType("Orientation::FaceOperators::"+name,
                                   numFaces, numOrts, 2);
    }
    
    // NOTE: the OrientationOperators within operatorsHost contain *device* views.
    auto operatorsHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), operators);

    if (matDim2 == 0)
    {
      // no non-trivial face operators (no face dofs)
      return operators;
    }
    
    // determine if there are edge dofs
    ordinal_type existEdgeDofs = (matDim1 > 0) ? 1 : 0;
    
    if ((basis->getFunctionSpace() != FUNCTION_SPACE_HDIV) || (cellTopo.getDimension() == 3)) // no face orientations for H(div) except in 3D
    {
      auto ordinalToTag = basis->getAllDofTags();
      auto tagToOrdinal = basis->getAllDofOrdinal();
      
      const ordinal_type numEdges = cellTopo.getSubcellCount(1);
      const ordinal_type numFaces = cellTopo.getSubcellCount(2);
      
      if (numFaces > 0)
      {
        for (ordinal_type faceId=0;faceId<numFaces;++faceId)
        {
          const ordinal_type ordFace = (2 < tagToOrdinal.extent(0) ? (static_cast<size_type>(faceId) < tagToOrdinal.extent(1) ? tagToOrdinal(FACE_DIM, faceId, 0) : -1) : -1);
          
          const ordinal_type numFaceOrts = 2*ordinal_type(cellTopo.getSideCount(FACE_DIM,faceId));
          // 2*(#face edges): a formula that happens to work for triangles and quads: 6 triangle orientations, 8 quad orientations.
          for (ordinal_type faceOrt=0; faceOrt<numFaceOrts; faceOrt++)
          {
            std::vector<ordinal_type> nonIdentityDofs;
            std::vector<ordinal_type> rowOffsets; // within the column storage
            std::vector<ordinal_type> colIDs;
            std::vector<double> weights;
            
            ordinal_type rowOffset = 0;
            
            if (ordFace != -1) {
              const ordinal_type ndofFace = ordinalToTag(ordFace, 3);
              const auto mat = Kokkos::subview(matDataHost,
                                               numEdges*existEdgeDofs+faceId, faceOrt,
                                               Kokkos::ALL(), Kokkos::ALL());
              for (ordinal_type i=0;i<ndofFace;++i) {
                
                const ordinal_type ii = tagToOrdinal(FACE_DIM, faceId, i);
                
                // first pass for ii:
                // check whether this is different from the identity
                // count number of nonzeros
                int nnz = 0;
                bool deviatesFromIdentity = false;
                for (ordinal_type l=0;l<ndofFace;++l) {
                  const ordinal_type ll = tagToOrdinal(FACE_DIM, faceId, l);
                  auto & mat_il = mat(i,l);
                  if (mat_il != 0.0)
                  {
                    nnz++;
                    if ((mat_il != 1.0) || (ii != ll))
                    {
                      deviatesFromIdentity = true;
                    }
                  }
                  else if (ii == ll)
                  {
                    // zero entry on the diagonal is also a deviation from the identity
                    deviatesFromIdentity = true;
                  }
                } // column
                INTREPID2_TEST_FOR_EXCEPTION(nnz == 0, std::invalid_argument, "Each dof should have *some* nonzero weight");
                if (deviatesFromIdentity)
                {
                  // then we store the nonzeros for ii
                  nonIdentityDofs.push_back(ii);
                  rowOffsets.push_back(rowOffset);
                  rowOffset += nnz;
                  
                  for (ordinal_type l=0;l<ndofFace;++l)
                  {
                    const ordinal_type ll = tagToOrdinal(FACE_DIM, faceId, l);
                    auto & mat_il = mat(i,l);
                    if (mat_il != 0.0)
                    {
                      colIDs.push_back(ll);
                      weights.push_back(mat_il);
                    }
                  }
                } // if (deviatesFromIdentity)
              } // row
              rowOffsets.push_back(rowOffset);
              std::vector<bool> transposeVector {false, true};
              for (const bool transpose : transposeVector)
              {
                ordinal_type transposeInt = transpose ? 1 : 0;
                OrientationOperator<DT> orientationOperator = constructOrientationOperatorInternal(nonIdentityDofs, rowOffsets, colIDs, weights, transpose);
                operatorsHost(faceId, faceOrt, transposeInt) = orientationOperator;
              }
            } // if (ordFace != -1)
          }
        }
      }
    }
    
    Kokkos::deep_copy(operators, operatorsHost);
    
    return operators;
  }

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
        if(cellBasis->getDofCount(1, edgeId) < 1) continue;
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
    static bool hookRegistered = false;
    if (!hookRegistered)
    {
      Kokkos::push_finalize_hook( [=] {
        ortCoeffData=std::map<std::pair<std::string,ordinal_type>, typename OrientationTools<DT>::CoeffMatrixDataViewType>();
      });
      hookRegistered = true;
    }

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
    static bool hookRegistered = false;
    if (!hookRegistered)
    {
      Kokkos::push_finalize_hook( [=] {
        ortInvCoeffData=std::map<std::pair<std::string,ordinal_type>, typename OrientationTools<DT>::CoeffMatrixDataViewType>();
      });
      hookRegistered = true;
    }

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
    static bool hookRegistered = false;
    if (!hookRegistered)
    {
      Kokkos::push_finalize_hook( [=] {
        edgeOperatorData.clear();
        faceOperatorData.clear();
      });
      hookRegistered = true;
    }

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
    static bool hookRegistered = false;
    if (!hookRegistered)
    {
      Kokkos::push_finalize_hook( [=] {
        invEdgeOperatorData.clear();
        invFaceOperatorData.clear();
      });
      hookRegistered = true;
    }

    const std::pair<std::string,ordinal_type> key(basis->getName(), basis->getDegree());
    const auto edgeFound = invEdgeOperatorData.find(key);
    const auto faceFound = invFaceOperatorData.find(key);
    
    OperatorViewType edgeOperator, faceOperator;
    if (edgeFound == invEdgeOperatorData.end()) {
      {
        auto basis_host = basis->getHostBasis();
        auto matData = createInvCoeffMatrix(basis);
        
        edgeOperator = createEdgeOperatorsInternal(basis_host.get(), matData);
        invEdgeOperatorData.insert(std::make_pair(key, edgeOperator));
        
        faceOperator = createFaceOperatorsInternal(basis_host.get(), matData);
        invFaceOperatorData.insert(std::make_pair(key, faceOperator));
      }
    } else {
      edgeOperator = edgeFound->second;
      INTREPID2_TEST_FOR_EXCEPTION(faceFound == invFaceOperatorData.end(), std::invalid_argument, "edgeOperator found while faceOperator was not");
      faceOperator = faceFound->second;
    }
      
    return std::make_tuple(edgeOperator,faceOperator);
  }

  template<typename DT>
  void OrientationTools<DT>::clearCoeffMatrix() {
    ortCoeffData.clear();
    ortInvCoeffData.clear();
  }

  template<typename DT>
  OrientationOperator<DT> OrientationTools<DT>::constructOrientationOperatorInternal(const std::vector<ordinal_type> &nonIdentityDofs,
                                                                                     const std::vector<ordinal_type> &rowOffsets,
                                                                                     const std::vector<ordinal_type> &colIDs,
                                                                                     const std::vector<double> &weights, const bool transpose)
  {
    static bool hookRegistered = false;
    if (!hookRegistered)
    {
      Kokkos::push_finalize_hook( [=] {
        doubleViewAllocations.clear();
        ordinalViewAllocations.clear();
      });
      hookRegistered = true;
    }
    
    const int numRows = static_cast<int>(nonIdentityDofs.size());
    if (numRows > 0)
    {
      if (!transpose)
      {
        const int numWeights = rowOffsets[numRows];
        
        const bool isWeightedPermutation = (nonIdentityDofs.size() == weights.size());
        const int numOffsetsToStore = isWeightedPermutation ? 0 : numRows+1; // convenient to be able to check the "next" row offset to get a column count, even when on the last row
        
        Kokkos::View<ordinal_type*,DT> rowIndices("OrientationOperator: rowIndices", numRows);
        Kokkos::View<ordinal_type*,DT> offsetForRowOrdinal("OrientationOperator: rowOffsets", numOffsetsToStore);
        Kokkos::View<ordinal_type*,DT> packedColumnIndices("OrientationOperator: packedColumnIndices", numWeights);
        Kokkos::View<double*,      DT> packedWeights("OrientationOperator: packedWeights", numWeights);
        
        auto rowIndicesHost          = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), rowIndices);
        auto offsetForRowOrdinalHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), offsetForRowOrdinal);
        auto packedColumnIndicesHost = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), packedColumnIndices);
        auto packedWeightsHost       = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), packedWeights);
        
        for (int rowOrdinal=0; rowOrdinal<numRows; rowOrdinal++)
        {
          rowIndicesHost(rowOrdinal)          = nonIdentityDofs[rowOrdinal];
          int thisRowOffset                   = rowOffsets[rowOrdinal];
          int nextRowOffset                   = rowOffsets[rowOrdinal+1];
          if (isWeightedPermutation)
          {
            INTREPID2_TEST_FOR_EXCEPTION(nextRowOffset - thisRowOffset != 1, std::invalid_argument, "Invalid orientation operator arguments: if number of weights is equal to the number of row indices, then row offsets should exactly match row ordinals (i.e., there should be exactly one column for every row entry).");
          }
          else
          {
            offsetForRowOrdinalHost(rowOrdinal) = thisRowOffset;
          }
          for (int i=thisRowOffset; i<nextRowOffset; i++)
          {
            packedColumnIndicesHost(i) = colIDs[i];
            packedWeightsHost(i)       = weights[i];
          }
        }
        if (!isWeightedPermutation)
        {
          offsetForRowOrdinalHost(numRows) = numWeights;
        }
        
        Kokkos::deep_copy(rowIndices, rowIndicesHost);
        Kokkos::deep_copy(offsetForRowOrdinal, offsetForRowOrdinalHost);
        Kokkos::deep_copy(packedColumnIndices, packedColumnIndicesHost);
        Kokkos::deep_copy(packedWeights, packedWeightsHost);
        
        // statically store the managed views
        ordinalViewAllocations.push_back(rowIndices);
        ordinalViewAllocations.push_back(offsetForRowOrdinal);
        ordinalViewAllocations.push_back(packedColumnIndices);
        doubleViewAllocations.push_back(packedWeights);
        
        using UnmanagedDoubleView  = typename OrientationOperator<DT>::UnmanagedDoubleView;
        using UnmanagedOrdinalView = typename OrientationOperator<DT>::UnmanagedOrdinalView;
        
        UnmanagedOrdinalView rowIndicesUnmanaged         (rowIndices.data(),          rowIndices.extent(0));
        UnmanagedOrdinalView offsetForRowOrdinalUnmanaged(offsetForRowOrdinal.data(), offsetForRowOrdinal.extent(0));
        UnmanagedOrdinalView packedColumnIndicesUnmanaged(packedColumnIndices.data(), packedColumnIndices.extent(0));
        UnmanagedDoubleView  packedWeightsUnmanaged      (packedWeights.data(),       packedWeights.extent(0));
        
        if (!isWeightedPermutation)
        {
          OrientationOperator<DT> orientationOperator(rowIndicesUnmanaged, offsetForRowOrdinalUnmanaged,
                                                      packedColumnIndicesUnmanaged, packedWeightsUnmanaged);
          
          return orientationOperator;
        }
        else
        {
          OrientationOperator<DT> orientationOperator(rowIndicesUnmanaged,
                                                      packedColumnIndicesUnmanaged, packedWeightsUnmanaged);
          
          return orientationOperator;
        }
      }
      else
      {
        // for the transpose case, we construct the arguments for the non-transpose case,
        // and then call this method with transpose = false.
        std::map<ordinal_type, std::map<ordinal_type,double> > transposeOperator; // column to (row -> weight) lookup
        std::set<ordinal_type> opRows(nonIdentityDofs.begin(), nonIdentityDofs.end());
        
        const ordinal_type numRows = static_cast<ordinal_type>(nonIdentityDofs.size());
        for (ordinal_type rowOrdinal = 0; rowOrdinal < numRows; rowOrdinal++)
        {
          const ordinal_type & rowID = nonIdentityDofs[rowOrdinal];
          const ordinal_type & rowOffset = rowOffsets[rowOrdinal];
          
          const ordinal_type numCols = rowOffsets[rowOrdinal+1] - rowOffset;
          for (ordinal_type colOrdinal=0; colOrdinal<numCols; colOrdinal++)
          {
            const ordinal_type & colID = colIDs[rowOffset + colOrdinal];
            const double      & weight = weights[rowOffset + colOrdinal];
            transposeOperator[colID][rowID] = weight;
            if (opRows.find(colID) == opRows.end())
            {
              // then the original operator had an implicit identity row for colID; the transpose should have a 1 in diagonal
              transposeOperator[colID][colID] = 1.0;
            }
          }
        }
        
        std::vector<ordinal_type> nonIdentityColDofs;
        std::vector<ordinal_type> colOffsets;
        std::vector<ordinal_type> rowIDs;
        std::vector<double> weightsTranspose;
        
        colOffsets.push_back(0);
        for (const auto & entry : transposeOperator)
        {
          const ordinal_type & colID = entry.first;
          nonIdentityColDofs.push_back(colID);
          const auto & rowMap = entry.second;
          for (const auto & rowEntry : rowMap)
          {
            const ordinal_type & rowID = rowEntry.first;
            const double      & weight = rowEntry.second;
            
            rowIDs.push_back(rowID);
            weightsTranspose.push_back(weight);
          }
          colOffsets.push_back(static_cast<ordinal_type>(rowIDs.size()));
        }
        
        return constructOrientationOperatorInternal(nonIdentityColDofs, colOffsets, rowIDs, weightsTranspose, false);
      }
    }
    // identity; nothing to allocate or store
    return OrientationOperator<DT>();
  }
}

#endif

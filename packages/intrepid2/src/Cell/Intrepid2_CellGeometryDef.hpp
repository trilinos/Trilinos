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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov),
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_CellGeometry.hpp
    \brief  Implementation file for Intrepid2::CellGeometry class.

    \author Nathan V. Roberts
*/

#ifndef Intrepid2_CellGeometryDef_h
#define Intrepid2_CellGeometryDef_h

namespace Intrepid2
{

  namespace Impl
  {
    /** \brief Store host-only "members" of CellGeometry using a static map indexed on the CellGeometry pointer.  This allows us to avoid issues related to non-CUDA-aware members with a lambda capture of a CellGeometry object.
     */
    template<class PointScalar, int spaceDim, typename DeviceType>
    class CellGeometryHostMembers
    {
      using BasisPtr = Teuchos::RCP<Intrepid2::Basis<DeviceType,PointScalar,PointScalar> >;
      using CellGeometryType = CellGeometry<PointScalar,spaceDim,DeviceType>;
    public:
      // conceptually, these should be private members, but for the definition of these, we need them to be externally accessible.
      static std::map<const CellGeometryType *, shards::CellTopology> cellTopology_;
      static std::map<const CellGeometryType *, BasisPtr> basisForNodes_;
      
    public:
      static void constructorCalled(const CellGeometryType *cellGeometry, const shards::CellTopology &cellTopo, BasisPtr basisForNodes)
      {
        cellTopology_[cellGeometry]  = cellTopo;
        basisForNodes_[cellGeometry] = basisForNodes;
      }
      
      static void destructorCalled(const CellGeometryType *cellGeometry)
      {
        cellTopology_.erase(cellGeometry);
        basisForNodes_.erase(cellGeometry);
      }
      
      static BasisPtr getBasis(const CellGeometryType *cellGeometry)
      {
        return basisForNodes_[cellGeometry];
      }
      
      static const shards::CellTopology & getCellTopology(const CellGeometryType *cellGeometry)
      {
        return cellTopology_[cellGeometry];
      }
    };
  
    // member lookup map definitions for CellGeometryHostMembers:
    template< class PointScalar, int spaceDim, typename DeviceType > typename std::map<const CellGeometry<PointScalar,spaceDim,DeviceType> *, shards::CellTopology> CellGeometryHostMembers< PointScalar,spaceDim,DeviceType>::cellTopology_;
  
    template< class PointScalar, int spaceDim, typename DeviceType > typename std::map<const CellGeometry<PointScalar,spaceDim,DeviceType> *, Teuchos::RCP<Intrepid2::Basis<DeviceType,PointScalar,PointScalar> >> CellGeometryHostMembers< PointScalar,spaceDim,DeviceType>::basisForNodes_;
  }

  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  CellGeometry<PointScalar,spaceDim,DeviceType>::CellGeometry(const CellGeometry<PointScalar,spaceDim,DeviceType> &cellGeometry)
:
  nodeOrdering_(cellGeometry.nodeOrdering_),
  cellGeometryType_(cellGeometry.cellGeometryType_),
  subdivisionStrategy_(cellGeometry.subdivisionStrategy_),
  affine_(cellGeometry.affine_),
  orientations_(cellGeometry.orientations_),
  origin_(cellGeometry.origin_),
  domainExtents_(cellGeometry.domainExtents_),
  gridCellCounts_(cellGeometry.gridCellCounts_),
  tensorVertices_(cellGeometry.tensorVertices_),
  cellToNodes_(cellGeometry.cellToNodes_),
  nodes_(cellGeometry.nodes_),
  numCells_(cellGeometry.numCells_),
  numNodesPerCell_(cellGeometry.numNodesPerCell_)
  {
    // host-only registration with HostMemberLookup:
#ifndef __CUDA_ARCH__
    shards::CellTopology cellTopo = cellGeometry.cellTopology();
    BasisPtr basisForNodes = cellGeometry.basisForNodes();
    using HostMemberLookup = ::Intrepid2::Impl::CellGeometryHostMembers<PointScalar, spaceDim, DeviceType>;
    HostMemberLookup::constructorCalled(this, cellTopo, basisForNodes);
#endif
  }

  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  CellGeometry<PointScalar,spaceDim,DeviceType>::~CellGeometry()
  {
    // host-only deregistration with HostMemberLookup:
#ifndef __CUDA_ARCH__
    using HostMemberLookup = ::Intrepid2::Impl::CellGeometryHostMembers<PointScalar, spaceDim, DeviceType>;
    HostMemberLookup::destructorCalled(this);
#endif
  }

  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  int CellGeometry<PointScalar,spaceDim,DeviceType>::numCellsPerGridCell(SubdivisionStrategy subdivisionStrategy) const
  {
    switch (subdivisionStrategy) {
      case NO_SUBDIVISION:
        return 1;
      case TWO_TRIANGLES_LEFT:
      case TWO_TRIANGLES_RIGHT:
        return 2;
      case FOUR_TRIANGLES:
        return 4;
      case FIVE_TETRAHEDRA:
        return 5;
      case SIX_TETRAHEDRA:
        return 6;
    }
    return -1;
  }

  template<class PointScalar, int spaceDim, typename DeviceType>
  Data<PointScalar,DeviceType>
  CellGeometry<PointScalar,spaceDim,DeviceType>::allocateJacobianDataPrivate(const TensorPoints<PointScalar,DeviceType> &points, const int &pointsPerCell, const int startCell, const int endCell) const
  {
    ScalarView<PointScalar,DeviceType> data;
    const int rank = 4; // C,P,D,D
    const int CELL_DIM  = 0;
    const int POINT_DIM = 1;
    const int D1_DIM    = 2;
    const int D2_DIM    = 3;
    
    const int numCellsWorkset = (endCell == -1) ? (numCells_ - startCell) : (endCell - startCell);
    
    Kokkos::Array<int,7>               extents        { numCellsWorkset, pointsPerCell, spaceDim, spaceDim,        1,        1,        1 };
    Kokkos::Array<DataVariationType,7> variationType  {        CONSTANT,      CONSTANT, CONSTANT, CONSTANT, CONSTANT, CONSTANT, CONSTANT };
    
    int blockPlusDiagonalLastNonDiagonal = -1;
          
    if (cellGeometryType_ == UNIFORM_GRID)
    {
      if (uniformJacobianModulus() != 1)
      {
        variationType[CELL_DIM]  = MODULAR;
        variationType[POINT_DIM] = CONSTANT;
        variationType[D1_DIM]    = GENERAL;
        variationType[D2_DIM]    = GENERAL;
        
        int cellTypeModulus = uniformJacobianModulus();
        
        data = getMatchingViewWithLabel(points.getTensorComponent(0), "CellGeometryProvider: Jacobian data", cellTypeModulus, spaceDim, spaceDim);
      }
      else
      {
        // diagonal Jacobian
        variationType[D1_DIM] = BLOCK_PLUS_DIAGONAL;
        variationType[D2_DIM] = BLOCK_PLUS_DIAGONAL;
        blockPlusDiagonalLastNonDiagonal = -1;
        
        data = getMatchingViewWithLabel(points.getTensorComponent(0), "CellGeometryProvider: Jacobian data", spaceDim);
      }
    }
    else if (cellGeometryType_ == TENSOR_GRID)
    {
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "tensor grid support not yet implemented");
    }
    else if (cellGeometryType_ == FIRST_ORDER)
    {
      const bool simplex = (spaceDim + 1 == cellToNodes_.extent_int(1));
      if (simplex)
      {
        variationType[CELL_DIM]  = GENERAL;
        variationType[POINT_DIM] = CONSTANT; // affine: no point variation
        variationType[D1_DIM]    = GENERAL;
        variationType[D2_DIM]    = GENERAL;
        
        data = getMatchingViewWithLabel(data, "CellGeometryProvider: Jacobian data", numCells_, spaceDim, spaceDim);
      }
      else
      {
        variationType[CELL_DIM]  = GENERAL;
        variationType[D1_DIM]    = GENERAL;
        variationType[D2_DIM]    = GENERAL;
        if (affine_)
        {
          // no point variation
          variationType[POINT_DIM] = CONSTANT;
          data = getMatchingViewWithLabel(data, "CellGeometryProvider: Jacobian data", numCellsWorkset, spaceDim, spaceDim);
        }
        else
        {
          variationType[POINT_DIM] = GENERAL;
          data = getMatchingViewWithLabel(data, "CellGeometryProvider: Jacobian data", numCellsWorkset, pointsPerCell, spaceDim, spaceDim);
        }
      }
    }
    else if (cellGeometryType_ == HIGHER_ORDER)
    {
      // most general case: varies in all 4 dimensions
      variationType[CELL_DIM]  = GENERAL;
      variationType[POINT_DIM] = GENERAL;
      variationType[D1_DIM]    = GENERAL;
      variationType[D2_DIM]    = GENERAL;
      data = getMatchingViewWithLabel(data, "CellGeometryProvider: Jacobian data", numCellsWorkset, pointsPerCell, spaceDim, spaceDim);
    }
    else
    {
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "support for this CellGeometryType is not yet implemented");
    }
    
    Data<PointScalar,DeviceType> jacobianData(data,rank,extents,variationType,blockPlusDiagonalLastNonDiagonal);
    return jacobianData;
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  void CellGeometry<PointScalar,spaceDim,DeviceType>::setJacobianDataPrivate(Data<PointScalar,DeviceType> &jacobianData, const TensorPoints<PointScalar,DeviceType> &points,
                                                                                const int &pointsPerCell, const Data<PointScalar,DeviceType> &refData,
                                                                                const int startCell, const int endCell) const
  {
    const int numCellsWorkset = (endCell == -1) ? (numCells_ - startCell) : (endCell - startCell);
    
    if (cellGeometryType_ == UNIFORM_GRID)
    {
      if (uniformJacobianModulus() != 1)
      {
        int cellTypeModulus = uniformJacobianModulus();
        
        auto dataView3 = jacobianData.getUnderlyingView3(); // (cellTypeModulus, spaceDim, spaceDim) allocated in allocateJacobianDataPrivate()
        auto dataHost = Kokkos::create_mirror_view(dataView3);
        
        const int startCellType = startCell % cellTypeModulus;
        const int endCellType   = (numCellsWorkset >= cellTypeModulus) ? startCellType + cellTypeModulus : startCellType + numCellsWorkset;
        const int gridCellOrdinal = 0; // sample cell
        for (int cellType=startCellType; cellType<endCellType; cellType++)
        {
          const int subdivisionOrdinal = cellType % cellTypeModulus;
          const int nodeZero = 0;
          // simplex Jacobian formula is J_00 = x1 - x0, J_01 = x2 - x0, etc.
          for (int i=0; i<spaceDim; i++)
          {
            for (int j=0; j<spaceDim; j++)
            {
              const int node = j+1; // this is the only node other than the 0 node that has non-zero derivative in the j direction -- and this has unit derivative
              // nodeZero has derivative -1 in every dimension.
              const auto J_ij = subdivisionCoordinate(gridCellOrdinal, subdivisionOrdinal, node, i) - subdivisionCoordinate(gridCellOrdinal, subdivisionOrdinal, nodeZero, i);
              dataHost(cellType,i,j) = J_ij;
            }
          }
        }
        
        Kokkos::deep_copy(dataView3,dataHost);
      }
      else
      {
        // diagonal Jacobian
        auto dataView1 = jacobianData.getUnderlyingView1(); // (spaceDim) allocated in allocateJacobianDataPrivate()
        const auto domainExtents = domainExtents_;
        const auto gridCellCounts = gridCellCounts_;
        
        using ExecutionSpace = typename DeviceType::execution_space;
        auto policy = Kokkos::RangePolicy<>(ExecutionSpace(),0,spaceDim);
        Kokkos::parallel_for("fill jacobian", policy, KOKKOS_LAMBDA(const int d1)
        {
          // diagonal jacobian
          const double REF_SPACE_EXTENT = 2.0;
          dataView1(d1) = (domainExtents[d1] / REF_SPACE_EXTENT) / gridCellCounts[d1];
        });
        ExecutionSpace().fence();
      }
    }
    else if (cellGeometryType_ == TENSOR_GRID)
    {
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "tensor grid support not yet implemented");
    }
    else if ((cellGeometryType_ == FIRST_ORDER) || (cellGeometryType_ == HIGHER_ORDER))
    {
      const bool simplex = (spaceDim + 1 == cellToNodes_.extent_int(1));
      if (simplex)
      {
        auto dataView3 = jacobianData.getUnderlyingView3(); // (numCells_, spaceDim, spaceDim) allocated in allocateJacobianDataPrivate()
        
        // get local (shallow) copies to avoid implicit references to this
        auto cellToNodes = cellToNodes_;
        auto nodes       = nodes_;
        
        using ExecutionSpace = typename DeviceType::execution_space;
        auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({startCell,0,0},{numCellsWorkset,spaceDim,spaceDim});
        
        Kokkos::parallel_for("compute first-order simplex Jacobians", policy,
        KOKKOS_LAMBDA (const int &cellOrdinal, const int &d1, const int &d2) {
          const int nodeZero = 0;    // nodeZero has derivative -1 in every dimension.
          const int node     = d2+1; // this is the only node other than the 0 node that has non-zero derivative in the d2 direction -- and this has unit derivative (except in 1D, where each derivative is ±0.5)
          const auto & nodeCoord     = nodes(cellToNodes(cellOrdinal,node),     d1);
          const auto & nodeZeroCoord = nodes(cellToNodes(cellOrdinal,nodeZero), d1);
          const PointScalar J_ij = nodeCoord - nodeZeroCoord;
          dataView3(cellOrdinal,d1,d2) = (spaceDim != 1) ? J_ij : J_ij * 0.5;
        });
      }
      else
      {
        using CellTools = Intrepid2::CellTools<DeviceType>;
        auto basisForNodes = this->basisForNodes();
        
        if (affine_)
        {
          // no point variation
          auto dataView3 = jacobianData.getUnderlyingView3(); // (numCellsWorkset, spaceDim, spaceDim) allocated in allocateJacobianDataPrivate()

          // TODO: find an allocation-free way to do this… (consider modifying CellTools::setJacobian() to support affine case.)
          const int onePoint = 1;
          auto testPointView = getMatchingViewWithLabel(dataView3, "CellGeometryProvider: test point", onePoint, spaceDim);
          auto tempData      = getMatchingViewWithLabel(dataView3, "CellGeometryProvider: temporary Jacobian data", numCellsWorkset, onePoint, spaceDim, spaceDim);
          
          Kokkos::deep_copy(testPointView, 0.0);
          
          CellTools::setJacobian(tempData, testPointView, *this, basisForNodes, startCell, endCell);
          
          auto tempDataSubview = Kokkos::subview(tempData, Kokkos::ALL(), 0, Kokkos::ALL(), Kokkos::ALL());
          Kokkos::deep_copy(dataView3, tempDataSubview);
        }
        else
        {
          auto dataView = jacobianData.getUnderlyingView(); // (numCellsWorkset, pointsPerCell, spaceDim, spaceDim) allocated in allocateJacobianDataPrivate()
          
          // refData should contain the basis gradients; shape is (F,P,D)
          INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(!refData.isValid(), std::invalid_argument, "refData should be a valid container for cases with non-affine geometry");
          INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(refData.rank() != 3, std::invalid_argument, "refData should have shape (F,P,D)");
          INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(refData.extent_int(0) != basisForNodes->getCardinality(), std::invalid_argument, "refData should have shape (F,P,D)");
          INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(refData.extent_int(1) != points.extent_int(0), std::invalid_argument, "refData should have shape (F,P,D)");
          INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(refData.extent_int(2) != spaceDim, std::invalid_argument, "refData should have shape (F,P,D)");
          
          CellTools::setJacobian(dataView, *this, refData, startCell, endCell);
        }
      }
    }
    else
    {
      // TODO: handle the other cases
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "support for this CellGeometryType is not yet implemented");
    }
  }

  // Uniform grid constructor, with optional subdivision into simplices
  template<class PointScalar, int spaceDim, typename DeviceType>
  CellGeometry<PointScalar,spaceDim,DeviceType>::CellGeometry(const Kokkos::Array<PointScalar,spaceDim> &origin,
                                                              const Kokkos::Array<PointScalar,spaceDim> &domainExtents,
                                                              const Kokkos::Array<int,spaceDim> &gridCellCounts,
                                                              SubdivisionStrategy subdivisionStrategy,
                                                              HypercubeNodeOrdering nodeOrdering)
  :
  nodeOrdering_(nodeOrdering),
  cellGeometryType_(UNIFORM_GRID),
  subdivisionStrategy_(subdivisionStrategy),
  affine_(true),
  origin_(origin),
  domainExtents_(domainExtents),
  gridCellCounts_(gridCellCounts)
  {
    numCells_ = 1;
    for (int d=0; d<spaceDim; d++)
    {
      numCells_ *= gridCellCounts_[d];
    }
    numCells_ *= numCellsPerGridCell(subdivisionStrategy_);
    
    shards::CellTopology cellTopo; // will register with HostMemberLookup below
    if (subdivisionStrategy_ == NO_SUBDIVISION)
    {
      // hypercube
      numNodesPerCell_ = 1 << spaceDim; // 2^D vertices in a D-dimensional hypercube
      
      if (spaceDim == 1)
      {
        cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >());
      }
      else if (spaceDim == 2)
      {
        cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<> >());
      }
      else if (spaceDim == 3)
      {
        cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<> >());
      }
      else
      {
        // TODO: Once shards supports higher-dimensional hypercubes, initialize cellTopo accordingly
      }
    }
    else
    {
      // simplex
      numNodesPerCell_ = spaceDim + 1; // D+1 vertices in a D-dimensional simplex
      if (spaceDim == 2)
      {
        cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<> >());
      }
      else if (spaceDim == 3)
      {
        cellTopo = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<> >());
      }
      else
      {
        // TODO: Once shards supports higher-dimensional simplices, initialize cellTopo_ accordingly
      }
    }
    
    using BasisFamily = DerivedNodalBasisFamily<DeviceType,PointScalar,PointScalar>;
    const int linearPolyOrder = 1;
    BasisPtr basisForNodes = getBasis<BasisFamily>(cellTopo, FUNCTION_SPACE_HGRAD, linearPolyOrder);
    
    if (nodeOrdering_ == HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS)
    {
      // override basisForNodes for quad, hexahedron.  Apparently the lowest-order bases below are *not* in the same order as their
      // arbitrary-polynomial-order counterparts; the latter do not match the order of the shards::CellTopology nodes.
      if (cellTopo.getKey() == shards::Quadrilateral<>::key)
      {
        basisForNodes = Teuchos::rcp( new Basis_HGRAD_QUAD_C1_FEM<DeviceType,PointScalar,PointScalar>() );
      }
      else if (cellTopo.getKey() == shards::Hexahedron<>::key)
      {
        basisForNodes = Teuchos::rcp( new Basis_HGRAD_HEX_C1_FEM<DeviceType,PointScalar,PointScalar>() );
      }
    }
    
    using HostMemberLookup = ::Intrepid2::Impl::CellGeometryHostMembers<PointScalar, spaceDim, DeviceType>;
    HostMemberLookup::constructorCalled(this, cellTopo, basisForNodes);
  }
    
  // Node-based constructor for straight-edged cell geometry.
  // If claimAffine is true, we assume (without checking) that the mapping from reference space is affine.
  // (If claimAffine is false, we check whether the topology is simplicial; if so, we conclude that the mapping must be affine.)
  template<class PointScalar, int spaceDim, typename DeviceType>
  CellGeometry<PointScalar,spaceDim,DeviceType>::CellGeometry(const shards::CellTopology &cellTopo,
                                                              ScalarView<int,DeviceType> cellToNodes,
                                                              ScalarView<PointScalar,DeviceType> nodes,
                                                              const bool claimAffine,
                                                              const HypercubeNodeOrdering nodeOrdering)
  :
  nodeOrdering_(nodeOrdering),
  cellGeometryType_(FIRST_ORDER),
  cellToNodes_(cellToNodes),
  nodes_(nodes)
  {
    numCells_ = cellToNodes.extent_int(0);
    numNodesPerCell_ = cellToNodes.extent_int(1);
    INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(numNodesPerCell_ != cellTopo.getNodeCount(), std::invalid_argument, "cellToNodes.extent(1) does not match the cell topology node count");
  
    if (!claimAffine)
    {
      // if cellTopo is simplicial, then since the geometry is straight-edged, it is also affine
      const bool simplicialTopo = (cellTopo.getNodeCount() == cellTopo.getDimension() + 1);
      affine_ = simplicialTopo;
    }
    else
    {
      affine_ = true;
    }
    
    using BasisFamily = DerivedNodalBasisFamily<DeviceType,PointScalar,PointScalar>;
    const int linearPolyOrder = 1;
    BasisPtr basisForNodes = getBasis<BasisFamily>(cellTopo, FUNCTION_SPACE_HGRAD, linearPolyOrder);
    
    if (nodeOrdering_ == HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS)
    {
      // override basisForNodes for quad, hexahedron.  Apparently the lowest-order bases below are *not* in the same order as their
      // arbitrary-polynomial-order counterparts; the latter do not match the order of the shards::CellTopology nodes.
      if (cellTopo.getKey() == shards::Quadrilateral<>::key)
      {
        basisForNodes = Teuchos::rcp( new Basis_HGRAD_QUAD_C1_FEM<DeviceType,PointScalar,PointScalar>() );
      }
      else if (cellTopo.getKey() == shards::Hexahedron<>::key)
      {
        basisForNodes = Teuchos::rcp( new Basis_HGRAD_HEX_C1_FEM<DeviceType,PointScalar,PointScalar>() );
      }
    }
    
    using HostMemberLookup = ::Intrepid2::Impl::CellGeometryHostMembers<PointScalar, spaceDim, DeviceType>;
    HostMemberLookup::constructorCalled(this, cellTopo, basisForNodes);
  }

  // Constructor for higher-order geometry
  template<class PointScalar, int spaceDim, typename DeviceType>
  CellGeometry<PointScalar,spaceDim,DeviceType>::CellGeometry(Teuchos::RCP<Intrepid2::Basis<DeviceType,PointScalar,PointScalar> > basisForNodes,
                                                              ScalarView<PointScalar,DeviceType> cellNodes)
  :
  nodeOrdering_(HYPERCUBE_NODE_ORDER_TENSOR),
  cellGeometryType_(HIGHER_ORDER),
  nodes_(cellNodes)
  {
    numCells_ = cellNodes.extent_int(0);
    numNodesPerCell_ = cellNodes.extent_int(1);
    
    // if basis degree is 1, mark as first-order geometry
    const bool firstOrderGeometry = (basisForNodes->getDegree() == 1);
    cellGeometryType_ = firstOrderGeometry ? FIRST_ORDER : HIGHER_ORDER;
    
    shards::CellTopology cellTopo = basisForNodes->getBaseCellTopology();
    
    if (firstOrderGeometry && (cellTopo.getNodeCount() == spaceDim + 1)) // lowest-order and simplicial
    {
      affine_ = true;
    }
    else
    {
      affine_ = false;
    }
    using HostMemberLookup = ::Intrepid2::Impl::CellGeometryHostMembers<PointScalar, spaceDim, DeviceType>;
    HostMemberLookup::constructorCalled(this, cellTopo, basisForNodes);
  }
    
  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  bool CellGeometry<PointScalar,spaceDim,DeviceType>::affine() const
  {
    return affine_;
  }
    
  template<class PointScalar, int spaceDim, typename DeviceType>
  TensorData<PointScalar,DeviceType>
  CellGeometry<PointScalar,spaceDim,DeviceType>::allocateCellMeasure( const Data<PointScalar,DeviceType> & jacobianDet,
                                                                      const TensorData<PointScalar,DeviceType> & cubatureWeights ) const
  {
    // Output possibilities for a cubatureWeights with N components:
    // 1. For AFFINE elements (jacobianDet cell-wise constant), returns a container with N+1 tensorial components; the first component corresponds to cells
    // 2. Otherwise, returns a container with 1 tensorial component
    
    INTREPID2_TEST_FOR_EXCEPTION(cubatureWeights.rank() != 1, std::invalid_argument, "cubatureWeights container must have shape (P)");
    
    const int numTensorComponents = affine_ ? cubatureWeights.numTensorComponents() + 1 : 1;
    std::vector< Data<PointScalar,DeviceType> > tensorComponents(numTensorComponents);
    
    if (affine_)
    {
      const int cellExtent = jacobianDet.extent_int(0);
      Kokkos::Array<DataVariationType,7> cellVariationTypes {jacobianDet.getVariationTypes()[0], CONSTANT, CONSTANT, CONSTANT, CONSTANT, CONSTANT, CONSTANT};
      const int cellDataDim = jacobianDet.getDataExtent(0);
      Kokkos::Array<int,7> cellExtents{cellExtent,1,1,1,1,1,1};
      
      ScalarView<PointScalar,DeviceType> detDataView ("cell relative volumes", cellDataDim);
      tensorComponents[0] = Data<PointScalar,DeviceType>(detDataView,1,cellExtents,cellVariationTypes);
      
      for (int cubTensorComponent=0; cubTensorComponent<numTensorComponents-1; cubTensorComponent++)
      {
        auto cubatureComponent = cubatureWeights.getTensorComponent(cubTensorComponent);
        const auto cubatureExtents        = cubatureComponent.getExtents();
        const auto cubatureVariationTypes = cubatureComponent.getVariationTypes();
        const int numPoints               = cubatureComponent.getDataExtent(0);
        ScalarView<PointScalar,DeviceType> cubatureWeightView ("cubature component weights", numPoints);
        const int pointComponentRank = 1;
        tensorComponents[cubTensorComponent+1] = Data<PointScalar,DeviceType>(cubatureWeightView,pointComponentRank,cubatureExtents,cubatureVariationTypes);
      }
    }
    else
    {
      const int cellExtent = jacobianDet.extent_int(0);
      Kokkos::Array<DataVariationType,7> variationTypes {jacobianDet.getVariationTypes()[0], GENERAL, CONSTANT, CONSTANT, CONSTANT, CONSTANT, CONSTANT};
      const int cellDataDim = jacobianDet.getDataExtent(0);
      
      const int numPoints = cubatureWeights.extent_int(0);
      Kokkos::Array<int,7> extents{cellExtent,numPoints,1,1,1,1,1};
      
      ScalarView<PointScalar,DeviceType> cubatureWeightView;
      if (variationTypes[0] != CONSTANT)
      {
        cubatureWeightView = ScalarView<PointScalar,DeviceType>("cell measure", cellDataDim, numPoints);
      }
      else
      {
        cubatureWeightView = ScalarView<PointScalar,DeviceType>("cell measure", numPoints);
      }
      const int cellMeasureRank = 2;
      tensorComponents[0] = Data<PointScalar,DeviceType>(cubatureWeightView,cellMeasureRank,extents,variationTypes);
    }
    const bool separateFirstComponent = (numTensorComponents > 1);
    return TensorData<PointScalar,DeviceType>(tensorComponents, separateFirstComponent);
  }
    
  template<class PointScalar, int spaceDim, typename DeviceType>
  void CellGeometry<PointScalar,spaceDim,DeviceType>::computeCellMeasure( TensorData<PointScalar,DeviceType> &cellMeasure,
                                                                          const Data<PointScalar,DeviceType> & jacobianDet,
                                                                          const TensorData<PointScalar,DeviceType> & cubatureWeights ) const
  {
    // Output possibilities for a cubatureWeights with N components:
    // 1. For AFFINE elements (jacobianDet constant on each cell), returns a container with N+1 tensorial components; the first component corresponds to cells
    // 2. Otherwise, returns a container with 1 tensorial component
    
    INTREPID2_TEST_FOR_EXCEPTION((cellMeasure.numTensorComponents() != cubatureWeights.numTensorComponents() + 1) && (cellMeasure.numTensorComponents() != 1), std::invalid_argument,
                                 "cellMeasure must either have a tensor component count of 1 or a tensor component count that is one higher than that of cubatureWeights");
    
    INTREPID2_TEST_FOR_EXCEPTION(cubatureWeights.rank() != 1, std::invalid_argument, "cubatureWeights container must have shape (P)");
    
    if (cellMeasure.numTensorComponents() == cubatureWeights.numTensorComponents() + 1)
    {
      // affine case; the first component should contain the cell volume divided by ref cell volume; this should be stored in jacobianDet
      Kokkos::deep_copy(cellMeasure.getTensorComponent(0).getUnderlyingView1(), jacobianDet.getUnderlyingView1()); // copy point-invariant data from jacobianDet to the first tensor component of cell measure container
      const int numTensorDimensions = cubatureWeights.numTensorComponents();
      for (int i=1; i<numTensorDimensions+1; i++)
      {
        Kokkos::deep_copy(cellMeasure.getTensorComponent(i).getUnderlyingView1(), cubatureWeights.getTensorComponent(i-1).getUnderlyingView1());
      }
    }
    else
    {
      auto detVaries  = jacobianDet.getVariationTypes();
      
      const bool detCellVaries  = detVaries[0] != CONSTANT;
      const bool detPointVaries = detVaries[1] != CONSTANT;
      
      if (detCellVaries && detPointVaries)
      {
        auto cellMeasureData = cellMeasure.getTensorComponent(0).getUnderlyingView2();
        auto detData = jacobianDet.getUnderlyingView2();
        const int numCells = detData.extent_int(0);
        const int numPoints = detData.extent_int(1);
        // This lambda rewritten using single-dimension RangePolicy instead of MDRangePolicy to work around an apparent CUDA 10.1.243 compiler bug
        Kokkos::parallel_for( numCells * numPoints,
        KOKKOS_LAMBDA (const int &cellPointOrdinal) {
          const int cellOrdinal  = cellPointOrdinal / numPoints;
          const int pointOrdinal = cellPointOrdinal % numPoints;
          cellMeasureData(cellOrdinal,pointOrdinal) = detData(cellOrdinal,pointOrdinal) * cubatureWeights(pointOrdinal);
        });
      }
      else if (detCellVaries && !detPointVaries)
      {
        auto cellMeasureData = cellMeasure.getTensorComponent(0).getUnderlyingView2();
        auto detData = jacobianDet.getUnderlyingView1();
        using ExecutionSpace = typename DeviceType::execution_space;
        Kokkos::parallel_for(
        Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{detData.extent_int(0),cubatureWeights.extent_int(0)}),
        KOKKOS_LAMBDA (int cellOrdinal, int pointOrdinal) {
          cellMeasureData(cellOrdinal,pointOrdinal) = detData(cellOrdinal) * cubatureWeights(pointOrdinal);
        });
      }
      else
      {
        // constant jacobian det case
        // cell measure data has shape (P)
        auto cellMeasureData = cellMeasure.getTensorComponent(0).getUnderlyingView1();
        auto detData = jacobianDet.getUnderlyingView1();
        using ExecutionSpace = typename DeviceType::execution_space;
        Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace>(0,cellMeasureData.extent_int(0)),
        KOKKOS_LAMBDA (const int &pointOrdinal) {
          cellMeasureData(pointOrdinal) = detData(0) * cubatureWeights(pointOrdinal);
        });
      }
    }
  }
    
  template<class PointScalar, int spaceDim, typename DeviceType>
  typename CellGeometry<PointScalar,spaceDim,DeviceType>::BasisPtr
  CellGeometry<PointScalar,spaceDim,DeviceType>::basisForNodes() const
  {
    using HostMemberLookup = ::Intrepid2::Impl::CellGeometryHostMembers<PointScalar, spaceDim, DeviceType>;
    return HostMemberLookup::getBasis(this);
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  const shards::CellTopology & CellGeometry<PointScalar,spaceDim,DeviceType>::cellTopology() const
  {
    using HostMemberLookup = ::Intrepid2::Impl::CellGeometryHostMembers<PointScalar, spaceDim, DeviceType>;
    return HostMemberLookup::getCellTopology(this);
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  DataVariationType CellGeometry<PointScalar,spaceDim,DeviceType>::cellVariationType() const
  {
    if (cellGeometryType_ == UNIFORM_GRID)
    {
      const int numSubdivisions = numCellsPerGridCell(subdivisionStrategy_);
      if (numSubdivisions == 1)
      {
        return CONSTANT;
      }
      else
      {
        return MODULAR;
      }
    }
    else return GENERAL;
  }

  template<class PointScalar, int spaceDim, typename DeviceType>
  Data<PointScalar,DeviceType> CellGeometry<PointScalar,spaceDim,DeviceType>::getJacobianRefData(const TensorPoints<PointScalar,DeviceType> &points) const
  {
    Data<PointScalar,DeviceType> emptyRefData;
    if (cellGeometryType_ == UNIFORM_GRID)
    {
      // no need for basis computations
      return emptyRefData;
    }
    else if (cellGeometryType_ == TENSOR_GRID)
    {
      // no need for basis values
      return emptyRefData;
    }
    else if ((cellGeometryType_ == FIRST_ORDER) || (cellGeometryType_ == HIGHER_ORDER))
    {
      const bool simplex = (spaceDim + 1 == cellToNodes_.extent_int(1));
      if (simplex)
      {
        // no need for precomputed basis values
        return emptyRefData;
      }
      else
      {
        auto basisForNodes = this->basisForNodes();
        
        if (affine_)
        {
          // no need for precomputed basis values
          return emptyRefData;
        }
        else
        {
          auto basisGradients = basisForNodes->allocateBasisValues(points, OPERATOR_GRAD);
          basisForNodes->getValues(basisGradients, points, OPERATOR_GRAD);
          
          int numPoints = points.extent_int(0);
          int numFields = basisForNodes->getCardinality();
          
          // At some (likely relatively small) memory cost, we copy the BasisGradients into an explicit (F,P,D) container.
          // Given that we expect to reuse this for a non-trivial number of cell in the common use case, the extra memory
          // cost is likely worth the increased flop count, etc.  (One might want to revisit this in cases of high spaceDim
          // and/or very high polynomial order.)
          
          auto firstPointComponentView = points.getTensorComponent(0); // (P,D0)
          auto basisGradientsView = getMatchingViewWithLabel(firstPointComponentView, "CellGeometryProvider: temporary basisGradients", numFields, numPoints, spaceDim);
          
          using ExecutionSpace = typename DeviceType::execution_space;
          auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{numFields,numPoints,spaceDim});
          
          Kokkos::parallel_for("copy basis gradients", policy,
          KOKKOS_LAMBDA (const int &fieldOrdinal, const int &pointOrdinal, const int &d) {
            basisGradientsView(fieldOrdinal,pointOrdinal,d) = basisGradients(fieldOrdinal,pointOrdinal,d);
          });
          ExecutionSpace().fence();
          
          Data<PointScalar,DeviceType> basisRefData(basisGradientsView);
          return basisRefData;
        }
      }
    }
    else
    {
      // TODO: handle the other cases
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "support for this CellGeometryType is not yet implemented");
    }
    return emptyRefData;
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  int CellGeometry<PointScalar,spaceDim,DeviceType>::hypercubeComponentNodeNumber(int hypercubeNodeNumber, int d) const
  {
    if (nodeOrdering_ == HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS)
    {
      // note that Shards numbers nodes for quad counter-clockwise
      // cube is tensor-product of the (x,y) quad with a z line segment
      if (d==0)
      {
        if ((hypercubeNodeNumber % 4 == 1) || (hypercubeNodeNumber % 4 == 2))
          return 1;
        else
          return 0;
      }
      else if (d==1)
      {
        if ((hypercubeNodeNumber % 4 == 2) || (hypercubeNodeNumber % 4 == 3))
          return 1;
        else
          return 0;
      }
    }
    // tensor formula coincides with shards formula for d ≥ 2
    const int nodesForPriorDimensions = 1 << d;
    if ((hypercubeNodeNumber / nodesForPriorDimensions) % 2 == 1)
      return 1;
    else
      return 0;
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  void CellGeometry<PointScalar,spaceDim,DeviceType>::initializeOrientations()
  {
    using HostExecSpace = typename Kokkos::Impl::is_space<DeviceType>::host_mirror_space::execution_space ;
    
    const bool isGridType = (cellGeometryType_ == TENSOR_GRID) || (cellGeometryType_ == UNIFORM_GRID);
    const int numOrientations = isGridType ? numCellsPerGridCell(subdivisionStrategy_) : numCells();
    
    const int nodesPerCell = numNodesPerCell();
    
    ScalarView<Orientation, DeviceType> orientationsView("orientations", numOrientations);
    auto orientationsHost = Kokkos::create_mirror_view(typename HostExecSpace::memory_space(), orientationsView);
    
    ScalarView<PointScalar, HostExecSpace> cellNodesHost("cellNodesHost",numOrientations,nodesPerCell); // (C,N) -- where C = numOrientations
    
    DataVariationType cellVariationType;
    
    if (isGridType)
    {
      // then there are as many distinct orientations possible as there are there are cells per grid cell
      // fill cellNodesHost with sample nodes from grid cell 0
      
      const int numSubdivisions = numCellsPerGridCell(subdivisionStrategy_); // can be up to 6
      const int gridCellOrdinal = 0;
      auto hostPolicy = Kokkos::MDRangePolicy<HostExecSpace,Kokkos::Rank<2>>({0,0},{numSubdivisions,nodesPerCell});
      Kokkos::parallel_for("fill cellNodesHost", hostPolicy,
      KOKKOS_LAMBDA (const int &subdivisionOrdinal, const int &nodeInCell) {
        auto node = gridCellNodeForSubdivisionNode(gridCellOrdinal, subdivisionOrdinal, nodeInCell);
        cellNodesHost(subdivisionOrdinal,nodeInCell) = node;
      });
      
      cellVariationType = (numSubdivisions == 1) ? CONSTANT : MODULAR;
    }
    else
    {
      cellVariationType = GENERAL;
      auto cellToNodesHost = Kokkos::create_mirror_view_and_copy(typename HostExecSpace::memory_space(), cellToNodes_);
    }
    
    OrientationTools<HostExecSpace>::getOrientation(orientationsHost,cellNodesHost,this->cellTopology());
    Kokkos::deep_copy(orientationsView,orientationsHost);
    
    const int orientationsRank = 1; // shape (C)
    const Kokkos::Array<int,7>               orientationExtents {static_cast<int>(numCells_),1,1,1,1,1,1};
    const Kokkos::Array<DataVariationType,7> orientationVariationTypes { cellVariationType, CONSTANT, CONSTANT, CONSTANT, CONSTANT, CONSTANT, CONSTANT};
    orientations_ = Data<Orientation,DeviceType>(orientationsView, orientationsRank, orientationExtents, orientationVariationTypes);
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  size_t CellGeometry<PointScalar,spaceDim,DeviceType>::extent(const int& r) const {
    if (r == 0)
    {
      return numCells_;
    }
    else if (r == 1)
    {
      return numNodesPerCell_;
    }
    else if (r == 2)
    {
      return spaceDim;
    }
    else
    {
      return 1;
    }
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  template <typename iType>
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if<std::is_integral<iType>::value, int>::type
  CellGeometry<PointScalar,spaceDim,DeviceType>::extent_int(const iType& r) const
  {
    return static_cast<int>(extent(r));
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  typename CellGeometry<PointScalar,spaceDim,DeviceType>::HypercubeNodeOrdering
  CellGeometry<PointScalar,spaceDim,DeviceType>::nodeOrderingForHypercubes() const
  {
    return nodeOrdering_;
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  int CellGeometry<PointScalar,spaceDim,DeviceType>::numCells() const
  {
    return numCells_;
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  int CellGeometry<PointScalar,spaceDim,DeviceType>::numCellsInDimension(const int &dim) const
  {
    if (cellGeometryType_ == UNIFORM_GRID)
    {
      return gridCellCounts_[dim];
    }
    else if (cellGeometryType_ == TENSOR_GRID)
    {
      return tensorVertices_.extent_int(dim);
    }
    else
    {
      return -1; // not valid for this cell geometry type
    }
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  int CellGeometry<PointScalar,spaceDim,DeviceType>::numNodesPerCell() const
  {
    return numNodesPerCell_;
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  Orientation CellGeometry<PointScalar,spaceDim,DeviceType>::getOrientation(int &cellNumber) const
  {
    INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(!orientations_.isValid(), std::invalid_argument, "orientations_ not initialized; call initializeOrientations() first");
    return orientations_(cellNumber);
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  Data<Orientation,DeviceType> CellGeometry<PointScalar,spaceDim,DeviceType>::getOrientations()
  {
    if (!orientations_.isValid())
    {
      initializeOrientations();
    }
    return orientations_;
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  PointScalar CellGeometry<PointScalar,spaceDim,DeviceType>::gridCellCoordinate(const int &gridCellOrdinal, const int &localNodeNumber, const int &dim) const
  {
    const int componentNode = hypercubeComponentNodeNumber(localNodeNumber, dim);
    int cellCountForPriorDimensions = 1;
    for (int d=0; d<dim; d++)
    {
      cellCountForPriorDimensions *= numCellsInDimension(d);
    }
    const int componentGridCellOrdinal = (gridCellOrdinal / cellCountForPriorDimensions) % numCellsInDimension(dim);
    const int vertexOrdinal = componentGridCellOrdinal + componentNode;
    if (cellGeometryType_ == UNIFORM_GRID)
    {
      return origin_[dim] + (vertexOrdinal * domainExtents_[dim]) / gridCellCounts_[dim];
    }
    else if (cellGeometryType_ == TENSOR_GRID)
    {
      Kokkos::Array<int,spaceDim> pointOrdinalComponents;
      for (int d=0; d<spaceDim; d++)
      {
        pointOrdinalComponents[d] = 0;
      }
      pointOrdinalComponents[dim] = vertexOrdinal;
      return tensorVertices_(pointOrdinalComponents,dim);
    }
    else
    {
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "Unsupported geometry type");
      return 0; // unreachable; here to avoid compiler warnings
    }
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  unsigned CellGeometry<PointScalar,spaceDim,DeviceType>::rank() const
  {
    return 3; // (C,N,D)
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  int CellGeometry<PointScalar,spaceDim,DeviceType>::gridCellNodeForSubdivisionNode(const int &gridCellOrdinal, const int &subdivisionOrdinal,
                                                                                       const int &subdivisionNodeNumber) const
  {
    // TODO: do something to reuse the nodeLookup containers
    switch (subdivisionStrategy_)
    {
      case NO_SUBDIVISION:
        return subdivisionNodeNumber;
      case TWO_TRIANGLES_RIGHT:
      case TWO_TRIANGLES_LEFT:
      case FOUR_TRIANGLES:
      {
        Kokkos::Array<int,3> nodeLookup;
        if (subdivisionStrategy_ == TWO_TRIANGLES_RIGHT)
        {
          if (subdivisionOrdinal == 0)
          {
            // bottom-right cell: node numbers coincide with quad node numbers
            nodeLookup = {0,1,2};
          }
          else if (subdivisionOrdinal == 1)
          {
            // node 0 --> node 2
            // node 1 --> node 3
            // node 2 --> node 0
            nodeLookup = {2,3,0};
          }
          else
          {
            INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "Unsupported subdivision ordinal");
          }
        }
        else if (subdivisionStrategy_ == TWO_TRIANGLES_LEFT)
        {
          if (subdivisionOrdinal == 0)
          {
            // bottom-left cell:
            // node 0 --> node 3
            // node 1 --> node 0
            // node 2 --> node 1
            nodeLookup = {3,0,1};
          }
          else if (subdivisionOrdinal == 1)
          {
            // node 0 --> node 2
            // node 1 --> node 3
            // node 2 --> node 0
            nodeLookup = {2,3,0};
          }
          else
          {
            INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "Unsupported subdivision ordinal");
          }
        }
        else // FOUR_TRIANGLES
        {
          // counter-clockwise, bottom triangle first
          // bottom triangle goes:
          // 0 --> 1
          // 1 --> center
          // 2 --> 0
          // and similarly for the other triangles, proceeding counter-clockwise
          // node 1 always being the center, we special-case that
          if (subdivisionNodeNumber == 1)
          {
            // center coordinate: we call this node 4 in the quadrilateral
            return 4;
          }
          else
          {
            nodeLookup = {(subdivisionOrdinal + 1) % 4, -1, subdivisionOrdinal};
          }
        }
        const int gridCellNodeNumber = nodeLookup[subdivisionNodeNumber];
        return gridCellNodeNumber;
      }
      case FIVE_TETRAHEDRA:
      case SIX_TETRAHEDRA:
      {
        Kokkos::Array<int,4> nodeLookup;
        if (subdivisionStrategy_ == FIVE_TETRAHEDRA)
        {
          /*
                  // to discretize a unit cube into 5 tetrahedra, we can take the four vertices
                  // (1,1,1)
                  // (0,0,1)
                  // (0,1,0)
                  // (1,0,0)
                  // as an interior tetrahedron.  Call this cell 0.  The remaining 4 cells can be determined
                  // by selecting three of the above points (there are exactly 4 such combinations) and then selecting
                  // from the remaining four vertices of the cube the one nearest the plane defined by those three points.
                  // The remaining four vertices are:
                  // (0,0,0)
                  // (1,1,0)
                  // (1,0,1)
                  // (0,1,1)
                  // For each of these four, we pick one, and then take the three nearest vertices from cell 0 to form a new tetrahedron.
                  // We enumerate as follows:
                  // cell 0: (1,1,1), (0,0,1), (0,1,0), (1,0,0)
                  // cell 1: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
                  // cell 2: (1,1,0), (1,1,1), (0,1,0), (1,0,0)
                  // cell 3: (1,0,1), (1,1,1), (0,0,1), (1,0,0)
                  // cell 4: (0,1,1), (1,1,1), (0,0,1), (0,1,0)
                  */
          // tetrahedra are as follows:
          // 0: {1,3,4,6}
          // 1: {0,1,3,4}
          // 2: {1,2,3,6}
          // 3: {1,4,5,6}
          // 4: {3,4,6,7}
          switch (subdivisionOrdinal) {
            case 0:
              nodeLookup = {1,3,4,6};
              break;
            case 1:
              nodeLookup = {0,1,3,4};
              break;
            case 2:
              nodeLookup = {1,2,3,6};
              break;
            case 3:
              nodeLookup = {1,4,5,6};
              break;
            case 4:
              nodeLookup = {3,4,6,7};
              break;
            default:
              INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "invalid subdivisionOrdinal");
              break;
          }
        }
        else if (subdivisionStrategy_ == SIX_TETRAHEDRA)
        {
          INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "support for SIX_TETRAHEDRA not yet implemented");
        }
        const int gridCellNodeNumber = nodeLookup[subdivisionNodeNumber];
        return gridCellNodeNumber;
      }
      default:
        INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(true, std::invalid_argument, "Subdivision strategy not yet implemented!");
        // some compilers complain about missing return
        return 0; // statement should be unreachable...
    }
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  PointScalar CellGeometry<PointScalar,spaceDim,DeviceType>::subdivisionCoordinate(const int &gridCellOrdinal, const int &subdivisionOrdinal,
                                                                                      const int &subdivisionNodeNumber, const int &d) const
  {
    int gridCellNode = gridCellNodeForSubdivisionNode(gridCellOrdinal, subdivisionOrdinal, subdivisionNodeNumber);
    
    if (subdivisionStrategy_ == FOUR_TRIANGLES)
    {
      // this is the one case in which the gridCellNode may not actually be a node in the grid cell
      if (gridCellNode == 4) // center vertex
      {
        // d == 0 means quad vertices 0 and 1 suffice;
        // d == 1 means quad vertices 0 and 3 suffice
        const int gridVertex0 = 0;
        const int gridVertex1 = (d == 0) ? 1 : 3;
        return 0.5 * (gridCellCoordinate(gridCellOrdinal, gridVertex0, d) + gridCellCoordinate(gridCellOrdinal, gridVertex1, d));
      }
    }
    return gridCellCoordinate(gridCellOrdinal, gridCellNode, d);
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  PointScalar
  CellGeometry<PointScalar,spaceDim,DeviceType>::operator()(const int& cell, const int& node, const int& dim) const {
    if ((cellGeometryType_ == UNIFORM_GRID) || (cellGeometryType_ == TENSOR_GRID))
    {
      const int numSubdivisions = numCellsPerGridCell(subdivisionStrategy_);
      if (numSubdivisions == 1)
      {
        // hypercube
        return gridCellCoordinate(cell, node, dim);
      }
      else
      {
        const int subdivisionOrdinal = cell % numSubdivisions;
        const int gridCellOrdinal    = cell / numSubdivisions;
        return subdivisionCoordinate(gridCellOrdinal, subdivisionOrdinal, node, dim);
      }
    }
    else
    {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE((cell < 0),                                     std::invalid_argument, "cell out of bounds");
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(static_cast<unsigned>(cell) > numCells_,        std::invalid_argument, "cell out of bounds");
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE((node < 0),                                     std::invalid_argument, "node out of bounds");
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(static_cast<unsigned>(node) > numNodesPerCell_, std::invalid_argument, "node out of bounds");
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE((dim < 0),                                      std::invalid_argument, "dim out of bounds" );
      INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(dim > spaceDim,                                 std::invalid_argument, "dim out of bounds" );
#endif
      if (cellToNodes_.is_allocated())
      {
        const int nodeNumber = cellToNodes_(cell,node);
        return nodes_(nodeNumber,dim);
      }
      else
      {
        return nodes_(cell,node,dim);
      }
    }
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  KOKKOS_INLINE_FUNCTION
  int CellGeometry<PointScalar,spaceDim,DeviceType>::uniformJacobianModulus() const
  {
    if (cellGeometryType_ == UNIFORM_GRID)
    {
      return numCellsPerGridCell(subdivisionStrategy_);
    }
    else
    {
      return numCells_;
    }
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  Data<PointScalar,DeviceType> CellGeometry<PointScalar,spaceDim,DeviceType>::allocateJacobianData(const TensorPoints<PointScalar,DeviceType> &points, const int startCell, const int endCell) const
  {
    const int pointsPerCell = points.extent_int(0);
    return allocateJacobianDataPrivate(points,pointsPerCell,startCell,endCell);
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  Data<PointScalar,DeviceType> CellGeometry<PointScalar,spaceDim,DeviceType>::allocateJacobianData(const ScalarView<PointScalar,DeviceType> &points, const int startCell, const int endCell) const
  {
    // if points is rank 3, it has shape (C,P,D).  If it's rank 2, (P,D).
    const int pointDimension = (points.rank() == 3) ? 1 : 0;
    const int pointsPerCell = points.extent_int(pointDimension);
    TensorPoints<PointScalar,DeviceType> tensorPoints(points);
    return allocateJacobianDataPrivate(tensorPoints,pointsPerCell,startCell,endCell);
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  Data<PointScalar,DeviceType> CellGeometry<PointScalar,spaceDim,DeviceType>::allocateJacobianData(const int &numPoints, const int startCell, const int endCell) const
  {
    INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(!affine_, std::invalid_argument, "this version of allocateJacobianData() is only supported for affine CellGeometry");
    
    TensorPoints<PointScalar,DeviceType> emptyPoints;
    return allocateJacobianDataPrivate(emptyPoints,numPoints,startCell,endCell);
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  void CellGeometry<PointScalar,spaceDim,DeviceType>::setJacobian(Data<PointScalar,DeviceType> &jacobianData, const TensorPoints<PointScalar,DeviceType> &points,
                                                                  const Data<PointScalar,DeviceType> &refData, const int startCell, const int endCell) const
  {
    const int pointsPerCell = points.extent_int(0);
    setJacobianDataPrivate(jacobianData,points,pointsPerCell,refData,startCell,endCell);
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  void CellGeometry<PointScalar,spaceDim,DeviceType>::setJacobian(Data<PointScalar,DeviceType> &jacobianData, const ScalarView<PointScalar,DeviceType> &points,
                                                                  const Data<PointScalar,DeviceType> &refData, const int startCell, const int endCell) const
  {
    // if points is rank 3, it has shape (C,P,D).  If it's rank 2, (P,D).
    const int pointDimension = (points.rank() == 3) ? 1 : 0;
    const int pointsPerCell = points.extent_int(pointDimension);
    TensorPoints<PointScalar,DeviceType> tensorPoints(points);
    setJacobianDataPrivate(jacobianData,tensorPoints,pointsPerCell,refData,startCell,endCell);
  }
  
  template<class PointScalar, int spaceDim, typename DeviceType>
  void CellGeometry<PointScalar,spaceDim,DeviceType>::setJacobian(Data<PointScalar,DeviceType> &jacobianData, const int &numPoints, const int startCell, const int endCell) const
  {
    INTREPID2_TEST_FOR_EXCEPTION_DEVICE_SAFE(!affine_, std::invalid_argument, "this version of setJacobian() is only supported for affine CellGeometry");
    
    TensorPoints<PointScalar,DeviceType> emptyPoints;
    Data<PointScalar,DeviceType> emptyRefData;
    setJacobianDataPrivate(jacobianData,emptyPoints,numPoints,emptyRefData,startCell,endCell);
  }
} // namespace Intrepid2

#endif /* Intrepid2_CellGeometryDef_h */

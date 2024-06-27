// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CellGeometry.hpp
    \brief  Allows definition of cell geometry information, including uniform and curvilinear mesh definition, in a manner that exposes several types of mesh regularity, which in turn enables various optimizations.

    \author Nathan V. Roberts
*/

#ifndef Intrepid2_CellGeometry_h
#define Intrepid2_CellGeometry_h

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_Data.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"
#include "Intrepid2_Orientation.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_TensorData.hpp"
#include "Intrepid2_TensorPoints.hpp"
#include "Intrepid2_TransformedBasisValues.hpp"
#include "Intrepid2_Utils.hpp"
#include "Intrepid2_VectorData.hpp"

#include "Intrepid2_ScalarView.hpp"

namespace Intrepid2
{
  /**
   \class  Intrepid2::CellGeometry
   \brief CellGeometry provides the nodes for a set of cells; has options that support efficient definition of uniform grids as well as options for arbitrary geometry, including curvilinear.
   
   CellGeometry implements allocate and set methods for cell Jacobians, their inverses, and their determinants.  These are computed and stored in a fashion that avoids redundant computation and storage when CellGeometry parameters imply structure, either because the mapping is affine, or because there is axis-aligned tensorial structure that implies a diagonal Jacobian.
   
   CellGeometry also handles orientations for the cells, accessible through the getOrientation() and getOrientations().  These are backed by lazily-initialized storage of orientations for all cells.
      
   \note Conceptually, it would make sense to use class inheritance and have different member data for each type of geometry supported.  We instead glom all the options together into one multi-modal class; this is basically to avoid certain difficulties with vtables under CUDA.
   */
  template<class PointScalar, int spaceDim, typename DeviceType>
  class CellGeometry
  {
    public:
    /*! Supported types for CellGeometry */
    enum CellGeometryType {
        UNIFORM_GRID   /// each grid division has the same dimensions
      , TENSOR_GRID    /// grid expressed as a Cartesian product of 1D grids  (could be a Shishkin mesh, e.g.)
      , EXTRUDED_GRID  /// lower-dimensional geometry that is orthogonally extruded in higher dimensions
      , FIRST_ORDER    /// geometry expressible in terms of vertices of the cell
      , HIGHER_ORDER   /// geometry expressible in terms of a higher-order basis (must be specified)
    };
    
    /*! Options for uniform, tensor grids subdivided into simplices */
    enum SubdivisionStrategy
    {
        NO_SUBDIVISION      /// no subdivision
      , TWO_TRIANGLES_RIGHT /// square --> two triangles, with a hypotenuse of slope  1
      , TWO_TRIANGLES_LEFT  /// square --> two triangles, with a hypotenuse of slope -1
      , FOUR_TRIANGLES      /// square --> four triangles, with a new vertex at center
      , FIVE_TETRAHEDRA     /// cube   --> five tetrahedra
      , SIX_TETRAHEDRA      /// cube   --> six tetrahedra
      , SIX_PYRAMIDS        /// cube   --> six pyramids
    };
    
    /*! Distinguish between "classic" (counterclockwise in 2D, counterclockwise bottom followed by counterclockwise top in 3D) Shards ordering and a more natural tensor ordering.  The tensor ordering is such that the lowest-order bit in the composite vertex number corresponds to the x dimension, the second-lowest-order bit corresponds to the y dimension, etc.  (There is not yet a non-"classic" Shards ordering, but we hope to add one.) */
    enum HypercubeNodeOrdering
    {
        HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS /// classic shards ordering
      , HYPERCUBE_NODE_ORDER_TENSOR         /// a more natural tensor ordering
    };
    
    /** \brief  Helper method that returns the number of cells into which each grid cell will be subdivided based on the input subdivision strategy.
       \param [in] subdivisionStrategy - the subdivision strategy.
       \return the number of subdivisions specified by the subdivision strategy; -1 for unsupported subdivisionStrategy values.
    */
    KOKKOS_INLINE_FUNCTION
    int numCellsPerGridCell(SubdivisionStrategy subdivisionStrategy) const;
    
  public:
    /** \brief  Notionally-private method that provides a common interface for multiple public-facing allocateJacobianData() methods.  (Marked as public due to compiler constraints.)
       \param [in] pointComponentView - typically either the first component of a TensorPoints container, or a View containing all the points, but may be empty.  Only used for scalar-type-size-matched construction of new views (important for views of Fad types with hidden dimensions).
       \param [in] pointsPerCell - the number of points at which the Jacobian will be evaluated in each cell.  If points is a valid container, pointsPerCell must match its first dimension.
       \param [in] startCell - the first cell ordinal for which the Jacobian will be evaluated (used to define worksets smaller than the whole CellGeometry).
       \param [in] endCell - the first cell ordinal for which the Jacobian will not be evaluated (used to define worksets smaller than the whole CellGeometry); use -1 to indicate that all cells from startCell on should be evaluated.
       \return a Data object appropriately sized to accomodate the specified Jacobian values.
    */
    Data<PointScalar,DeviceType> allocateJacobianDataPrivate(const ScalarView<PointScalar,DeviceType> &pointComponentView, const int &pointsPerCell, const int startCell, const int endCell) const;
    
    /** \brief  Notionally-private method that provides a common interface for multiple public-facing setJacobianData() methods.  (Marked as public due to compiler constraints.)
       \param [out] jacobianData - a container, allocated by allocateJacobianData(), into which the evaluated Jacobians will be placed.
       \param [in] pointsPerCell - the number of points at which the Jacobian will be evaluated in each cell.  If points is a valid container, pointsPerCell must match its first dimension.
       \param [in] refData - the return from getJacobianRefData(); may be an empty container, depending on details of CellGeometry (e.g. if it is affine)
       \param [in] startCell - the first cell ordinal for which the Jacobian will be evaluated (used to define worksets smaller than the whole CellGeometry).
       \param [in] endCell - the first cell ordinal for which the Jacobian will not be evaluated (used to define worksets smaller than the whole CellGeometry); use -1 to indicate that all cells from startCell on should be evaluated.
    */
    void setJacobianDataPrivate(Data<PointScalar,DeviceType> &jacobianData, const int &pointsPerCell, const Data<PointScalar,DeviceType> &refData, const int startCell, const int endCell) const;
  protected:
    HypercubeNodeOrdering nodeOrdering_;
    CellGeometryType    cellGeometryType_;
    SubdivisionStrategy subdivisionStrategy_ = NO_SUBDIVISION;
    bool affine_; // if true, each cell has constant Jacobian across the cell
    Data<Orientation, DeviceType> orientations_; // for grid types, this could have either a single entry or one matching numCellsPerGridCell().  For other types, it has as many entries as there are cells.
    
    // uniform grid data -- used for UNIFORM_GRID type
    Kokkos::Array<PointScalar,spaceDim> origin_;         // point specifying a corner of the mesh
    Kokkos::Array<PointScalar,spaceDim> domainExtents_;  // how far the domain extends in each dimension
    Kokkos::Array<int,spaceDim>         gridCellCounts_; // how many grid cells wide the mesh is in each dimension
    
    // tensor grid data -- only used for TENSOR_GRID type
    TensorPoints<PointScalar, DeviceType> tensorVertices_;
    
    // arbitrary cell node data, used for both higher-order and first-order
    // (here, nodes are understood as geometry degrees of freedom)
    ScalarView<int,DeviceType>         cellToNodes_; // (C,N) -- N is the number of nodes per cell; values are global node ordinals
    ScalarView<PointScalar,DeviceType> nodes_;       // (GN,D) or (C,N,D) -- GN is the number of global nodes; (C,N,D) used only if cellToNodes_ is empty.
    using BasisPtr = Teuchos::RCP< Basis<DeviceType,PointScalar,PointScalar> >;
    
    unsigned numCells_        = 0;
    unsigned numNodesPerCell_ = 0;
  public:
    /** \brief  Uniform grid constructor, with optional subdivision into simplices.
        \param [in] origin - location of the "origin" of the mesh, coordinates start here; node locations are based on the sum of origin and domain extents, subdivided according to cell counts in each dimension.
        \param [in] domainExtents - the size of the domain in each coordinate dimension.
        \param [in] gridCellCounts - the number of grid cells in each coordinate dimension.
        \param [in] subdivisionStrategy - whether (and how) to subdivide the (hypercube) grid elements into simplices.
        \param [in] nodeOrdering - applicable for hypercube cell topologies; specifies whether to use the order used by Shards (and lowest-order Intrepid2 bases), or the one used by higher-order Intrepid2 bases.
    */
    CellGeometry(const Kokkos::Array<PointScalar,spaceDim> &origin,
                 const Kokkos::Array<PointScalar,spaceDim> &domainExtents,
                 const Kokkos::Array<int,spaceDim> &gridCellCounts,
                 SubdivisionStrategy subdivisionStrategy = NO_SUBDIVISION,
                 HypercubeNodeOrdering nodeOrdering = HYPERCUBE_NODE_ORDER_TENSOR);
    
    /** \brief  Node-based constructor for straight-edged geometry. 
        \param [in] cellTopo - the cell topology for all cells.
        \param [in] cellToNodes - (C,LN) container specifying the global index of the local node.
        \param [in] nodes - (GN,D) container specifying the coordinate weight for the global node in the specified dimension; if cellToNodes is not allocated, this must be a (C,N,D) container
        \param [in] claimAffine - whether to assume (without checking) that the mapping from reference space is affine.  (If claimAffine is false, we check whether the cell topology is simplicial, and if so, set the affine_ member variable to true.)
        \param [in] nodeOrdering - applicable for hypercube cell topologies; specifies whether to use the order used by Shards (and lowest-order Intrepid2 bases), or the one used by higher-order Intrepid2 bases.
    */
    CellGeometry(const shards::CellTopology &cellTopo,
                 ScalarView<int,DeviceType> cellToNodes,
                 ScalarView<PointScalar,DeviceType> nodes,
                 const bool claimAffine = false,
                 const HypercubeNodeOrdering nodeOrdering = HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS);

    /** \brief  Node-based constructor for curvilinear geometry.
        \param [in] basisForNodes - the basis in terms of which the reference-to-physical transformation is expressed.
        \param [in] cellNodes - (C,F,D) container specifying the coordinate weight of the node; F is the (local) field ordinal for the basis in terms of which the reference-to-physical transformation is expressed.
     */
    CellGeometry(Teuchos::RCP<Intrepid2::Basis<DeviceType,PointScalar,PointScalar> > basisForNodes,
                 ScalarView<PointScalar,DeviceType> cellNodes);
    
    /** \brief  Copy constructor.
        \param [in] cellGeometry - the object being copied.
     */
    KOKKOS_INLINE_FUNCTION CellGeometry(const CellGeometry &cellGeometry);
    
    /** \brief  Destructor.
    */
    KOKKOS_INLINE_FUNCTION ~CellGeometry();
    
    //! Returns true if Jacobian is constant within each cell
    KOKKOS_INLINE_FUNCTION
    bool affine() const;
    
    /** \brief  Allocate a TensorData object appropriate for passing to computeCellMeasure().
       \param [in] jacobianDet - a container approriately sized to store the determinants of the Jacobians of the reference-to-physical mapping.
       \param [in] cubatureWeights - a container approriately sized to store the quadrature weights.
       \return TensorData object sized to accept the result of computeCellMeasure() when called with the provided jacobianDet and cubatureWeights arguments.
    */
    TensorData<PointScalar,DeviceType> allocateCellMeasure( const Data<PointScalar,DeviceType> & jacobianDet, const TensorData<PointScalar,DeviceType> & cubatureWeights ) const;
    
    /** \brief  Compute cell measures that correspond to provided Jacobian determinants and
       \param [out] cellMeasure - a container, usually provided by allocateCellMeasure(), to store the cell measures.
       \param [in] jacobianDet -  the determinants of the Jacobians of the reference-to-physical mapping.
       \param [in] cubatureWeights - the quadrature weights.
    */
    void computeCellMeasure( TensorData<PointScalar,DeviceType> &cellMeasure, const Data<PointScalar,DeviceType> & jacobianDet, const TensorData<PointScalar,DeviceType> & cubatureWeights ) const;
    
    //! H^1 Basis used in the reference-to-physical transformation.  Linear for straight-edged geometry; higher-order for curvilinear.
    BasisPtr basisForNodes() const;
    
    //! The shards CellTopology for each cell within the CellGeometry object.  Note that this is always a lowest-order CellTopology, even for higher-order curvilinear geometry.  Higher-order geometry for CellGeometry is expressed in terms of the basis returned by basisForNodes().
    const shards::CellTopology & cellTopology() const;
    
    //! Returns a global node number corresponding to the specified (cell, node) combination.  If uniform grid (possibly subdivided), this number is computed dynamically; for more general meshes, this simply returns the result of a lookup in the cellToNodes container provided at construction.
    KOKKOS_INLINE_FUNCTION
    ordinal_type cellToNode(ordinal_type cellIndex, ordinal_type cellNodeNumber) const;
    
    //! Indicates the type of geometric variation from one cell to the next.  If all cells are known to have the same geometry, returns CONSTANT.
    //! Returns CONSTANT for uniform grids with no subdivisions, MODULAR for uniform grids with subdivisions, GENERAL for all others.
    KOKKOS_INLINE_FUNCTION
    DataVariationType cellVariationType() const;
    
    //! For hypercube vertex number hypercubeNodeNumber, returns the component node number in specified dimension d.  This is either 0 or 1.  (The mapping depends on the node ordering specified in nodeOrdering_.)
    KOKKOS_INLINE_FUNCTION
    int hypercubeComponentNodeNumber(int hypercubeNodeNumber, int d) const;
    
    //! Initialize the internal orientations_ member with the orientations of each member cell.  These are used for projected geometry, and have a logical shape (C).
    void initializeOrientations();
    
    //! Returns the logical extent of the container in the specified dimension; the shape of CellGeometry is always (C,N,D), where C is the number of cells, N is the number of nodes per cell (may be more than the number of vertices, in the case of curvilinear geometry), and D is the spatial dimension.
    KOKKOS_INLINE_FUNCTION
    size_t extent(const int& r) const;
    
    //! Returns the logical extent of the container in the specified dimension as an int; the shape of CellGeometry is always (C,N,D), where C is the number of cells, N is the number of nodes per cell (may be more than the number of vertices, in the case of curvilinear geometry), and D is the spatial dimension.
    template <typename iType>
    KOKKOS_INLINE_FUNCTION
    typename std::enable_if<std::is_integral<iType>::value, int>::type
    extent_int(const iType& r) const;
    
    //! Returns the node ordering used for hypercubes.
    KOKKOS_INLINE_FUNCTION
    HypercubeNodeOrdering nodeOrderingForHypercubes() const;
    
    //! Returns the number of cells.
    KOKKOS_INLINE_FUNCTION
    int numCells() const;
    
    //! For uniform grid and tensor grid CellGeometry, returns the number of cells in the specified component dimension.  For other CellGeometry types, returns -1.
    KOKKOS_INLINE_FUNCTION
    int numCellsInDimension(const int &dim) const;
    
    //! Returns the number of nodes per cell; may be more than the number of vertices in the corresponding CellTopology, in the case of higher-order (curvilinear) geometry.
    KOKKOS_INLINE_FUNCTION
    int numNodesPerCell() const;
    
    //! Returns the orientation for the specified cell.  Requires that initializeOrientations() has been called previously.
    KOKKOS_INLINE_FUNCTION
    Orientation getOrientation(int &cellNumber) const;
    
    //! Returns the orientations for all cells.  Calls initializeOrientations() if it has not previously been called.
    Data<Orientation,DeviceType> getOrientations();
    
    //! \brief Fills the provided container with the orientations for the specified cell range.  Calls getOrientations() and copies the orientations from the Data container into the ScalarView container.
    //! \param [out] orientationsView - the container that will be filled.
    //! \param [in] startCell - the first cell ordinal whose orientation will be copied.
    //! \param [in] endCell - the first cell ordinal whose orientation will not be copied; use -1 to indicate that orientations for all cells from startCell on should be copied.
    void orientations(ScalarView<Orientation,DeviceType> orientationsView, const int &startCell = 0, const int &endCell = -1);
    
    //! returns coordinate in dimension dim of the indicated node in the indicated grid cell
    KOKKOS_INLINE_FUNCTION
    PointScalar gridCellCoordinate(const int &gridCellOrdinal, const int &localNodeNumber, const int &dim) const;
    
    //! Returns the logical rank of this container.  This is always 3.
    KOKKOS_INLINE_FUNCTION
    unsigned rank() const;
    
    //! returns coordinate in dimension d for the indicated subdivision of the indicated grid cell
    KOKKOS_INLINE_FUNCTION
    int gridCellNodeForSubdivisionNode(const int &gridCellOrdinal, const int &subdivisionOrdinal,
                                       const int &subdivisionNodeNumber) const;
    
    //! returns coordinate in dimension d for the indicated subdivision of the indicated grid cell
    KOKKOS_INLINE_FUNCTION
    PointScalar subdivisionCoordinate(const int &gridCellOrdinal, const int &subdivisionOrdinal,
                                      const int &subdivisionNodeNumber, const int &d) const;
    
    //! Return the coordinate (weight) of the specified node.  For straight-edged geometry, this is simply the physical coordinate of the vertex.  For all geometries, this can be understood as a weight on the corresponding H^1 basis function used in the reference-to-physical map.
    KOKKOS_INLINE_FUNCTION
    PointScalar
    operator()(const int& cell, const int& node, const int& dim) const;
    
    //! Returns an integer indicating the number of distinct cell types vis-a-vis Jacobians
    KOKKOS_INLINE_FUNCTION
    int uniformJacobianModulus() const;
    
    /** \brief  Allocate a container into which Jacobians of the reference-to-physical mapping can be placed.
       \param [in] points - the points at which the Jacobian will be evaluated.
       \param [in] startCell - the first cell ordinal for which the Jacobian will be evaluated (used to define worksets smaller than the whole CellGeometry).
       \param [in] endCell - the first cell ordinal for which the Jacobian will not be evaluated (used to define worksets smaller than the whole CellGeometry); use -1 to indicate that all cells from startCell on should be evaluated.
       \return a Data object appropriately sized to accomodate the specified Jacobian values.
    */
    Data<PointScalar,DeviceType> allocateJacobianData(const TensorPoints<PointScalar,DeviceType> &points, const int startCell=0, const int endCell=-1) const;
    
    /** \brief  Allocate a container into which Jacobians of the reference-to-physical mapping can be placed.
       \param [in] points - the points at which the Jacobian will be evaluated.
       \param [in] startCell - the first cell ordinal for which the Jacobian will be evaluated (used to define worksets smaller than the whole CellGeometry).
       \param [in] endCell - the first cell ordinal for which the Jacobian will not be evaluated (used to define worksets smaller than the whole CellGeometry); use -1 to indicate that all cells from startCell on should be evaluated.
       \return a Data object appropriately sized to accomodate the specified Jacobian values.
    */
    Data<PointScalar,DeviceType> allocateJacobianData(const ScalarView<PointScalar,DeviceType> &points, const int startCell=0, const int endCell=-1) const;
    
    /** \brief  Allocate a container into which Jacobians of the reference-to-physical mapping can be placed (variant for affine geometry).
       \param [in] numPoints - the number of points at which the Jacobian will be defined.
       \param [in] startCell - the first cell ordinal for which the Jacobian will be evaluated (used to define worksets smaller than the whole CellGeometry).
       \param [in] endCell - the first cell ordinal for which the Jacobian will not be evaluated (used to define worksets smaller than the whole CellGeometry); use -1 to indicate that all cells from startCell on should be evaluated.
       \return a Data object appropriately sized to accomodate the specified Jacobian values.
    */
    Data<PointScalar,DeviceType> allocateJacobianData(const int &numPoints, const int startCell=0, const int endCell=-1) const;
    
    /** \brief Computes reference-space data for the specified points, to be used in setJacobian().
       \param [in] points - the points at which the Jacobian will be evaluated; shape (P,D) or (C,P,D).
       \return a Data object with any reference-space data required.  This may be empty, if no reference-space data is required in setJacobian().  If filled, will have shape (F,P,D) or (C,F,P,D).
    */
    Data<PointScalar,DeviceType> getJacobianRefData(const ScalarView<PointScalar,DeviceType> &points) const;
    
    /** \brief Computes reference-space data for the specified points, to be used in setJacobian().
       \param [in] points - the points at which the Jacobian will be evaluated, with shape (P,D).
       \return a Data object with any reference-space data required.  This may be empty, if no reference-space data is required in setJacobian().  If filled, will have shape (F,P,D).
    */
    Data<PointScalar,DeviceType> getJacobianRefData(const TensorPoints<PointScalar,DeviceType> &points) const;
    
    /** \brief Compute Jacobian values for the reference-to-physical transformation, and place them in the provided container.
       \param [out] jacobianData - a container, allocated by allocateJacobianData(), into which the evaluated Jacobians will be placed.
       \param [in] points - the points at which the Jacobian will be evaluated.
       \param [in] refData - the return from getJacobianRefData(); may be an empty container, depending on details of CellGeometry (e.g. if it is affine)
       \param [in] startCell - the first cell ordinal for which the Jacobian will be evaluated (used to define worksets smaller than the whole CellGeometry).
       \param [in] endCell - the first cell ordinal for which the Jacobian will not be evaluated (used to define worksets smaller than the whole CellGeometry); use -1 to indicate that all cells from startCell on should be evaluated.
    */
    void setJacobian(Data<PointScalar,DeviceType> &jacobianData, const TensorPoints<PointScalar,DeviceType> &points, const Data<PointScalar,DeviceType> &refData,
                     const int startCell=0, const int endCell=-1) const;
    
    /** \brief Compute Jacobian values for the reference-to-physical transformation, and place them in the provided container.
       \param [out] jacobianData - a container, allocated by allocateJacobianData(), into which the evaluated Jacobians will be placed.
       \param [in] points - the points at which the Jacobian will be evaluated.
       \param [in] refData - the return from getJacobianRefData(); may be an empty container, depending on details of CellGeometry (e.g. if it is affine)
       \param [in] startCell - the first cell ordinal for which the Jacobian will be evaluated (used to define worksets smaller than the whole CellGeometry).
       \param [in] endCell - the first cell ordinal for which the Jacobian will not be evaluated (used to define worksets smaller than the whole CellGeometry); use -1 to indicate that all cells from startCell on should be evaluated.
    */
    void setJacobian(Data<PointScalar,DeviceType> &jacobianData, const ScalarView<PointScalar,DeviceType> &points, const Data<PointScalar,DeviceType> &refData,
                     const int startCell=0, const int endCell=-1) const;
    
    /** \brief Compute Jacobian values for the reference-to-physical transformation, and place them in the provided container (variant for affine geometry).
       \param [out] jacobianData - a container, allocated by allocateJacobianData(), into which the evaluated Jacobians will be placed.
       \param [in] numPoints - the number of points at which the Jacobian will be defined.
       \param [in] startCell - the first cell ordinal for which the Jacobian will be evaluated (used to define worksets smaller than the whole CellGeometry).
       \param [in] endCell - the first cell ordinal for which the Jacobian will not be evaluated (used to define worksets smaller than the whole CellGeometry); use -1 to indicate that all cells from startCell on should be evaluated.
    */
    void setJacobian(Data<PointScalar,DeviceType> &jacobianData, const int &numPoints, const int startCell=0, const int endCell=-1) const;
  };
} // namespace Intrepid2

#include <Intrepid2_CellGeometryDef.hpp>

#endif /* Intrepid2_CellGeometry_h */

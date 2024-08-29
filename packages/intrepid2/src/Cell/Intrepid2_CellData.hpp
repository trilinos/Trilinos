// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CellData.hpp
    \brief  Header file for the classes:
    Intrepid2::RefSubcellParametrization,
    Intrepid2::RefCellNodes,
    Intrepid2::RefCellCenter.
    \author Kyungjoo Kim
    \author Mauro Perego
 */

#ifndef __INTREPID2_CELLDATA_HPP__
#define __INTREPID2_CELLDATA_HPP__

#include "Intrepid2_ConfigDefs.hpp"

#include "Shards_CellTopology.hpp"

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"
#include "Intrepid2_Kernels.hpp"

namespace Intrepid2 {

//============================================================================================//
//                                                                                            //
//                                RefSubcellParametrization                                   //
//                                                                                            //
//============================================================================================//

/** \class  Intrepid2::RefSubcellParametrization
      \brief  This class defines the parametrizations of edges and faces of
      supported reference cells.
      The parametrization mappings are stored in static Kokkos views.
      The class is templated on the Kokkos::Device Type which is used to determine layout
      and memory space of the views.

      Given an edge {V0, V1} of some reference cell, its parametrization is a mapping from
      [-1,1] onto the edge. Parametrization of a triangular face {V0,V1,V2} is mapping from
      the standard 2-simplex {(0,0,0), (1,0,0), (0,1,0)}, embedded in 3D onto that face.
      Parametrization of a quadrilateral face {V0,V1,V2,V3} is mapping from the standard
      2-cube {(-1,-1,0),(1,-1,0),(1,1,0),(-1,1,0)}, embedded in 3D, onto that face.

      All mappings are affine and orientation-preserving, i.e., they preserve the tangent
      and normal directions implied by the vertex order of the edge or the face relative to
      the reference cell:

      \li     the tangent on [-1,1] from -1 in the direction of 1 is mapped to a tangent on edge {V0,V1}
      from V0 in the direction of V1  (the forward direction of the edge determined by its
      start and end vertices)

      \li     the normal in the direction of (0,0,1) to the standard 2-simplex {(0,0,0),(1,0,0),(0,1,0)}
      and the standard 2-cube {(-1,-1,0),(1,-1,0),(1,1,0),(-1,1,0)} is mapped to a normal
      on {V0,V1,V2} and {V0,V1,V2,V3}, determined according to the right-hand rule
      (see http://mathworld.wolfram.com/Right-HandRule.html for definition of right-hand rule
      and Section \ref Section sec_cell_topology_subcell_map for further details).

      Because faces of all reference cells supported in Intrepid2 are affine images of either
      the standard 2-simplex or the standard 2-cube, the coordinate functions of the respective
      parmetrization maps are linear polynomials in the parameter variables (u,v), i.e., they
      are of the form \c F_i(u,v)=C_0(i)+C_1(i)u+C_2(i)v;  \c 0<=i<3 (face parametrizations
      are supported only for 3D cells, thus parametrization maps have 3 coordinate functions).
      As a result, application of these maps is independent of the face type which is convenient
      for cells such as Wedge or Pyramid that have both types of faces. Also, coefficients of
      coordinate functions for all faces can be stored together in the same array.
 */


template<typename DeviceType>
class RefSubcellParametrization {
public:
  using ConstViewType = Kokkos::DynRankView<const double,DeviceType>;

  /** \brief  Default constructor.
   */
  RefSubcellParametrization() = default;

  /** \brief  Destructor
   */
  ~RefSubcellParametrization() = default;


  /** \brief  Checks if a cell topology has a reference parametrization
      \param  cellTopoKey   [in]  - key of the cell topology
      \return true if the cell topology has a cell parametrization, false otherwise
   */
  static inline bool
  isSupported( const unsigned cellTopoKey );


  /** \brief  Returns a Kokkos view with the coefficients of the parametrization maps
      for the edges or faces of a reference cell topology.

      See Intrepid2::RefSubcellParametrization class description for details
      about parametrization maps.


      \param  subcellDim        [in]  - dimension of subcells whose parametrization map is returned
      \param  parentCellKey     [in]  - topology key of the reference cell owning the subcells

      \return a rank-3 Kokkos::View containing coefficients of the parameterization map for
              all subcells of the specified dimension.
              Returned view has dimensions (\# subcells, parentCellDim, subcellDim+1)

   */
  static inline
  ConstViewType
  get( const ordinal_type          subcellDim,
      const unsigned              parentCellKey );

private:

  /** \brief  Computes and stores static views containing the parametrizations maps of
      edges and faces of all reference cells.
   */
  static void set();

  /** \brief  Sets parametrizations of reference edges and faces of cell
      topologies with reference cells. Used to populate Host views from array data.

      See Intrepid2::RefSubcellParametrization::set() for more info about parametrization maps.

      \param  subcellParam           [out]  - host view with the coefficients of the parametrization map
      \param  subcellDim             [in]   - dimension of the subcells being parametrized (1 or 2)
      \param  parentCell             [in]   - topology of the parent cell owning the subcells.
   */
  template <typename HostViewType>
  static void
  set(       HostViewType          subcellParam,
      const ordinal_type          subcellDim,
      const shards::CellTopology  parentCell );

  //! static views containing the parametrization maps, allocated on DeviceType::memory_space
  using ViewType = Kokkos::DynRankView<double,DeviceType>;
  static ViewType lineEdgesParam;  // edge maps for 2d non-standard cells; shell line and beam
  static ViewType triEdgesParam, quadEdgesParam; // edge maps for 2d standard cells
  static ViewType shellTriEdgesParam, shellQuadEdgesParam; // edge maps for 3d non-standard cells; shell tri and quad
  static ViewType tetEdgesParam, hexEdgesParam, pyrEdgesParam, wedgeEdgesParam; // edge maps for 3d standard cells
  static ViewType shellTriFacesParam, shellQuadFacesParam; // face maps for 3d non-standard cells
  static ViewType tetFacesParam, hexFacesParam, pyrFacesParam, wedgeFacesParam; // face maps for 3d standard cells

  //! whether the parametrizations have been already computed using the method set()
  static bool isSubcellParametrizationSet_;


};



//============================================================================================//
//                                                                                            //
//                                        RefCellNodes                                        //
//                                                                                            //
//============================================================================================//


/** \class  Intrepid2::RefCellNodes
    \brief  This class defines the coordinates of the nodes of reference cells according for
    supported cell topologies.
    The node coordinates are stored in static views.
    The class is templated on the Kokkos::Device Type which is used to
    determine layout and memory space of the views.
 */
template<typename DeviceType>
class RefCellNodes {
public:
  using ConstViewType = Kokkos::DynRankView<const double,DeviceType>;

  /** \brief  Default constructor.
   */
  RefCellNodes() = default;

  /** \brief  Destructor
   */
  ~RefCellNodes() = default;


  /** \brief  Retrieves the Cartesian coordinates of reference cell nodes.

    Returns a Kokkos view containing the coordinates of reference cell nodes.
    Requires the key of the reference cell topology.
    Node coordinates are always returned as an (x,y,z)-triple
    regardless of the actual topological cell dimension. The unused coordinates are
    set to zero, e.g., node 0 of Line<2> is returned as {-1,0,0}.

    \param  cellTopoKey              [in]  - key of the cell topology
    \return a rank-2 Kokkos::View containing the coordinates of the cell nodes
            The returned view has dimensions (\# nodes, 3)

   */
  static inline
  ConstViewType
  get(const unsigned      cellTopoKey);

private:
  /** \brief Set reference nodes coordinates for supported topologies.
   */
  static void set();

  //! static views containing the node coordinates allocated on DeviceType::memory_space
  using ViewType = Kokkos::DynRankView<double,DeviceType>;
  static ViewType lineNodes, line3Nodes;
  static ViewType triangleNodes, triangle4Nodes, triangle6Nodes;
  static ViewType quadrilateralNodes, quadrilateral8Nodes, quadrilateral9Nodes;
  static ViewType tetrahedronNodes, tetrahedron8Nodes, tetrahedron10Nodes, tetrahedron11Nodes;
  static ViewType hexahedronNodes, hexahedron20Nodes, hexahedron27Nodes;
  static ViewType pyramidNodes, pyramid13Nodes, pyramid14Nodes;
  static ViewType wedgeNodes, wedge15Nodes, wedge18Nodes;


  /** \struct Intrepid2::RefCellNodes::ReferenceNodeDataStatic
    \brief Reference node containers for each supported topology
   */
  struct ReferenceNodeDataStatic {
    double line[2][3], line_3[3][3];
    double triangle[3][3], triangle_4[4][3], triangle_6[6][3];
    double quadrilateral[4][3], quadrilateral_8[8][3], quadrilateral_9[9][3];
    double tetrahedron[4][3], tetrahedron_8[8][3], tetrahedron_10[10][3], tetrahedron_11[10][3];
    double hexahedron[8][3], hexahedron_20[20][3], hexahedron_27[27][3];
    double pyramid[5][3], pyramid_13[13][3], pyramid_14[14][3];
    double wedge[6][3], wedge_15[15][3], wedge_18[18][3];
  };

  //! static struct containing the nodes coordinates on host
  static const ReferenceNodeDataStatic refNodeDataStatic_;

  //! whether the nodes coordinates have been already set using the method set()
  static bool isReferenceNodeDataSet_;

};

//============================================================================================//
//                                                                                            //
//                                       RefCellCenter                                        //
//                                                                                            //
//============================================================================================//


/** \class  Intrepid2::RefCellCenter
    \brief  This class defines the coordinates of the barycenter of the supported reference cells.
    The barycenter coordinates are stored in static views.
    The class is templated on the Kokkos::Device Type which is used to
    determine layout and memory space of the views.
 */

template<typename DeviceType>
class RefCellCenter {
public:
  using ConstViewType = Kokkos::DynRankView<const double,DeviceType>;


  /** \brief  Default constructor.
   */
  RefCellCenter() = default;

  /** \brief  Destructor
   */
  ~RefCellCenter() = default;


  /** \brief  Retrieves the Cartesian coordinates of a reference cell barycenter.

    Returns Cartesian coordinates of a reference cell barycenter. Requires cell topology
    with a reference cell. Barycenter coordinates are always returned as an (x,y,z)-triple
    regardless of the actual topological cell dimension. The unused coordinates are
    set to zero, e.g., center of Line<2> is returned as {0,0,0}.

    \param  cell              [in]  - key of the cell topology

    \return a rank-1 Kokkos::View containing the coordinates of the cell nodes
            The returned view has dimension 3.
   */
  static inline
  ConstViewType
  get(const unsigned      cellTopoKey);

private:
  /** \brief Set center coordinates of reference cell for supported topologies.
   */
  static void set();

  //! static views containing the center coordinates allocated on DeviceType::memory_space
  using ViewType = Kokkos::DynRankView<double,DeviceType>;
  static ViewType lineCenter;
  static ViewType triangleCenter;
  static ViewType quadrilateralCenter;
  static ViewType tetrahedronCenter;
  static ViewType hexahedronCenter;
  static ViewType pyramidCenter;
  static ViewType wedgeCenter;

  /** \struct Intrepid2::RefCellCenter::ReferenceNodeDataStatic
    \brief Reference node containers for each supported topology
   */
  struct ReferenceCenterDataStatic {
    double line[3];
    double triangle[3];
    double quadrilateral[3];
    double tetrahedron[3];
    double hexahedron[3];
    double pyramid[3];
    double wedge[3];
  };

  //! static struct containing the nodes coordinates on host
  static  const ReferenceCenterDataStatic refCenterDataStatic_;

  //! whether the center coordinates have been already set using the method set()
  static bool isReferenceCellCenterDataSet_;
};

//============================================================================================//
//                                                                                            //
//                                       PointInclusion                                       //
//                                                                                            //
//============================================================================================//


/** \class  Intrepid2::PointInclusion
    \brief  This class implements a <var>check</var> function that determines whether a given point is 
            inside or outside the reference cell for a specific topology.
            The class has a template argument for the key of the shards topology
 */

template<unsigned CellTopologyKey>
  struct PointInclusion;
  
  /** 
   \brief Line topology
  */
  template<>
  struct PointInclusion<shards::Line<>::key> {
    template<typename PointViewType, typename ScalarType>
    KOKKOS_INLINE_FUNCTION
    static bool
    check(const PointViewType &point, const ScalarType threshold);   
  };
  
  /** 
   \brief Triangle topology
  */
  template<>
  struct PointInclusion<shards::Triangle<>::key> {
    template<typename PointViewType, typename ScalarType>
    KOKKOS_INLINE_FUNCTION
    static bool
    check(const PointViewType &point, const ScalarType threshold);
  };
  
  /** 
   \brief Quadrilateral topology
  */
  template<>
  struct PointInclusion<shards::Quadrilateral<>::key> {

    template<typename PointViewType, typename ScalarType>
    KOKKOS_INLINE_FUNCTION
    static bool
    check(const PointViewType &point, const ScalarType threshold);
  };
    
  /** 
   \brief Tetrahedron topology
  */
  template<>
  struct PointInclusion<shards::Tetrahedron<>::key> {
    template<typename PointViewType, typename ScalarType>
    KOKKOS_INLINE_FUNCTION
    static bool
    check(const PointViewType &point, const ScalarType threshold);
  };

  /** 
   \brief Hexahedron topology
  */
  template<>
  struct PointInclusion<shards::Hexahedron<>::key> {
    template<typename PointViewType, typename ScalarType>
    KOKKOS_INLINE_FUNCTION
    static bool
    check(const PointViewType &point, const ScalarType threshold);
  };
  
  /** 
   \brief Pyramid topology
  */
  template<>
  struct PointInclusion<shards::Pyramid<>::key> {
    template<typename PointViewType, typename ScalarType>
    KOKKOS_INLINE_FUNCTION
    static bool
    check(const PointViewType &point, const ScalarType threshold);
  };

  /** 
   \brief Wedge topology
  */
  template<>
  struct PointInclusion<shards::Wedge<>::key> {
    template<typename PointViewType, typename ScalarType>
    KOKKOS_INLINE_FUNCTION
    static bool
    check(const PointViewType &point, const ScalarType threshold);
  };

  const CellTopologyData* getCellTopologyData(const unsigned& cellTopologyKey);
  
}

#include "Intrepid2_CellDataDef.hpp"

#endif


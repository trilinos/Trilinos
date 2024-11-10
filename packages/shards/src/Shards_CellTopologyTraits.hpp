// @HEADER
// *****************************************************************************
//                Shards : Shared Discretization Tools
//
// Copyright 2008-2011 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef Shards_CellTopologyTraits_hpp
#define Shards_CellTopologyTraits_hpp

#include <Shards_TypeList.hpp>
#include <Shards_IndexList.hpp>
#include <Shards_CellTopologyData.h>

namespace shards {

/** \addtogroup shards_package_cell_topology
 *  \{
 */

/**
 *  \brief  Return a CellTopology singleton for the given cell topology traits.
 */
template< class Traits >
const CellTopologyData * getCellTopologyData();

template< unsigned Dimension ,
          unsigned Number_Vertex ,
          unsigned Number_Node ,
          class    EdgeList = TypeListEnd ,
          class    EdgeMaps = TypeListEnd ,
          class    FaceList = TypeListEnd ,
          class    FaceMaps = TypeListEnd ,
          class    PermutationMaps = TypeListEnd ,
          class    PermutationPolarity = IndexList<> >
struct CellTopologyTraits ;

struct Node ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Implementation details for a much of this file ...

#ifndef DOXYGEN_COMPILE

template< class CellTop , class CellMap , unsigned Index , bool Good >
struct SubcellNodeIndex ;

template< class CellTop , class CellMap , unsigned Index >
struct SubcellNodeIndex< CellTop , CellMap , Index , false >
{ enum { value = ~0u }; };

template< class CellTop , class CellMap , unsigned Index >
struct SubcellNodeIndex< CellTop , CellMap , Index , true >
{
private:
  typedef typename CellTop::template subcell<0> subcell_node ;
public:
  enum { value = Index < subcell_node::count
               ? IndexListAt< CellMap , Index >::value : ~0u };
};

//----------------------------------------------------------------------

template< unsigned SubcellDim , unsigned SubcellOrd , unsigned NodeIndex ,
          unsigned Dimension , unsigned Number_Vertex , unsigned Number_Node ,
          class EdgeList , class EdgeMaps ,
          class FaceList , class FaceMaps ,
          class PermMaps , class Pol>
struct SubcellTopologyTraits ;

template< class ListType > struct TypeListHomogeneous ;

//----------------------------------------------------------------------
// Self-subcell reference

// Node
template<>
struct SubcellTopologyTraits<0,0,0,0,0,0,TypeListEnd,TypeListEnd,
                                         TypeListEnd,TypeListEnd,
                                         TypeListEnd,IndexList<> >
{
  typedef CellTopologyTraits<0,0,0> topology ;
  enum { count = 1 };
  enum { node = ~0u };
  enum { homogeneity = true };
};

// Particle
template<>
struct SubcellTopologyTraits<0,0,0,0,1,1,TypeListEnd,TypeListEnd,
                                         TypeListEnd,TypeListEnd,
                                         TypeListEnd,IndexList<> >
{
  typedef CellTopologyTraits<0,1,1> topology ;
  enum { count = 1 };
  enum { node = 0 };  // A Particle has 1 node, and NodeIndex (3rd tmpl arg) is 0, so it's valid
  enum { homogeneity = true };
};

// Line
template< unsigned NodeIndex ,
          unsigned NV , unsigned NN ,
          class EList , class EMaps ,
          class FList , class FMaps ,
          class PMaps , class Pol>
struct SubcellTopologyTraits<1,0,NodeIndex, 1,NV,NN,EList,EMaps,FList,FMaps,PMaps,Pol>
{
  typedef CellTopologyTraits<1,NV,NN,EList,EMaps,FList,FMaps,PMaps,Pol> topology ;
  enum { count = 1 };
  enum { node = NodeIndex < NN ? NodeIndex : ~0u };
  enum { homogeneity = true };
};

// Face
template< unsigned NodeIndex ,
          unsigned NV , unsigned NN ,
          class EList , class EMaps ,
          class FList , class FMaps ,
          class PMaps , class Pol>
struct SubcellTopologyTraits<2,0,NodeIndex, 2,NV,NN,EList,EMaps,FList,FMaps,PMaps,Pol>
{
  typedef CellTopologyTraits<2,NV,NN,EList,EMaps,FList,FMaps,PMaps,Pol> topology ;
  enum { count = 1 };
  enum { node = NodeIndex < NN ? NodeIndex : ~0u };
  enum { homogeneity = true };
};

// Volume
template< unsigned NodeIndex ,
          unsigned NV , unsigned NN ,
          class EList , class EMaps ,
          class FList , class FMaps ,
          class PMaps , class Pol>
struct SubcellTopologyTraits<3,0,NodeIndex, 3,NV,NN,EList,EMaps,FList,FMaps,PMaps,Pol>
{
  typedef CellTopologyTraits<3,NV,NN,EList,EMaps,FList,FMaps,PMaps,Pol> topology ;
  enum { count = 1 };
  enum { node = NodeIndex < NN ? NodeIndex : ~0u };
  enum { homogeneity = true };
};

//----------------------------------------------------------------------
// Node-subcell reference:

template< unsigned SubcellOrd ,
          unsigned D , unsigned NV , unsigned NN ,
          class EList , class EMaps ,
          class FList , class FMaps ,
          class PMaps , class Pol>
struct SubcellTopologyTraits<0,SubcellOrd,0, D,NV,NN,EList,EMaps,FList,FMaps,PMaps,Pol>
{
  typedef CellTopologyTraits<0,0,0> topology ;
  enum { count = NN };
  enum { node = SubcellOrd < NN ? SubcellOrd : ~0u };
  enum { homogeneity = true };
};

// Edge-subcell reference:

template< unsigned SubcellOrd , unsigned NodeIndex ,
          unsigned D , unsigned NV , unsigned NN ,
          class EList , class EMaps ,
          class FList , class FMaps ,
          class PMaps , class Pol>
struct SubcellTopologyTraits<1,SubcellOrd,NodeIndex,
                             D,NV,NN,EList,EMaps,FList,FMaps,PMaps,Pol>
{
private:
  typedef typename TypeListAt<EMaps,SubcellOrd>::type node_map ;
public:

  typedef typename TypeListAt<EList,SubcellOrd>::type topology ;

  enum { count = TypeListLength<EList>::value };

  enum { node = SubcellNodeIndex< topology , node_map , NodeIndex ,
                                  SubcellOrd < count >::value };

  enum { homogeneity = TypeListHomogeneous<EList>::value };
};

// Face-subcell reference:

template< unsigned SubcellOrd , unsigned NodeIndex ,
          unsigned D , unsigned NV , unsigned NN ,
          class EList , class EMaps ,
          class FList , class FMaps ,
          class PMaps , class Pol>
struct SubcellTopologyTraits< 2, SubcellOrd, NodeIndex,
                              D,NV,NN,EList,EMaps,FList,FMaps,PMaps,Pol >
{
private:
  typedef typename TypeListAt<FMaps,SubcellOrd>::type node_map ;
public:

  typedef typename TypeListAt<FList,SubcellOrd>::type topology ;

  enum { count = TypeListLength<FList>::value };

  enum { node = SubcellNodeIndex< topology , node_map , NodeIndex ,
                                  SubcellOrd < count >::value };

  enum { homogeneity = TypeListHomogeneous<FList>::value };
};

//----------------------------------------------------------------------
// Only partially specialized subcell references are valid.

template< unsigned SubcellDim , unsigned SubcellOrd , unsigned NodeIndex ,
          unsigned Dimension ,
          unsigned Number_Vertex , unsigned Number_Node ,
          class EdgeList , class EdgeMaps ,
          class FaceList , class FaceMaps ,
          class PermMaps , class Pol >
struct SubcellTopologyTraits
{
  typedef void topology ;
  enum { count = 0 };
  enum { node = ~0u };
  enum { homogeneity = false };
};

//----------------------------------------------------------------------

template<>
struct TypeListHomogeneous<TypeListEnd> {
  enum { value = true };
};

template< class T >
struct TypeListHomogeneous< TypeList<T,TypeListEnd> > {
  enum { value = true };
};

template< class T , class Tail >
struct TypeListHomogeneous< TypeList< T, TypeList< T , Tail > > > {
  enum { value = TypeListHomogeneous< TypeList<T,Tail> >::value };
};

template< class ListType >
struct TypeListHomogeneous 
{
  enum { value = false };
};

//----------------------------------------------------------------------

template< unsigned I , unsigned J > struct AssertEqual ;

template< unsigned I > struct AssertEqual<I,I> { enum { OK = true }; };

#endif /* DOXYGEN_COMPILE */

//----------------------------------------------------------------------
/**  \brief  Compile-time traits for a cell topology. */
template< unsigned Dimension , unsigned Number_Vertex , unsigned Number_Node ,
          class EdgeList , class EdgeMaps ,
          class FaceList , class FaceMaps ,
          class PermutationMaps ,
          class PermutationPolarity >
struct CellTopologyTraits
{
  /** \brief  The <em> self </em> type for the traits */
  typedef CellTopologyTraits< Dimension, Number_Vertex, Number_Node,
                              EdgeList, EdgeMaps, FaceList, FaceMaps,
                              PermutationMaps, PermutationPolarity > Traits ;

  enum {
    /** \brief  Topological dimension */
    dimension    = Dimension ,

    /** \brief  Number of vertices */
    vertex_count = Number_Vertex ,

    /** \brief  Number of nodes (a.k.a. Cell^0 subcells). */
    node_count   = Number_Node ,

    /** \brief  Number of edges (a.k.a. Cell^1 subcells). */
    edge_count   = TypeListLength<EdgeList>::value ,

#ifndef DOXYGEN_COMPILE
    face_count   = TypeListLength<FaceList>::value ,
#endif

    /** \brief  Number of sides (a.k.a. Cell^(D-1) subcells). */
    side_count   = Dimension == 3 ? face_count : (
                   Dimension == 2 ? edge_count : 0 ),

    /** \brief  Unique key for this topology.
     *
     *  Uniqueness assumes that extended topology nodes (non-vertex nodes)
     *  are placed regularly throughout the cell topology.  For example,
     *  if any edge has an interior node then all edges have an interior node.
     *  If this assumption is violated then the key cannot guarantee uniqueness.
     */
    key  = ( dimension    << 28 /*  4 bits, max    7 */ ) |
           ( face_count   << 22 /*  6 bits, max   63 */ ) |
           ( edge_count   << 16 /*  6 bits, max   63 */ ) |
           ( vertex_count << 10 /*  6 bits, max   63 */ ) |
           ( node_count         /* 10 bits, max 1023 */ ) };

  /** \brief Subcell information
   *
   *  - <b> subcell<Dim>::count        </b> Number of subcells of this dimension
   *  - <b> subcell<Dim>::homogeneity  </b> Homogeneous subcells of this dim
   *  - <b> subcell<Dim,Ord>::topology </b> topology of the subcell
   *  - <b> subcell<Dim,Ord,J>::node   </b> node ordinal of subcell's node J
   */
  template< unsigned Dim, unsigned Ord = 0, unsigned J = 0 >
  struct subcell :
    public SubcellTopologyTraits< Dim , Ord , J ,
                                  dimension , vertex_count , node_count ,
                                  EdgeList , EdgeMaps ,
                                  FaceList , FaceMaps ,
                                  PermutationMaps, PermutationPolarity > {};

  /** \brief Side subcell information
   *
   *  - <b> side<>::count       </b> Number of sides
   *  - <b> side<>::homogeneity </b> Homogeneous sides
   *  - <b> side<Ord>::topology </b> topology of the side
   *  - <b> side<Ord,J>::node   </b> node ordinal of side's node J
   */
  template< unsigned Ord = 0 , unsigned J = 0 >
  struct side :
    public SubcellTopologyTraits< ( 1 < dimension ? dimension - 1 : 4 ) ,
                                  Ord , J ,
                                  dimension , vertex_count , node_count ,
                                  EdgeList , EdgeMaps ,
                                  FaceList , FaceMaps ,
                                  TypeListEnd , IndexList<> > {};

  /** \brief Edge subcell information
   *
   *  - <b> edge<>::count       </b> Number of edge
   *  - <b> edge<>::homogeneity </b> Homogeneous edges
   *  - <b> edge<Ord>::topology </b> topology of the edge
   *  - <b> edge<Ord,J>::node   </b> node ordinal of edge's node J
   */
  template< unsigned Ord = 0 , unsigned J = 0 >
  struct edge :
    public SubcellTopologyTraits< ( 1 < dimension ? 1 : 4 ) , Ord , J ,
                                  dimension , vertex_count , node_count ,
                                  EdgeList , EdgeMaps ,
                                  TypeListEnd , TypeListEnd ,
                                  TypeListEnd , IndexList<> > {};

  //--------------------------------------------------------------------
  /** \brief  Node permutations for proper subcells.
   *
   *  ParentCell and SubCell are connected if every node of SubCell
   *  is also a node of ParentCell.  However, the connection may be
   *  permuted.
   *
   *  Let ParentCell be dimension D and SubCell be dimension dim < D.
   *  Let SubCell be connected as subcell Ord with permutation P.
   *
   *  Then <b> ParentCell.node(K) == SubCell.node(J) </b> where:
   *  -  SubCellTopology == ParentCellTopology::subcell<dim,Ord>::topology
   *  -  K  = ParentCellTopology::subcell<dim,Ord,JP>::node
   *  -  JP = SubCellTopology::permutation<P,J>::node
   *  -  J  = SubCellTopology::permutation_inverse<P,JP>::node
   *
   *  The permutation map for P == 0 is required to be identity.
   */
  template< unsigned Perm , unsigned J = 0 >
  struct permutation {
  private:
    typedef typename TypeListAt< PermutationMaps , Perm >::type node_map ;
  public:
    enum { node = J < node_count ? IndexListAt< node_map , J >::value : ~0u };
    enum { polarity = IndexListAt< PermutationPolarity , Perm >::value };
  };

  template< unsigned Perm , unsigned J = 0 >
  struct permutation_inverse {
  private:
    typedef typename TypeListAt< PermutationMaps , Perm >::type forward_map ;
    typedef typename IndexListInverse< forward_map >::type node_map ;
  public:
    enum { node = J < node_count ? IndexListAt< node_map , J >::value : ~0u };
    enum { polarity = IndexListAt< PermutationPolarity , Perm >::value };
  };

  enum { permutation_count = TypeListLength< PermutationMaps >::value };

  //--------------------------------------------------------------------

private:

#ifndef DOXYGEN_COMPILE

  enum { nedge_map = TypeListLength<EdgeMaps>::value ,
         nface_map = TypeListLength<FaceMaps>::value ,
         polarity_count = IndexListLength< PermutationPolarity >::value };

  enum { OK_edge  = AssertEqual< edge_count , nedge_map >::OK };
  enum { OK_face  = AssertEqual< face_count , nface_map >::OK };
  enum { OK_dimen = AssertEqual< 0 , (dimension    >>  3) >::OK };
  enum { OK_faceN = AssertEqual< 0 , (face_count   >>  6) >::OK };
  enum { OK_edgeN = AssertEqual< 0 , (edge_count   >>  6) >::OK };
  enum { OK_vertN = AssertEqual< 0 , (vertex_count >>  6) >::OK };
  enum { OK_nodeN = AssertEqual< 0 , (node_count   >> 10) >::OK };
  enum { OK_permN = AssertEqual< permutation_count, polarity_count >::OK };

#endif

};

/** \} */

} // namespace shards

#endif // Shards_CellTopologyTraits_hpp


/*
//@HEADER
// ************************************************************************
//
//                Shards : Shared Discretization Tools
//                 Copyright 2008 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Carter Edwards (hcedwar@sandia.gov),
//                    Pavel Bochev (pbboche@sandia.gov), or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
//@HEADER
*/

#ifndef Shards_CellTopologyData_h
#define Shards_CellTopologyData_h

#if defined( __cplusplus )
extern "C" {
#endif

/** \addtogroup shards_package_cell_topology
 *  \{
 */

/*----------------------------------------------------------------------*/
 
struct CellTopologyData ;
struct CellTopologyData_Subcell ;
struct CellTopologyData_Permutation ;

/** \brief A simple 'C' struct of cell topology attributes.
 *
 *  The topology may be extended such that the number of nodes
 *  (subcells of dimension zero) is greater than the number of
 *  vertices.  In this case the vertices must be ordered first.
 *
 * Nodes, edges, and sides are subcells with a particular dimension.
 * A cell has edges and sides only if its dimension is greater than one.
 *  - node has Dim == 0
 *  - edge has Dim == 1
 *  - side has Dim == dimension - 1.
 */
struct CellTopologyData {
  /** \brief  Base, a.k.a. not-extended, version of this topology
   *          where vertex_count == node_count.
   */
  const struct CellTopologyData * base ;

  /** \brief Intuitive name for this topology */
  const char * name ;

  /** \brief Unique key for this topology */
  unsigned key ;

  /** \brief Topological dimension */
  unsigned dimension ;

  /** \brief Number of vertices. */
  unsigned vertex_count ;

  /** \brief Number of nodes (a.k.a. \f$ {Cell}^{0} \f$ subcells).
   *
   *  A topology is <em> extended </em> if node_count > vertex_count
   */
  unsigned node_count ;

  /** \brief Number of edges (a.k.a. \f$ {Cell}^{1} \f$ boundary subcells). */
  unsigned edge_count ;

  /** \brief Number of sides (a.k.a. \f$ {Cell}^{D-1} \f$ boundary subcells). */
  unsigned side_count ;

  /** \brief Number of defined permutations */
  unsigned permutation_count ;

  /** \brief Flag if the subcells of a given dimension are homogeneous */
  unsigned subcell_homogeneity[4] ;

  /** \brief Number of subcells of each dimension. */
  unsigned subcell_count[4] ;

  /** \brief  Array of subcells of each dimension
   *
   *  The length of each subcell array is subcell_count[Dim]
   *  - <b> subcell[Dim][Ord].topology </b> topology of the subcell
   *  - <b> subcell[Dim][Ord].node[I]  </b> node ordinal of the subcell's node I
   */
  const struct CellTopologyData_Subcell * subcell[4] ;

  /** \brief  Array of side subcells of length side_count
   *
   *  The length of the side array is side_count
   *  - <b> side[Ord].topology </b> topology of the side
   *  - <b> side[Ord].node[I]  </b> node ordinal of the side's node I
   */
  const struct CellTopologyData_Subcell * side ;

  /** \brief  Array of edges subcells of length edge_count
   *
   *  The length of the edge array is edge_count
   *  - <b> edge[Ord].topology </b> topology of the edge
   *  - <b> edge[Ord].node[I]  </b> node ordinal of the edge's node I
   */
  const struct CellTopologyData_Subcell * edge ;

  /** \brief  Array of node permutations.
   *
   *  - required: 0 <= P < permutation_count
   *  - required: 0 <= I < node_count
   *
   *  Let ParentCell be dimension D and SubCell be dimension dim < D. 
   *  Let SubCell be connected as subcell Ord with permutation P. 
   * 
   *  Then <b> ParentCell.node(K) == SubCell.node(I) </b> where: 
   *  -  SubCellTopology == ParentCellTopology->subcell[dim][Ord].topology
   *  -  K  = ParentCellTopology->subcell[dim][Ord].node[IP]
   *  -  IP = SubCellTopology->permutation[P].node[I]
   *  -  I  = SubCellTopology->permutation_inverse[P].node[IP]
   *
   *  The permutation map for P == 0 is required to be identity. 
   */
  const struct CellTopologyData_Permutation * permutation ;
  const struct CellTopologyData_Permutation * permutation_inverse ;
};

/** \brief Subcell information.
 *
 *  - required: 0 <= Dim <= 3
 *  - required: 0 <= Ord <= subcell_count[Dim]
 *  - required: 0 <= J   <  subcell[Dim][Ord]->subcell_count[0]
 *  - subcell[Dim][Ord].topology
 *  - subcell[Dim][Ord].node[J]
 */
struct CellTopologyData_Subcell {
  /** \brief Subcell topology */
  const struct CellTopologyData * topology ;

  /** \brief Subcell indexing of \f$ {Cell}^{0} \f$
   *         with respect to parent cell. */
  const unsigned * node ;
};

/** \brief  Self-typedef */
typedef struct CellTopologyData  CellTopologyData ;

/** \brief  Array of node permutations.
 *
 *  - required: 0 <= P < permutation_count
 *  - required: 0 <= I < node_count
 *
 *  Let ParentCell be dimension D and SubCell be dimension dim < D. 
 *  Let SubCell be connected as subcell Ord with permutation P. 
 * 
 *  Then <b> ParentCell.node(K) == SubCell.node(I) </b> where: 
 *  -  SubCellTopology == ParentCellTopology->subcell[dim][Ord].topology
 *  -  K  = ParentCellTopology->subcell[dim][Ord].node[IP]
 *  -  IP = SubCellTopology->permutation[P].node[I]
 *  -  I  = SubCellTopology->permutation_inverse[P].node[IP]
 *
 *  The permutation map for P == 0 is required to be identity. 
 */
struct CellTopologyData_Permutation {
  const unsigned * node ;
  unsigned         polarity ;
};

/** \brief  Values for the CellTopologyData_Permutation polarity */
enum {
  CELL_PERMUTATION_POLARITY_IRRELEVANT = 0 ,
  CELL_PERMUTATION_POLARITY_POSITIVE   = 1 ,
  CELL_PERMUTATION_POLARITY_NEGATIVE   = 2
};

/** \brief  Map a cell->face->edge ordinal to the cell->edge ordinal.
 *          Return -1 for erroneous input.
 */
extern
int mapCellFaceEdge( const CellTopologyData * cell_topology ,
                     unsigned face_ordinal ,
                     unsigned face_edge_ordinal );

/** \} */

#if defined( __cplusplus )
} /* extern "C" */
#endif

#endif /* Shards_CellTopologyData_h */


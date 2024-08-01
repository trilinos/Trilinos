// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_IntrepidFieldPattern_hpp__
#define __Panzer_IntrepidFieldPattern_hpp__

#include "Panzer_FieldPattern.hpp"

// Trilinos includes
#include "Kokkos_Core.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Intrepid2_Basis.hpp"

#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Teuchos_RCP.hpp"
#include <set>

namespace panzer {

  /** This is a derived class that specializes based
   * on a single intrepid basis function.
   */
  class Intrepid2FieldPattern : public FieldPattern {
  public:
    Intrepid2FieldPattern(const Teuchos::RCP< Intrepid2::Basis<PHX::Device,double,double> > &intrepidBasis);

    virtual int getSubcellCount(int dim) const;
    virtual const std::vector<int> & getSubcellIndices(int dim, int cellIndex) const;
    virtual int getDimension() const;
    virtual shards::CellTopology getCellTopology() const;

    virtual void getSubcellClosureIndices(int dim,int cellIndex,std::vector<int> & indices) const;

    // static functions for examining shards objects

    /** For a given sub cell find the set of sub cells at all dimensions contained
     * internally.  This is inclusive, so that (dim,subCell) will be in the set.
     *
     * \param[in] cellTopo Parent cell topology being used.
     * \param[in] dim      Dimension of sub cell
     * \param[in] subCell  Ordinal of sub cell at specified dimension
     * \param[in,out] closure Set of sub cells associated with specified sub cell.
     *
     * \note Sub cell dimension and ordinals are inserted into <code>closure</code>.
     *       Previous information will not be removed.
     */
    static void buildSubcellClosure(const shards::CellTopology & cellTopo,unsigned dim,unsigned subCell,
                                    std::set<std::pair<unsigned,unsigned> > & closure);

    /** Search a cell topology for sub cells containing a specfic set of nodes.
     * This is a downward search (inclusive) from a user specified dimension.
     *
     * \param[in] cellTopo Parent cell topology being used.
     * \param[in] dim      Dimension of sub cell
     * \param[in] nodes    Nodes forming the super set
     * \param[in,out] subCells Specific sub cells containing the nodes.
     *
     * \note Sub cell dimension and ordinals are inserted into <code>subCells</code>.
     *       Previous information will not be removed.
     */
    static void findContainedSubcells(const shards::CellTopology & cellTopo,unsigned dim,
                                      const std::vector<unsigned> & nodes,
                                      std::set<std::pair<unsigned,unsigned> > & subCells);

    /** Get the set of nodes making up the user specified sub cells.
     *
     * \param[in] cellTopo Parent cell topology being used.
     * \param[in] dim      Dimension of sub cell
     * \param[in] subCell  Ordinal of sub cell at specified dimension
     * \param[in,out] nodes Nodes associated with sub cell.
     */
    static void getSubcellNodes(const shards::CellTopology & cellTopo,unsigned dim,unsigned subCell,
                                std::vector<unsigned> & nodes);

    /** \brief Does this field pattern support interpolatory coordinates?
     *
     * If this method returns true then <code>getInterpolatoryCoordinates</code> will
     * succeed, otherwise it will throw.
     *
     * \returns True if this pattern supports interpolatory coordinates.
     */
    bool supportsInterpolatoryCoordinates() const;

    /** Get the local coordinates for this field. This is independent of element
     * locations.
     *
     * \param[in,out] coords   Coordinates associated with this field type.
     */
    void getInterpolatoryCoordinates(Kokkos::DynRankView<double,PHX::Device> & coords) const;

    /** Get the local coordinates for this field.
     *
     * \param[in] cellNodes   Coordinates of the cell nodes.
     * \param[in,out] coords   Coordinates associated with this field type.
     * \param[in] meshCellTopology Mesh cell topology
     * 
     * \note If no cell topology is supplied, it will be inferred from the intrepid basis.
     * This will result in incorrect behavior for second order or higher meshes.
     */
    void getInterpolatoryCoordinates(const Kokkos::DynRankView<double,PHX::Device> & cellNodes,
                                     Kokkos::DynRankView<double,PHX::Device> & coords,
                                     Teuchos::RCP<const shards::CellTopology> meshCellTopology=Teuchos::null) const;

    /// Returns the underlying Intrepid2::Basis object
    Teuchos::RCP< Intrepid2::Basis<PHX::Device,double,double> > getIntrepidBasis() const;

  protected:
    Teuchos::RCP< Intrepid2::Basis<PHX::Device,double,double> > intrepidBasis_;

    //mutable std::vector<int> subcellIndices_;
    mutable std::vector<std::vector<std::vector<int> > > subcellIndicies_;
    std::vector<int> empty_;
  };

}

#endif

// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_NodalFieldPattern_hpp__
#define __Panzer_NodalFieldPattern_hpp__

#include "PanzerDofMgr_config.hpp"

#include <ostream>
#include <vector>

#include "Panzer_FieldPattern.hpp"

#include "Shards_CellTopology.hpp"

namespace panzer {

/** A simple field pattern that is useful for guranteeing
  * you can compute orientations.
  */
class NodalFieldPattern : public FieldPattern {
public:

   NodalFieldPattern();

   NodalFieldPattern(const shards::CellTopology & ct);

   //! Do nothing destructor
   virtual ~NodalFieldPattern() {}

   //! Set the cell topology for this field pattern
   void setCellTopology(const shards::CellTopology & ct);

   /** How many sub cells of a particular dimension are there?
     * For instance a "quad" element as four sub cells of dimension
     * 0 (node type),four sub cells of dimension one (edge type), and
     * one sub cell of dimension two (area type).
     *
     * \param[in] dim Dimension of the sub cell of interest
     *
     * \returns Number of sub cells of dimension <code>dim</code>
     */
   virtual int getSubcellCount(int dim) const;
 
   /** Get the local indices associated with a particular sub cell.
     * The sub cell is specified through its dimension and cell index.
     * A vector is returned that gives the indices.
     * 
     * \param[in] dim Dimension of the sub cell of interest
     * \param[in] cellIndex Index of the sub cell (must be lest than
     *                      <code>getSubcellCount(dim)</code> and greater
     *                      than or equal to 0.
     *
     * \returns A vector of indices for the specified sub cell.
     */
   virtual const std::vector<int> & getSubcellIndices(int dim,int cellIndex) const;

   /** Get the set of indices that are on the sub cell. The ``closure'' means
     * that this set includes indices of all the child sub cells as well.
     *
     * \param[in] dim Dimension of the sub cell of interest
     * \param[in] cellIndex Index of the sub cell (must be lest than
     *                      <code>getSubcellCount(dim)</code> and greater
     *                      than or equal to 0.
     * \param[out] indices Vector, to be overwritten, containing the closure
     *                     indices on return.
     */
   virtual void getSubcellClosureIndices(int dim,int cellIndex,std::vector<int> & indices) const;

   /** What is the dimension of this pattern. For
     * instance a "quad" element is two dimensional. This
     * is a purely geometric quantity.
     */
   virtual int getDimension() const;

   /** Get the cell topology associated with this field pattern.
     */
   virtual shards::CellTopology getCellTopology() const
   { return cellTopo_; }

public:
   shards::CellTopology cellTopo_;
   std::vector<std::vector<int> > nodeIndices_;
   std::vector<int> empty_;
};

}

#endif

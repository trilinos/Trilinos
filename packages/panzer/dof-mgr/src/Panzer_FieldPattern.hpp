// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_FieldPattern_hpp__
#define __Panzer_FieldPattern_hpp__

#include <ostream>
#include <vector>

#include "Shards_CellTopology.hpp"

namespace panzer {

class FieldPattern {
public:
   //! Do nothing destructor
   virtual ~FieldPattern() = 0;

   /** How many sub cells of a particular dimension are there?
     * For instance a "quad" element as four sub cells of dimension
     * 0 (node type),four sub cells of dimension one (edge type), and
     * one sub cell of dimension two (area type).
     *
     * \param[in] dim Dimension of the sub cell of interest
     *
     * \returns Number of sub cells of dimension <code>dim</code>
     */
   virtual int getSubcellCount(int dim) const = 0;
 
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
   virtual const std::vector<int> & getSubcellIndices(int dim,int cellIndex) const = 0;

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
   virtual void getSubcellClosureIndices(int dim,int cellIndex,std::vector<int> & indices) const = 0;

   /** What is the dimension of this pattern. For
     * instance a "quad" element is two dimensional. This
     * is a purely geometric quantity.
     */
   virtual int getDimension() const = 0;

   /** Get the cell topology associated with this field pattern.
     */
   virtual shards::CellTopology getCellTopology() const = 0;
   
   /** How many degrees of freedom are associated with
     * this pattern. In this terminology an ID is a single
     * degree of freedom defined at the element level. By default
     * this is the sum of the size of all the sub cells. Where
     * the sub cells are defined by the <code>getSubcellIndices</code>.
     *
     * \returns The number of IDs associated with this pattern.
     */
   virtual int numberIds() const;

   /** Do these patterns share the same geometry? Meaning 
     * are these two element shapes (on the canonical domain)
     * the same?
     *
     * \returns True if both geometries are the same, false otherwise.
     */
   virtual bool sameGeometry(const FieldPattern & fp) const;

   /** Are the subcells consistent for each dimension. A set of sub
     * cells is consistent if they have the same number of IDs.
     * The default behavior calls <code>getSubcellIndices</code> for
     * each dimension and cell index, and gurantees the number of IDs
     * are equal.
     *
     * \returns True if the pattern has consistent sub cells.
     */
   virtual bool consistentSubcells() const;

   /** Print this pattern to a stream.
     *
     * \param[in] os Desired output stream.
     */
   virtual void print(std::ostream & os) const;

   /* Compare two field patterns, a call to <code>sameGeometry()</code>
    * must return true and calls to <code>getSubcellIndices()</code> must
    * return the equivalent vectors
    */
   virtual bool equals(const FieldPattern & fp) const;
};

//! Stream IO
std::ostream & operator<<(std::ostream & os,const FieldPattern & fp);

}

#endif

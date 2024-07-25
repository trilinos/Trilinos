// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_GeometricAggFieldPattern_hpp__
#define __Panzer_GeometricAggFieldPattern_hpp__

#include "Panzer_FieldPattern.hpp"
#include "Panzer_FieldType.hpp"
#include "Teuchos_RCP.hpp"

namespace panzer {

/** This class takes a vector of field patterns and
  * finds the union of all the patterns. That is the
  * set of degrees of freedom needed to be specified
  * on the mesh to encapsulate all the patterns.
  *
  * For instance if one field pattern specifies one
  * degrees of freedom per node and another one specifies
  * one degree of freedom per edge, this pattern will
  * specify one degree of freedom per node and one per edge.
  *
  * Essentially if you took all of the FieldPatterns and 
  * stacked them on top of each other then the 
  * shadow cast would be the geometric field pattern.
  *
  * \note Because the <code>GeometricAggFieldPattern</code> is passed
  *       to the <code>ConnManager</code> this gives an implicit
  *       ordering to the mesh IDs that is associated
  *       with the logical progression of the sub cells.
  */
class GeometricAggFieldPattern : public FieldPattern {
public:
   GeometricAggFieldPattern();

   /** Construct the pattern, this automatically calls 
     * <code>buildPattern()</code> and is equivalent to calling the
     * default constructor and then <code>buildPattern()</code>.
     */
  GeometricAggFieldPattern(std::vector<std::pair<FieldType,Teuchos::RCP<const FieldPattern>>> & patterns);

   /** Construct the geometric version of this pattern, this automatically calls 
     * <code>buildPattern()</code> and is equivalent to calling the
     * default constructor and then <code>buildPattern()</code>.
     */
  GeometricAggFieldPattern(const FieldType& fieldType,
                           const Teuchos::RCP<const FieldPattern> & pattern);

   virtual ~GeometricAggFieldPattern() {}

   /** Construct the underlying data for this object.
     * Before calling this function getDimension(), 
     * getSubcellCount() and getSubcellIndices methods
     * will throw an exception.
     *
     * \note This can be called multiple times, however
     *       any previous state is not saved and the new
     *       state is constructed.
     */
   virtual void buildPattern(const std::vector<std::pair<FieldType,Teuchos::RCP<const FieldPattern>>> & patterns);
   
   /** Convenience function, simply builds a vector and calls the vectorized
     * form of <code>buildPattern</code>.
     */
   virtual void buildPattern(const FieldType& fieldType,
                             const Teuchos::RCP<const FieldPattern> & pattern);

   /** Returns the sub cell count (see <code>FieldPattern</code>) if 
     * <code>buildPattern</code> has been called. Otherwise it will throw an exception.
     */
   virtual int getSubcellCount(int dim) const;

   /** Returns the sub cell indices (see <code>FieldPattern</code>) if 
     * <code>buildPattern</code> has been called. Otherwise it will throw an exception.
     */
   virtual const std::vector<int> & getSubcellIndices(int dim,int cellIndex) const;

   /* This function has no functionality in this case.
    * If called it will throw an assertion failure
    */
   virtual void getSubcellClosureIndices(int /* dim */, int /* cellIndex */, std::vector<int>& /* indices */) const
   { TEUCHOS_ASSERT(false); }

   /** Returns the dimension (see <code>FieldPattern</code>) if 
     * <code>buildPattern</code> has been called. Otherwise it will throw an exception.
     */
   virtual int getDimension() const;

   /** Returns cell topology from parent field parents
     */
   virtual shards::CellTopology getCellTopology() const;

protected:
   bool patternBuilt_;
   std::size_t dimension_;
   std::vector<std::vector<std::vector<int> > > patternData_;
   shards::CellTopology cellTopo_;
};

}

#endif

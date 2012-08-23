// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_FieldAggPattern_hpp__
#define __Panzer_FieldAggPattern_hpp__

#include "Panzer_FieldPattern.hpp"

#include <map>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_Tuple.hpp"

namespace panzer {

/** This class performs the interlacing between all the
  * passed in field pattern. It is meant to be used internally
  * and is specialized to return the FEI arrays required.
  * Hence it requires a map containing both the desired field
  * pattern objects and the integer used by FEI to specify 
  * the field.
  */
class FieldAggPattern : public FieldPattern {
public:
   // The FEI format for prescribing fields seems quite general,
   // so it is my belief that an abstract can be written to it.

   FieldAggPattern();

   /** Automatically calls <code>buildPattern</code>
     */
   FieldAggPattern(std::vector<std::pair<int,Teuchos::RCP<const FieldPattern> > > & patterns,
                   const Teuchos::RCP<const FieldPattern> & geomAggPattern=Teuchos::null);

   virtual ~FieldAggPattern() {}

   /** Return the geometric field pattern associated with the set of field patterns.
     */
   virtual Teuchos::RCP<const FieldPattern> getGeometricAggFieldPattern() const;

   /** Build the pattern associated with a set of patterns with their field IDs.
     * Patterns must be on the same geometry, otherwise this function throws a
     * std::logic_error
     */
   virtual void buildPattern(const std::vector<std::pair<int,Teuchos::RCP<const FieldPattern> > > & patterns,
                             const Teuchos::RCP<const FieldPattern> & geomAggPattern=Teuchos::null);

   // functions from the abstract FieldPattern
   virtual int getDimension() const;
   virtual int getSubcellCount(int dimension) const;
   virtual const std::vector<int> & getSubcellIndices(int dimension, int subcell) const;
   virtual void getSubcellClosureIndices(int, int, std::vector<int> &) const;
   virtual shards::CellTopology getCellTopology() const;

   //! Print this pattern
   virtual void print(std::ostream & os) const;

   /** Extract the field pattern associated with an argument
     *
     * \param[in] fieldId ID for field to lookup
     *
     * \returns Field pattern associated with a particular ID
     */
   virtual Teuchos::RCP<const FieldPattern> getFieldPattern(int fieldId) const;

   //! \defgroup FEISupport Support functions for building FEI patterns
   //@{ 

   //! Lenght of vector is number of Ids, value is how many ids per field
   const std::vector<int> & numFieldsPerId() const { return numFields_; }

   /** Length of vector is sum(numFieldsPerId()).  Value is field associated with
     * the ID.
     */
   const std::vector<int> & fieldIds() const { return fieldIds_; }

   //@}

   //! \defgroup LOCALOFFSET Functions for indexing into a GID array ordered according to this pattern
   //@{ 

   /** This function produces a map between the ordering specified by
     * this pattern and an ordering required by the pattern associated with
     * <code>fieldId</code>. 
     *
     * For instance if you have a vector called <code>GIDs</code>
     * of length <code>this->numberIds()</code> and you want the GIDs
     * associated with <code>fieldId = 0 </code>. Simply take the <code>offsets</code>
     * vector and index it into the <code>GIDs</code> vector:
     * <code>GIDs[offsets[*]]</code>.
     *
     * \param[in] fieldId Field to look up
     *
     * \returns offsets Offsets into global IDs vector. The order of this vector 
     *                  is defined by the underlying FieldPattern defining the requested
     *                  field. For the IntrepidFieldPattern this will correspond to the 
     *                  required order of the basis functions
     */
   const std::vector<int> & localOffsets(int fieldId) const;
     
   /** Returns a pair of vectors. The first is the offsets into the global ID 
     * array for the field and subcell specified. This will be the
     * indexing in the <code>GIDs</code> vector (see <code>localOffsets</code>). Notice
     * that this returns the "closure".  Meaning that all offsets for a particular subcell
     * and the IDs that are on lower dimensional sub cells will be returned. This function
     * has no (currently) specified order. The second vector is the index of the field's
     * basis functions for the closure indices. The first and second vectors should
     * have the same size.
     *
     * \note You cannot depend on the order of the IDs will remain consistent in future
     *       versions.
     */
   const std::pair<std::vector<int>,std::vector<int> > &
   localOffsets_closure(int fieldId,int subcellDim,int subcellId) const;

   //@}

protected:
   typedef Teuchos::RCP<const FieldPattern> FPPtr; // for convenience

   /** Build the map that takes a field ID and maps it to a particular
     * index in the patterns vector
     */
   void buildFieldIdToPatternIdx();

   /** merge the current set of patterns in the patterns_ vector
     * into the numFields_ and fieldIds_ vectors. This looks 
     * at the patterns_ vector and has side effects, wiping out
     * any previous state stored in those vectors.
     */
   void buildFieldIdsVector();

   /** merge the current set of patterns in the patterns_ vector
     * into the numFields_ and fieldIds_ vectors. This looks 
     * at the patterns_ vector and has side effects.
     */
   void mergeFieldPatterns(int dim,int subcell);

   /** Adds all the pattern's sub cell indices to the vector
     * specified.  The sub cell (dim and sub cell ordinal) is also specified.
     */ 
   void addAllPatternSubcellIndices(int dim,int sc,std::vector<int> & indices);

   /** Build Pattern data vectors.
     */
   void buildFieldPatternData();

   /** Helper function that builds the local offsets for this field.
     */
   void localOffsets_build(int fieldId,std::vector<int> & offsets) const;
   
   std::size_t dimension_;
   Teuchos::RCP<const FieldPattern> geomAggPattern_;

   // field pattern data
   std::vector<std::vector<std::vector<int> > > patternData_;

   // FEI patterns
   std::vector<int> numFields_; // number of fields at each geometric ID
   std::vector<int> fieldIds_;  // field IDs at each geometric ID

   // storage for patterns
   std::vector<std::pair<int,Teuchos::RCP<const FieldPattern> > > patterns_;
   std::map<int,int> fieldIdToPatternIdx_;

   mutable std::map<int, std::vector<int> > fieldOffsets_;

   struct LessThan  
   { bool operator()(const Teuchos::Tuple<int,3> & a,const Teuchos::Tuple<int,3> & b) const; };
   mutable std::map<Teuchos::Tuple<int,3>, std::pair<std::vector<int>,std::vector<int> >,LessThan>
         fieldSubcellOffsets_closure_;
};

}

#endif

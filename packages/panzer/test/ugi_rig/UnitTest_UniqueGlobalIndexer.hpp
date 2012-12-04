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

#ifndef __UnitTest_UniqueGlobalIndexer_hpp__
#define __UnitTest_UniqueGlobalIndexer_hpp__

#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Panzer_UniqueGlobalIndexer.hpp"

#ifdef HAVE_MPI
   #include "Teuchos_DefaultMpiComm.hpp"
   #include "mpi.h"
#else
   #include "Teuchos_DefaultSerialComm.hpp"
#endif

#include "Intrepid_FieldContainer.hpp"

namespace panzer {
namespace unit_test {

/** This a dummy global indexer used for testing. It has two fields
  *    "TEMPERATURE" and "ION_TEMPERATURE"
  * There is one element block called "block_0" with two elements. Each element
  * resides on a different processor.
  */
class UniqueGlobalIndexer : public virtual panzer::UniqueGlobalIndexer<short,int> {
public:
   UniqueGlobalIndexer(int rank,int procCount);

   ~UniqueGlobalIndexer() {}

   /** \brief Get the number used for access to this
     *        field
     *
     * Get the number used for access to this
     * field. This is used as the input parameter
     * to the other functions that provide access
     * to the global unknowns.
     *
     * \param[in] str Human readable name of the field
     *
     * \returns A unique integer associated with the
     *          field if the field exisits. Otherwise
     *          a -1 is returned.
     */
   virtual int getFieldNum(const std::string & str) const;

   virtual int getNumFields() const { return 2; }
   virtual void getFieldOrder(std::vector<std::string> & order) const 
   { order.push_back("U"); order.push_back("T"); }

   virtual const std::string & getFieldString(int field) const;

   virtual Teuchos::RCP<Teuchos::Comm<int> > getComm() const
   { 
      #ifdef HAVE_MPI
         return Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
      #else
         return Teuchos::rcp(Teuchos::SerialComm<int>());
      #endif
   }

   /** What are the blockIds included in this connection manager?
     */
   virtual void getElementBlockIds(std::vector<std::string> & elementBlockIds) const; 

   /** Is the specified field in the element block? 
     */
   virtual bool fieldInBlock(const std::string & field, const std::string & block) const;

   /** Get the local element IDs for a paricular element
     * block.
     *
     * \param[in] blockId Block ID
     *
     * \returns Vector of local element IDs.
     */
   virtual const std::vector<short> & getElementBlock(const std::string & blockId) const;

   /** Get field numbers associated with a particular element block.
     */
   virtual const std::vector<int> & getBlockFieldNumbers(const std::string & blockId) const;

   /** \brief Get the global IDs for a particular element. This function
     * overwrites the <code>gids</code> variable.
     */
   virtual void getElementGIDs(short localElmtId,std::vector<int> & gids,const std::string & blockId="") const;

   virtual void getElementOrientation(short localElmtId,std::vector<double> & gidsOrientation) const
   { TEUCHOS_ASSERT(false); }

   /** \brief Use the field pattern so that you can find a particular
     *        field in the GIDs array.
     */
   virtual const std::vector<int> & getGIDFieldOffsets(const std::string & blockId,int fieldNum) const;

   /** \brief Use the field pattern so that you can find a particular
     *        field in the GIDs array. This version lets you specify the sub
     *        cell you are interested in and gets the closure. Meaning all the
     *        IDs of equal or lesser sub cell dimension that are contained within
     *        the specified sub cell. For instance for an edge, this function would
     *        return offsets for the edge and the nodes on that edge.
     *
     * \param[in] blockId
     * \param[in] fieldNum
     * \param[in] subcellDim
     * \param[in] subcellId
     */
   // virtual const std::vector<int> & 
   virtual const std::pair<std::vector<int>,std::vector<int> > & 
   getGIDFieldOffsets_closure(const std::string & blockId, int fieldNum,
                                                               int subcellDim,int subcellId) const;

   /** Get set of indices owned by this processor
     */
   virtual void getOwnedIndices(std::vector<int> & indices) const;

   /** Get set of indices owned and shared by this processor.
     * This can be thought of as the ``ghosted'' indices.
     */
   virtual void getOwnedAndSharedIndices(std::vector<int> & indices) const;

   /** Get a yes/no on ownership for each index in a vector
     */
   virtual void ownedIndices(const std::vector<int> & indices,std::vector<bool> & isOwned) const;

   void getCoordinates(short localElementId,Intrepid::FieldContainer<double> & points);

private:
   int procRank_;
   mutable Teuchos::RCP<std::vector<short> > elements_; // local element IDs
   mutable Teuchos::RCP<std::vector<int> > field0Offset_; // local element IDs
   mutable Teuchos::RCP<std::vector<int> > field1Offset_; // local element IDs

   mutable Teuchos::RCP<std::map<std::pair<int,int>,std::pair<std::vector<int>,std::vector<int> > > > field0_closures_;
   mutable Teuchos::RCP<std::map<std::pair<int,int>,std::pair<std::vector<int>,std::vector<int> > > > field1_closures_;
};

} // end unit test
} // end panzer

#endif

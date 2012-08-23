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

#ifndef __Panzer_DOFManager_decl_hpp__
#define __Panzer_DOFManager_decl_hpp__

#include <map>

// FEI includes
#include "fei_base.hpp"
#include "fei_Factory.hpp"

#ifdef HAVE_MPI
   #include "mpi.h"
#endif

#include "Panzer_config.hpp"
#include "Panzer_FieldPattern.hpp"
#include "Panzer_FieldAggPattern.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"

#include "Teuchos_RCP.hpp"

#include <boost/unordered_set.hpp>

namespace panzer {

template <typename LocalOrdinalT,typename GlobalOrdinalT>
class DOFManager : public UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> {
public:
   typedef GlobalOrdinalT GlobalOrdinal;
   typedef LocalOrdinalT LocalOrdinal;
   typedef std::map<int,std::string>::const_iterator const_field_iterator;

   virtual ~DOFManager() {}

   DOFManager();

   /** Constructor that sets the connection manager and communicator
     * objects. This is equivalent to calling the default constructor and
     * then "setConnManager" routine.
     */
   DOFManager(const Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > & connMngr,MPI_Comm mpiComm);

   /** \brief Set the connection manager and MPI_Comm objects.
     *
     * Set the connection manager and MPI_Comm objects. If this method
     * is called more than once, the behavior is to reset the indices in
     * the DOF manager.  However, the fields will be the same (this assumes
     * that the element blocks are consistent with the fields). The indices
     * will need to be rebuilt by calling <code>buildGlobalUnknowns</code>.
     *
     * \param[in] connMngr Connection manager to use.
     * \param[in] mpiComm  Communicator to use.
     */
   void setConnManager(const Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > & connMngr,MPI_Comm mpiComm);

   /** Get communicator associated with this manager.
     */
   virtual Teuchos::RCP<Teuchos::Comm<int> > getComm() const
   { return communicator_; }

   /** Get the FieldPattern describing the geometry used for this problem.
     * If it has not been constructed then null is returned.
     */
   Teuchos::RCP<const FieldPattern> getGeometricFieldPattern() const
   { return geomPattern_; }
   

   /** \brief Reset the indicies for this DOF manager.
     *
     * This method resets the indices and wipes out internal state. This method
     * does preserve the fields and the patterns added to the object. Also the
     * old connection manager is returned.
     *
     * \returns Old connection manager.
     */
   Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > resetIndices();

   /** \brief Add a field to the DOF manager.
     *
     * Add a field to the DOF manager. Immediately after
     * adding the field the field number and field size
     * will be available for a user to access
     *
     * \param[in] str Human readable name of the field
     * \param[in] pattern Pattern defining the basis function to be used
     *
     * \note <code>addField</code> cannot be called after <code>buildGlobalUnknowns</code> 
     *       or <code>registerFields</code>.
     */
   void addField(const std::string & str,const Teuchos::RCP<const FieldPattern> & pattern);

   void addField(const std::string & blockId,const std::string & str,const Teuchos::RCP<const FieldPattern> & pattern);

   /** Set the ordering of the fields to be used internally.  This controls
     * to some extent the local ordering (on a node or edge) of the individual fields.
     *
     * \param[in] fieldOrder Vector of field IDs order in the correct way
     *
     * \note If no ordering is set then the default ordering is alphabetical on 
     *       the field names (as dictated by <code>std::map<std::string,*></code>).
     */
   void setFieldOrder(const std::vector<int> & fieldOrder);

   /** Set the ordering of the fields to be used internally.  This controls
     * to some extent the local ordering (on a node or edge) of the individual fields.
     *
     * \param[in] fieldOrder Vector of field IDs order in the correct way
     *
     * \note If no ordering is set then the default ordering is alphabetical on 
     *       the field names (as dictated by <code>std::map<std::string,*></code>).
     */
   void setFieldOrder(const std::vector<std::string> & fieldOrder);

   /** Get the field order used. Return the field IDs.
     */
   void getFieldOrder(std::vector<int> & fieldOrder) const;

   /** Get the field order used. Return the field strings.
     */
   void getFieldOrder(std::vector<std::string> & fieldOrder) const;

   /** \brief Find a field pattern stored for a particular block and field number. This will
     *        retrive the pattern added with <code>addField(blockId,fieldNum)</code>.
     *
     * Find a field pattern stored for a particular block and field number. This will
     * retrive the pattern added with <code>addField(blockId,fieldNum)</code>. If no pattern
     * is found this function returns <code>Teuchos::null</code>.
     *
     * \param[in] blockId Element block id
     * \param[in] fieldNum Field integer identifier
     *
     * \returns Pointer to <code>FieldPattern</code> requested if the field exists,
     *          otherwise <code>Teuchos::null</code> is returned.
     */
   Teuchos::RCP<const FieldPattern> getFieldPattern(const std::string & blockId, int fieldNum) const;

   /** \brief Find a field pattern stored for a particular block and field number. This will
     *        retrive the pattern added with <code>addField(blockId,fieldNum)</code>.
     *
     * Find a field pattern stored for a particular block and field number. This will
     * retrive the pattern added with <code>addField(blockId,fieldNum)</code>. If no pattern
     * is found this function returns <code>Teuchos::null</code>.
     *
     * \param[in] blockId Element block id
     * \param[in] fieldName Field string identifier
     *
     * \returns Pointer to <code>FieldPattern</code> requested if the field exists,
     *          otherwise <code>Teuchos::null</code> is returned.
     */
   Teuchos::RCP<const FieldPattern> getFieldPattern(const std::string & blockId, const std::string & fieldName) const;
   
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
   int getFieldNum(const std::string & str) const;

   /** \brief Reverse lookup of the field string from
     *        a field number.
     *
     * \param[in] num Field number. Assumed to be 
     *                a valid field number.  Computed
     *                from <code>getFieldNum</code>.
     *
     * \returns Field name. 
     */
   const std::string & getFieldString(int num) const
   { return intToFieldStr_.find(num)->second; }
 
   /** \brief How many fields are handled by this manager.
     *
     * How many fields are handled by this manager. 
     *
     * \returns The number of fields used by this
     *          manager.
     */
   int getNumFields() const;
   
   /**  Returns the connection manager current being used.
     */
   Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > getConnManager() const 
   { return connMngr_; } 

   /** build the global unknown numberings
     *   1. this builds the pattens
     *   2. Build a default geometric pattern to pass to the connection manager
     *   3. initializes the connectivity
     *   4. calls initComplete
     */
   virtual void buildGlobalUnknowns();

   /** build the global unknown numberings
     *   1. this builds the pattens
     *   2. calls initComplete
     *
     * This method allows a different geometric
     * field pattern to used. It does not call the
     * ConnManger::buildConnectivity, and just
     * uses the provided field pattern as a the
     * geometric pattern. Note this requires that
     * ConnManager::buildConnectivity has already
     * been called.
     *
     * \note It might be useful (and fun!) to add a listener on to this function
     *       to notify interested parties of possible changes to the unknown structure
     *       and CRS matrix graph.
     */
   virtual void buildGlobalUnknowns(const Teuchos::RCP<const FieldPattern> & geomPattern);

   /** Builds the orientations for each unknown. Notice that unknowns
     * will be either 1 or -1. If the orientation is not required then
     * the rule is for the orientation is to be 1.
     */
   virtual void buildUnknownsOrientation();

   /** Prints to an output stream the information about
     * the aggregated field.
     */
   void printFieldInformation(std::ostream & os) const;

   //! \defgroup FieldAssembly_Indices Methods to access the global indices
   //{@

   /** What are the blockIds included in this connection manager?
     */
   virtual void getElementBlockIds(std::vector<std::string> & elementBlockIds) const
   { getConnManager()->getElementBlockIds(elementBlockIds); }

   /** Get field numbers associated with a particular element block.
     */
   virtual const std::vector<int> & getBlockFieldNumbers(const std::string & block) const
   { return fieldAggPattern_.find(block)->second->fieldIds(); }

   /** Is the specified field in the element block? 
     */
   virtual bool fieldInBlock(const std::string & field, const std::string & block) const
   { return fieldStringToPattern_.find(std::make_pair(block,field))->second != Teuchos::null; }

   /** Get the local element IDs for a paricular element
     * block.
     *
     * \param[in] blockId Block ID
     *
     * \returns Vector of local element IDs.
     */
   virtual const std::vector<LocalOrdinal> & getElementBlock(const std::string & blockId) const
   { return getConnManager()->getElementBlock(blockId); }

   /** \brief Get the global IDs for a particular element. This function
     * overwrites the <code>gids</code> variable.
     */
   void getElementGIDs(LocalOrdinalT localElmtId,std::vector<GlobalOrdinalT> & gids,const std::string & blockIdHint="") const;

   /** \brief Get a vector containg the orientation of the GIDs relative to the neighbors.
     */
   virtual void getElementOrientation(LocalOrdinalT localElmtId,std::vector<double> & gidsOrientation) const;

   /** \brief Use the field pattern so that you can find a particular
     *        field in the GIDs array.
     */
   const std::vector<int> & getGIDFieldOffsets(const std::string & blockId,int fieldNum) const
   { return fieldAggPattern_.find(blockId)->second->localOffsets(fieldNum); }

   /** \brief Use the field pattern so that you can find a particular
     *        field in the GIDs array. This version lets you specify the sub
     *        cell you are interested in and gets the closure. Meaning all the
     *        IDs of equal or lesser sub cell dimension that are contained within
     *        the specified sub cell. For instance for an edge, this function would
     *        return offsets for the edge and the nodes on that edge. The first
     *        vector returned contains the index into the GIDs array. The second vector
     *        specifies the basis function IDs.
     *
     * \param[in] blockId
     * \param[in] fieldNum
     * \param[in] subcellDim
     * \param[in] subcellId
     */
   const std::pair<std::vector<int>,std::vector<int> > & 
   getGIDFieldOffsets_closure(const std::string & blockId,int fieldNum,int subcellDim,int subcellId) const
   { return fieldAggPattern_.find(blockId)->second->localOffsets_closure(fieldNum,subcellDim,subcellId); }

   //@}

   /** Return an iterator that iterates over the 
     * <code>std::pair<int,std::string></code> that defines
     * a field.
     */
   const_field_iterator beginFieldIter() const
   { return intToFieldStr_.begin(); }

   /** Return an end iterator that signals the termination of 
     * <code>std::pair<int,std::string></code> that defines
     * a field. (Ends <code>beginFieldIter</code>)
     */
   const_field_iterator endFieldIter() const
   { return intToFieldStr_.end(); }

   /** Get set of indices owned by this processor
     */
   virtual void getOwnedIndices(std::vector<GlobalOrdinalT> & indices) const;

   /** Get set of indices owned and shared by this processor.
     * This can be thought of as the ``ghosted'' indices.
     */
   virtual void getOwnedAndSharedIndices(std::vector<GlobalOrdinalT> & indices) const;

   /** Get a yes/no on ownership for each index in a vector
     */
   virtual void ownedIndices(const std::vector<GlobalOrdinalT> & indices,std::vector<bool> & isOwned) const;

   /** Get a set of field numbers associated with a particular element block.
     */
   const std::set<int> & getFields(const std::string & blockId) const
   { return blockToField_.find(blockId)->second; }

   /** Check the validity of a field order. This is used internally
     * as a sanity check. Checks for no repeats, bogus fields, and all fields
     * being included.
     *
     * \param[in] fieldOrder_ut Field order vector under test (ut).
     *
     * \returns true if the vector is valid, false otherwise.
     */
   bool validFieldOrder(const std::vector<std::string> & fieldOrder_ut,const std::set<std::string> & fields) const;

   /** This builds all numbers for the fields as well as
     * constructing a default field orderand validating the user specified field order.
     */
   void registerFields();

   /** How any GIDs are associate with a particular element block
     */
   inline int getElementBlockGIDCount(const std::string & blockId) const
   { return getElementBlockGIDCount(blockIdToIndex(blockId)); }

   /** How any GIDs are associate with a particular element block
     */
   inline int getElementBlockGIDCount(std::size_t blockIndex) const
   { int cnt = matrixGraph_->getConnectivityNumIndices(blockIndex); 
     if(cnt<0)
        return 0;
     return cnt;
   }

   /** Return if the orientations have been set to required.
     */
   bool getOrientationsRequired() const
   { return requireOrientations_; }

   /** Enable computation of the orientations.
     */
   void setOrientationsRequired(bool ro) 
   { requireOrientations_ = ro; }

protected:
   
   /** Get ordered field IDs associated with a particular element
     * block.
     */
   void getOrderedBlock(const std::vector<std::string> & fieldOrder,
                        const std::string & blockId,
                        std::vector<int> & orderedBlock) const;

   /** Using the natural ordering associated with the std::vector
     * retrieved from the connection manager
     */
   std::size_t blockIdToIndex(const std::string & blockId) const;

   /** Access the block Id to index vector directly.
     */
   const std::map<std::string,std::size_t> & blockIdToIndexMap() const;

   //! build the pattern associated with this manager
   bool buildPattern(const std::vector<std::string> & fieldOrder,
                     const std::string & blockId);

   // computes connectivity
   Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > connMngr_; 
   
   //! \defgroup MapFunctions Mapping objects
   //@{ 
   //! field string ==> field id
   std::map<std::string,int> fieldStrToInt_;
   std::map<int,std::string> intToFieldStr_;

   //! (block ID x field string) ==> pattern
   std::map<std::pair<std::string,std::string>,Teuchos::RCP<const FieldPattern> > fieldStringToPattern_;

   //! (block ID x field id) ==> pattern
   std::map<std::pair<std::string,int>,Teuchos::RCP<const FieldPattern> > fieldIntToPattern_;

   //! block id ==> Aggregate field pattern
   std::map<std::string,Teuchos::RCP<FieldAggPattern> > fieldAggPattern_;

   //! block id ==> set of field ids
   std::map<std::string,std::set<int> > blockToField_; // help define the pattern
   //@}

   // FEI based DOF management stuff
   Teuchos::RCP<fei::Factory> feiFactory_;
   fei::SharedPtr<fei::VectorSpace> vectorSpace_;
   fei::SharedPtr<fei::MatrixGraph> matrixGraph_;

   // map from a field to a vector of local element IDs
   std::map<int,std::vector<int> > field2ElmtIDs_;

   // storage for fast lookups of GID ownership
   boost::unordered_set<GlobalOrdinal> ownedGIDHashTable_;

   // maps blockIds to indices
   mutable Teuchos::RCP<std::map<std::string,std::size_t> > blockIdToIndex_;

   std::vector<std::string> fieldOrder_;

   // counters
   int nodeType_;
   int edgeType_;
   int numFields_;
   std::vector<int> patternNum_;

   bool fieldsRegistered_;

   Teuchos::RCP<const FieldPattern> geomPattern_;
   Teuchos::RCP<Teuchos::Comm<int> > communicator_;

   bool requireOrientations_;
   // this vector will be # of local elements, by number of GIDs on element block
   std::vector<std::vector<char> > orientation_; // we are using chars here
                                                 // to minimize storage and also
                                                 // we need only to store +/-1
};

}

#endif

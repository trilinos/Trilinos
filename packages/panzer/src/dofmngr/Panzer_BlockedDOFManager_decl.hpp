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

#ifndef __Panzer_BlockedDOFManager_decl_hpp__
#define __Panzer_BlockedDOFManager_decl_hpp__

#include <map>

#ifdef HAVE_MPI
   #include "mpi.h"
#endif

#include "Panzer_config.hpp"
#include "Panzer_FieldPattern.hpp"
#include "Panzer_FieldAggPattern.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_DOFManagerFEI.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultMpiComm.hpp"

#include <boost/unordered_set.hpp>

namespace panzer {

template <typename LocalOrdinalT,typename GlobalOrdinalT>
class BlockedDOFManager : public UniqueGlobalIndexer<LocalOrdinalT,std::pair<int,GlobalOrdinalT> > {
public:
   typedef std::pair<int,GlobalOrdinalT> GlobalOrdinal;
   typedef LocalOrdinalT LocalOrdinal;
   typedef std::map<int,std::string>::const_iterator const_field_iterator;

   virtual ~BlockedDOFManager() {}

   BlockedDOFManager();

   /** Constructor that sets the connection manager and communicator
     * objects. This is equivalent to calling the default constructor and
     * then "setConnManager" routine.
     */
   BlockedDOFManager(const Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > & connMngr,MPI_Comm mpiComm);

   ////////////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////////////
   //! \defgroup UniqueGlobalIndexer_methods Methods required by the <code>UniqueGlobalIndexer</code> interface
   //@{

   /** Get communicator associated with this manager.
     */
   virtual Teuchos::RCP<Teuchos::Comm<int> > getComm() const
   { return communicator_; }

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
   int getFieldNum(const std::string & str) const; // ?

   /** \brief Get the string name associated with a field number.
     *
     * Get the string used for access to this
     * field. 
     *
     * \param[in] int A unique integer associated with the
     *                field.
     * 
     * \returns Human readable name of the field
     *
     * \note This method will throw if invalid field number is 
     *       passed in as an argument.
     */
   const std::string & getFieldString(int num) const;

   /** What are the blockIds included in this connection manager?
     */
   virtual void getElementBlockIds(std::vector<std::string> & elementBlockIds) const
   { getConnManager()->getElementBlockIds(elementBlockIds); }

   /** Is the specified field in the element block? 
     */
   virtual bool fieldInBlock(const std::string & field, const std::string & block) const; // ?

   /** Get the local element IDs for a paricular element
     * block.
     *
     * \param[in] blockId Block ID
     *
     * \returns Vector of local element IDs.
     */
   virtual const std::vector<LocalOrdinal> & getElementBlock(const std::string & blockId) const // ?
   { return getConnManager()->getElementBlock(blockId); }

   /** Get field numbers associated with a particular element block.
     */
   virtual const std::vector<int> & getBlockFieldNumbers(const std::string & block) const; // ?

   /** \brief Get the global IDs for a particular element. This function
     * overwrites the <code>gids</code> variable.
     */
   void getElementGIDs(LocalOrdinalT localElmtId,std::vector<GlobalOrdinal> & gids,const std::string & blockIdHint="") const; // ?

   /** \brief Get a vector containg the orientation of the GIDs relative to the neighbors.
     */
   virtual void getElementOrientation(LocalOrdinalT localElmtId,std::vector<double> & gidsOrientation) const; // ?

   /** \brief Use the field pattern so that you can find a particular
     *        field in the GIDs array.
     */
   virtual const std::vector<int> & getGIDFieldOffsets(const std::string & blockId,int fieldNum) const; // ?

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
   virtual const std::pair<std::vector<int>,std::vector<int> > & 
   getGIDFieldOffsets_closure(const std::string & blockId,int fieldNum,int subcellDim,int subcellId) const; // ?

   /** Get set of indices owned by this processor
     */
   virtual void getOwnedIndices(std::vector<GlobalOrdinal> & indices) const; // ?

   /** Get set of indices owned and shared by this processor.
     * This can be thought of as the ``ghosted'' indices.
     */
   virtual void getOwnedAndSharedIndices(std::vector<GlobalOrdinal> & indices) const; // ?

   /** Get a yes/no on ownership for each index in a vector
     */
   virtual void ownedIndices(const std::vector<GlobalOrdinal> & indices,std::vector<bool> & isOwned) const; // ?

   //@}
   ////////////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////////////

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


   /** Get the FieldPattern describing the geometry used for this problem.
     * If it has not been constructed then null is returned.
     */
   Teuchos::RCP<const FieldPattern> getGeometricFieldPattern() const // ?
   { return geomPattern_; }
   
   /** \brief Reset the indicies for this DOF manager.
     *
     * This method resets the indices and wipes out internal state. This method
     * does preserve the fields and the patterns added to the object. Also the
     * old connection manager is returned.
     *
     * \returns Old connection manager.
     */
   Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > resetIndices(); // ?

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
   void setFieldOrder(const std::vector<std::vector<std::string> > & fieldOrder);

   /** Return the number of field blocks associated with this manager. If the field
     * order has not been set then this method returns 1.
     */
   int getNumFieldBlocks() const;

   /** Get the field order used. Return the field strings.
     */
   void getFieldOrder(std::vector<std::vector<std::string> > & fieldOrder) const;

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
   Teuchos::RCP<const FieldPattern> getFieldPattern(const std::string & blockId, const std::string & fieldName) const; // ?
 
   /** \brief How many fields are handled by this manager.
     *
     * How many fields are handled by this manager. 
     *
     * \returns The number of fields used by this
     *          manager.
     */
   int getNumFields() const; // ?

   /**  Returns the connection manager current being used.
     */
   Teuchos::RCP<const ConnManager<LocalOrdinalT,GlobalOrdinalT> > getConnManager() const 
   { return connMngr_; } 

   /**  Returns the connection manager current being used.
     */
   Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > getConnManager() 
   { return connMngr_; } 

   /** build the global unknown numberings
     *   1. this builds the pattens
     *   2. Build a default geometric pattern to pass to the connection manager
     *   3. initializes the connectivity
     *   4. calls initComplete
     */
   virtual void buildGlobalUnknowns(); // ?

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
   virtual void buildGlobalUnknowns(const Teuchos::RCP<const FieldPattern> & geomPattern); // ?

   /** Prints to an output stream the information about
     * the aggregated field.
     */
   void printFieldInformation(std::ostream & os) const; // ?

   /** Check the validity of a field order. This is used internally
     * as a sanity check. Checks for no repeats, bogus fields, and all fields
     * being included.
     *
     * \param[in] fieldOrder_ut Field order vector under test (ut).
     *
     * \returns true if the vector is valid, false otherwise.
     */
   bool validFieldOrder(const std::vector<std::vector<std::string> > & fieldOrder_ut,const std::set<std::string> & fields) const;

   /** This builds all numbers for the fields as well as
     * constructing a default field orderand validating the user specified field order.
     */
   void registerFields(); // ?

   /** Has the method <code>registerFields</code> been called?
     */
   bool fieldsRegistered() const 
   { return fieldsRegistered_; }

   /** Extract the field DOFManagers used underneath to define the
     * global unknowns.
     */ 
   const std::vector<Teuchos::RCP<DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT> > > &
   getFieldDOFManagers() const
   { return fieldBlockManagers_; }

   /** Return the maximum field number returned by any sub DOFManager.
     * Mostly exposed for testing purposes.
     */ 
   inline int getMaxSubFieldNumber() const
   { return maxSubFieldNum_; }
 
   /** Get field block associated with this field number.
     *
     * \note No bounds checking is performed in this method
     */
   int getFieldBlock(int fieldNum) const
   { return fieldNumToFieldBlk_.find(fieldNum)->second; }

   /** Get GID offset for a particular field block.
     *
     * \note No bounds checking is performed in this method
     */
   int getBlockGIDOffset(const std::string & elementBlock,int fieldBlock) const
   { 
      std::map<std::pair<std::string,int>,int>::const_iterator itr = 
            blockGIDOffset_.find(std::make_pair(elementBlock,fieldBlock));

      if(itr==blockGIDOffset_.end())
         return -1;
      else
         return itr->second;
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
  
   /** This routine calls the <code>addField</code> method on the fieldBlockManager adding all
     * the fields it is supposed to control, and then calls registerFields.
     *
     * This method assumes that the activeFields are a legitimate ordering for the local field block.
     */
   void addFieldsToFieldBlockManager(const std::vector<std::string> & activeFields,
                                     DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT> & fieldBlockManager) const;


   // computes connectivity
   Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > connMngr_; 
   
   //! \defgroup MapFunctions Mapping objects
   //@{ 
   //! field string ==> field number
   std::map<std::string,int> fieldStrToNum_;

   //! field number ==> field string
   std::map<int,std::string> fieldNumToStr_;

   //! field number ==> field block
   std::map<int,int> fieldNumToFieldBlk_;

   //! (block ID x field string) ==> pattern
   std::map<std::pair<std::string,std::string>,Teuchos::RCP<const FieldPattern> > fieldStringToPattern_;

   //! block ID ==> field strings
   std::map<std::string,std::set<std::string> > blockIdToFieldStrings_;

   //! block ID ==> field numbers
   std::map<std::string,std::vector<int> > blockIdToFieldNumbers_;

   //! (element block,field block) ==> gid offset
   std::map<std::pair<std::string,int>,int> blockGIDOffset_;

   //@}

   // storage for fast lookups of GID ownership
   boost::unordered_set<GlobalOrdinal> ownedGIDHashTable_;

   std::vector<std::vector<std::string> > fieldOrder_;

   bool fieldsRegistered_;

   Teuchos::RCP<const FieldPattern> geomPattern_;
   Teuchos::RCP<Teuchos::MpiComm<int> > communicator_;

   std::vector<Teuchos::RCP<DOFManagerFEI<LocalOrdinalT,GlobalOrdinalT> > > fieldBlockManagers_;

   MPI_Comm mpiComm_;
   int maxSubFieldNum_;

   /** Maps: elem block ids ==> (fieldNum ==> gidFieldOffsets vector) 
     * This uses lazy evaluation for construction.
     */
   mutable std::map<std::string,std::map<int,std::vector<int> > > gidFieldOffsets_;

   struct LessThan
   { bool operator()(const Teuchos::Tuple<int,3> & a,const Teuchos::Tuple<int,3> & b) const; };
   typedef std::map<Teuchos::Tuple<int,3>, std::pair<std::vector<int>,std::vector<int> >,LessThan> TupleToVectorPairMap;

   /** Maps: elem block ids ==> ( (fieldNum,subcellDim,subcellId) ==> closure vector pair)
     * This uses lazy evaluation for construction.
     */
   mutable std::map<std::string,TupleToVectorPairMap> gidFieldOffsets_closure_;

   bool requireOrientations_;
};

}

#endif

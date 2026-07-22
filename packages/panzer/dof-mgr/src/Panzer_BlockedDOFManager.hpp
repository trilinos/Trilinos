// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_BlockedDOFManager_hpp__
#define __Panzer_BlockedDOFManager_hpp__

#include <map>
#include <set>

#ifdef HAVE_MPI
   #include <mpi.h>
#endif

#include "PanzerDofMgr_config.hpp"
#include "Panzer_FieldPattern.hpp"
#include "Panzer_FieldAggPattern.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_NodeType.hpp"
#include "Panzer_HashUtils.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultMpiComm.hpp"

#include <unordered_set>

namespace panzer {

class BlockedDOFManager : public GlobalIndexer {
public:
   // typedef std::pair<int,GlobalOrdinalT> GlobalOrdinal;
  using GlobalOrdinal = panzer::GlobalOrdinal;
  using LocalOrdinal = panzer::LocalOrdinal;
  using const_field_iterator = std::map<int,std::string>::const_iterator;

   virtual ~BlockedDOFManager() {}

   BlockedDOFManager();

   /** Constructor that sets the connection manager and communicator
     * objects. This is equivalent to calling the default constructor and
     * then "setConnManager" routine.
     */
   BlockedDOFManager(const Teuchos::RCP<ConnManager> & connMngr,MPI_Comm mpiComm);

   ////////////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////////////
   //! \defgroup GlobalIndexer_methods Methods required by the <code>GlobalIndexer</code> interface
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
   void getElementGIDs(panzer::LocalOrdinal localElmtId,std::vector<GlobalOrdinal> & gids,const std::string & blockIdHint="") const; // ?

   /** \brief Get the global IDs for a particular element. This function
     * overwrites the <code>gids</code> variable.
     *
     * The pair consists of the field block index (pair.first) and the
     * corresponding gid (pair.second) within that field block. The
     * field block is used to access the sub-DOFManager that the field
     * is assigned to.
     *
     * NOTE: This function is temporary and is only used by the Epetra
     * Blocked Gather/Scatters. It is an inefficient path to access
     * data and has been deprecated. When Epetra support is dropped,
     * this function will be removed.
     */
  void getElementGIDsPair(panzer::LocalOrdinal localElmtId,std::vector<std::pair<int,GlobalOrdinal>> & gids,const std::string & blockIdHint="") const;

   /** \brief Get a vector containg the orientation of the GIDs relative to the neighbors.
     */
   virtual void getElementOrientation(panzer::LocalOrdinal localElmtId,std::vector<double> & gidsOrientation) const; // ?

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

   /**
    *  \brief Get the set of indices owned by this processor.
    *
    *  \param[out] indices A `vector` that will be filled with the indices
    *              owned by this processor.
    */
   virtual void getOwnedIndices(std::vector<GlobalOrdinal>& indices) const;

   /**
    *  \brief Get the set of indices ghosted for this processor.
    *
    *  \param[out] indices A `vector` that will be filled with the indices
    *              ghosted for this processor.
    */
   virtual void
   getGhostedIndices(std::vector<GlobalOrdinal>& indices) const;

   /**
    *  \brief Get the set of owned and ghosted indices for this processor.
    *
    *  \param[out] indices A `vector` that will be filled with the owned and
    *              ghosted indices for this processor.
    */
   virtual void
   getOwnedAndGhostedIndices(std::vector<GlobalOrdinal>& indices) const;
  
   // For backwards compatibility with Epetra. Will be deprecated.
   void getElementGIDsAsInt(panzer::LocalOrdinal localElmtId,std::vector<int> & gids,const std::string & blockIdHint="") const;
   virtual void getOwnedIndicesAsInt(std::vector<int>& indices) const;
   virtual void getGhostedIndicesAsInt(std::vector<int>& indices) const;
   virtual void getOwnedAndGhostedIndicesAsInt(std::vector<int>& indices) const;

   /**
    *  \brief Get the number of indices owned by this processor.
    *
    *  \returns The number of indices owned by this processor.
    */
   virtual int
   getNumOwned() const;

   /**
    *  \brief Get the number of indices ghosted for this processor.
    *
    *  \returns The number of indices ghosted for this processor.
    */
   virtual int
   getNumGhosted() const;

   /**
    *  \brief Get the number of owned and ghosted indices for this processor.
    *
    *  \returns The number of owned and ghosted indices for this processor.
    */
   virtual int
   getNumOwnedAndGhosted() const;

   /** Get a yes/no on ownership for each index in a vector
     */
   virtual void ownedIndices(const std::vector<GlobalOrdinal> & indices,std::vector<bool> & isOwned) const; // ?

   //@}
   ////////////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////////////

   /** \brief Use the FEI DOF manager internally, or the standard version.
     */
   void setUseDOFManagerFEI(bool useFEI)
   { useDOFManagerFEI_ = useFEI; }

   /** \brief which DOF Manager is used internally?
     */
   bool getUseDOFManagerFEI() const
   {
     return false;
   }

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
   void setConnManager(const Teuchos::RCP<ConnManager> & connMngr,MPI_Comm mpiComm);


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
   Teuchos::RCP<ConnManager> resetIndices(); // ?

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

   void getFieldOrder(std::vector<std::string>& /* fieldOrder */) const { TEUCHOS_ASSERT(false); } // what???

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
   Teuchos::RCP<const ConnManager> getConnManager() const
   { return connMngr_; }

   /**  Returns the connection manager current being used.
     */
   Teuchos::RCP<ConnManager> getConnManager()
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

   /** This method simply builds the global unknowns by using the passed in global indexers.
     * The internal connection manager must use the underlying connection manager for all
     * the global indexers. Finally only global indexers of type
     * <code>DOFManager</code> can be used at the moment.
     *
     * \note The type of global indexer, and agreement with the geometric field pattern are all checked.
     */
   virtual void buildGlobalUnknowns(const std::vector<Teuchos::RCP<GlobalIndexer>> & fieldBlockManagers);

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
   void registerFields(bool buildSubUGIs);

   /** Has the method <code>registerFields</code> been called?
     */
   bool fieldsRegistered() const
   { return fieldsRegistered_; }

   /** Extract the field DOFManagers used underneath to define the
     * global unknowns.
     */
   const std::vector<Teuchos::RCP<GlobalIndexer>> &
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

   /** Enable TieBreak in sub dofmanger
     */
   void enableTieBreak(bool useTieBreak)
   { useTieBreak_ = useTieBreak; }

   /** \brief How any GIDs are associate with a particular element block
     *
     * This is a per-element count. If you have a quad element with two
     * piecewise bi-linear fields this method returns 8.
     */
   virtual int getElementBlockGIDCount(const std::string & blockId) const;

   /** \brief How any GIDs are associate with a particular element block.
     *
     * This is a per-element count. If you have a quad element with two
     * piecewise bi-linear fields this method returns 8.
     */
   virtual int getElementBlockGIDCount(const std::size_t & blockIndex) const;

protected:

   /** Build a new indexer. The concrete type is specified internally by this object (FEI version standard)
     */
   Teuchos::RCP<GlobalIndexer> buildNewIndexer(const Teuchos::RCP<ConnManager> & connManager,
                                                     MPI_Comm mpiComm) const;

   /** Do appropriate casting below and set orientations for a particular indexer. (handles FEI versus standard DOFManager)
     */
   void setOrientationsRequired(const Teuchos::RCP<GlobalIndexer> & indexer,bool required) const;

   /** Do appropriate casting below and call buildGlobalUnknowns for a particular indexer. (handles FEI versus standard DOFManager)
     */
   void buildGlobalUnknowns(const Teuchos::RCP<GlobalIndexer> & indexer,const Teuchos::RCP<const FieldPattern> & geomPattern) const;

   /** Do appropriate casting below and call getElementBlockGIDCount for a particular indexer. (handles FEI versus standard DOFManager)
     */
   int getElementBlockGIDCount(const Teuchos::RCP<GlobalIndexer> & indexer,const std::string & elementBlock) const;

   /** Do appropriate casting below and call getElementBlockGIDCount for a particular indexer. (handles FEI versus standard DOFManager)
     */
   int getElementBlockGIDCount(const Teuchos::RCP<GlobalIndexer> & indexer,const std::size_t & blockIndex) const;

   /** Do appropriate casting below and call printFieldInformation for a particular indexer. (handles FEI versus standard DOFManager)
     */
   void printFieldInformation(const Teuchos::RCP<GlobalIndexer> & indexer,std::ostream & os) const;

   /** This routine calls the <code>addField</code> method on the fieldBlockManager adding all
     * the fields it is supposed to control, and then calls registerFields.
     *
     * This method assumes that the activeFields are a legitimate ordering for the local field block.
     */
   void addFieldsToFieldBlockManager(const std::vector<std::string> & activeFields,
                                     GlobalIndexer & fieldBlockManager) const;

   /** This routine calls the <code>addField</code> method on the fieldBlockManager adding all
     * the fields it is supposed to control, and then calls registerFields.
     *
     * This method assumes that the activeFields are a legitimate ordering for the local field block.
     */
   void addFieldsToFieldBlockManager(const std::vector<std::string> & activeFields,
                                     DOFManager & fieldBlockManager) const;

   // computes connectivity
   Teuchos::RCP<ConnManager> connMngr_;

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
  //std::unordered_set<GlobalOrdinal,panzer::pair_hash> ownedGIDHashTable_;

   std::vector<std::vector<std::string> > fieldOrder_;

   bool fieldsRegistered_;

   Teuchos::RCP<const FieldPattern> geomPattern_;
   Teuchos::RCP<Teuchos::MpiComm<int> > communicator_;

   std::vector<Teuchos::RCP<GlobalIndexer>> fieldBlockManagers_;

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

   bool useDOFManagerFEI_;
   bool useTieBreak_;
};

}

#endif

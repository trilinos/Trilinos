#ifndef __Panzer_DOFManager_hpp__
#define __Panzer_DOFManager_hpp__

#include <map>

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

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
   Teuchos::RCP<const FieldPattern> getFieldPattern(const std::string & blockId, const std::string & fieldName) const
   { return getFieldPattern(blockId,getFieldNum(fieldName)); }
   
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
   Teuchos::RCP<const ConnManager<LocalOrdinalT,GlobalOrdinalT> > getConnManager() const 
   { return connMngr_; } 

   /** build the global unknown numberings
     *   1. this builds the pattens
     *   2. initializes the connectivity
     *   3. calls initComplete
     */
   virtual void buildGlobalUnknowns();

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
   void getElementGIDs(LocalOrdinalT localElmtId,std::vector<GlobalOrdinalT> & gids) const;

   /** \brief Use the field pattern so that you can find a particular
     *        field in the GIDs array.
     */
   const std::vector<int> & getGIDFieldOffsets(const std::string & blockId,int fieldNum) const
   { return fieldAggPattern_.find(blockId)->second->localOffsets(fieldNum); }

   /** \brief Use the field pattern so that you can find a particular
     *        field in the GIDs array. This version lets you specify the sub
     *        cell you are interested in. 
     *
     * \param[in] blockId
     * \param[in] fieldNum
     * \param[in] subCellDim
     * \param[in] subCellId
     */
   const std::vector<int> & getGIDFieldOffsets(const std::string & blockId,int fieldNum,int subCellDim,int subCellId) const
   { TEUCHOS_ASSERT(false); } 

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
   const std::vector<int> & getGIDFieldOffsets_closure(const std::string & blockId,int fieldNum,int subcellDim,int subcellId) const
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
 

   //! \defgroup LA_DS_Info Methods that are specific to the linear algebra
   //@{

   //! get the map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getMap() const;

   //! get the overlapped map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> getOverlapMap() const;

   //! get the graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> getGraph() const;

   //! get the overlapped graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> getOverlapGraph() const;

   //! get importer for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Import> getOverlapImport() const;

   //! get exporter for converting an overalapped object to a "normal" object
   virtual const Teuchos::RCP<Epetra_Export> getOverlapExport() const;

   //@}

protected:
   /** Build the default field ordering: simply uses ordering
     * imposed by <code>fieldStrToInt_</code> (alphabetical on field name)
     */
   void buildDefaultFieldOrder();
   
   std::vector<int> getOrderedBlock(const std::string & blockId);

   /** Using the natural ordering associated with the std::vector
     * retrieved from the connection manager
     */
   std::size_t blockIdToIndex(const std::string & blockId) const;

   //! build the pattern associated with this manager
   void buildPattern(const std::string & blockId,const Teuchos::RCP<FieldPattern> & geomPattern);

   // get the map from the matrix
   virtual const Teuchos::RCP<Epetra_Map> buildMap() const;
   virtual const Teuchos::RCP<Epetra_Map> buildOverlapMap() const;

   // get the graph of the crs matrix
   virtual const Teuchos::RCP<Epetra_CrsGraph> buildGraph() const;
   virtual const Teuchos::RCP<Epetra_CrsGraph> buildOverlapGraph() const;

   // computes connectivity
   Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > connMngr_; 
   
   //! \defgroup MapFunctions Mapping objects
   //@{ 
   //! field string ==> field id
   std::map<std::string,int> fieldStrToInt_;
   std::map<int,std::string> intToFieldStr_;

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

   // storage for Epetra graphs and maps
   Teuchos::RCP<Epetra_Comm> comm_;
   mutable Teuchos::RCP<Epetra_Map> map_;
   mutable Teuchos::RCP<Epetra_Map> overlappedMap_;
   mutable Teuchos::RCP<Epetra_CrsGraph> graph_;
   mutable Teuchos::RCP<Epetra_CrsGraph> overlappedGraph_;

   // maps blockIds to indices
   mutable Teuchos::RCP<std::map<std::string,std::size_t> > blockIdToIndex_;

   std::vector<int> fieldOrder_;

   // counters
   int nodeType_;
   int edgeType_;
   int numFields_;
   std::vector<int> patternNum_;
};

}

#include "Panzer_DOFManagerT.hpp"

#endif

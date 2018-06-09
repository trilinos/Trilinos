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

#ifndef __Panzer_DOFManager_hpp__
#define __Panzer_DOFManager_hpp__
#include <map>

#include <mpi.h>

#include "PanzerDofMgr_config.hpp"
#include "Panzer_FieldPattern.hpp"
#include "Panzer_FieldAggPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_NodeType.hpp"
#include "Panzer_FieldType.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Teuchos_RCP.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"

namespace panzer {


template <typename LocalOrdinalT, typename GlobalOrdinalT>
class DOFManager : public UniqueGlobalIndexer<LocalOrdinalT, GlobalOrdinalT> {
public:
  typedef GlobalOrdinalT GO;
  typedef LocalOrdinalT LO;

  virtual ~DOFManager() {}

  DOFManager();

  /** Constructor that sets the connection manager and communicator
    * objects. This is equivalent to calling the default constructor and
    * then "setConnManager(...)" routine.
    */
  DOFManager(const Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > & connMngr,MPI_Comm mpiComm);

  //! Adds a Connection Manager that will be associated with this DOFManager.
  void setConnManager(const Teuchos::RCP<ConnManager<LO,GO> > & connMngr, MPI_Comm mpiComm);

  Teuchos::RCP<ConnManager<LO,GO> > getConnManager() const
  { return connMngr_; }

  virtual Teuchos::RCP<const ConnManagerBase<LocalOrdinalT> > getConnManagerBase() const
  { return getConnManager(); }

  /** \brief Add a field to the DOF manager.
    *
    * Add a field to the DOF manager. Immediately after
    * adding the field the field number and field size
    * will be available for a user to access
    *
    * \param[in] str Human readable name of the field
    * \param[in] pattern Pattern defining the basis function to be used
    * \param[in] type Type of the Field (CG/DG) for generating GIDs
    *
    * \note <code>addField</code> cannot be called after <code>buildGlobalUnknowns</code> 
    *       or <code>registerFields</code>.
    */
  int addField(const std::string & str, const Teuchos::RCP<const FieldPattern> & pattern,
               const panzer::FieldType& type = panzer::FieldType::CG);

  /** \brief Add a field with an option for specifying the block.
    *
    * Add a field to the DOF manager. Immediately after
    * adding the field the field number and field size
    * will be available for a user to access
    *
    * \param[in] blockID Name of the element block that this field should be added to
    * \param[in] str Human readable name of the field
    * \param[in] pattern Pattern defining the basis function to be used
    * \param[in] type Type of the Field (CG/DG) for generating GIDs
    *
    * \note <code>addField</code> cannot be called after <code>buildGlobalUnknowns</code> 
    *       or <code>registerFields</code>.
    */
  int addField(const std::string & blockID, const std::string & str,
               const Teuchos::RCP<const FieldPattern> & pattern,
               const panzer::FieldType& type = panzer::FieldType::CG);

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
  Teuchos::RCP<const FieldPattern> getFieldPattern(const std::string & name) const;

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

  /**
   *  \brief Get the set of indices owned by this processor.
   *
   *  \param[out] indices A `vector` that will be fille with the indices owned
   *                      by this processor.
   */
  void
  getOwnedIndices(
    std::vector<GlobalOrdinalT>& indices) const;

  /**
   *  \brief Get the set of indices ghosted for this processor.
   *
   *  \param[out] indices A `vector` that will be fille with the indices
   *                      ghosted for this processor.
   */
  void
  getGhostedIndices(
    std::vector<GlobalOrdinalT>& indices) const;

  /**
   *  \brief Get the set of owned and ghosted indices for this processor.
   *
   *  \param[out] indices A `vector` that will be fille with the owned and
   *                      ghosted indices for this processor.
   */
  void
  getOwnedAndGhostedIndices(
    std::vector<GlobalOrdinalT>& indices) const;

  /**
   *  \brief Get the number of indices owned by this processor.
   *
   *  \returns The number of indices owned by this processor.
   */
  int
  getNumOwned() const;

  /**
   *  \brief Get the number of indices ghosted for this processor.
   *
   *  \returns The number of indices ghosted for this processor.
   */
  int
  getNumGhosted() const;

  /**
   *  \brief Get the number of owned and ghosted indices for this processor.
   *
   *  \returns The number of owned and ghosted indices for this processor.
   */
  int
  getNumOwnedAndGhosted() const;

  //! gets the number of fields
  int getNumFields() const;

  /** gets the field pattern so you can find a particular
    * field in the GIDs array.
    */
  const std::vector<int> & getGIDFieldOffsets(const std::string & blockID, int fieldNum) const;

  /** gets the field pattern so you can find a particular
    * field in the GIDs array.
    */
  const Kokkos::View<const int*,PHX::Device> getGIDFieldOffsetsKokkos(const std::string & blockID, int fieldNum) const;

  //! get associated GIDs for a given local element
  void getElementGIDs(LO localElementID, std::vector<GO> & gids, const std::string & blockIdHint="") const;

  //! builds the global unknowns array
  void buildGlobalUnknowns();

  //! builds the global unknowns array
  void buildGlobalUnknowns(const Teuchos::RCP<const FieldPattern> & geomPattern);
  
  int getFieldNum(const std::string & string) const;

  Teuchos::RCP<Teuchos::Comm<int> > getComm() const
  { return communicator_; }

  Teuchos::RCP<const FieldPattern> getGeometricFieldPattern() const
  { return ga_fp_; }

  void getElementBlockIds(std::vector<std::string> & elementBlockIds) const
  { connMngr_->getElementBlockIds(elementBlockIds); }
  
  bool fieldInBlock(const std::string & field, const std::string & block) const;

  const std::vector<int> & getBlockFieldNumbers(const std::string & blockId) const;

//************************************************************************

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
  const std::pair<std::vector<int>,std::vector<int> > & 
  getGIDFieldOffsets_closure(const std::string & blockId, int fieldNum, int subcellDim,int subcellId) const;

  //! Get the owned element block
  const std::vector<LocalOrdinalT> & getElementBlock(const std::string & blockId) const
  { return connMngr_->getElementBlock(blockId); }

  void ownedIndices(const std::vector<GlobalOrdinalT> & indices,std::vector<bool> & isOwned) const;

  void setFieldOrder(const std::vector<std::string> & fieldOrder );

  void getFieldOrder(std::vector<std::string> & fieldOrder) const;

  bool validFieldOrder(const std::vector<std::string> & proposed_fieldOrder);

  //TODO:this
  void buildUnknownsOrientation();

  bool getOrientationsRequired() const
  { return requireOrientations_; }

  void setOrientationsRequired(bool ro)
  { requireOrientations_ = ro; }

   /** \brief Get a vector containg the orientation of the GIDs relative to the neighbors.
     */
  void getElementOrientation(LocalOrdinalT localElmtId,std::vector<double> & gidsOrientation) const;

  const std::string & getFieldString(int num) const;

  /** \brief Reset the indices for this DOF manager.
    *
    * This method resets the indices and wipes out internal state. This method
    * does preserve the fields and the patterns added to the object. Also the
    * old connection manager is returned.
    *
    * \returns Old connection manager.
    */
  Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > resetIndices();

  /** \brief How any GIDs are associate with a particular element block
    *
    * This is a per-element count. If you have a quad element with two
    * piecewise bi-linear fields this method returns 8.
    */
  virtual int getElementBlockGIDCount(const std::string & blockId) const
  { return getElementBlockGIDCount(blockIdToIndex(blockId)); }

  /** \brief How any GIDs are associate with a particular element block.
    *
    * This is a per-element count. If you have a quad element with two
    * piecewise bi-linear fields this method returns 8.
    */
  virtual int getElementBlockGIDCount(const std::size_t & blockIndex) const
  { return elementBlockGIDCount_[blockIndex]; }

  /** Prints to an output stream the information about
    * the aggregated field.
    */
  void printFieldInformation(std::ostream & os) const;

  /** Turn on/off the use of a tie break object in the
    * createOneToOne algorithm. Turning this one gives 
    * better load balancing.
    */
  void enableTieBreak(bool enable)   
  { useTieBreak_ = enable; }

  /** Turn on/off the use of neighbor elements in the construction of the
    * global ids. If on, then the ghosted GIDs will include GIDs from neighbor
    * elements, and you will be able to call getElement(G/L)IDs for elements in
    * the one ring of this processor.
    */
  void useNeighbors(bool flag)   
  { useNeighbors_ = flag; }

  // These functions are primarily for testing purposes
  // they are not intended to be useful otherwise (thus they are not 
  // documented in the Doxygen style

  // Return the number of elemnts as measured by the count of GID arrays.
  // note that this will include ghosted elements!
  std::size_t getNumberElementGIDArrays() const
  { return elementGIDs_.size(); }

protected:

  /** Use Zoltan2 to locally reorder with RCM.
    */
  Teuchos::RCP<const Tpetra::Map<LO,GO,panzer::TpetraNodeType> >
  runLocalRCMReordering(const Teuchos::RCP<const Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,panzer::TpetraNodeType> > &);

  /** Using the natural ordering associated with the std::vector
    * retrieved from the connection manager
    */
  std::size_t blockIdToIndex(const std::string & blockId) const;

  /** This small struct is a utility meant to unify access
    * to elements and allow better code reuse. Basically
    * it provides a switch between the neighbor element blocks
    * and the owned element blocks.
    */
  class ElementBlockAccess {
    bool useOwned_;
    Teuchos::RCP<const ConnManager<LO,GO> > connMngr_;
  public:
    ElementBlockAccess(bool owned,const Teuchos::RCP<const ConnManager<LO,GO> > & connMngr) 
      : useOwned_(owned), connMngr_(connMngr) {}
    
    const std::vector<LO> & getElementBlock(const std::string & eBlock) const 
    {
      if(useOwned_==true)
        return connMngr_->getElementBlock(eBlock);
      else
        return connMngr_->getNeighborElementBlock(eBlock);
    }
  };

  /** Build the overlapped communication map given an element access object.
    * This map is used to construct the GIDs, and also to communicate the used
    * GIDs. (this is steps 1 and 2)
    */
  Teuchos::RCP<const Tpetra::Map<LO,GO,panzer::TpetraNodeType> >
  buildOverlapMapFromElements(const ElementBlockAccess & access) const; 

  /** Build a tagged multivector (as defined in GUN paper) to use in global unknown numbering algorithm.
    * Note that this is non-const. It does modify the <code>elementBlockGIDCount</code> member variable.
    */
  Teuchos::RCP<Tpetra::MultiVector<GO,LO,GO,panzer::TpetraNodeType> >
  buildTaggedMultiVector(const ElementBlockAccess & access);

  /** Build global unknowns using the algorithm in the Global Unknowns Numbering paper (GUN). This  
    * returns a non-overlapped multi-vector with the unique global IDs as owned by this processor. The input
    * tagged overlapped multi-vector (<code>overlap_mv</code>) is overwritten with the global IDs. Note
    * fields on geometric entities that are not assigned a global ID are given an entry of -1.
    */
  std::pair<Teuchos::RCP<Tpetra::MultiVector<GO,LO,GO,panzer::TpetraNodeType> >,
            Teuchos::RCP<Tpetra::MultiVector<GO,LO,GO,panzer::TpetraNodeType> > >
  buildGlobalUnknowns_GUN(const Tpetra::MultiVector<GO,LO,GO,panzer::TpetraNodeType> & tagged_overlap_mv,
                          Tpetra::MultiVector<GO,LO,GO,panzer::TpetraNodeType> & overlap_mv) const;

  void fillGIDsFromOverlappedMV(const ElementBlockAccess & access,
                                std::vector<std::vector< GO > > & elementGIDs,
                                const Tpetra::Map<LO,GO,panzer::TpetraNodeType> & overlapmap,
                                const Tpetra::MultiVector<GO,LO,GO,panzer::TpetraNodeType> & overlap_mv) const;
  void buildLocalIdsFromOwnedAndGhostedElements();
  
  Teuchos::RCP<ConnManager<LO,GO> > connMngr_;
  Teuchos::RCP<Teuchos::Comm<int> > communicator_;

  //Please note: AID=absolute ID. This is an attempt to remember that
  // fieldPatterns_ is unchanging storage for FPs.
  std::vector<Teuchos::RCP<const FieldPattern> > fieldPatterns_;
  std::vector<FieldType> fieldTypes_; // FieldType for a Field Pattern. Use AID to access just like fieldPatterns_.
  std::map<std::string,int> fieldNameToAID_;

  std::vector<std::string> blockOrder_; // To be got from the ConnManager.
  std::map<std::string,int> blockNameToID_; // I'm not sure the above vector is needed, this might suffice.
  std::vector<std::vector<int> > blockToAssociatedFP_; // each sub-vector is associated by
  // a block, with ordering given in blockOrder_. ints refer to the order in fieldPatterns_;
  std::vector<std::string> fieldStringOrder_;
  std::vector<int> fieldAIDOrder_; // Both of these must be updated and edited together.
  // The AID offers a simpler way to manage FPs internally.

  Teuchos::RCP<const panzer::FieldPattern> ga_fp_; // geometric aggregate field pattern
  std::vector<Teuchos::RCP<panzer::FieldAggPattern> > fa_fps_; //Ordered by blockOrder_;

  std::vector<GO> owned_;
  std::vector<GO> ghosted_;

  // Element GIDS ordered by LID.
  std::vector<std::vector< GO > > elementGIDs_;

  // Mimics the functionality of the getElemenentBlockGIDCount in
  // the original DOFManager. Indexed according to blockOrder_.
  std::vector<int> elementBlockGIDCount_;

  int numFields_;

  bool buildConnectivityRun_;

  bool requireOrientations_;
  std::vector<std::vector<signed char> > orientation_;

  bool useTieBreak_;
  bool useNeighbors_;
};

}

#endif

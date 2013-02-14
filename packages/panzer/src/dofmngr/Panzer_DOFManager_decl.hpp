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

#include "mpi.h"

#include "Panzer_config.hpp"
#include "Panzer_FieldPattern.hpp"
#include "Panzer_FieldAggPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"

#include "Teuchos_RCP.hpp"

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
    * then "setConnManager" routine.
    */
  DOFManager(const Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > & connMngr,MPI_Comm mpiComm);

  //! Adds a Connection Manager that will be associated with this DOFManager.
  void setConnManager(const Teuchos::RCP<ConnManager<LO,GO> > & connMngr, MPI_Comm mpiComm);

  Teuchos::RCP<ConnManager<LO,GO> > getConnManager() const
  { return connMngr_; }

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
  int addField(const std::string & str, const Teuchos::RCP<const FieldPattern> & pattern);

  //! Adds a field with an option for specifying the block.
  int addField(const std::string & blockID, const std::string & str, const Teuchos::RCP<const FieldPattern> & pattern);


  /** Returns the fieldpattern of the given name
    * This could also be done using the number you'd get from getFieldNum which
    * isn't yet included.
    */
  Teuchos::RCP<const FieldPattern> getFieldPattern(const std::string & name) const;

  void getOwnedIndices(std::vector<GlobalOrdinalT> & indices) const;

  void getOwnedAndSharedIndices(std::vector<GlobalOrdinalT> & indices) const;

  //! gets the number of fields
  int getNumFields() const;

  /** gets the field pattern so you can find a particular
    * field in the GIDs aray.
    */
  const std::vector<int> & getGIDFieldOffsets(const std::string & blockID, int fieldNum) const;

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
  
  /** Because we only have one element block, if the field is present, and
    * the block is valid. It works
    */
  bool fieldInBlock(const std::string & field, const std::string & block) const;

  /** Because we only have one element block, we are guarenteed 
    * that all fields are going to be in the one block we have
    */
  const std::vector<int> & getBlockFieldNumbers(const std::string & blockId) const;

//************************************************************************

  //TODO:this.
  const std::pair<std::vector<int>,std::vector<int> > & 
  getGIDFieldOffsets_closure(const std::string & blockId, int fieldNum, int subcellDim,int subcellId) const;

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

  //TODO:this
  void getElementOrientation(LocalOrdinalT localElmtId,std::vector<double> & gidsOrientation) const;

  const std::string & getFieldString(int num) const;

  /** \brief Reset the indicies for this DOF manager.
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
  int getElementBlockGIDCount(const std::string & blockId) const
  { return getElementBlockGIDCount(blockIdToIndex(blockId)); }

  /** \brief How any GIDs are associate with a particular element block.
    *
    * This is a per-element count. If you have a quad element with two
    * piecewise bi-linear fields this method returns 8.
    */
  int getElementBlockGIDCount(const std::size_t & blockIndex) const
  { return elementBlockGIDCount_[blockIndex]; }

  /** Prints to an output stream the information about
    * the aggregated field.
    */
  void printFieldInformation(std::ostream & os) const;

protected:

   /** Using the natural ordering associated with the std::vector
     * retrieved from the connection manager
     */
   std::size_t blockIdToIndex(const std::string & blockId) const;
  
  Teuchos::RCP<ConnManager<LO,GO> > connMngr_;
  Teuchos::RCP<Teuchos::Comm<int> > communicator_;

  //Please note: AID=absolute ID. This is an attempt to remember that
  //fieldPatterns_ is unchanging storage for FPs.
  std::vector<Teuchos::RCP<const FieldPattern> > fieldPatterns_;
  std::map<std::string,int> fieldNameToAID_;
  std::vector<std::string> blockOrder_; //To be got from the ConnManager.
  std::map<std::string,int> blockNameToID_; //I'm not sure the above vector is needed, this might suffice.
  std::vector<std::vector<int> > blockToAssociatedFP_; //each sub-vector is associated by
  //a block, with ordering given in blockOrder_. ints refer to the order in fieldPatterns_;
  std::vector<int> FPsInAll_;
  std::vector<std::string> fieldStringOrder_;
  std::vector<int> fieldAIDOrder_; //Both of these must be updated and edited together.
  //The AID offers a simpler way to manage FPs internally.

  Teuchos::RCP<const panzer::FieldPattern> ga_fp_; // geometric aggregate field pattern
  std::vector<Teuchos::RCP<panzer::FieldAggPattern> > fa_fps_; //Ordered by blockOrder_;

  std::vector<GO> owned_;
  std::vector<GO> owned_and_ghosted_;

  //Element GIDS ordered by LID.
  std::vector<std::vector< GO > > elementGIDs_;

  //Mimics the functionality of the getElemenentBlockGIDCount in
  //the original DOFManager. Indexed according to blockOrder_.
  std::vector<int> elementBlockGIDCount_;

  int numFields_;

  bool buildConnectivityRun_;

  bool requireOrientations_;
  std::vector<std::vector<char> > orientation_;
};

}

#endif

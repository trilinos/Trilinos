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

#ifndef PANZER_DOF_MANAGER2_IMPL_HPP
#define PANZER_DOF_MANAGER2_IMPL_HPP

#include <map>

#include "mpi.h"

#include "PanzerDofMgr_config.hpp"
#include "Panzer_FieldPattern.hpp"
#include "Panzer_FieldAggPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_UniqueGlobalIndexer_Utilities.hpp"
#include "Panzer_DOFManager_Functors.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Panzer_NodalFieldPattern.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"

#include <unordered_set> // a hash table

/*
#define HAVE_ZOLTAN2
#ifdef HAVE_ZOLTAN2
#include "Zoltan2_XpetraCrsGraphInput.hpp"
#include "Zoltan2_OrderingProblem.hpp"
#endif
*/

namespace panzer {

namespace {
template <typename LocalOrdinal,typename GlobalOrdinal>
class HashTieBreak : public Tpetra::Details::TieBreak<LocalOrdinal,GlobalOrdinal> {
  const unsigned int seed_;

public:
  HashTieBreak(const unsigned int seed = (2654435761U))
    : seed_(seed) { }

  virtual std::size_t selectedIndex(GlobalOrdinal GID,
      const std::vector<std::pair<int,LocalOrdinal> > & pid_and_lid) const
  {
    // this is Epetra's hash/Tpetra's default hash: See Tpetra HashTable
    int intkey = (int) ((GID & 0x000000007fffffffLL) +
       ((GID & 0x7fffffff80000000LL) >> 31));
    return std::size_t((seed_ ^ intkey) % pid_and_lid.size());
  }
};

}

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayRCP;
using Teuchos::Array;
using Teuchos::ArrayView;

template <typename LO, typename GO>
DOFManager<LO,GO>::DOFManager()
  : numFields_(0),buildConnectivityRun_(false),requireOrientations_(false), useTieBreak_(false), buildGhosted_(false)
{ }

template <typename LO, typename GO>
DOFManager<LO,GO>::DOFManager(const Teuchos::RCP<ConnManager<LO,GO> > & connMngr,MPI_Comm mpiComm)
  : numFields_(0),buildConnectivityRun_(false),requireOrientations_(false), useTieBreak_(false), buildGhosted_(false)
{
  setConnManager(connMngr,mpiComm);
}

template <typename LO, typename GO>
void DOFManager<LO,GO>::setConnManager(const Teuchos::RCP<ConnManager<LO,GO> > & connMngr, MPI_Comm mpiComm)
{
  TEUCHOS_TEST_FOR_EXCEPTION(buildConnectivityRun_,std::logic_error,
                      "DOFManager::setConnManager: setConnManager cannot be called after "
                      "buildGlobalUnknowns has been called");
  connMngr_ = connMngr;
  //We acquire the block ordering from the connmanager.
  connMngr_->getElementBlockIds(blockOrder_);
  for (size_t i = 0; i < blockOrder_.size(); ++i) {
    blockNameToID_.insert(std::map<std::string,int>::value_type(blockOrder_[i],i));
    //We must also initialize vectors for FP associations.
  }
  blockToAssociatedFP_.resize(blockOrder_.size());
  communicator_ = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(mpiComm)));
}


//Adds a field to be used in creating the Global Numbering
//Returns the index for the field pattern
template <typename LO, typename GO>
int DOFManager<LO,GO>::addField(const std::string & str, const Teuchos::RCP<const FieldPattern> & pattern)
{
  TEUCHOS_TEST_FOR_EXCEPTION(buildConnectivityRun_,std::logic_error,
                      "DOFManager::addField: addField cannot be called after "
                      "buildGlobalUnknowns has been called");

  fieldPatterns_.push_back(pattern);
  fieldNameToAID_.insert(std::map<std::string,int>::value_type(str, numFields_));

  //The default values for IDs are the sequential order they are added in.
  fieldStringOrder_.push_back(str);
  fieldAIDOrder_.push_back(numFields_);

  for(std::size_t i=0;i<blockOrder_.size();i++) {
    blockToAssociatedFP_[i].push_back(numFields_);
  }

  ++numFields_;
  return numFields_-1;
}

template <typename LO, typename GO>
int DOFManager<LO,GO>::addField(const std::string & blockID, const std::string & str, const Teuchos::RCP<const FieldPattern> & pattern)
{
  TEUCHOS_TEST_FOR_EXCEPTION(buildConnectivityRun_,std::logic_error,
                      "DOFManager::addField: addField cannot be called after "
                      "buildGlobalUnknowns has been called");
  TEUCHOS_TEST_FOR_EXCEPTION((connMngr_==Teuchos::null),std::logic_error,
                             "DOFManager::addField: you must add a ConnManager before"
                             "you can associate a FP with a given block.")
  bool found=false;
  size_t blocknum=0;
  while(!found && blocknum<blockOrder_.size()){
    if(blockOrder_[blocknum]==blockID){
      found=true;
      break;
    }
    blocknum++;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(!found,std::logic_error, "DOFManager::addField: Invalid block name.");

  //This will be different if the FieldPattern is already present.
  //We need to check for that.
  found=false;
  std::map<std::string,int>::const_iterator fpIter = fieldNameToAID_.find(str);
  if(fpIter!=fieldNameToAID_.end())
    found=true;

  if(!found){
    fieldPatterns_.push_back(pattern);
    fieldNameToAID_.insert(std::map<std::string,int>::value_type(str, numFields_));
    //The default values for IDs are the sequential order they are added in.
    fieldStringOrder_.push_back(str);
    fieldAIDOrder_.push_back(numFields_);

    //This is going to be associated with blocknum.
    blockToAssociatedFP_[blocknum].push_back(numFields_);
    ++numFields_;
    return numFields_-1;
  }
  else{
    blockToAssociatedFP_[blocknum].push_back(fpIter->second);
    return numFields_;
  }
}



  //Returns the fieldpattern of the given name
template <typename LO, typename GO>
Teuchos::RCP<const FieldPattern> DOFManager<LO,GO>::getFieldPattern(const std::string & name) const
{
  std::map<std::string,int>::const_iterator fitr = fieldNameToAID_.find(name);
  if(fitr==fieldNameToAID_.end())
    return Teuchos::null;

  if(fitr->second<int(fieldPatterns_.size()))
    return fieldPatterns_[fitr->second];

  return Teuchos::null;
}

//Returns the fieldpattern of the given name
template <typename LO, typename GO>
Teuchos::RCP<const FieldPattern> DOFManager<LO,GO>::getFieldPattern(const std::string & blockId, const std::string & fieldName) const
{
  std::map<std::string,int>::const_iterator fitr = fieldNameToAID_.find(fieldName);
  if(fitr==fieldNameToAID_.end())
    return Teuchos::null;

  bool found=false;
  size_t blocknum=0;
  while(!found && blocknum<blockOrder_.size()){
    if(blockOrder_[blocknum]==blockId){
      found=true;
      break;
    }
    blocknum++;
  }

  if(!found)
    return Teuchos::null;

  std::vector<int>::const_iterator itr = std::find(blockToAssociatedFP_[blocknum].begin(),
                                                   blockToAssociatedFP_[blocknum].end(),fitr->second);
  if(itr!=blockToAssociatedFP_[blocknum].end()) {
    if(fitr->second<int(fieldPatterns_.size()))
      return fieldPatterns_[fitr->second];
  }

  return Teuchos::null;
}

template <typename LO, typename GO>
void DOFManager<LO,GO>::getOwnedIndices(std::vector<GO> & indicies) const
{
  indicies.resize(owned_.size());
  for (size_t i = 0; i < owned_.size(); ++i) {
    indicies[i]=owned_[i];
  }
}

template <typename LO, typename GO>
void DOFManager<LO,GO>::getOwnedAndSharedIndices(std::vector<GO> & indicies) const
{
  indicies.resize(owned_and_ghosted_.size());
  for (size_t i = 0; i < owned_and_ghosted_.size(); ++i) {
    indicies[i]=owned_and_ghosted_[i];
  }
}

  //gets the number of fields
template <typename LO, typename GO>
int DOFManager<LO,GO>::getNumFields() const
{
  return numFields_;
}

template <typename LO, typename GO>
const std::vector<int> & DOFManager<LO,GO>::getGIDFieldOffsets(const std::string & blockID, int fieldNum) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!buildConnectivityRun_,std::logic_error, "DOFManager::getGIDFieldOffsets: cannot be called before buildGlobalUnknowns has been called");
  std::map<std::string,int>::const_iterator bitr = blockNameToID_.find(blockID);
  if(bitr==blockNameToID_.end())
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"DOFManager::fieldInBlock: invalid block name");
  int bid=bitr->second;
  if(fa_fps_[bid]!=Teuchos::null)
    return fa_fps_[bid]->localOffsets(fieldNum);

  static const std::vector<int> empty;
  return empty;
}

template <typename LO, typename GO>
void DOFManager<LO,GO>::getElementGIDs(LO localElementID, std::vector<GO> & gids, const std::string & blockIdHint) const
{
  gids = elementGIDs_[localElementID];
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::buildGlobalUnknowns()
{
  /* STEPS.
   * 1.  Build GA_FP and all block's FA_FP's and place into respective data structures.
   */
  if(requireOrientations_){
    fieldPatterns_.push_back(Teuchos::rcp(new NodalFieldPattern(fieldPatterns_[0]->getCellTopology())));
  }
  RCP<GeometricAggFieldPattern> aggFieldPattern = Teuchos::rcp(new GeometricAggFieldPattern);;
  aggFieldPattern = Teuchos::rcp(new GeometricAggFieldPattern(fieldPatterns_));

  connMngr_->buildConnectivity(*aggFieldPattern);

  // using new geometric pattern, build global unknowns
  buildGlobalUnknowns(aggFieldPattern);
}

//builds the global unknowns array
template <typename LO, typename GO>
void DOFManager<LO,GO>::buildGlobalUnknowns(const Teuchos::RCP<const FieldPattern> & geomPattern)
{
  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setOutputToRootOnly(-1);
  out.setShowProcRank(true);

  TEUCHOS_TEST_FOR_EXCEPTION(buildConnectivityRun_,std::logic_error,
                      "DOFManager::buildGlobalUnknowns: buildGlobalUnknowns cannot be called again "
                      "after buildGlobalUnknowns has been called");
  //Some stuff for the Map.
  typedef panzer::TpetraNodeType Node;
  typedef Tpetra::Map<LO, GO, Node> Map;

  typedef Tpetra::Export<LO,GO,Node> Export;
  typedef Tpetra::Import<LO,GO,Node> Import;

  //the GIDs are of type GO.
  typedef Tpetra::MultiVector<GO,LO,GO,Node> MultiVector;

  // this is a safety check to make sure that nodes are included
  // in the geometric field pattern when orientations are required
  if(getOrientationsRequired()) {
    std::size_t sz = geomPattern->getSubcellIndices(0,0).size();

    TEUCHOS_TEST_FOR_EXCEPTION(sz==0,std::logic_error,
                               "DOFManager::buildGlobalUnknowns requires a geometric pattern including "
                                 "the nodes when orientations are needed!");
  }

  RCP<const Teuchos::Comm<int> > comm = communicator_;

  /* STEPS.
   * 1.  Build all block's FA_FP's and place into respective data structures.
   */
  ga_fp_ = geomPattern;

  //We will iterate through all of the blocks, building a FieldAggPattern for
  //each of them.

  for (size_t b = 0; b < blockOrder_.size(); ++b) {
    std::vector<std::pair< int, RCP<const panzer::FieldPattern> > > faConstruct;
    //The ID is going to be the AID, and then everything will work.
    //The ID should not be the AID, it should be the ID it has in the ordering.

    for (size_t i = 0; i < fieldAIDOrder_.size(); ++i) {
      int looking = fieldAIDOrder_[i];

      //Check if in b's fp list
      std::vector<int>::const_iterator reu = std::find(blockToAssociatedFP_[b].begin(), blockToAssociatedFP_[b].end(), looking);
      if(!(reu==blockToAssociatedFP_[b].end())){
        faConstruct.push_back(std::make_pair(i, fieldPatterns_[fieldAIDOrder_[i]]));
      }

    }

    if(faConstruct.size()>0) {
      fa_fps_.push_back(rcp(new FieldAggPattern(faConstruct, ga_fp_)));

      // how many global IDs are in this element block?
      int gidsInBlock = fa_fps_[fa_fps_.size()-1]->numberIds();
      elementBlockGIDCount_.push_back(gidsInBlock);
    }
    else {
      fa_fps_.push_back(Teuchos::null);
      elementBlockGIDCount_.push_back(0);
    }
  }

  // given a set of elements over each element block build an overlap
  // map that will provide the required element entities for the
  // set of elements requested.
  ElementBlockAccess ownedAccess(true,connMngr_);

  /* Steps 2. and 3.
   */
  RCP<const Map> overlapmap       = buildOverlapMapFromElements(ownedAccess);

 /* 4.  Create an overlapped multivector from the overlap map.
   */
  Teuchos::RCP<MultiVector> overlap_mv;
  overlap_mv = Tpetra::createMultiVector<GO>(overlapmap,(size_t)numFields_);

 /* 5.  Iterate through all local elements again, checking with the FP
   *     information. Mark up the overlap map accordingly.
   */


  ArrayRCP<ArrayRCP<GO> > edittwoview = overlap_mv->get2dViewNonConst();
  for (size_t b = 0; b < blockOrder_.size(); ++b) {
    // there has to be a field pattern assocaited with the block
    if(fa_fps_[b]==Teuchos::null)
      continue;

    const std::vector<LO> & numFields= fa_fps_[b]->numFieldsPerId();
    const std::vector<LO> & fieldIds= fa_fps_[b]->fieldIds();
    const std::vector<LO> & myElements = connMngr_->getElementBlock(blockOrder_[b]);
    for (size_t l = 0; l < myElements.size(); ++l) {
      LO connSize = connMngr_->getConnectivitySize(myElements[l]);
      const GO * elmtConn = connMngr_->getConnectivity(myElements[l]);
      int offset=0;
      for (int c = 0; c < connSize; ++c) {
        size_t lid = overlapmap->getLocalElement(elmtConn[c]);
        for (int n = 0; n < numFields[c]; ++n) {
          int whichField = fieldIds[offset];
          //Row will be lid. column will be whichField.
          //Shove onto local ordering
          edittwoview[whichField][lid]=1;
          offset++;
        }
      }
    }
  }

 /* 6.  Create a OneToOne map from the overlap map.
   */

  RCP<const Map> non_overlap_map;
  if(!useTieBreak_) {
    non_overlap_map = Tpetra::createOneToOne<LO,GO,Node>(overlapmap);
  }
  else {
    // use a hash tie break to get better load balancing from create one to one
    HashTieBreak<LO,GO> tie_break;
    non_overlap_map = Tpetra::createOneToOne<LO,GO,Node>(overlapmap,tie_break);
  }

 /* 7.  Create a non-overlapped multivector from OneToOne map.
   */
  Teuchos::RCP<MultiVector> non_overlap_mv;
  non_overlap_mv = Tpetra::createMultiVector<GO>(non_overlap_map,(size_t)numFields_);

 /* 8.  Create an export between the two maps.
   */
  Export e(overlapmap,non_overlap_map);

  // Note:  ETP 04/26/16  Temporarily create an importer for all of the
  // doImport() calls below.  This works around mysterious failures when
  // using the exporter for Cuda builds.
  Import imp(non_overlap_map,overlapmap);

 /* 9.  Export data using ABSMAX.
   */
  non_overlap_mv->doExport(*overlap_mv,e,Tpetra::ABSMAX);

 /* 10. Compute the local sum using Kokkos.
   */
  GO localsum=0;
  {
    typedef typename Tpetra::MultiVector<GO,Node> MV;
    typedef typename MV::dual_view_type::t_dev KV;
    typedef typename MV::dual_view_type::t_dev::memory_space DMS;
    KV values = non_overlap_mv->template getLocalView<DMS>();
    auto mv_size = values.dimension_0();
    Kokkos::parallel_reduce(mv_size,panzer::dof_functors::SumRank2<GO,KV>(values),localsum);
  }

 /* 11. Create a map using local sums to generate final GIDs.
   */
  RCP<const Map> gid_map =
    rcp (new Map (Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid (),
                  static_cast<size_t> (localsum), static_cast<GO> (0), comm));
  // mfh 28 Apr 2015: This doesn't work because createContigMap
  // assumes the default Node type, but Panzer might use a different
  // Node type.  Just call the Map constructor; don't call those
  // nonmember "constructors."
  //RCP<const Map> gid_map = Tpetra::createContigMap<LO,GO>(-1,localsum, comm);

 /* 12. Iterate through the non-overlapping MV and assign GIDs to
   *     the necessary points. (Assign a -1 elsewhere.)
   */
  ArrayView<const GO> owned_ids = gid_map->getNodeElementList();
  int which_id=0;
  ArrayRCP<ArrayRCP<GO> > editnonoverlap = non_overlap_mv->get2dViewNonConst();
  for(size_t i=0; i<non_overlap_mv->getLocalLength(); ++i){
    for(int j=0; j<numFields_; ++j){
      if(editnonoverlap[j][i]!=0){
        editnonoverlap[j][i]=owned_ids[which_id];
        which_id++;
      }
      else{
        editnonoverlap[j][i]=-1;
      }

    }
  }

 /* 13. Import data back to the overlap MV using REPLACE.
   */
  overlap_mv->doImport(*non_overlap_mv,imp,Tpetra::REPLACE);

 /* 14. Cross reference local element connectivity and overlap map to
   *     create final GID vectors.
   */

  // this bit of code takes the uniquely assigned GIDs and spreads them
  // out for processing by local element ID
  fillGIDsFromOverlappedMV(ownedAccess,elementGIDs_,*overlapmap,*overlap_mv);

  // if neighbor unknowns are required, then make sure they are included
  // in the elementGIDs_
  if (buildGhosted_) { // enabling this turns on GID construction for
                       // neighbor processors
    ElementBlockAccess naborAccess(false,connMngr_);
    RCP<const Map> overlapmap_nabor = buildOverlapMapFromElements(naborAccess);

    // Export e(overlapmap_nabor,non_overlap_map);
    Import imp_nabor(non_overlap_map,overlapmap_nabor);

    Teuchos::RCP<MultiVector> overlap_mv_nabor
        = Tpetra::createMultiVector<GO>(overlapmap_nabor,(size_t)numFields_);

    // get all neighbor information
    overlap_mv_nabor->doImport(*non_overlap_mv,imp_nabor,Tpetra::REPLACE);

    fillGIDsFromOverlappedMV(naborAccess,elementGIDs_,*overlapmap_nabor,*overlap_mv_nabor);
  }

  //////////////////////////////////////////////////////////////////
  // this is where the code is modified to artificially induce GIDs
  // over 2 Billion unknowns
  //////////////////////////////////////////////////////////////////
  #if 0
  {
    panzer::Ordinal64 offset = 0xFFFFFFFFLL;

    for (size_t b = 0; b < blockOrder_.size(); ++b) {
      const std::vector<LO> & myElements = connMngr_->getElementBlock(blockOrder_[b]);
      for (std::size_t l = 0; l < myElements.size(); ++l) {
        std::vector<GO> & localGIDs = elementGIDs_[myElements[l]];
        for(std::size_t c=0;c<localGIDs.size();c++)
          localGIDs[c] += offset;
      }
    }

    Teuchos::ArrayRCP<GO> nvals = non_overlap_mv->get1dViewNonConst();
    for (int j = 0; j < nvals.size(); ++j)
      nvals[j] += offset;
  }
  #endif

  // build owned vector
  {
    typedef std::unordered_set<GO> HashTable;
    HashTable isOwned, remainingOwned;
    //owned_ is made up of owned_ids.
    Teuchos::ArrayRCP<const GO> nvals = non_overlap_mv->get1dView();
    for (int j = 0; j < nvals.size(); ++j) {
      if(nvals[j]!=-1) {
        isOwned.insert(nvals[j]);
        remainingOwned.insert(nvals[j]);
      }
    }

    HashTable hashTable; // use to detect if global ID has been added to owned_and_ghosted_
    for (size_t b = 0; b < blockOrder_.size(); ++b) {

      if(fa_fps_[b]==Teuchos::null)
        continue;

      const std::vector<LO> & myElements = connMngr_->getElementBlock(blockOrder_[b]);

      for (size_t l = 0; l < myElements.size(); ++l) {
        const std::vector<GO> & localOrdering = elementGIDs_[myElements[l]];

        // add "novel" global ids that are also owned to owned array
        for(std::size_t i=0;i<localOrdering.size();i++) {
          // don't try to add if ID is not owned
          if(isOwned.find(localOrdering[i])==isOwned.end())
            continue;

          // only add a global ID if it has not been added
          std::pair<typename HashTable::iterator,bool> insertResult = hashTable.insert(localOrdering[i]);
          if(insertResult.second) {
            owned_.push_back(localOrdering[i]);
            remainingOwned.erase(localOrdering[i]);
          }
        }
      }
    }

    // add any other owned GIDs that were not already included.
    // these are ones that are owned locally but not required by any
    // element owned by this processor (uggh!)
    for(typename HashTable::const_iterator itr=remainingOwned.begin();itr!=remainingOwned.end();itr++)
      owned_.push_back(*itr);

    if(owned_.size()!=isOwned.size()) {
      out << "I'm about to hang because of unknown numbering failure ... sorry! (line = " << __LINE__ << ")" << std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(owned_.size()!=isOwned.size(),std::logic_error,
                                 "DOFManager::buildGlobalUnkonwns: Failure because not all owned unknowns have been accounted for.");
    }

  }

  // build owned and ghosted array: The old simple way led to slow
  // Jacobian assembly, the new way speeds up Jacobian assembly
  {
    // loop over all elements. do greedy ordering of local values over elements for
    // building owned_and_ghosted, hopefully this gives a better layout
    // for element ordered assembly
    typedef std::unordered_set<GO> HashTable;
    HashTable hashTable; // use to detect if global ID has been added to owned_and_ghosted_

    for(std::size_t i=0;i<owned_.size();i++) {
      hashTable.insert(owned_[i]);
      owned_and_ghosted_.push_back(owned_[i]);
    }

    // this cute trick (blah) of constructing a accessor vector is to eliminate a copy
    // of the block of code computing shared/owned DOFs
    std::vector<ElementBlockAccess> blockAccessVec;
    blockAccessVec.push_back(ElementBlockAccess(true,connMngr_));
    if(buildGhosted_)
      blockAccessVec.push_back(ElementBlockAccess(false,connMngr_));
    // all owned will be processed first followed by those that are
    // optionally ghosted

    for(std::size_t a=0;a < blockAccessVec.size(); a++) {
      // get access type (owned or neighbor)
      const ElementBlockAccess & access = blockAccessVec[a];

      for (size_t b = 0; b < blockOrder_.size(); ++b) {
        if(fa_fps_[b]==Teuchos::null)
          continue;

        const std::vector<LO> & myElements = access.getElementBlock(blockOrder_[b]);

        for (size_t l = 0; l < myElements.size(); ++l) {
          const std::vector<GO> & localOrdering = elementGIDs_[myElements[l]];

          // add "novel" global ids into owned_and_ghosted_ vector.
          for(std::size_t i=0;i<localOrdering.size();i++) {
            std::pair<typename HashTable::iterator,bool> insertResult = hashTable.insert(localOrdering[i]);

            // if insertion succeeds, then this is "novel" to owned_and_ghosted_
            // vector so include it
            if(insertResult.second)
              owned_and_ghosted_.push_back(localOrdering[i]);
          }
        }
      }
    }
  }

  buildConnectivityRun_ = true;

  // build orientations if required
  if(requireOrientations_)
    buildUnknownsOrientation();

  // allocate the local IDs
  if (buildGhosted_)
    this->buildLocalIdsFromOwnedAndSharedElements();
  else
    this->buildLocalIds();
}

template <typename LO, typename GO>
int DOFManager<LO,GO>::getFieldNum(const std::string & string) const
{
  int ind=0;
  bool found=false;
  while(!found && (size_t)ind<fieldStringOrder_.size()){
    if(fieldStringOrder_[ind]==string)
      found=true;
    else
      ind++;
  }

  if(found)
    return ind;

  // didn't find anything...return -1
  return -1;
}

template <typename LO, typename GO>
void DOFManager<LO,GO>::getFieldOrder(std::vector<std::string> & fieldOrder) const
{
  fieldOrder.resize(fieldStringOrder_.size());
  for (size_t i = 0; i < fieldStringOrder_.size(); ++i)
    fieldOrder[i]=fieldStringOrder_[i];
}

template <typename LO, typename GO>
bool DOFManager<LO,GO>::fieldInBlock(const std::string & field, const std::string & block) const
{
  std::map<std::string,int>::const_iterator fitr = fieldNameToAID_.find(field);
  if(fitr==fieldNameToAID_.end()) {
    std::stringstream ss;
    printFieldInformation(ss);
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"DOFManager::fieldInBlock: invalid field name. DOF information is:\n"+ss.str());
  }
  std::map<std::string,int>::const_iterator bitr = blockNameToID_.find(block);
  if(bitr==blockNameToID_.end()) {
    std::stringstream ss;
    printFieldInformation(ss);
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"DOFManager::fieldInBlock: invalid block name. DOF information is:\n"+ss.str());
  }
  int fid=fitr->second;
  int bid=bitr->second;

  bool found=false;
  for (size_t i = 0; i < blockToAssociatedFP_[bid].size(); ++i) {
    if(blockToAssociatedFP_[bid][i]==fid){
      found=true;
      break;
    }
  }

  return found;
}

template <typename LO, typename GO>
const std::vector<int> & DOFManager<LO,GO>::getBlockFieldNumbers(const std::string & blockId) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!buildConnectivityRun_,std::logic_error,"DOFManager::getBlockFieldNumbers: BuildConnectivity must be run first.");

  std::map<std::string,int>::const_iterator bitr = blockNameToID_.find(blockId);
  if(bitr==blockNameToID_.end())
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"DOFManager::fieldInBlock: invalid block name");
  int bid=bitr->second;

  // there has to be a field pattern assocaited with the block
  if(fa_fps_[bid]!=Teuchos::null)
    return fa_fps_[bid]->fieldIds();

  // nothing to return
  static std::vector<int> empty;
  return empty;
}

template <typename LO, typename GO>
const std::pair<std::vector<int>,std::vector<int> > &
DOFManager<LO,GO>::getGIDFieldOffsets_closure(const std::string & blockId, int fieldNum, int subcellDim,int subcellId) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!buildConnectivityRun_,std::logic_error,"DOFManager::getGIDFieldOffsets_closure: BuildConnectivity must be run first.");
  std::map<std::string,int>::const_iterator bitr = blockNameToID_.find(blockId);
  if(bitr==blockNameToID_.end())
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error, "DOFManager::getGIDFieldOffsets_closure: invalid block name.");

  // there has to be a field pattern assocaited with the block
  if(fa_fps_[bitr->second]!=Teuchos::null)
    return fa_fps_[bitr->second]->localOffsets_closure(fieldNum, subcellDim, subcellId);

  static std::pair<std::vector<int>,std::vector<int> > empty;
  return empty;
}

template <typename LO, typename GO>
void DOFManager<LO,GO>::ownedIndices(const std::vector<GO> & indices,std::vector<bool> & isOwned) const
{
  //Resizes the isOwned array.
  if(indices.size()!=isOwned.size())
    isOwned.resize(indices.size(),false);
  typename std::vector<GO>::const_iterator endOf = owned_.end();
  for (std::size_t i = 0; i < indices.size(); ++i) {
    isOwned[i] = ( std::find(owned_.begin(), owned_.end(), indices[i])!=endOf );
  }
}

template <typename LO, typename GO>
void DOFManager<LO,GO>::setFieldOrder(const std::vector<std::string> & fieldOrder)
{
  TEUCHOS_TEST_FOR_EXCEPTION(buildConnectivityRun_,std::logic_error,
                      "DOFManager::setFieldOrder: setFieldOrder cannot be called after "
                      "buildGlobalUnknowns has been called");
  if(validFieldOrder(fieldOrder)){
    fieldStringOrder_=fieldOrder;
    //We also need to reassign the fieldAIDOrder_.
    for (size_t i = 0; i < fieldStringOrder_.size(); ++i) {
      fieldAIDOrder_[i]=fieldNameToAID_[fieldStringOrder_[i]];
    }
  }
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"DOFManager::setFieldOrder: Invalid Field Ordering!");

}


template <typename LO, typename GO>
bool DOFManager<LO,GO>::validFieldOrder(const std::vector<std::string> & proposed_fieldOrder)
{
  if(fieldStringOrder_.size()!=proposed_fieldOrder.size())
    return false;
  //I'm using a not very efficient way of doing this, but there should never
  //be that many fields, so it isn't that big of a deal.
  for (size_t i = 0; i < fieldStringOrder_.size(); ++i) {
    bool found=false;
    for (size_t j = 0; j < proposed_fieldOrder.size(); ++j) {
      if(fieldStringOrder_[i]==proposed_fieldOrder[j])
        found=true;
    }
    if(!found)
      return false;
  }
  return true;
}

template <typename LO, typename GO>
const std::string & DOFManager<LO,GO>::getFieldString(int num) const
{
  //A reverse lookup through fieldStringOrder_.
  if(num>=(int)fieldStringOrder_.size())
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "DOFManager::getFieldString: invalid number");
  return fieldStringOrder_[num];
}

//Everything associated with orientation is not yeat built, but this
//is the method as found in Panzer_DOFManager_impl.hpp. There are
//going to need to be some substantial changes to the code as it applies
//to this DOFManager, but the basic ideas and format should be similar.
//
template <typename LO, typename GO>
void DOFManager<LO,GO>::buildUnknownsOrientation()
{
  orientation_.clear(); // clean up previous work

  std::vector<std::string> elementBlockIds;
  connMngr_->getElementBlockIds(elementBlockIds);

  // figure out how many total elements are owned by this processor (why is this so hard!)
  std::size_t myElementCount = 0;
  for(std::vector<std::string>::const_iterator blockItr=elementBlockIds.begin(); blockItr!=elementBlockIds.end();++blockItr)
    myElementCount += connMngr_->getElementBlock(*blockItr).size();

  // allocate for each block
  orientation_.resize(myElementCount);

  // loop over all element blocks
  for(std::vector<std::string>::const_iterator blockItr=elementBlockIds.begin();
    blockItr!=elementBlockIds.end();++blockItr) {
      const std::string & blockName = *blockItr;

     // this block has no unknowns (or elements)
    std::map<std::string,int>::const_iterator fap = blockNameToID_.find(blockName);
    if(fap==blockNameToID_.end()) {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"DOFManager::buildUnknownsOrientation: invalid block name");
    }

    int bid=fap->second;

    if(fa_fps_[bid]==Teuchos::null)
      continue;

     // grab field patterns, will be necessary to compute orientations
    const FieldPattern & fieldPattern = *fa_fps_[bid];

    //Should be ga_fp_ (geometric aggregate field pattern)
    std::vector<std::pair<int,int> > topEdgeIndices;
    orientation_helpers::computePatternEdgeIndices(*ga_fp_,topEdgeIndices);

    // grab face orientations if 3D
    std::vector<std::vector<int> > topFaceIndices;
    if(ga_fp_->getDimension()==3)
      orientation_helpers::computePatternFaceIndices(*ga_fp_,topFaceIndices);

    //How many GIDs are associated with a particular element bloc
    const std::vector<LO> & elmts = connMngr_->getElementBlock(blockName);
    for(std::size_t e=0;e<elmts.size();e++) {
       // this is the vector of orientations to fill: initialize it correctly
      std::vector<char> & eOrientation = orientation_[elmts[e]];

      // This resize seems to be the same as fieldPattern.numberIDs().
      // When computer edge orientations is called, that is the assert.
      // There should be no reason to make it anymore complicated.
      eOrientation.resize(fieldPattern.numberIds());
      for(std::size_t s=0;s<eOrientation.size();s++)
        eOrientation[s] = 1; // put in 1 by default

      // get geometry ids
      LO connSz = connMngr_->getConnectivitySize(elmts[e]);
      const GO * connPtr = connMngr_->getConnectivity(elmts[e]);
      const std::vector<GO> connectivity(connPtr,connPtr+connSz);

      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity, fieldPattern, eOrientation);

      // compute face orientations in 3D
      if(ga_fp_->getDimension()==3)
        orientation_helpers::computeCellFaceOrientations(topFaceIndices, connectivity, fieldPattern, eOrientation);
    }
  }
}

template <typename LO, typename GO>
void DOFManager<LO,GO>::getElementOrientation(LO localElmtId,std::vector<double> & gidsOrientation) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(orientation_.size()==0,std::logic_error,
                              "DOFManager::getElementOrientations: Orientations were not constructed!");

   const std::vector<char> & local_o = orientation_[localElmtId];
   gidsOrientation.resize(local_o.size());
   for(std::size_t i=0;i<local_o.size();i++) {
      gidsOrientation[i] = double(local_o[i]);
   }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > DOFManager<LocalOrdinalT,GlobalOrdinalT>::resetIndices()
{
   Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > connMngr = connMngr_;

   connMngr_ = Teuchos::null;

   // wipe out FEI objects
   orientation_.clear(); // clean up previous work
   fa_fps_.clear();
   elementGIDs_.clear();
   owned_.clear();
   owned_and_ghosted_.clear();
   elementBlockGIDCount_.clear();

   return connMngr;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
std::size_t DOFManager<LocalOrdinalT,GlobalOrdinalT>::blockIdToIndex(const std::string & blockId) const
{
  std::map<std::string,int>::const_iterator bitr = blockNameToID_.find(blockId);
  if(bitr==blockNameToID_.end())
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"DOFManager::fieldInBlock: invalid block name");
  return bitr->second;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::printFieldInformation(std::ostream & os) const
{
  os << "DOFManager Field Information: " << std::endl;

  // sanity check
  TEUCHOS_ASSERT(blockOrder_.size()==blockToAssociatedFP_.size());

  for(std::size_t i=0;i<blockOrder_.size();i++) {
    os << "   Element Block = " << blockOrder_[i] << std::endl;

    // output field information
    const std::vector<int> & fieldIds = blockToAssociatedFP_[i];
    for(std::size_t f=0;f<fieldIds.size();f++)
      os << "      \"" << getFieldString(fieldIds[f]) << "\" is field ID " << fieldIds[f] << std::endl;
  }
}

template <typename LO,typename GO>
Teuchos::RCP<const Tpetra::Map<LO,GO,panzer::TpetraNodeType> >
DOFManager<LO,GO>::
buildOverlapMapFromElements(const ElementBlockAccess & access) const
{
  /*
   * 2.  Iterate through all local elements and create the overlapVector
   *     of concerned elements.
   */

  std::set<GO> overlapset;
  for (size_t i = 0; i < blockOrder_.size(); ++i) {
    const std::vector<LO> & myElements = access.getElementBlock(blockOrder_[i]);
    for (size_t e = 0; e < myElements.size(); ++e) {
      LO connSize = connMngr_->getConnectivitySize(myElements[e]);
      const GO * elmtConn = connMngr_->getConnectivity(myElements[e]);
      for (int k = 0; k < connSize; ++k) {
        overlapset.insert(elmtConn[k]);
      }
    }
  }

  Array<GO> overlapVector;
  for (typename std::set<GO>::const_iterator itr = overlapset.begin(); itr!=overlapset.end(); ++itr) {
    overlapVector.push_back(*itr);
  }

  /* 3.  Construct an overlap map from this structure.
   */
  return Tpetra::createNonContigMap<LO,GO>(overlapVector,communicator_);
}

template <typename LO,typename GO>
void DOFManager<LO,GO>::
fillGIDsFromOverlappedMV(const ElementBlockAccess & access,
                         std::vector<std::vector< GO > > & elementGIDs,
                         const Tpetra::Map<LO,GO,panzer::TpetraNodeType> & overlapmap,
                         const Tpetra::MultiVector<GO,LO,GO,panzer::TpetraNodeType> & overlap_mv) const
{
  using Teuchos::ArrayRCP;

  //To generate elementGIDs we need to go through all of the local elements.
  ArrayRCP<ArrayRCP<const GO> > twoview = overlap_mv.get2dView();

  //And for each of the things in fa_fp.fieldIds we go to that column. To the the row,
  //we move from globalID to localID in the map and use our local value for something.
  for (size_t b = 0; b < blockOrder_.size(); ++b) {
    const std::vector<LO> & myElements = access.getElementBlock(blockOrder_[b]);

    if(fa_fps_[b]==Teuchos::null) {
      // fill elements that are not used with empty vectors
      for (size_t l = 0; l < myElements.size(); ++l) {
        LO thisID=myElements[l];
        if(elementGIDs.size()<=(size_t)thisID)
          elementGIDs.resize(thisID+1);
      }
      continue;
    }

    const std::vector<int> & numFields= fa_fps_[b]->numFieldsPerId();
    const std::vector<int> & fieldIds= fa_fps_[b]->fieldIds();
    //
    //
    for (size_t l = 0; l < myElements.size(); ++l) {
      LO connSize = connMngr_->getConnectivitySize(myElements[l]);
      const GO * elmtConn = connMngr_->getConnectivity(myElements[l]);
      std::vector<GO> localOrdering;
      int offset=0;
      for (int c = 0; c < connSize; ++c) {
        size_t lid = overlapmap.getLocalElement(elmtConn[c]);
        for (int n = 0; n < numFields[c]; ++n) {
          int whichField = fieldIds[offset];
          offset++;
          //Row will be lid. column will be whichField.
          //Shove onto local ordering
          localOrdering.push_back(twoview[whichField][lid]);
        }
      }
      LO thisID=myElements[l];
      if(elementGIDs.size()<=(size_t)thisID){
        elementGIDs.resize(thisID+1);
      }
      elementGIDs[thisID]=localOrdering;
    }
  }
}

template <typename LO, typename GO>
void DOFManager<LO,GO>::buildLocalIdsFromOwnedAndSharedElements()
{
  std::vector<std::vector<LO> > elementLIDs(elementGIDs_.size());

  std::vector<GO> ownedAndShared;
  this->getOwnedAndSharedIndices(ownedAndShared);

  // build global to local hash map (temporary and used only once)
  std::unordered_map<GO,LO> hashMap;
  for(std::size_t i = 0; i < ownedAndShared.size(); ++i)
    hashMap[ownedAndShared[i]] = i;

  for (std::size_t i = 0; i < elementGIDs_.size(); ++i) {
    const std::vector<GO>& gids = elementGIDs_[i];
    std::vector<LO>& lids = elementLIDs[i];
    lids.resize(gids.size());
    for (std::size_t g = 0; g < gids.size(); ++g)
      lids[g] = hashMap[gids[g]];
  }

  this->setLocalIds(elementLIDs);
}

/*
template <typename LO,typename GO>
Teuchos::RCP<const Tpetra::Map<LO,GO,panzer::TpetraNodeType> >
DOFManager<LO,GO>::runLocalRCMReordering(const Teuchos::RCP<const Tpetra::Map<LO,GO,panzer::TpetraNodeType> > & map)
{
  typedef panzer::TpetraNodeType Node;
  typedef Tpetra::Map<LO, GO, Node> Map;
  typedef Tpetra::CrsGraph<LO, GO, Node> Graph;

  Teuchos::RCP<Graph> graph = Teuchos::rcp(new Graph(map,0));

  // build Crs Graph from the mesh
  for (size_t b = 0; b < blockOrder_.size(); ++b) {
    if(fa_fps_[b]==Teuchos::null)
      continue;

    const std::vector<LO> & myElements = connMngr_->getElementBlock(blockOrder_[b]);
    for (size_t l = 0; l < myElements.size(); ++l) {
      LO connSize = connMngr_->getConnectivitySize(myElements[l]);
      const GO * elmtConn = connMngr_->getConnectivity(myElements[l]);
      for (int c = 0; c < connSize; ++c) {
        LO lid = map->getLocalElement(elmtConn[c]);
        if(Teuchos::OrdinalTraits<LO>::invalid()!=lid)
          continue;

        graph->insertGlobalIndices(elmtConn[c],Teuchos::arrayView(elmtConn,connSize));
      }
    }
  }

  graph->fillComplete();

  std::vector<GO> newOrder(map->getNodeNumElements());
  {
    // graph is constructed, now run RCM using zoltan2
    typedef Zoltan2::XpetraCrsGraphInput<Graph> SparseGraphAdapter;

    Teuchos::ParameterList params;
    params.set("order_method", "rcm");
    SparseGraphAdapter adapter(graph);

    Zoltan2::OrderingProblem<SparseGraphAdapter> problem(&adapter, &params,MPI_COMM_SELF);
    problem.solve();

    // build a new global ording array using permutation
    Zoltan2::OrderingSolution<GO,LO> * soln = problem.getSolution();

    size_t dummy;
    size_t checkLength = soln->getPermutationSize();
    LO * checkPerm = soln->getPermutation(&dummy);

    Teuchos::ArrayView<const GO > oldOrder = map->getNodeElementList();
    TEUCHOS_ASSERT(checkLength==oldOrder.size());
    TEUCHOS_ASSERT(checkLength==newOrder.size());

    for(size_t i=0;i<checkLength;i++)
      newOrder[checkPerm[i]] = oldOrder[i];
  }

  return Tpetra::createNonContigMap<LO,GO>(newOrder,communicator_);
}
*/

} /*panzer*/

#endif

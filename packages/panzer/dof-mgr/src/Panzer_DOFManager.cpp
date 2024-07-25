// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_DOF_MANAGER2_IMPL_HPP
#define PANZER_DOF_MANAGER2_IMPL_HPP

#include <map>
#include <set>

#include <mpi.h>

#include "PanzerDofMgr_config.hpp"
#include "Panzer_FieldPattern.hpp"
#include "Panzer_FieldAggPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_GlobalIndexer_Utilities.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_DOFManager_Functors.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Panzer_NodalFieldPattern.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"

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

///////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////
namespace {
template <typename LocalOrdinal,typename GlobalOrdinal>
class GreedyTieBreak : public Tpetra::Details::TieBreak<LocalOrdinal,GlobalOrdinal> {

public:
  GreedyTieBreak() { }

  virtual bool mayHaveSideEffects() const {
    return true;
  }

  virtual std::size_t selectedIndex(GlobalOrdinal /* GID */,
                                    const std::vector<std::pair<int,LocalOrdinal> > & pid_and_lid) const
  {
    // always choose index of pair with smallest pid
    const std::size_t numLids = pid_and_lid.size();
    std::size_t idx = 0;
    int minpid = pid_and_lid[0].first;
    std::size_t minidx = 0;
    for (idx = 0; idx < numLids; ++idx) {
      if (pid_and_lid[idx].first < minpid) {
        minpid = pid_and_lid[idx].first;
        minidx = idx;
      }
    }
    return minidx;
  }
};

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayRCP;
using Teuchos::Array;
using Teuchos::ArrayView;

///////////////////////////////////////////////////////////////////////////////
DOFManager::DOFManager()
  : numFields_(0),buildConnectivityRun_(false),requireOrientations_(false), useTieBreak_(false), useNeighbors_(false)
{ }

///////////////////////////////////////////////////////////////////////////////
DOFManager::DOFManager(const Teuchos::RCP<ConnManager> & connMngr,MPI_Comm mpiComm)
  : numFields_(0),buildConnectivityRun_(false),requireOrientations_(false), useTieBreak_(false), useNeighbors_(false)
{
  setConnManager(connMngr,mpiComm);
}

///////////////////////////////////////////////////////////////////////////////
void DOFManager::setConnManager(const Teuchos::RCP<ConnManager> & connMngr, MPI_Comm mpiComm)
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

///////////////////////////////////////////////////////////////////////////////
//Adds a field to be used in creating the Global Numbering
//Returns the index for the field pattern
int DOFManager::addField(const std::string & str, const Teuchos::RCP<const FieldPattern> & pattern,
                         const panzer::FieldType& type)
{
  TEUCHOS_TEST_FOR_EXCEPTION(buildConnectivityRun_,std::logic_error,
                      "DOFManager::addField: addField cannot be called after "
                      "buildGlobalUnknowns has been called");

  fieldPatterns_.push_back(pattern);
  fieldTypes_.push_back(type);
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

///////////////////////////////////////////////////////////////////////////////
int DOFManager::addField(const std::string & blockID, const std::string & str, const Teuchos::RCP<const FieldPattern> & pattern,
                         const panzer::FieldType& type)
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
    fieldTypes_.push_back(type);
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

///////////////////////////////////////////////////////////////////////////////
Teuchos::RCP<const FieldPattern> DOFManager::getFieldPattern(const std::string & name) const
{
  std::map<std::string,int>::const_iterator fitr = fieldNameToAID_.find(name);
  if(fitr==fieldNameToAID_.end())
    return Teuchos::null;

  if(fitr->second<int(fieldPatterns_.size()))
    return fieldPatterns_[fitr->second];

  return Teuchos::null;
}

///////////////////////////////////////////////////////////////////////////////
Teuchos::RCP<const FieldPattern> DOFManager::getFieldPattern(const std::string & blockId, const std::string & fieldName) const
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

///////////////////////////////////////////////////////////////////////////////
void DOFManager::getOwnedIndices(std::vector<panzer::GlobalOrdinal>& indices) const
{
  indices = owned_;
}

///////////////////////////////////////////////////////////////////////////////
void DOFManager::getGhostedIndices(std::vector<panzer::GlobalOrdinal>& indices) const
{
  indices = ghosted_;
}

///////////////////////////////////////////////////////////////////////////////
void DOFManager::getOwnedAndGhostedIndices(std::vector<panzer::GlobalOrdinal>& indices) const
{
  using std::size_t;
  indices.resize(owned_.size() + ghosted_.size());
  for (size_t i(0); i < owned_.size(); ++i)
    indices[i] = owned_[i];
  for (size_t i(0); i < ghosted_.size(); ++i)
    indices[owned_.size() + i] = ghosted_[i];
}

///////////////////////////////////////////////////////////////////////////////
void DOFManager::getElementGIDsAsInt(panzer::LocalOrdinal localElementID, std::vector<int>& gids, const std::string& /* blockIdHint */) const
{
  const auto & element_ids = elementGIDs_[localElementID];
  gids.resize(element_ids.size());
  for (std::size_t i=0; i < gids.size(); ++i)
    gids[i] = element_ids[i];
}

///////////////////////////////////////////////////////////////////////////////
void DOFManager::getOwnedIndicesAsInt(std::vector<int>& indices) const
{
  indices.resize(owned_.size());
  for (std::size_t i=0; i < owned_.size(); ++i)
    indices[i] = owned_[i];
}

///////////////////////////////////////////////////////////////////////////////
void DOFManager::getGhostedIndicesAsInt(std::vector<int>& indices) const
{
  indices.resize(ghosted_.size());
  for (std::size_t i=0; i < ghosted_.size(); ++i)
    indices[i] = ghosted_[i];
}

///////////////////////////////////////////////////////////////////////////////
void DOFManager::getOwnedAndGhostedIndicesAsInt(std::vector<int>& indices) const
{
  indices.resize(owned_.size() + ghosted_.size());
  for (std::size_t i=0; i < owned_.size(); ++i)
    indices[i] = owned_[i];
  for (std::size_t i=0; i < ghosted_.size(); ++i)
    indices[owned_.size() + i] = ghosted_[i];
}

///////////////////////////////////////////////////////////////////////////////
int DOFManager::getNumOwned() const
{
  return owned_.size();
}

///////////////////////////////////////////////////////////////////////////////
int DOFManager::getNumGhosted() const
{
  return ghosted_.size();
}

///////////////////////////////////////////////////////////////////////////////
int DOFManager::getNumOwnedAndGhosted() const
{
  return owned_.size() + ghosted_.size();
}

///////////////////////////////////////////////////////////////////////////////
int DOFManager::getNumFields() const
{
  return numFields_;
}

///////////////////////////////////////////////////////////////////////////////
const std::vector<int> &
DOFManager::getGIDFieldOffsets(const std::string & blockID, int fieldNum) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!buildConnectivityRun_,std::logic_error, "DOFManager::getGIDFieldOffsets: cannot be called before "
                                                                      "buildGlobalUnknowns has been called");
  std::map<std::string,int>::const_iterator bitr = blockNameToID_.find(blockID);
  if(bitr==blockNameToID_.end())
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"DOFManager::fieldInBlock: invalid block name");
  int bid=bitr->second;
  if(fa_fps_[bid]!=Teuchos::null)
    return fa_fps_[bid]->localOffsets(fieldNum);

  static const std::vector<int> empty;
  return empty;
}

///////////////////////////////////////////////////////////////////////////////
const PHX::View<const int*>
DOFManager::getGIDFieldOffsetsKokkos(const std::string & blockID, int fieldNum) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!buildConnectivityRun_,std::logic_error, "DOFManager::getGIDFieldOffsets: cannot be called before "
                                                                      "buildGlobalUnknowns has been called");
  std::map<std::string,int>::const_iterator bitr = blockNameToID_.find(blockID);
  if(bitr==blockNameToID_.end()) {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"DOFManager::fieldInBlock: invalid block name");
  }

  int bid=bitr->second;
  if(fa_fps_[bid]!=Teuchos::null)
    return fa_fps_[bid]->localOffsetsKokkos(fieldNum);

  static const PHX::View<int*> empty("panzer::DOFManager::getGIDFieldOffsetsKokkos() empty",0);
  return empty;
}

///////////////////////////////////////////////////////////////////////////////
const PHX::ViewOfViews<1,PHX::View<const int*>>
DOFManager::getGIDFieldOffsetsKokkos(const std::string & blockID, const std::vector<int> & fieldNums) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!buildConnectivityRun_,std::logic_error, "DOFManager::getGIDFieldOffsets: cannot be called before "
                                                                      "buildGlobalUnknowns has been called");
  std::map<std::string,int>::const_iterator bitr = blockNameToID_.find(blockID);
  if(bitr==blockNameToID_.end()) {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"DOFManager::fieldInBlock: invalid block name");
  }

  PHX::ViewOfViews<1,PHX::View<const int*>> vov("panzer::getGIDFieldOffsetsKokkos vector version",fieldNums.size());
  vov.disableSafetyCheck(); // Its going to be moved/copied

  int bid=bitr->second;

  for (size_t i=0; i < fieldNums.size(); ++i) {
    if(fa_fps_[bid]!=Teuchos::null) {
      vov.addView(fa_fps_[bid]->localOffsetsKokkos(fieldNums[i]),i);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"panzer::DOFManager::getGIDFieldOffsets() - no field pattern exists in block " 
                                 << blockID << "for field number " << fieldNums[i] << " exists!");
      // vov.addView(PHX::View<int*>("panzer::DOFManager::getGIDFieldOffsetsKokkos() empty",0),i);
    }
  }

  vov.syncHostToDevice();

  return vov;
}

///////////////////////////////////////////////////////////////////////////////
void DOFManager::getElementGIDs(panzer::LocalOrdinal localElementID, std::vector<panzer::GlobalOrdinal>& gids, const std::string& /* blockIdHint */) const
{
  gids = elementGIDs_[localElementID];
}

///////////////////////////////////////////////////////////////////////////////
void DOFManager::buildGlobalUnknowns()
{
  /* STEPS.
   * 1.  Build GA_FP and all block's FA_FP's and place into respective data structures.
   */
  if(requireOrientations_){
    fieldPatterns_.push_back(Teuchos::rcp(new NodalFieldPattern(fieldPatterns_[0]->getCellTopology())));
    fieldTypes_.push_back(FieldType::CG);
  }

  TEUCHOS_ASSERT(fieldPatterns_.size() == fieldTypes_.size());
  std::vector<std::pair<FieldType,Teuchos::RCP<const FieldPattern>>> tmp;
  for (std::size_t i=0; i < fieldPatterns_.size(); ++i)
    tmp.push_back(std::make_pair(fieldTypes_[i],fieldPatterns_[i]));

  RCP<GeometricAggFieldPattern> aggFieldPattern = Teuchos::rcp(new GeometricAggFieldPattern(tmp));

  connMngr_->buildConnectivity(*aggFieldPattern);

  // using new geometric pattern, build global unknowns
  this->buildGlobalUnknowns(aggFieldPattern);
}

///////////////////////////////////////////////////////////////////////////////
void DOFManager::buildGlobalUnknowns(const Teuchos::RCP<const FieldPattern> & geomPattern)
{
  // some typedefs
  typedef panzer::TpetraNodeType Node;
  typedef Tpetra::Map<panzer::LocalOrdinal, panzer::GlobalOrdinal, Node> Map;

  typedef Tpetra::Import<panzer::LocalOrdinal,panzer::GlobalOrdinal,Node> Import;

  //the GIDs are of type panzer::GlobalOrdinal.
  typedef Tpetra::MultiVector<panzer::GlobalOrdinal,panzer::LocalOrdinal,panzer::GlobalOrdinal,Node> MultiVector;

  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManager::buildGlobalUnknowns",buildGlobalUnknowns);

  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setOutputToRootOnly(-1);
  out.setShowProcRank(true);

  TEUCHOS_TEST_FOR_EXCEPTION(buildConnectivityRun_,std::logic_error,
                      "DOFManager::buildGlobalUnknowns: buildGlobalUnknowns cannot be called again "
                      "after buildGlobalUnknowns has been called");

  // this is a safety check to make sure that nodes are included
  // in the geometric field pattern when orientations are required
  if(getOrientationsRequired()) {
    std::size_t sz = geomPattern->getSubcellIndices(0,0).size();

    TEUCHOS_TEST_FOR_EXCEPTION(sz==0,std::logic_error,
                               "DOFManager::buildGlobalUnknowns requires a geometric pattern including "
                                 "the nodes when orientations are needed!");
  }

  /* STEPS.
   * 1.  Build all block's FA_FP's and place into respective data structures.
   */
  ga_fp_ = geomPattern;

  // given a set of elements over each element block build an overlap
  // map that will provide the required element entities for the
  // set of elements requested.
  ElementBlockAccess ownedAccess(true,connMngr_);

  // INPUT: To the algorithm in the GUN paper
  RCP<MultiVector> tagged_overlap_mv = this->buildTaggedMultiVector(ownedAccess);
  RCP<const Map> overlap_map   = tagged_overlap_mv->getMap();

  RCP<MultiVector> overlap_mv = Tpetra::createMultiVector<panzer::GlobalOrdinal>(overlap_map,(size_t)numFields_);

  // call the GUN paper algorithm
  auto non_overlap_pair = this->buildGlobalUnknowns_GUN(*tagged_overlap_mv,*overlap_mv);
  RCP<MultiVector> non_overlap_mv = non_overlap_pair.first;
  RCP<MultiVector> tagged_non_overlap_mv = non_overlap_pair.second;
  RCP<const Map> non_overlap_map = non_overlap_mv->getMap();

 /* 14. Cross reference local element connectivity and overlap map to
   *     create final GID vectors.
   */

  // this bit of code takes the uniquely assigned GIDs and spreads them
  // out for processing by local element ID
  this->fillGIDsFromOverlappedMV(ownedAccess,elementGIDs_,*overlap_map,*overlap_mv);

  // if neighbor unknowns are required, then make sure they are included
  // in the elementGIDs_
  if (useNeighbors_) { // enabling this turns on GID construction for
                       // neighbor processors
    ElementBlockAccess neighborAccess(false,connMngr_);
    RCP<const Map> overlap_map_neighbor =
      this->buildOverlapMapFromElements(neighborAccess);

    // Export e(overlap_map_neighbor,non_overlap_map);
    Import imp_neighbor(non_overlap_map,overlap_map_neighbor);

    Teuchos::RCP<MultiVector> overlap_mv_neighbor =
      Tpetra::createMultiVector<panzer::GlobalOrdinal>(overlap_map_neighbor, (size_t)numFields_);

    // get all neighbor information
    overlap_mv_neighbor->doImport(*non_overlap_mv, imp_neighbor,
      Tpetra::REPLACE);

    this->fillGIDsFromOverlappedMV(neighborAccess, elementGIDs_,
      *overlap_map_neighbor, *overlap_mv_neighbor);
  }

  //////////////////////////////////////////////////////////////////
  // this is where the code is modified to artificially induce GIDs
  // over 2 Billion unknowns
  //////////////////////////////////////////////////////////////////
  #if 0
  {
    panzer::GlobalOrdinal offset = 0xFFFFFFFFLL;

    for (size_t b = 0; b < blockOrder_.size(); ++b) {
      const std::vector<panzer::LocalOrdinal> & myElements = connMngr_->getElementBlock(blockOrder_[b]);
      for (std::size_t l = 0; l < myElements.size(); ++l) {
        std::vector<panzer::GlobalOrdinal> & localGIDs = elementGIDs_[myElements[l]];
        for(std::size_t c=0;c<localGIDs.size();c++)
          localGIDs[c] += offset;
      }
    }

    Teuchos::ArrayRCP<panzer::GlobalOrdinal> nvals = non_overlap_mv->get1dViewNonConst();
    for (int j = 0; j < nvals.size(); ++j)
      nvals[j] += offset;
  }
  #endif

  // build owned vector
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManager::buildGlobalUnknowns::build_owned_vector",build_owned_vector);

    typedef std::unordered_set<panzer::GlobalOrdinal> HashTable;
    HashTable isOwned, remainingOwned;

    // owned_ is made up of owned_ids.: This doesn't work for high order
    auto nvals = non_overlap_mv->getLocalViewHost(Tpetra::Access::ReadOnly);
    auto tagged_vals = tagged_non_overlap_mv->getLocalViewHost(Tpetra::Access::ReadOnly);
    TEUCHOS_ASSERT(nvals.size()==tagged_vals.size());
    for (size_t i = 0; i < nvals.extent(1); ++i) {
      for (size_t j = 0; j < nvals.extent(0); ++j) {
	if (nvals(j,i) != -1) {
	  for(panzer::GlobalOrdinal offset=0;offset<tagged_vals(j,i);++offset)
	    isOwned.insert(nvals(j,i)+offset);
	}
	else {
	  // sanity check
	  TEUCHOS_ASSERT(tagged_vals(j,i)==0);
	}
      }
    }
    remainingOwned = isOwned;

    HashTable hashTable; // use to detect if global ID has been added to owned_
    for (size_t b = 0; b < blockOrder_.size(); ++b) {

      if(fa_fps_[b]==Teuchos::null)
        continue;

      const std::vector<panzer::LocalOrdinal> & myElements = connMngr_->getElementBlock(blockOrder_[b]);

      for (size_t l = 0; l < myElements.size(); ++l) {
        const std::vector<panzer::GlobalOrdinal> & localOrdering = elementGIDs_[myElements[l]];

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

  // Build the ghosted_ array.  The old simple way led to slow Jacobian
  // assembly; the new way speeds up Jacobian assembly.
  {
    // Loop over all the elements and do a greedy ordering of local values over
    // the elements for building ghosted_.  Hopefully this gives a better
    // layout for an element-ordered assembly.
    PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManager::buildGlobalUnknowns::build_ghosted_array",BGA);

    // Use a hash table to detect if global IDs have been added to owned_.
    typedef std::unordered_set<panzer::GlobalOrdinal> HashTable;
    HashTable hashTable;
    for (std::size_t i = 0; i < owned_.size(); i++)
      hashTable.insert(owned_[i]);

    // Here we construct an accessor vector, such that we first process
    // everything in the current element block, optionally followed by
    // everything in the neighbor element block.
    std::vector<ElementBlockAccess> blockAccessVec;
    blockAccessVec.push_back(ElementBlockAccess(true,connMngr_));
    if(useNeighbors_)
      blockAccessVec.push_back(ElementBlockAccess(false,connMngr_));
    for (std::size_t a = 0; a < blockAccessVec.size(); ++a)
    {
      // Get the access type (owned or neighbor).
      const ElementBlockAccess& access = blockAccessVec[a];
      for (size_t b = 0; b < blockOrder_.size(); ++b)
      {
        if (fa_fps_[b] == Teuchos::null)
          continue;
        const std::vector<panzer::LocalOrdinal>& myElements =
          access.getElementBlock(blockOrder_[b]);
        for (size_t l = 0; l < myElements.size(); ++l)
        {
          const std::vector<panzer::GlobalOrdinal>& localOrdering = elementGIDs_[myElements[l]];

          // Add "novel" global IDs into the ghosted_ vector.
          for (std::size_t i = 0; i < localOrdering.size(); ++i)
          {
            std::pair<typename HashTable::iterator, bool> insertResult =
              hashTable.insert(localOrdering[i]);

            // If the insertion succeeds, then this is "novel" to the owned_
            // and ghosted_ vectors, so include it in ghosted_.
            if(insertResult.second)
              ghosted_.push_back(localOrdering[i]);
          }
        }
      }
    }
  }

  buildConnectivityRun_ = true;

  // build orientations if required
  if(requireOrientations_) {
    PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManager::buildGlobalUnknowns::build_orientation",BO);
    buildUnknownsOrientation();
  }

  // allocate the local IDs
  if (useNeighbors_) {
    PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManager::buildGlobalUnknowns::build_local_ids_from_owned_and_ghosted",BLOFOG);
    this->buildLocalIdsFromOwnedAndGhostedElements();
  }
  else {
    PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManager::buildGlobalUnknowns::build_local_ids",BLI);
    this->buildLocalIds();
  }
}

///////////////////////////////////////////////////////////////////////////////
std::pair<Teuchos::RCP<Tpetra::MultiVector<panzer::GlobalOrdinal,panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> >,
          Teuchos::RCP<Tpetra::MultiVector<panzer::GlobalOrdinal,panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> > >
DOFManager::buildGlobalUnknowns_GUN(const Tpetra::MultiVector<panzer::GlobalOrdinal,panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> & tagged_overlap_mv,
                                    Tpetra::MultiVector<panzer::GlobalOrdinal,panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> & overlap_mv) const
{
  // some typedefs
  typedef panzer::TpetraNodeType Node;
  typedef Tpetra::Map<panzer::LocalOrdinal, panzer::GlobalOrdinal, Node> Map;

  typedef Tpetra::Export<panzer::LocalOrdinal,panzer::GlobalOrdinal,Node> Export;
  typedef Tpetra::Import<panzer::LocalOrdinal,panzer::GlobalOrdinal,Node> Import;

  //the GIDs are of type panzer::GlobalOrdinal.
  typedef Tpetra::MultiVector<panzer::GlobalOrdinal,panzer::LocalOrdinal,panzer::GlobalOrdinal,Node> MultiVector;

  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManager::buildGlobalUnknowns_GUN",BGU_GUN);

  // LINE 2: In the GUN paper
  RCP<const Map> overlap_map   = tagged_overlap_mv.getMap();

 /* 6.  Create a OneToOne map from the overlap map.
   */

  // LINE 4: In the GUN paper

  RCP<const Map> non_overlap_map;
  if(!useTieBreak_) {
    PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManager::buildGlobalUnknowns_GUN::line_04 createOneToOne",GUN04);

    GreedyTieBreak<panzer::LocalOrdinal,panzer::GlobalOrdinal> greedy_tie_break;
    non_overlap_map = Tpetra::createOneToOne<panzer::LocalOrdinal,panzer::GlobalOrdinal,Node>(overlap_map, greedy_tie_break);
  }
  else {
    // use a hash tie break to get better load balancing from create one to one
    // Aug. 4, 2016...this is a bad idea and doesn't work
    HashTieBreak<panzer::LocalOrdinal,panzer::GlobalOrdinal> tie_break;
    non_overlap_map = Tpetra::createOneToOne<panzer::LocalOrdinal,panzer::GlobalOrdinal,Node>(overlap_map,tie_break);
  }

 /* 7.  Create a non-overlapped multivector from OneToOne map.
   */

  // LINE 5: In the GUN paper

  Teuchos::RCP<MultiVector> tagged_non_overlap_mv;
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManager::buildGlobalUnknowns_GUN::line_05 alloc_unique_mv",GUN05);

    tagged_non_overlap_mv = Tpetra::createMultiVector<panzer::GlobalOrdinal>(non_overlap_map,(size_t)numFields_);
  }

 /* 8.  Create an export between the two maps.
   */

  // LINE 6: In the GUN paper
  RCP<Export> exp;
  RCP<Import> imp;
  RCP<MultiVector> non_overlap_mv;
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManager::buildGlobalUnknowns_GUN::line_06 export",GUN06);

    exp = rcp(new Export(overlap_map,non_overlap_map));

    /* 9.  Export data using ABSMAX.
      */
    tagged_non_overlap_mv->doExport(tagged_overlap_mv,*exp,Tpetra::ABSMAX);

    // copy the tagged one, so as to preserve the tagged MV so we can overwrite
    // the non_overlap_mv
    non_overlap_mv = rcp(new MultiVector(*tagged_non_overlap_mv,Teuchos::Copy));
  }

 /* 10. Compute the local sum using Kokkos.
   */

  // LINES 7-9: In the GUN paper

  panzer::GlobalOrdinal localsum=0;
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManager::buildGlobalUnknowns_GUN::line_07-09 local_count",GUN07_09);
    auto values = non_overlap_mv->getLocalViewDevice(Tpetra::Access::ReadOnly);
    panzer::dof_functors::SumRank2<panzer::GlobalOrdinal, decltype(values)>{values}.apply(localsum);
  }

 /* 11. Create a map using local sums to generate final GIDs.
   */

  // LINE 10: In the GUN paper

  panzer::GlobalOrdinal myOffset = -1;
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManager::buildGlobalUnknowns_GUN::line_10 prefix_sum",GUN_10);

    // do a prefix sum
    panzer::GlobalOrdinal scanResult = 0;
    Teuchos::scan<int, panzer::GlobalOrdinal> (*getComm(), Teuchos::REDUCE_SUM, static_cast<panzer::GlobalOrdinal> (localsum), Teuchos::outArg (scanResult));
    myOffset = scanResult - localsum;
  }

  // LINE 11 and 12: In the GUN paper, these steps are eliminated because
  // the non_overlap_mv is reused

 /* 12. Iterate through the non-overlapping MV and assign GIDs to
   *     the necessary points. (Assign a -1 elsewhere.)
   */

  // LINES 13-21: In the GUN paper

  {
    PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManager::buildGlobalUnknowns_GUN::line_13-21 gid_assignment",GUN13_21);
    int which_id=0;
    auto editnonoverlap = non_overlap_mv->getLocalViewHost(Tpetra::Access::ReadWrite);
    for(size_t i=0; i<non_overlap_mv->getLocalLength(); ++i){
      for(int j=0; j<numFields_; ++j){
        if(editnonoverlap(i,j)!=0){
          // editnonoverlap[j][i]=myOffset+which_id;
          int ndof = Teuchos::as<int>(editnonoverlap(i,j));
          editnonoverlap(i,j)=myOffset+which_id;
          which_id+=ndof;
        }
        else{
          editnonoverlap(i,j)=-1;
        }

      }
    }
  }

  // LINE 22: In the GUN paper. Were performed above, and the overlaped_mv is
  //          abused to handle input tagging.

 /* 13. Import data back to the overlap MV using REPLACE.
   */

  // LINE 23: In the GUN paper

  {
    PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManager::buildGlobalUnknowns_GUN::line_23 final_import",GUN23);

    // use exporter to save on communication setup costs
    overlap_mv.doImport(*non_overlap_mv,*exp,Tpetra::REPLACE);
  }

  //std::cout << Teuchos::describe(*non_overlap_mv,Teuchos::VERB_EXTREME)  << std::endl;

  // return non_overlap_mv;
  return std::make_pair(non_overlap_mv,tagged_non_overlap_mv);
}

///////////////////////////////////////////////////////////////////////////////
Teuchos::RCP<Tpetra::MultiVector<panzer::GlobalOrdinal,panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> >
DOFManager::buildTaggedMultiVector(const ElementBlockAccess & ownedAccess)
{
  // some typedefs
  typedef panzer::TpetraNodeType Node;
  typedef Tpetra::Map<panzer::LocalOrdinal, panzer::GlobalOrdinal, Node> Map;
  typedef Tpetra::MultiVector<panzer::GlobalOrdinal,panzer::LocalOrdinal,panzer::GlobalOrdinal,Node> MultiVector;

  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManager::buildTaggedMultiVector",BTMV);

  //We will iterate through all of the blocks, building a FieldAggPattern for
  //each of them.

  for (size_t b = 0; b < blockOrder_.size(); ++b) {
    std::vector<std::tuple< int, panzer::FieldType, RCP<const panzer::FieldPattern> > > faConstruct;
    //The ID is going to be the AID, and then everything will work.
    //The ID should not be the AID, it should be the ID it has in the ordering.

    for (size_t i = 0; i < fieldAIDOrder_.size(); ++i) {
      int looking = fieldAIDOrder_[i];

      //Check if in b's fp list
      std::vector<int>::const_iterator reu = std::find(blockToAssociatedFP_[b].begin(), blockToAssociatedFP_[b].end(), looking);
      if(!(reu==blockToAssociatedFP_[b].end())){
        faConstruct.push_back(std::make_tuple(i, fieldTypes_[fieldAIDOrder_[i]], fieldPatterns_[fieldAIDOrder_[i]]));
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

  RCP<const Map> overlapmap       = this->buildOverlapMapFromElements(ownedAccess);

  // LINE 22: In the GUN paper...the overlap_mv is reused for the tagged multivector.
  //          This is a bit of a practical abuse of the algorithm presented in the paper.

  Teuchos::RCP<MultiVector> overlap_mv;
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManager::buildTaggedMultiVector::allocate_tagged_multivector",ATMV);

    overlap_mv = Tpetra::createMultiVector<panzer::GlobalOrdinal>(overlapmap,(size_t)numFields_);
    overlap_mv->putScalar(0); // if tpetra is not initialized with zeros
  }

  /* 5.  Iterate through all local elements again, checking with the FP
   *     information. Mark up the overlap map accordingly.
   */

  {
    PANZER_FUNC_TIME_MONITOR_DIFF("panzer::DOFManager::buildTaggedMultiVector::fill_tagged_multivector",FTMV);

    // temporary working vector to fill each row in tagged array
    std::vector<int> working(overlap_mv->getNumVectors());
    auto edittwoview_host = overlap_mv->getLocalViewHost(Tpetra::Access::ReadWrite);
    for (size_t b = 0; b < blockOrder_.size(); ++b) {
      // there has to be a field pattern associated with the block
      if(fa_fps_[b]==Teuchos::null)
        continue;

      const std::vector<panzer::LocalOrdinal> & numFields= fa_fps_[b]->numFieldsPerId();
      const std::vector<panzer::LocalOrdinal> & fieldIds= fa_fps_[b]->fieldIds();
      const std::vector<panzer::LocalOrdinal> & myElements = connMngr_->getElementBlock(blockOrder_[b]);
      for (size_t l = 0; l < myElements.size(); ++l) {
        auto connSize = connMngr_->getConnectivitySize(myElements[l]);
        const auto * elmtConn = connMngr_->getConnectivity(myElements[l]);
        int offset=0;
        for (int c = 0; c < connSize; ++c) {
          size_t lid = overlapmap->getLocalElement(elmtConn[c]);
          for(std::size_t i=0;i<working.size();i++)
            working[i] = 0;
          for (int n = 0; n < numFields[c]; ++n) {
            int whichField = fieldIds[offset];
            //Row will be lid. column will be whichField.
            //Shove onto local ordering
            working[whichField]++;
            offset++;
          }
          for(std::size_t i=0;i<working.size();i++) {
            auto current = edittwoview_host(lid,i);
            edittwoview_host(lid,i) = (current > working[i]) ? current : working[i];
          }

        }
      }
    }
    // // verbose output for inspecting overlap_mv
    // for(int i=0;i<overlap_mv->getLocalLength(); i++) {
    //   for(int j=0;j<overlap_mv->getNumVectors() ; j++)
    //     std::cout << edittwoview[j][i] << " ";
    //   std::cout << std::endl;
    // }
  }

  return overlap_mv;
}

///////////////////////////////////////////////////////////////////////////////
int DOFManager::getFieldNum(const std::string & string) const
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

///////////////////////////////////////////////////////////////////////////////
void DOFManager::getFieldOrder(std::vector<std::string> & fieldOrder) const
{
  fieldOrder.resize(fieldStringOrder_.size());
  for (size_t i = 0; i < fieldStringOrder_.size(); ++i)
    fieldOrder[i]=fieldStringOrder_[i];
}

///////////////////////////////////////////////////////////////////////////////
bool DOFManager::fieldInBlock(const std::string & field, const std::string & block) const
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

///////////////////////////////////////////////////////////////////////////////
const std::vector<int> & DOFManager::getBlockFieldNumbers(const std::string & blockId) const
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

///////////////////////////////////////////////////////////////////////////////
const std::pair<std::vector<int>,std::vector<int> > &
DOFManager::getGIDFieldOffsets_closure(const std::string & blockId, int fieldNum, int subcellDim,int subcellId) const
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

///////////////////////////////////////////////////////////////////////////////
void DOFManager::ownedIndices(const std::vector<panzer::GlobalOrdinal> & indices,std::vector<bool> & isOwned) const
{
  //Resizes the isOwned array.
  if(indices.size()!=isOwned.size())
    isOwned.resize(indices.size(),false);
  typename std::vector<panzer::GlobalOrdinal>::const_iterator endOf = owned_.end();
  for (std::size_t i = 0; i < indices.size(); ++i) {
    isOwned[i] = ( std::find(owned_.begin(), owned_.end(), indices[i])!=endOf );
  }
}

///////////////////////////////////////////////////////////////////////////////
void DOFManager::setFieldOrder(const std::vector<std::string> & fieldOrder)
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
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"DOFManager::setFieldOrder: Invalid Field Ordering!");
  }
}

///////////////////////////////////////////////////////////////////////////////
bool DOFManager::validFieldOrder(const std::vector<std::string> & proposed_fieldOrder)
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

///////////////////////////////////////////////////////////////////////////////
const std::string & DOFManager::getFieldString(int num) const
{
  //A reverse lookup through fieldStringOrder_.
  if(num>=(int)fieldStringOrder_.size())
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "DOFManager::getFieldString: invalid number");
  return fieldStringOrder_[num];
}

///////////////////////////////////////////////////////////////////////////////
// Everything associated with orientation is not yet built, but this
// is the method as found in Panzer_DOFManager_impl.hpp. There are
// going to need to be some substantial changes to the code as it applies
// to this DOFManager, but the basic ideas and format should be similar.
//
void DOFManager::buildUnknownsOrientation()
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
    const std::vector<panzer::LocalOrdinal> & elmts = connMngr_->getElementBlock(blockName);
    for(std::size_t e=0;e<elmts.size();e++) {
       // this is the vector of orientations to fill: initialize it correctly
      std::vector<signed char> & eOrientation = orientation_[elmts[e]];

      // This resize seems to be the same as fieldPattern.numberIDs().
      // When computer edge orientations is called, that is the assert.
      // There should be no reason to make it anymore complicated.
      eOrientation.resize(fieldPattern.numberIds());
      for(std::size_t s=0;s<eOrientation.size();s++)
        eOrientation[s] = 1; // put in 1 by default

      // get geometry ids
      auto connSz = connMngr_->getConnectivitySize(elmts[e]);
      const panzer::GlobalOrdinal * connPtr = connMngr_->getConnectivity(elmts[e]);
      const std::vector<panzer::GlobalOrdinal> connectivity(connPtr,connPtr+connSz);

      orientation_helpers::computeCellEdgeOrientations(topEdgeIndices, connectivity, fieldPattern, eOrientation);

      // compute face orientations in 3D
      if(ga_fp_->getDimension()==3)
        orientation_helpers::computeCellFaceOrientations(topFaceIndices, connectivity, fieldPattern, eOrientation);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
void DOFManager::getElementOrientation(panzer::LocalOrdinal localElmtId,std::vector<double> & gidsOrientation) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(orientation_.size()==0,std::logic_error,
                              "DOFManager::getElementOrientations: Orientations were not constructed!");

   const std::vector<signed char> & local_o = orientation_[localElmtId];
   gidsOrientation.resize(local_o.size());
   for(std::size_t i=0;i<local_o.size();i++) {
      gidsOrientation[i] = double(local_o[i]);
   }
}

Teuchos::RCP<ConnManager> DOFManager::resetIndices()
{
   Teuchos::RCP<ConnManager> connMngr = connMngr_;

   connMngr_ = Teuchos::null;

   // wipe out FEI objects
   orientation_.clear(); // clean up previous work
   fa_fps_.clear();
   elementGIDs_.clear();
   owned_.clear();
   ghosted_.clear();
   elementBlockGIDCount_.clear();

   return connMngr;
}

///////////////////////////////////////////////////////////////////////////////
std::size_t DOFManager::blockIdToIndex(const std::string & blockId) const
{
  std::map<std::string,int>::const_iterator bitr = blockNameToID_.find(blockId);
  if(bitr==blockNameToID_.end())
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"DOFManager::fieldInBlock: invalid block name");
  return bitr->second;
}

///////////////////////////////////////////////////////////////////////////////
void DOFManager::printFieldInformation(std::ostream & os) const
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

///////////////////////////////////////////////////////////////////////////////
Teuchos::RCP<const Tpetra::Map<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> >
DOFManager::buildOverlapMapFromElements(const ElementBlockAccess & access) const
{
  PANZER_FUNC_TIME_MONITOR("panzer::DOFManager::builderOverlapMapFromElements");

  /*
   * 2.  Iterate through all local elements and create the overlapVector
   *     of concerned elements.
   */

  std::set<panzer::GlobalOrdinal> overlapset;
  for (size_t i = 0; i < blockOrder_.size(); ++i) {
    const std::vector<panzer::LocalOrdinal> & myElements = access.getElementBlock(blockOrder_[i]);
    for (size_t e = 0; e < myElements.size(); ++e) {
      auto connSize = connMngr_->getConnectivitySize(myElements[e]);
      const panzer::GlobalOrdinal * elmtConn = connMngr_->getConnectivity(myElements[e]);
      for (int k = 0; k < connSize; ++k) {
        overlapset.insert(elmtConn[k]);
      }
    }
  }

  Array<panzer::GlobalOrdinal> overlapVector;
  for (typename std::set<panzer::GlobalOrdinal>::const_iterator itr = overlapset.begin(); itr!=overlapset.end(); ++itr) {
    overlapVector.push_back(*itr);
  }

  /* 3.  Construct an overlap map from this structure.
   */
  return Tpetra::createNonContigMapWithNode<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType>(overlapVector,getComm());
}

///////////////////////////////////////////////////////////////////////////////
void DOFManager::
fillGIDsFromOverlappedMV(const ElementBlockAccess & access,
                         std::vector<std::vector< panzer::GlobalOrdinal > > & elementGIDs,
                         const Tpetra::Map<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> & overlapmap,
                         const Tpetra::MultiVector<panzer::GlobalOrdinal,panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> & const_overlap_mv) const
{
  using Teuchos::ArrayRCP;

  //To generate elementGIDs we need to go through all of the local elements.
  auto overlap_mv = const_cast<Tpetra::MultiVector<panzer::GlobalOrdinal,panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType>&>(const_overlap_mv);
  const auto twoview_host = overlap_mv.getLocalViewHost(Tpetra::Access::ReadOnly);

  //And for each of the things in fa_fp.fieldIds we go to that column. To the the row,
  //we move from globalID to localID in the map and use our local value for something.
  for (size_t b = 0; b < blockOrder_.size(); ++b) {
    const std::vector<panzer::LocalOrdinal> & myElements = access.getElementBlock(blockOrder_[b]);

    if(fa_fps_[b]==Teuchos::null) {
      // fill elements that are not used with empty vectors
      for (size_t l = 0; l < myElements.size(); ++l) {
        panzer::LocalOrdinal thisID=myElements[l];
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
      auto connSize = connMngr_->getConnectivitySize(myElements[l]);
      const panzer::GlobalOrdinal * elmtConn = connMngr_->getConnectivity(myElements[l]);
      std::vector<panzer::GlobalOrdinal> localOrdering;
      int offset=0;
      for (int c = 0; c < connSize; ++c) {
        size_t lid = overlapmap.getLocalElement(elmtConn[c]);
        std::vector<int> dofsPerField(numFields_,0);
        for (int n = 0; n < numFields[c]; ++n) {
          int whichField = fieldIds[offset];
          offset++;
          //Row will be lid. column will be whichField.
          //Shove onto local ordering
          localOrdering.push_back(twoview_host(lid,whichField)+dofsPerField[whichField]);

          dofsPerField[whichField]++;
        }
      }
      panzer::LocalOrdinal thisID=myElements[l];
      if(elementGIDs.size()<=(size_t)thisID){
        elementGIDs.resize(thisID+1);
      }
      elementGIDs[thisID]=localOrdering;
    }
  }
}

void DOFManager::buildLocalIdsFromOwnedAndGhostedElements()
{
  std::vector<std::vector<panzer::LocalOrdinal> > elementLIDs(elementGIDs_.size());

  std::vector<panzer::GlobalOrdinal> ownedAndGhosted;
  this->getOwnedAndGhostedIndices(ownedAndGhosted);

  // build global to local hash map (temporary and used only once)
  std::unordered_map<panzer::GlobalOrdinal,panzer::LocalOrdinal> hashMap;
  for(std::size_t i = 0; i < ownedAndGhosted.size(); ++i)
    hashMap[ownedAndGhosted[i]] = i;

  for (std::size_t i = 0; i < elementGIDs_.size(); ++i) {
    const std::vector<panzer::GlobalOrdinal>& gids = elementGIDs_[i];
    std::vector<panzer::LocalOrdinal>& lids = elementLIDs[i];
    lids.resize(gids.size());
    for (std::size_t g = 0; g < gids.size(); ++g)
      lids[g] = hashMap[gids[g]];
  }

  this->setLocalIds(elementLIDs);
}

/*
template <typename panzer::LocalOrdinal,typename panzer::GlobalOrdinal>
Teuchos::RCP<const Tpetra::Map<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> >
DOFManager<panzer::LocalOrdinal,panzer::GlobalOrdinal>::runLocalRCMReordering(const Teuchos::RCP<const Tpetra::Map<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> > & map)
{
  typedef panzer::TpetraNodeType Node;
  typedef Tpetra::Map<panzer::LocalOrdinal, panzer::GlobalOrdinal, Node> Map;
  typedef Tpetra::CrsGraph<panzer::LocalOrdinal, panzer::GlobalOrdinal, Node> Graph;

  Teuchos::RCP<Graph> graph = Teuchos::rcp(new Graph(map,0));

  // build Crs Graph from the mesh
  for (size_t b = 0; b < blockOrder_.size(); ++b) {
    if(fa_fps_[b]==Teuchos::null)
      continue;

    const std::vector<panzer::LocalOrdinal> & myElements = connMngr_->getElementBlock(blockOrder_[b]);
    for (size_t l = 0; l < myElements.size(); ++l) {
      panzer::LocalOrdinal connSize = connMngr_->getConnectivitySize(myElements[l]);
      const panzer::GlobalOrdinal * elmtConn = connMngr_->getConnectivity(myElements[l]);
      for (int c = 0; c < connSize; ++c) {
        panzer::LocalOrdinal lid = map->getLocalElement(elmtConn[c]);
        if(Teuchos::OrdinalTraits<panzer::LocalOrdinal>::invalid()!=lid)
          continue;

        graph->insertGlobalIndices(elmtConn[c],Teuchos::arrayView(elmtConn,connSize));
      }
    }
  }

  graph->fillComplete();

  std::vector<panzer::GlobalOrdinal> newOrder(map->getLocalNumElements());
  {
    // graph is constructed, now run RCM using zoltan2
    typedef Zoltan2::XpetraCrsGraphInput<Graph> SparseGraphAdapter;

    Teuchos::ParameterList params;
    params.set("order_method", "rcm");
    SparseGraphAdapter adapter(graph);

    Zoltan2::OrderingProblem<SparseGraphAdapter> problem(&adapter, &params,MPI_COMM_SELF);
    problem.solve();

    // build a new global ording array using permutation
    Zoltan2::OrderingSolution<panzer::GlobalOrdinal,panzer::LocalOrdinal> * soln = problem.getSolution();

    size_t dummy;
    size_t checkLength = soln->getPermutationSize();
    panzer::LocalOrdinal * checkPerm = soln->getPermutation(&dummy);

    Teuchos::ArrayView<const panzer::GlobalOrdinal > oldOrder = map->getLocalElementList();
    TEUCHOS_ASSERT(checkLength==oldOrder.size());
    TEUCHOS_ASSERT(checkLength==newOrder.size());

    for(size_t i=0;i<checkLength;i++)
      newOrder[checkPerm[i]] = oldOrder[i];
  }

  return Tpetra::createNonContigMap<panzer::LocalOrdinal,panzer::GlobalOrdinal>(newOrder,communicator_);
}
*/

} /*panzer*/

#endif

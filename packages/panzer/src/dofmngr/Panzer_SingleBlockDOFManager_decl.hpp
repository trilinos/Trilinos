#ifndef __Panzer_SingleBlockDOFManager_decl_hpp__
#define __Panzer_SingleBlockDOFManager_decl_hpp__

#include <map>

#include "mpi.h"

#include "Panzer_config.hpp"
#include "Panzer_FieldPattern.hpp"
#include "Panzer_FieldAggPattern.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"

#include "Teuchos_RCP.hpp"

namespace panzer {


template <typename LocalOrdinalT, typename GlobalOrdinalT>
class SingleBlockDOFManager : public UniqueGlobalIndexer<LocalOrdinalT, GlobalOrdinalT> {
public:
  typedef GlobalOrdinalT GO;
  typedef LocalOrdinalT LO;

  virtual ~SingleBlockDOFManager() {}

  SingleBlockDOFManager();

  //TO-WRITE: Cosntructor that takes a connection manager.
  
  //Adds a Connection Manager that will be associated with this SingleBlockDOFManager.
  void setConnManager(const Teuchos::RCP<ConnManager<LO,GO> > & connMngr, MPI_Comm mpiComm);

  Teuchos::RCP<const ConnManager<LO,GO> > getConnManager() const
  { return connMngr_; }
  //Adds a field to be used in creating the Global Numbering
  int addField(const std::string & str, const Teuchos::RCP<const FieldPattern> & pattern);


  //TO-WRITE: A method that differentiates between elementblocks
  //TO-WRITE: Field Ordering Method...I'm not sure why this is that important.

  //Returns the fieldpattern of the given name
  //This could also be done using the number you'd get from getFieldNum which
  //isn't yet included.
  Teuchos::RCP<const FieldPattern> getFieldPattern(const std::string & name);

  void getOwnedIndices(std::vector<GlobalOrdinalT> & indices) const;

  void getOwnedAndSharedIndices(std::vector<GlobalOrdinalT> & indices) const;
  //gets the number of fields
  int getNumFields() const;

  //gets the field pattern so you can find a particular
  //field in the GIDs aray.
  const std::vector<int> & getGIDFieldOffsets(const std::string & blockID, int fieldNum) const;


  //get associated GIDs for a given local element
  void getElementGIDs(LO localElementID, std::vector<GO> & gids, const std::string & blockIdHint="") const;

  //builds the global unknowns array
  void buildGlobalUnknowns();
  
  int getFieldNum(const std::string & string) const;

  Teuchos::RCP<Teuchos::Comm<int> > getComm() const
  { return communicator_; }

  void getElementBlockIds(std::vector<std::string> & elementBlockIds) const
  { connMngr_->getElementBlockIds(elementBlockIds); }
  
  //Because we only have one element block, if the field is present, and
  //the block is valid. It works
  bool fieldInBlock(const std::string & field, const std::string & block) const;

  //Because we only have one element block, we are guarenteed 
  //that all fields are going to be in the one block we have
  const std::vector<int> & getBlockFieldNumbers(const std::string & blockId) const;

//************************************************************************

  //TODO:this.
  const std::pair<std::vector<int>,std::vector<int> > & 
  getGIDFieldOffsets_closure(const std::string & blockId, int fieldNum,
                              int subcellDim,int subcellId) const
  { return fa_fp->localOffsets_closure(fieldNum,subcellDim,subcellId); }


  void getElementOrientation(LocalOrdinalT localElmtId,std::vector<double> & gidsOrientation) const
  { TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"SingleBlockDOFManager::getElementOrientation not implemented yet!"); }

  const std::vector<LocalOrdinalT> & getElementBlock(const std::string & blockId) const
  { return connMngr_->getElementBlock(blockId); }

  void ownedIndices(const std::vector<GlobalOrdinalT> & indices,std::vector<bool> & isOwned) const;

  void setFieldOrder(const std::vector<std::string> & fieldOrder );

  void getFieldOrder(std::vector<std::string> & fieldOrder)
  { fieldOrder=fieldOrder_;}

  bool validFieldOrder(const std::vector<std::string> & proposed_fieldOrder);

  virtual void buildUnknownsOrientation();

  bool getOrientationsRequired() const
  { return requireOrientations_; }

  void setOrientationsRequired(bool ro)
  { requireOrientations_ = ro; }

protected:
  Teuchos::RCP<ConnManager<LO,GO> > connMngr_;

  std::map<std::string,int > fieldStringToIndex_;
  std::map<std::string,int > fieldStringToID_;
  std::vector<Teuchos::RCP<const FieldPattern> > fieldPatterns_;

  std::vector<std::vector< GO > > elementGIDs_;

  std::vector<std::string> fieldOrder_;

  //Teuchos::RCP<const FieldPattern> geomPattern_;
  Teuchos::RCP<Teuchos::Comm<int> > communicator_;

  Teuchos::RCP<panzer::FieldAggPattern> fa_fp;
  Teuchos::RCP<panzer::GeometricAggFieldPattern> ga_fp;

  std::vector<GO> owned_;
  std::vector<GO> owned_and_ghosted_;

  int num_dof_;
  int numFields_;
  
  bool buildConnectivityRun_;

  bool requireOrientations_;
  std::vector<std::vector<char> > orientation_;

};



}

#endif

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_sstream.hpp>

#include <limits>
#include <cmath>

#include <fei_MatrixGraph_Impl2.hpp>

#include <fei_utils.hpp>
#include "fei_TemplateUtils.hpp"

#include <fei_Pattern.hpp>
#include <fei_LogManager.hpp>
#include <fei_TemplateUtils.hpp>
#include <fei_impl_utils.hpp>
#include <snl_fei_Utils.hpp>
#include <fei_FieldMask.hpp>
#include <fei_Record.hpp>
#include <snl_fei_RecordCollection.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_ParameterSet.hpp>
#include <fei_ostream_ops.hpp>
#include <fei_Reducer.hpp>
#include <fei_GraphReducer.hpp>
#include <fei_ConnectivityBlock.hpp>
#include <snl_fei_BlkSizeMsgHandler.hpp>
#include <fei_Graph_Impl.hpp>
#include <snl_fei_Constraint.hpp>

#include <fei_EqnBuffer.hpp>
#include <fei_EqnCommMgr.hpp>
#include <SNL_FEI_Structure.hpp>

#undef fei_file
#define fei_file "fei_MatrixGraph.cpp"

#include <fei_ErrMacros.hpp>

namespace snl_fei {
static unsigned getFieldSize(int fieldID,
                      fei::VectorSpace* space1,
                      fei::VectorSpace* space2)
{
  unsigned fieldsize = 0;
  bool foundfield = false;
  if (space1 != NULL) {
    try {
      fieldsize = space1->getFieldSize(fieldID);
      foundfield = true;
    }
    catch (...) {
      foundfield = false;
    }
  }

  if (!foundfield && space2 != NULL) {
    try {
      fieldsize = space2->getFieldSize(fieldID);
    }
    catch (std::runtime_error& exc) {
      std::string msg("snl_fei::getFieldSize: ");
      msg += exc.what();
      throw std::runtime_error(msg);
    }
  }

  return(fieldsize);
}

}//namespace snl_fei

//------------------------------------------------------------------------------
fei::SharedPtr<fei::MatrixGraph>
fei::MatrixGraph_Impl2::Factory::createMatrixGraph(fei::SharedPtr<fei::VectorSpace> rowSpace,
                                             fei::SharedPtr<fei::VectorSpace> columnSpace,
                                             const char* name)
{
  fei::SharedPtr<fei::MatrixGraph> mgptr(new fei::MatrixGraph_Impl2(rowSpace,
                                                              columnSpace,
                                                              name));

  return(mgptr);
}

//=====Constructor==============================================================
fei::MatrixGraph_Impl2::MatrixGraph_Impl2(fei::SharedPtr<fei::VectorSpace> rowSpace,
                              fei::SharedPtr<fei::VectorSpace> colSpace,
                              const char* name)
 : comm_(rowSpace->getCommunicator()),
   rowSpace_(rowSpace),
   colSpace_(colSpace),
   haveRowSpace_(false),
   haveColSpace_(false),
   symmetric_(false),
   remotelyOwnedGraphRows_(NULL),
   simpleProblem_(false),
   blockEntryGraph_(false),
   patterns_(),
   connectivityBlocks_(),
   arbitraryBlockCounter_(1),
   sparseBlocks_(),
   lagrangeConstraints_(),
   penaltyConstraints_(),
   slaveConstraints_(),
   ptEqualBlk_(false),
   newSlaveData_(false),
   localNumSlaves_(0),
   globalNumSlaves_(0),
   D_(NULL),
   g_(),
   g_nonzero_(false),
   reducer_(),
   name_(),
   dbgprefix_("MatGrph: "),
   tmpIntArray1_(),
   tmpIntArray2_()
{
  localProc_ = fei::localProc(comm_);
  numProcs_ = fei::numProcs(comm_);

  if (rowSpace_.get() == NULL) {
    voidERReturn;
  }
  else haveRowSpace_ = true;

  if (colSpace_.get() != NULL) haveColSpace_ = true;
  else colSpace_ = rowSpace_;

  setName(name);
}

//------------------------------------------------------------------------------
fei::MatrixGraph_Impl2::~MatrixGraph_Impl2()
{
  fei::destroyValues(patterns_);
  patterns_.clear();

  fei::destroyValues(connectivityBlocks_);
  connectivityBlocks_.clear();

  int i, len = sparseBlocks_.size();
  for(i=0; i<len; ++i) {
    delete sparseBlocks_[i];
  }

  fei::destroyValues(lagrangeConstraints_);
  lagrangeConstraints_.clear();

  fei::destroyValues(penaltyConstraints_);
  penaltyConstraints_.clear();

  fei::destroyValues(slaveConstraints_);
  slaveConstraints_.clear();
}

//------------------------------------------------------------------------------
void fei::MatrixGraph_Impl2::setParameters(const fei::ParameterSet& params)
{
  const fei::Param* param = 0;
  fei::Param::ParamType ptype = fei::Param::BAD_TYPE;

  param = params.get("FEI_OUTPUT_LEVEL");
  ptype = param != NULL ? param->getType() : fei::Param::BAD_TYPE;
  if (ptype == fei::Param::STRING) {
    fei::LogManager& log_manager = fei::LogManager::getLogManager();
    log_manager.setOutputLevel(param->getStringValue().c_str());
    setOutputLevel(fei::utils::string_to_output_level(param->getStringValue()));
  }

  param = params.get("FEI_LOG_EQN");
  ptype =  param != NULL ? param->getType() : fei::Param::BAD_TYPE;  
  if (ptype == fei::Param::INT) {
    addLogEqn(param->getIntValue());
  }

  param = params.get("FEI_LOG_ID");
  ptype =  param != NULL ? param->getType() : fei::Param::BAD_TYPE;  
  if (ptype == fei::Param::INT) {
    addLogID(param->getIntValue());
  }

  param = params.get("name");
  ptype = param != NULL ? param->getType() : fei::Param::BAD_TYPE;
  if (ptype == fei::Param::STRING) {
    setName(param->getStringValue().c_str());
  }

  param = params.get("BLOCK_GRAPH");
  ptype = param != NULL ? param->getType() : fei::Param::BAD_TYPE;
  if (ptype != fei::Param::BAD_TYPE) {
    blockEntryGraph_ = true;
  }

  param = params.get("BLOCK_MATRIX");
  ptype = param != NULL ? param->getType() : fei::Param::BAD_TYPE;
  if (ptype != fei::Param::BAD_TYPE) {
    if (ptype == fei::Param::BOOL) {
      blockEntryGraph_ = param->getBoolValue();
    }
    else {
      blockEntryGraph_ = true;
    }
  }
}

//----------------------------------------------------------------------------
void fei::MatrixGraph_Impl2::setRowSpace(fei::SharedPtr<fei::VectorSpace> rowSpace)
{
  rowSpace_ = rowSpace;
  haveRowSpace_ = true;
  if (!haveColSpace_) symmetric_ = true;
}

//----------------------------------------------------------------------------
void fei::MatrixGraph_Impl2::setColumnSpace(fei::SharedPtr<fei::VectorSpace> colSpace)
{
  colSpace_ = colSpace;
  haveColSpace_ = true;
  if (colSpace_ == rowSpace_) symmetric_ = true;
  else symmetric_ = false;
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::addPattern(fei::Pattern* pattern)
{
  std::map<int,fei::Pattern*>::iterator
    p_iter = patterns_.begin(), p_iter_end = patterns_.end();

  bool matches_existing_pattern = false;
  int pattern_id = -1;
  while(p_iter != p_iter_end && !matches_existing_pattern) {
    if (*pattern == *(p_iter->second)) {
      matches_existing_pattern = true;
      pattern_id = p_iter->first;
    }
    else ++p_iter;
  }

  if (!matches_existing_pattern) {
    pattern_id = patterns_.size();
    patterns_.insert(std::make_pair(pattern_id, pattern));
  }
  else {
    delete pattern;
  }

  return pattern_id;
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::definePattern(int numIDs,
                                    int idType)
{
  snl_fei::RecordCollection* rec_coll = NULL;
  rowSpace_->getRecordCollection(idType, rec_coll);

  fei::Pattern* pattern = new fei::Pattern(numIDs, idType, rec_coll);
  return addPattern(pattern);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::definePattern(int numIDs,
                                      int idType,
                                      int fieldID)
{
  unsigned fieldsize = 0;
  try {
    fieldsize = snl_fei::getFieldSize(fieldID, rowSpace_.get(), colSpace_.get());
  }
  catch (std::runtime_error& exc) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei::MatrixGraph_Impl2::definePattern caught error: "<<exc.what();
    throw std::runtime_error(osstr.str());
  }

  snl_fei::RecordCollection* rec_coll = NULL;
  rowSpace_->getRecordCollection(idType, rec_coll);

  fei::Pattern* pattern =
    new fei::Pattern(numIDs, idType, rec_coll, fieldID, fieldsize);
  return addPattern(pattern);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::definePattern(int numIDs,
                                     int idType,
                                     const int* numFieldsPerID,
                                     const int* fieldIDs)
{ 
  std::vector<int> fieldSizes;
  try {
    int offset = 0;
    for(int i=0; i<numIDs; ++i) {
      for(int j=0; j<numFieldsPerID[i]; ++j) {
        fieldSizes.push_back(snl_fei::getFieldSize(fieldIDs[offset++],
                                              rowSpace_.get(), colSpace_.get()));
      }
    }
  }
  catch (std::runtime_error& exc) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei::MatrixGraph_Impl2::definePattern caught error: "<<exc.what();
    throw std::runtime_error(osstr.str());
  }

  snl_fei::RecordCollection* rec_coll = NULL;
  rowSpace_->getRecordCollection(idType, rec_coll);

  fei::Pattern* pattern = new fei::Pattern(numIDs, idType, rec_coll,
                        numFieldsPerID, fieldIDs, &(fieldSizes[0]));
  return addPattern(pattern);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::definePattern(int numIDs,
                                      const int* idTypes,
                                      const int* numFieldsPerID,
                                      const int* fieldIDs)
{
  std::vector<int> fieldSizes;
  try {
    int offset = 0;
    for(int i=0; i<numIDs; ++i) {
      for(int j=0; j<numFieldsPerID[i]; ++j) {
        fieldSizes.push_back(snl_fei::getFieldSize(fieldIDs[offset++],
                                              rowSpace_.get(), colSpace_.get()));
      }
    }
  }
  catch (std::runtime_error& exc) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei::MatrixGraph_Impl2::definePattern caught error: "<<exc.what();
    throw std::runtime_error(osstr.str());
  }

  std::vector<snl_fei::RecordCollection*> recordCollections(numIDs);
  for(int i=0; i<numIDs; ++i) {
    rowSpace_->getRecordCollection(idTypes[i], recordCollections[i]);
  }

  fei::Pattern* pattern = new fei::Pattern(numIDs, idTypes, &recordCollections[0],
                        numFieldsPerID, fieldIDs, &(fieldSizes[0]));
  return addPattern(pattern);
}

//------------------------------------------------------------------------------
fei::Pattern* fei::MatrixGraph_Impl2::getPattern(int patternID)
{
  std::map<int,fei::Pattern*>::iterator
    p_iter = patterns_.find(patternID);

  if (p_iter == patterns_.end()) {
    return NULL;
  }

  fei::Pattern* pattern = (*p_iter).second;
  return pattern;
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::initConnectivityBlock(int blockID,
                                            int numConnectivityLists,
                                            int patternID,
                                            bool diagonal)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    (*output_stream_) <<dbgprefix_<< "initConnectivityBlock symm, blkID="
     << blockID<<", num="<<numConnectivityLists<<", patternID="
     <<patternID<<FEI_ENDL;
  }

  if (blockID < 0) {
    fei::console_out() << "fei::MatrixGraph_Impl2::initConnectivityBlock: blockID ("
         << blockID << ") must be non-negative." << FEI_ENDL;
    ERReturn(-1);
  }

  std::map<int,fei::Pattern*>::iterator
    p_iter = patterns_.find(patternID);

  if (p_iter == patterns_.end()) {
    ERReturn(-1);
  }

  fei::Pattern* pattern = (*p_iter).second;

  if (getConnectivityBlock(blockID) != NULL) ERReturn(-1);

  fei::ConnectivityBlock* cblock =
    new fei::ConnectivityBlock(blockID, pattern, numConnectivityLists);

  connectivityBlocks_.insert(std::pair<int,fei::ConnectivityBlock*>(blockID, cblock));
  if (diagonal) cblock->setIsDiagonal(diagonal);

  return(0);
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::initConnectivityBlock(int numConnectivityLists,
                                                  int patternID,
                                                  bool diagonal)
{
  int blockID = connectivityBlocks_.size();
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    (*output_stream_) <<dbgprefix_<< "initConnectivityBlock symm, blkID="
     << blockID<<", num="<<numConnectivityLists<<", patternID="
     <<patternID<<FEI_ENDL;
  }

  std::map<int,fei::Pattern*>::iterator
    p_iter = patterns_.find(patternID);

  if (p_iter == patterns_.end()) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei::MatrixGraph_Impl2::initConnectivityBlock: ERROR, patternID "
      << patternID << " not found.";
    throw std::runtime_error(osstr.str());
  }

  fei::Pattern* pattern = (*p_iter).second;

  fei::ConnectivityBlock* cblock =
    new fei::ConnectivityBlock(blockID, pattern, numConnectivityLists);

  connectivityBlocks_.insert(std::pair<int,fei::ConnectivityBlock*>(blockID, cblock));
  if (diagonal) cblock->setIsDiagonal(diagonal);

  return(blockID);
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::initConnectivityBlock(int blockID,
                                            int numConnectivityLists,
                                            int rowPatternID,
                                            int colPatternID)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    (*output_stream_) <<dbgprefix_<< "initConnectivityBlock nonsymm, blkID="
     << blockID<<", num="<<numConnectivityLists<<", row-patternID="
     <<rowPatternID<<", col-patternID="<<colPatternID<<FEI_ENDL;
  }

  if (blockID < 0) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei::MatrixGraph_Impl2::initConnectivityBlock: blockID ("
          << blockID << ") must be non-negative." << FEI_ENDL;
    throw std::runtime_error(osstr.str());
  }

  std::map<int,fei::Pattern*>::iterator
    rp_iter = patterns_.find(rowPatternID);
  if (rp_iter == patterns_.end()) ERReturn(-1);

  fei::Pattern* rowpattern = (*rp_iter).second;

  std::map<int,fei::Pattern*>::iterator
    cp_iter = patterns_.find(colPatternID);
  if (cp_iter == patterns_.end()) ERReturn(-1);

  fei::Pattern* colpattern = (*cp_iter).second;

  if (getConnectivityBlock(blockID) != NULL) ERReturn(-1);

  fei::ConnectivityBlock* cblock =
    new fei::ConnectivityBlock(blockID, rowpattern,
                                   colpattern, numConnectivityLists);

  connectivityBlocks_.insert(std::pair<int,fei::ConnectivityBlock*>(blockID, cblock));

  return(0);
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::initConnectivity(int blockID,
                                       int connectivityID,
                                       const int* connectedIdentifiers)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    (*output_stream_) <<dbgprefix_<< "initConnectivity blkID="
     <<blockID<<", connID="<<connectivityID<<FEI_ENDL;
  }

  //for each connected identifier, get a pointer to its record from the
  //solution-space.
  //store those pointers in the appropriate connectivity-table, mapped to
  //the user-provided connectivityID.

  fei::ConnectivityBlock* connblk = getConnectivityBlock(blockID);
  if (connblk == NULL) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei::MatrixGraph_Impl2::initConnectivity ERROR, blockID " << blockID
          << " doesn't correspond to an existing ConnectivityBlock.";
    throw std::runtime_error(osstr.str());
  }

  fei::Pattern* pattern = connblk->getRowPattern();
  int numIDs = pattern->getNumIDs();

  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    for(int ii=0; ii<numIDs; ++ii) {
      if (isLogID(connectedIdentifiers[ii])) {
        FEI_OSTREAM& os = *output_stream_;
        os << dbgprefix_<<"initConnectivity blockID " << blockID << ", connID "
           << connectivityID << ": ";
        for(int jj=0; jj<numIDs; ++jj) {
          os << connectedIdentifiers[jj] << " ";
        }
        os << FEI_ENDL;
      }
    }
  }

  if (rowSpace_.get() == NULL) ERReturn(-1);

  std::map<int,int>& connectivityIDs = connblk->getConnectivityIDs();

  int idOffset = -1;
  std::map<int,int>::iterator
   iter = connectivityIDs.find(connectivityID);
  if (iter == connectivityIDs.end()) {
    idOffset = connectivityIDs.size();
    connectivityIDs.insert(iter, std::make_pair(connectivityID,idOffset));
  }
  else {
    idOffset = iter->second;
  }

  idOffset *= pattern->getNumIDs();

  int* rlist = &(connblk->getRowConnectivities()[idOffset]);

  CHK_ERR( getConnectivityRecords(pattern, rowSpace_.get(),
                                  connectedIdentifiers, rlist) );

  for(int i=0; i<numIDs; ++i) {
    if (pattern->getNumFieldsPerID()[i] > 0) {
      pattern->getRecordCollections()[i]->getRecordWithLocalID(rlist[i])->isInLocalSubdomain_ = true;
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::initConnectivity(int idType,
                                           int numRows,
                                           const int* rowIDs,
                                           const int* rowOffsets,
                                           const int* packedColumnIDs)
{
  fei::ConnectivityBlock* block = new fei::ConnectivityBlock(numRows,
                                                   rowIDs, rowOffsets);

  int max_row_len = 0;
  for(int i=0; i<numRows; ++i) {
    int row_len = rowOffsets[i+1]-rowOffsets[i];
    if (row_len > max_row_len) max_row_len = row_len;
  }

  int patternID = definePattern(max_row_len, idType);
  fei::Pattern* pattern = getPattern(patternID);
  block->setRowPattern(pattern);

  sparseBlocks_.push_back(block);

  int* row_records = &(block->getRowConnectivities()[0]);

  CHK_ERR( getConnectivityRecords(rowSpace_.get(),
                                  idType, numRows, rowIDs, row_records) );

  fei::VectorSpace* colSpace = rowSpace_.get();
  if (colSpace_.get() != NULL) {
    colSpace = colSpace_.get();
  }

  int* col_records = &(block->getColConnectivities()[0]);

  CHK_ERR( getConnectivityRecords(colSpace, idType, rowOffsets[numRows],
                                  packedColumnIDs, col_records));

  return(0);
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::initConnectivity(int idType,
                                       int numRows,
                                       const int* rowIDs,
                                       const int* rowLengths,
                                       const int*const* columnIDs)
{
  fei::ConnectivityBlock* block = new fei::ConnectivityBlock(numRows,
                                                   rowIDs, rowLengths, true);

  int max_row_len = 0;
  for(int i=0; i<numRows; ++i) {
    int row_len = rowLengths[i];
    if (row_len > max_row_len) max_row_len = row_len;
  }

  int patternID = definePattern(max_row_len, idType);
  fei::Pattern* pattern = getPattern(patternID);
  block->setRowPattern(pattern);

  sparseBlocks_.push_back(block);

  int* row_records = &(block->getRowConnectivities()[0]);

  CHK_ERR( getConnectivityRecords(rowSpace_.get(),
                                  idType, numRows, rowIDs, row_records) );

  fei::VectorSpace* colSpace = rowSpace_.get();
  if (colSpace_.get() != NULL) {
    colSpace = colSpace_.get();
  }

  int* col_records = &(block->getColConnectivities()[0]);

  int offset = 0;
  for(int i=0; i<numRows; ++i) {
    CHK_ERR( getConnectivityRecords(colSpace, idType, rowLengths[i],
                                    columnIDs[i], &(col_records[offset])));
    offset += rowLengths[i];
  }

  return(0);
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::getConnectivityRecords(fei::VectorSpace* vecSpace,
                                                 int idType,
                                                 int numIDs,
                                                 const int* IDs,
                                                 int* records)
{
  snl_fei::RecordCollection* collection = NULL;
  CHK_ERR( vecSpace->getRecordCollection(idType, collection) );

  for(int i=0; i<numIDs; ++i) {
    records[i] = collection->getLocalID(IDs[i]);
    if (records[i] == -1) {
      CHK_ERR( vecSpace->addDOFs(idType, 1, &IDs[i], &records[i]) );
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::getConnectivityRecords(fei::VectorSpace* vecSpace,
                                                 int idType,
                                                 int fieldID,
                                                 int numIDs,
                                                 const int* IDs,
                                                 int* records)
{
  snl_fei::RecordCollection* collection = NULL;
  CHK_ERR( vecSpace->getRecordCollection(idType, collection) );

  for(int i=0; i<numIDs; ++i) {
    records[i] = collection->getLocalID(IDs[i]);
    if (records[i] == -1) {
      CHK_ERR( vecSpace->addDOFs(fieldID, idType, 1, &IDs[i], &records[i]));
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::getConnectivityRecords(fei::Pattern* pattern,
                                             fei::VectorSpace* vecSpace,
                                             const int* connectedIdentifiers,
                                             int* recordList)
{
  fei::Pattern::PatternType pType = pattern->getPatternType();
  int i, numIDs = pattern->getNumIDs();
  const int* numFieldsPerID = pattern->getNumFieldsPerID();
  const int* fieldIDs = pattern->getFieldIDs();

  if (pType == fei::Pattern::GENERAL) {
    const int* idTypes = pattern->getIDTypes();

    int fieldOffset = 0;
    for(i=0; i<numIDs; ++i) {
      int id = connectedIdentifiers[i];

      for(int nf=0; nf<numFieldsPerID[i]; ++nf) {
        CHK_ERR( vecSpace->addDOFs(fieldIDs[fieldOffset++],
                                               idTypes[i], 1, &id,
                                               &recordList[i]));
      }
    }
  }
  else if (pType == fei::Pattern::SINGLE_IDTYPE ||
           pType == fei::Pattern::SIMPLE ||
           pType == fei::Pattern::NO_FIELD) {

    int idType = pattern->getIDTypes()[0];

    switch(pType) {
    case fei::Pattern::SINGLE_IDTYPE:
      {
        int fieldOffset = 0;
        for(i=0; i<numIDs; ++i) {
          int id = connectedIdentifiers[i];
          for(int nf=0; nf<numFieldsPerID[i]; ++nf) {
            CHK_ERR( vecSpace->addDOFs(fieldIDs[fieldOffset++],
                                                   idType, 1, &id,
                                                   &(recordList[i])));
          }
        }
      }
      break;
    case fei::Pattern::SIMPLE:
      {
        CHK_ERR( vecSpace->addDOFs(fieldIDs[0], idType,
                                   numIDs, connectedIdentifiers,
                                   recordList) );
      }
      break;
    case fei::Pattern::NO_FIELD:
      CHK_ERR( vecSpace->addDOFs(idType, numIDs,
                                 connectedIdentifiers, recordList));
      break;
    case fei::Pattern::GENERAL:
      //Include a stub for this case to avoid compiler warning...
      std::abort(); break;
    }
  }
  else {
    ERReturn(-1);
  }

  return(0);
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::initConnectivity(int blockID,
                                       int connectivityID,
                                       const int* rowConnectedIdentifiers,
                                       const int* colConnectedIdentifiers)
{
  //for each connected identifier, get a pointer to its record from the
  //solution-space.
  //store those pointers in the appropriate connectivity-table, mapped to
  //the user-provided connectivityID.

  fei::ConnectivityBlock* connblk = getConnectivityBlock(blockID);
  if (connblk == NULL) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei::MatrixGraph_Impl2::initConnectivity ERROR, blockID " << blockID
          << " doesn't correspond to an existing ConnectivityBlock.";
    throw std::runtime_error(osstr.str());
  }

  fei::Pattern* pattern = connblk->getRowPattern();
  int numIDs = pattern->getNumIDs();
  fei::Pattern* colPattern = connblk->getColPattern();
  int numColIDs = colPattern->getNumIDs();
  if (rowSpace_.get() == NULL) {
    ERReturn(-1);
  }
  if (colSpace_.get() == NULL) {
    ERReturn(-1);
  }

  std::map<int,int>& connectivityIDs = connblk->getConnectivityIDs();

  int i, idOffset = -1;
  std::map<int,int>::iterator
    iter = connectivityIDs.find(connectivityID);
  if (iter == connectivityIDs.end()) {
    idOffset = connectivityIDs.size();
    connectivityIDs.insert(iter, std::make_pair(connectivityID,idOffset)); 
  }
  else {
    idOffset = iter->second;
  }
  int* row_rlist =
    &(connblk->getRowConnectivities()[idOffset*numIDs]);
  int* col_rlist =
    &(connblk->getColConnectivities()[idOffset*numColIDs]);

  CHK_ERR( getConnectivityRecords(pattern, rowSpace_.get(),
                                  rowConnectedIdentifiers, row_rlist) );

  for(i=0; i<numIDs; ++i)
    pattern->getRecordCollections()[i]->getRecordWithLocalID(row_rlist[i])->isInLocalSubdomain_ = true;
  CHK_ERR( getConnectivityRecords(colPattern, colSpace_.get(),
                                  colConnectedIdentifiers, col_rlist) );

  for(i=0; i<numColIDs; ++i)
    colPattern->getRecordCollections()[i]->getRecordWithLocalID(col_rlist[i])->isInLocalSubdomain_ = true;

  return(0);
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::initConnectivity(int patternID,
                                       const int* connectedIdentifiers)
{
  std::map<int,fei::Pattern*>::iterator
    p_iter = patterns_.find(patternID);
  if (p_iter == patterns_.end()) ERReturn(-1);

  fei::Pattern* pattern = p_iter->second;

  int blockID = -arbitraryBlockCounter_++;
  fei::ConnectivityBlock* cblock = new fei::ConnectivityBlock(blockID, pattern, 1);

  connectivityBlocks_.insert(std::pair<int,fei::ConnectivityBlock*>(blockID, cblock));

  CHK_ERR( initConnectivity(blockID, 0, connectedIdentifiers) );
  return(0);
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::initConnectivity(int rowPatternID,
                                       const int* rowConnectedIdentifiers,
                                       int colPatternID,
                                       const int* colConnectedIdentifiers)
{
  std::map<int,fei::Pattern*>::iterator
    rp_iter = patterns_.find(rowPatternID);
  if (rp_iter == patterns_.end()) ERReturn(-1);

  fei::Pattern* rowpattern = (*rp_iter).second;

  std::map<int,fei::Pattern*>::iterator
    cp_iter = patterns_.find(colPatternID);
  if (cp_iter == patterns_.end()) ERReturn(-1);

  fei::Pattern* colpattern = (*cp_iter).second;

  int blockID = -arbitraryBlockCounter_++;
  fei::ConnectivityBlock* cblock = new fei::ConnectivityBlock(blockID, rowpattern,
                                                    colpattern, 1);

  connectivityBlocks_.insert(std::pair<int,fei::ConnectivityBlock*>(blockID, cblock));

  CHK_ERR( initConnectivity(blockID, 0,
                            rowConnectedIdentifiers,
                            colConnectedIdentifiers) );
  return(0);
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::getPatternNumIndices(int patternID,
                                           int& numIndices)
{
  std::map<int,fei::Pattern*>::iterator
    p_iter = patterns_.find(patternID);
  if (p_iter == patterns_.end()) ERReturn(-1);
  fei::Pattern* pattern = (*p_iter).second;

  numIndices = pattern->getNumIndices();

  return(0);
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::getPatternIndices(int patternID,
                                        const int* IDs,
                                        std::vector<int>& indices)
{
  std::map<int,fei::Pattern*>::iterator
    p_iter = patterns_.find(patternID);
  if (p_iter == patterns_.end()) ERReturn(-1);

  fei::Pattern* pattern = (*p_iter).second;

  indices.resize(pattern->getNumIndices());

  int numIDs                = pattern->getNumIDs();
  const int* idTypes        = pattern->getIDTypes();
  const int* numFieldsPerID = pattern->getNumFieldsPerID();
  const int* fieldIDs       = pattern->getFieldIDs();

  int offset = 0, fieldOffset = 0;
  for(int i=0; i<numIDs; ++i) {
    for(int j=0; j<numFieldsPerID[i]; ++j) {
      CHK_ERR( rowSpace_->getGlobalIndices(1, &(IDs[i]), idTypes[i],
                                                  fieldIDs[fieldOffset],
                                                  &(indices[offset])) );
      offset += snl_fei::getFieldSize(fieldIDs[fieldOffset++],
                                        rowSpace_.get(), colSpace_.get());
    }
  }

  return(0);
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::initConnectivity(int idType,
                                       int fieldID,
                                       int numRows,
                                       const int* rowIDs,
                                       const int* rowOffsets,
                                       const int* packedColumnIDs)
{
  fei::ConnectivityBlock* block = new fei::ConnectivityBlock(fieldID, numRows,
                                                   rowIDs, rowOffsets);

  sparseBlocks_.push_back(block);

  int* row_records = &(block->getRowConnectivities()[0]);

  CHK_ERR( getConnectivityRecords(rowSpace_.get(), idType, fieldID,
                                  numRows, rowIDs, row_records) );

  fei::VectorSpace* colSpace = rowSpace_.get();
  if (colSpace_.get() != NULL) {
    colSpace = colSpace_.get();
  }

  int* col_records = &(block->getColConnectivities()[0]);

  CHK_ERR( getConnectivityRecords(colSpace, idType, fieldID, rowOffsets[numRows],
                                  packedColumnIDs, col_records));

  return(0);
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::initLagrangeConstraint(int constraintID,
                                               int constraintIDType,
                                               int numIDs,
                                               const int* idTypes,
                                               const int* IDs,
                                               const int* fieldIDs)
{
  if (rowSpace_.get() == NULL) ERReturn(-1);

  ConstraintType* constraint = getLagrangeConstraint(constraintID);

  CHK_ERR( rowSpace_->addDOFs(constraintIDType, 1, &constraintID) );

  if (haveColSpace_) {
    if (colSpace_.get() == NULL) {
      ERReturn(-1);
    }
    CHK_ERR( colSpace_->addDOFs(constraintIDType,
                                            1, &constraintID) );
  }

  ConstraintType* newconstraint =
    new ConstraintType(constraintID, constraintIDType,
                      false, //isSlave
                      false,  //isPenalty
                      numIDs, idTypes, IDs, fieldIDs,
                      0, //slaveOffset
                      0, //offsetIntoSlaveField
                      NULL, //weights
                      0.0, //rhsValue
                      rowSpace_.get());
 
  if (constraint != NULL) {
    if (*constraint != *newconstraint) {
      delete constraint;
    }
    else {
      delete newconstraint;
      return(0);
    }
  }

  lagrangeConstraints_[constraintID] = newconstraint;

  return(0);
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::initPenaltyConstraint(int constraintID,
                                            int constraintIDType,
                                            int numIDs,
                                            const int* idTypes,
                                            const int* IDs,
                                            const int* fieldIDs)
{
  if (rowSpace_.get() == NULL) ERReturn(-1);

  ConstraintType* constraint = getPenaltyConstraint(constraintID);

  ConstraintType* newconstraint =
    new ConstraintType(constraintID, constraintIDType,
                                              false, //isSlave
                                              true,  //isPenalty
                                              numIDs, idTypes, IDs, fieldIDs,
                                              0, //slaveOffset
                                              0, //offsetIntoSlaveField
                                              NULL, //weights
                                              0.0, //rhsValue
                                              rowSpace_.get());

  if (constraint != NULL) {
    if (*constraint != *newconstraint) {
      delete constraint;
    }
    else {
      delete newconstraint;
      return(0);
    }
  }

  penaltyConstraints_[constraintID] = newconstraint;

  return(0);
}

//------------------------------------------------------------------------------
bool fei::MatrixGraph_Impl2::hasSlaveDof(int ID, int idType)
{
  snl_fei::RecordCollection* collection = NULL;
  rowSpace_->getRecordCollection(idType, collection);
  if (collection == NULL) {
    throw std::runtime_error("fei::MatrixGraph_Impl2::hasSlaveDof: ERROR, unknown idType");
  }

  fei::Record<int>* rec = collection->getRecordWithID(ID);

  if (rec == NULL) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei::MatrixGraph_Impl2::hasSlaveDof: ERROR, specified ID ("
     << ID << ") not found.";
    throw std::runtime_error(osstr.str());
  }

  return( rec->hasSlaveDof() );
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::initSlaveConstraint(int numIDs,
                                          const int* idTypes,
                                          const int* IDs,
                                          const int* fieldIDs,
                                          int offsetOfSlave,
                                          int offsetIntoSlaveField,
                                          const double* weights,
                                          double rhsValue)
{
  if (rowSpace_.get() == NULL) ERReturn(-1);

  FEI_OSTRINGSTREAM idosstr;
  idosstr << IDs[offsetOfSlave] << fieldIDs[offsetOfSlave] << offsetIntoSlaveField;
  FEI_ISTRINGSTREAM idisstr(idosstr.str());
  int crID = IDs[offsetOfSlave];
  idisstr >> crID;

  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& dbgos = *output_stream_;

    dbgos << dbgprefix_<<"initSlaveConstraint crID=" << crID << ", slaveID="
       << IDs[offsetOfSlave];

    if (output_level_ >= fei::FULL_LOGS) {
      dbgos << " { ";
      for(int ii=0; ii<numIDs; ++ii) {
        dbgos << "("<<IDs[ii] << ","<<fieldIDs[ii]<<") ";
      }
      dbgos << "}"<<FEI_ENDL;
    }
    else dbgos << FEI_ENDL;
  }

  ConstraintType* constraint = getSlaveConstraint(crID);

  ConstraintType* newconstraint =
    new ConstraintType(crID, 0,
                       true, //isSlave
                       false,  //isPenalty
                       numIDs, idTypes, IDs, fieldIDs,
                       offsetOfSlave,
                       offsetIntoSlaveField,
                       weights,
                       rhsValue,
                       rowSpace_.get());
 
  if (constraint != NULL) {
    if (*constraint != *newconstraint) {
      if (!constraint->structurallySame(*newconstraint)) {
        FEI_OSTRINGSTREAM osstr;
        osstr << "fei::MatrixGraph_Impl2::initSlaveConstraint: slave ID "<<IDs[offsetOfSlave]
              << " is already constrained, with different connectivity. Changing the"
              << " the structure of an existing constraint is not allowed.";
        throw std::runtime_error(osstr.str());
      }
      newSlaveData_ = true;
      delete constraint;
    }
    else {
      delete newconstraint;
      return(0);
    }
  }
  else {
    newSlaveData_ = true;
  }

  slaveConstraints_[crID] = newconstraint;

  return(0);
}

//------------------------------------------------------------------------------
bool
fei::MatrixGraph_Impl2::newSlaveData()
{
  bool globalBool = false;
  fei::Allreduce(comm_, newSlaveData_, globalBool);
  newSlaveData_ = globalBool;

  return(newSlaveData_);
}

//------------------------------------------------------------------------------
fei::ConstraintType*
fei::MatrixGraph_Impl2::getLagrangeConstraint(int constraintID)
{
  std::map<int,ConstraintType*>::iterator
    c_iter = lagrangeConstraints_.find(constraintID);
  if (c_iter == lagrangeConstraints_.end()) return(NULL);

  return( (*c_iter).second );
}

//------------------------------------------------------------------------------
fei::ConstraintType*
fei::MatrixGraph_Impl2::getSlaveConstraint(int constraintID)
{
  std::map<int,ConstraintType*>::iterator
    c_iter = slaveConstraints_.find(constraintID);
  if (c_iter == slaveConstraints_.end()) return(NULL);

  return( (*c_iter).second );
}

//------------------------------------------------------------------------------
fei::ConstraintType*
fei::MatrixGraph_Impl2::getPenaltyConstraint(int constraintID)
{
  std::map<int,ConstraintType*>::iterator
    c_iter = penaltyConstraints_.find(constraintID);
  if (c_iter == penaltyConstraints_.end()) return(NULL);

  return( (*c_iter).second );
}

//------------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::initComplete()
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    (*output_stream_) <<dbgprefix_<< "initComplete" << FEI_ENDL;
  }
  if (rowSpace_.get() == NULL) {
    ERReturn(-1);
  }
  else {
    CHK_ERR( rowSpace_->initComplete() );
  }

  if (haveColSpace_ && colSpace_.get() == NULL) ERReturn(-1);
  if (haveColSpace_) {
    CHK_ERR( colSpace_->initComplete() );
  }

  if (rowSpace_->fieldMasks_.size() == 1 &&
      rowSpace_->fieldMasks_[0]->getNumFields() == 1) {
    simpleProblem_ = true;
  }

  std::vector<int>& eqnNums = rowSpace_->getEqnNumbers();
  vspcEqnPtr_ = eqnNums.size() > 0 ? &eqnNums[0] : NULL;

  //If there are slave constraints (on any processor), we need to create
  //a slave dependency matrix.
  localNumSlaves_ = slaveConstraints_.size();
  globalNumSlaves_ = 0;
  CHK_ERR( fei::GlobalSum(comm_, localNumSlaves_, globalNumSlaves_) );

  if (globalNumSlaves_ > 0) {
    //we need to allocate the slave dependency and reduction matrices too.
    CHK_ERR( createSlaveMatrices() );
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::createAlgebraicGraph(bool blockEntryGraph,
                                                 fei::Graph* graph,
                                                 bool gatherFromOverlap)
{
  //Now run the connectivity blocks, adding vertices (locations of nonzeros) to
  //the graph.

  bool wasAlreadyBlockEntry = blockEntryGraph_;

  blockEntryGraph_ = blockEntryGraph;

  for(unsigned i=0; i<sparseBlocks_.size(); ++i) {
    fei::ConnectivityBlock& cblock = *(sparseBlocks_[i]);
    CHK_ERR( addBlockToGraph_sparse(graph, &cblock) );
  }

  std::map<int,fei::ConnectivityBlock*>::const_iterator
    cb_iter = connectivityBlocks_.begin(),
    cb_end  = connectivityBlocks_.end();

  while(cb_iter != cb_end) {
    fei::ConnectivityBlock& cblock = *((*cb_iter).second);

    bool symmetricBlock = cblock.isSymmetric();

    if (symmetricBlock) {
      fei::Pattern* pattern = cblock.getRowPattern();
      if (pattern == NULL) {
        ERReturn(-1);
      }

      fei::Pattern::PatternType pType = pattern->getPatternType();
      if (pType == fei::Pattern::GENERAL || pType == fei::Pattern::SINGLE_IDTYPE) {
        CHK_ERR( addBlockToGraph_multiField_symmetric(graph, &cblock) );
      }
      else if (pType == fei::Pattern::SIMPLE) {
        CHK_ERR( addBlockToGraph_singleField_symmetric(graph, &cblock) );
      }
      else if (pType == fei::Pattern::NO_FIELD) {
        CHK_ERR( addBlockToGraph_noField_symmetric(graph, &cblock) );
      }
    }
    else {
      fei::Pattern* pattern = cblock.getRowPattern();
      fei::Pattern* colpattern = cblock.getColPattern();
      if (pattern == NULL || colpattern == NULL) {
        ERReturn(-1);
      }

      fei::Pattern::PatternType pType = pattern->getPatternType();
      fei::Pattern::PatternType colPType = colpattern->getPatternType();
      if (pType == fei::Pattern::SIMPLE && colPType == fei::Pattern::SIMPLE) {
        CHK_ERR( addBlockToGraph_singleField_nonsymmetric(graph, &cblock) );
      }
      else {
        CHK_ERR( addBlockToGraph_multiField_nonsymmetric(graph, &cblock) );
      }
    }
    ++cb_iter;
  }

  CHK_ERR( addLagrangeConstraintsToGraph(graph) );

  CHK_ERR( addPenaltyConstraintsToGraph(graph) );

  if (output_level_ >= fei::FULL_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"# before graph->gatherFromOverlap()" << FEI_ENDL;
    CHK_ERR( graph->writeLocalGraph(os) );
    CHK_ERR( graph->writeRemoteGraph(os) );
  }

  if (gatherFromOverlap) {
    CHK_ERR( graph->gatherFromOverlap() );
  }

  if (blockEntryGraph) {
    CHK_ERR( exchangeBlkEqnSizes(graph) );
  }

  if (output_level_ >= fei::FULL_LOGS && fei::numProcs(comm_)>1 &&
      output_stream_ != NULL && gatherFromOverlap) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<" after graph->gatherFromOverlap()" << FEI_ENDL;
    CHK_ERR( graph->writeLocalGraph(*output_stream_) );
  }

  if (!wasAlreadyBlockEntry) {
    setIndicesMode(POINT_ENTRY_GRAPH);
  }

  return(0);
}

//----------------------------------------------------------------------------
fei::SharedPtr<fei::SparseRowGraph>
fei::MatrixGraph_Impl2::createGraph(bool blockEntryGraph,
                                    bool localRowGraph_includeSharedRows)
{
  fei::SharedPtr<fei::SparseRowGraph> localRows;

  std::vector<int> globalOffsets;

  if (blockEntryGraph) {
    if (reducer_.get() != NULL) {
      throw std::runtime_error("fei::MatrixGraph_Impl2::createGraph ERROR, can't specify both block-entry assembly and slave-constraint reduction.");
    }

    rowSpace_->getGlobalBlkIndexOffsets(globalOffsets);
  }
  else {
    rowSpace_->getGlobalIndexOffsets(globalOffsets);
  }

  if ((int)globalOffsets.size() < localProc_+2) return localRows;

  int firstOffset = globalOffsets[localProc_];
  int lastOffset = globalOffsets[localProc_+1] - 1;

  if (reducer_.get() != NULL) {
    std::vector<int>& reduced_eqns = reducer_->getLocalReducedEqns();
    firstOffset = reduced_eqns[0];
    lastOffset = reduced_eqns[reduced_eqns.size()-1];
  }

  fei::SharedPtr<fei::Graph> inner_graph(new fei::Graph_Impl(comm_, firstOffset, lastOffset) );
  fei::SharedPtr<fei::Graph> graph;

  if (reducer_.get() == NULL) {
    graph = inner_graph;
  }
  else {
    graph.reset( new fei::GraphReducer(reducer_, inner_graph) );
  }

  bool gatherFromOverlap = !localRowGraph_includeSharedRows;
  int err = createAlgebraicGraph(blockEntryGraph, graph.get(),
                                 gatherFromOverlap);
  if (err != 0) {
    return(localRows);
  }

  localRows = fei::createSparseRowGraph(*(inner_graph->getLocalGraph()));

  remotelyOwnedGraphRows_ = fei::createSparseRowGraph(inner_graph->getRemoteGraph());

  remotelyOwnedGraphRows_->blockEntries = blockEntryGraph;

  if (localRowGraph_includeSharedRows &&
      remotelyOwnedGraphRows_->rowNumbers.size() > 0) {

    fei::SharedPtr<fei::SparseRowGraph> localPlusShared =
      snl_fei::mergeSparseRowGraphs(localRows.get(), remotelyOwnedGraphRows_.get());
    localRows = localPlusShared;
  }

  return(localRows);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::exchangeBlkEqnSizes(fei::Graph* graph)
{
  //If point-equals-block (meaning block-equations are of size 1) or
  //if blockEntryGraph_ is false (meaning we haven't constructed a block-entry-
  //graph) then there is nothing to do in this function.
  if ( rowSpace_->getPointBlockMap()->ptEqualBlk() ) {
    return(0);
  }

  snl_fei::BlkSizeMsgHandler blkHandler(rowSpace_.get(), graph, comm_);

  CHK_ERR( blkHandler.do_the_exchange() );

  return(0);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::createSlaveMatrices()
{
  //In this function we extract the dependency matrix of equation numbers 
  //from the slave-constraints locally on each processor, then gather those
  //slave-equations onto each processor (so that each processor holds ALL
  //slave-equations).
  //Then we remove any couplings that may exist among the slave-equations and
  //finally create a FillableMat object to hold the final dependency matrix D_.
  //

  if (!newSlaveData()) {
    return(0);
  }

  std::vector<int>& eqnNums = rowSpace_->getEqnNumbers();
  vspcEqnPtr_ = eqnNums.size() > 0 ? &eqnNums[0] : NULL;

  fei::SharedPtr<fei::FillableMat> local_D(new fei::FillableMat);
  fei::SharedPtr<fei::CSVec> local_g(new fei::CSVec);

  std::vector<int> masterEqns;
  std::vector<double> masterCoefs;

  std::map<int,ConstraintType*>::const_iterator
    cr_iter = slaveConstraints_.begin(),
    cr_end  = slaveConstraints_.end();

  for(; cr_iter != cr_end; ++cr_iter) {
    ConstraintType* cr = (*cr_iter).second;

    fei::Record<int>* slaveRecord = cr->getSlave();
    int slaveIDType = cr->getIDType();
    int slaveFieldID = cr->getSlaveFieldID();
    int offsetIntoSlaveField = cr->getOffsetIntoSlaveField();
    int slaveEqn = -1;
    CHK_ERR( rowSpace_->getGlobalIndex(slaveIDType, slaveRecord->getID(),
                                       slaveFieldID, 0, offsetIntoSlaveField,
                                       slaveEqn) );

    std::vector<int>& masterRecords_vec = cr->getMasters();
    int* masterRecords = &masterRecords_vec[0];
    std::vector<snl_fei::RecordCollection*>& masterRecColls = cr->getMasterRecordCollections();
    std::vector<int>& masterIDTypes = cr->getMasterIDTypes();
    std::vector<int>& masterFieldIDs = cr->getMasterFieldIDs();
    std::vector<double>& masterWeights = cr->getMasterWeights();
    double* masterWtPtr = &masterWeights[0];

    masterEqns.resize(masterWeights.size());
    masterCoefs.resize(masterWeights.size());

    int* masterEqnsPtr = &(masterEqns[0]);
    double* masterCoefsPtr = &(masterCoefs[0]);

    int offset = 0;
    for(size_t j=0; j<masterIDTypes.size(); ++j) {
      fei::Record<int>* masterRecord = masterRecColls[j]->getRecordWithLocalID(masterRecords[j]);
      int* eqnNumbers = vspcEqnPtr_+masterRecord->getOffsetIntoEqnNumbers();
      fei::FieldMask* mask = masterRecord->getFieldMask();
      int eqnOffset = 0;
      if (!simpleProblem_) {
        mask->getFieldEqnOffset(masterFieldIDs[j], eqnOffset);
      }

      unsigned fieldSize = rowSpace_->getFieldSize(masterFieldIDs[j]);

      int eqn = eqnNumbers[eqnOffset];
      for(unsigned k=0; k<fieldSize; ++k) {
        masterEqnsPtr[offset++] = eqn+k;
      }
    }

    double fei_eps = 1.e-49;

    offset = 0;
    for(size_t jj=0; jj<masterEqns.size(); ++jj) {
      if (std::abs(masterWtPtr[jj]) > fei_eps) {
        masterCoefsPtr[offset] = masterWtPtr[jj];
        masterEqnsPtr[offset++] = masterEqnsPtr[jj];
      }
    }

    if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
      bool log_eqn = isLogEqn(slaveEqn);
      if (!log_eqn) {
        for(unsigned ii=0; ii<masterEqns.size(); ++ii) {
          if (isLogEqn(masterEqnsPtr[ii])) {
            log_eqn = true;
            break;
          }
        }
      }

      if (log_eqn) {
        FEI_OSTREAM& os = *output_stream_;
        os << "createSlaveMatrices: " << slaveEqn << " = ";
        for(unsigned ii=0; ii<masterEqns.size(); ++ii) {
          if (ii!=0) os << "+ ";
          os << masterCoefsPtr[ii]<<"*"<<masterEqnsPtr[ii]<<" ";
        }
        os << FEI_ENDL;
      }
    }

    local_D->putRow(slaveEqn, masterEqnsPtr, masterCoefsPtr, offset);

    if (std::abs(cr->getRHSValue()) > fei_eps) {
      fei::put_entry(*local_g, slaveEqn, cr->getRHSValue());
    }
  }

  if (D_.get() == NULL) {
    D_.reset(new fei::FillableMat);
  }

  if (g_.get() == NULL) {
    g_.reset(new fei::CSVec);
  }

#ifndef FEI_SER
  if (numProcs_ > 1) {
    fei::impl_utils::global_union(comm_, *local_D, *D_);
    fei::impl_utils::global_union(comm_, *local_g, *g_);
  }
  else {
    *D_ = *local_D;
    *g_ = *local_g;
  }
#else
  *D_ = *local_D;
  *g_ = *local_g;
#endif

  if (output_level_ >= fei::FULL_LOGS && output_stream_ != NULL) {
    (*output_stream_) << "#  D_ (pre-removeCouplings):"<<FEI_ENDL;
    (*output_stream_) << *D_;
  }

  int levelsOfCoupling = fei::impl_utils::remove_couplings(*D_);
  (void)levelsOfCoupling;

  if (reducer_.get() == NULL) {
    reducer_.reset(new fei::Reducer(D_, g_, comm_));

    std::vector<int> indices;
    rowSpace_->getIndices_Owned(indices);

    reducer_->setLocalUnreducedEqns(indices);
  }
  else {
    reducer_->initialize();
  }

  if (output_level_ >= fei::FULL_LOGS && output_stream_ != NULL) {
    (*output_stream_) << "#  D_:"<<FEI_ENDL;
    (*output_stream_) << *D_;
  }

  double fei_eps = 1.e-49;

  g_nonzero_ = false;
  for(size_t j=0; j<g_->size(); ++j) {
    double coef = g_->coefs()[j];
    if (std::abs(coef) > fei_eps) {
      g_nonzero_ = true;
    }
  }

  newSlaveData_ = false;

  return(0);
}

//----------------------------------------------------------------------------
fei::SharedPtr<fei::Reducer>
fei::MatrixGraph_Impl2::getReducer()
{
  return( reducer_ );
}

//----------------------------------------------------------------------------
fei::SharedPtr<fei::SparseRowGraph>
fei::MatrixGraph_Impl2::getRemotelyOwnedGraphRows()
{
  return( remotelyOwnedGraphRows_ );
}

//----------------------------------------------------------------------------
void fei::MatrixGraph_Impl2::getConstrainedIndices(std::vector<int>& crindices) const
{
  if (constrained_indices_.empty()) {
    crindices.clear();
    return;
  }

  std::set<int>::const_iterator
    s_iter = constrained_indices_.begin(),
    s_end = constrained_indices_.end();

  crindices.resize(constrained_indices_.size());

  int offset = 0;
  for(; s_iter != s_end; ++s_iter) {
    crindices[offset++] = *s_iter;
  }
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::addLagrangeConstraintsToGraph(fei::Graph* graph)
{
  std::vector<int> indices;
  std::map<int,ConstraintType*>::const_iterator
    cr_iter = lagrangeConstraints_.begin(),
    cr_end  = lagrangeConstraints_.end();

  constrained_indices_.clear();

  while(cr_iter != cr_end) {
    ConstraintType* cr = (*cr_iter).second;
    int crID = cr->getConstraintID();

    CHK_ERR( getConstraintConnectivityIndices(cr, indices) );

    int numIndices = indices.size();
    int* indicesPtr = &(indices[0]);

    for(int i=0; i<numIndices; ++i) {
      constrained_indices_.insert(indicesPtr[i]);
    }

    int crEqnRow = -1, tmp;
    if (blockEntryGraph_) {
      CHK_ERR( rowSpace_->getGlobalBlkIndex(cr->getIDType(),
                                                   crID, crEqnRow) );
      cr->setBlkEqnNumber(crEqnRow);
      CHK_ERR( rowSpace_->getGlobalIndex(cr->getIDType(),
                                                crID, tmp) );
      cr->setEqnNumber(tmp);
    }
    else {
      CHK_ERR( rowSpace_->getGlobalIndex(cr->getIDType(),
                                                crID, crEqnRow) );
      cr->setEqnNumber(crEqnRow);
      CHK_ERR( rowSpace_->getGlobalBlkIndex(cr->getIDType(),
                                                   crID, tmp) );
      cr->setBlkEqnNumber(tmp);
    }

    //now add the row contribution
    CHK_ERR( graph->addIndices(crEqnRow, numIndices, &(indices[0])) );

    //Let's add a diagonal entry to the graph for this constraint-equation,
    //just in case we need to fiddle with this equation later (e.g. discard the
    //constraint equation and replace it with a dirichlet boundary condition...).
    CHK_ERR( graph->addIndices(crEqnRow, 1, &crEqnRow) );

    //and finally, add the column contribution (which is simply the transpose
    //of the row contribution)
    for(int k=0; k<numIndices; ++k) {
      CHK_ERR( graph->addIndices(indicesPtr[k], 1, &crEqnRow) );
    }

    ++cr_iter;
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::
getConstraintConnectivityIndices(ConstraintType* cr,
                                 std::vector<int>& globalIndices)
{
  std::vector<int>& fieldSizes = tmpIntArray1_;
  std::vector<int>& ones = tmpIntArray2_;

  std::vector<int>& constrainedRecords = cr->getMasters();
  std::vector<int>& constrainedFieldIDs = cr->getMasterFieldIDs();
  std::vector<snl_fei::RecordCollection*>& recordCollections = cr->getMasterRecordCollections();

  int len = constrainedRecords.size();
  fieldSizes.resize(len);

  ones.assign(len, 1);

  int numIndices = 0;
  if (blockEntryGraph_) {
    numIndices = len;
  }
  else {
    for(int j=0; j<len; ++j) {
      unsigned fieldSize = rowSpace_->getFieldSize(constrainedFieldIDs[j]);
      fieldSizes[j] = fieldSize;
      numIndices += fieldSize;
    }
  }

  globalIndices.resize(numIndices);

  int checkNum;
  CHK_ERR( getConnectivityIndices_multiField(&recordCollections[0],
                                             &constrainedRecords[0],
                                             len, &ones[0],
                                             &constrainedFieldIDs[0],
                                             &fieldSizes[0],
                                             numIndices, &globalIndices[0],
                                             checkNum) );
  if (numIndices != checkNum) {
    ERReturn(-1);
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::addPenaltyConstraintsToGraph(fei::Graph* graph)
{
  std::vector<int> indices;
  std::map<int,ConstraintType*>::const_iterator
    cr_iter = penaltyConstraints_.begin(),
    cr_end  = penaltyConstraints_.end();

  while(cr_iter != cr_end) {
    ConstraintType* cr = (*cr_iter).second;

    CHK_ERR( getConstraintConnectivityIndices(cr, indices) );

    int numIndices = indices.size();

    //now add the symmetric contributions to the graph
    CHK_ERR( graph->addSymmetricIndices(numIndices, &(indices[0])) );

    ++cr_iter;
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::compareStructure(const fei::MatrixGraph& matrixGraph,
                                       bool& equivalent) const
{
  equivalent = false;
  int localCode = 1;
  int globalCode = 0;

  int myNumBlocks = getNumConnectivityBlocks();
  int numBlocks = matrixGraph.getNumConnectivityBlocks();

  if (myNumBlocks != numBlocks) {
    CHK_ERR( fei::GlobalMax(comm_, localCode, globalCode) );
    equivalent = globalCode > 0 ? false : true;
    return(0);
  }

  if (numBlocks > 0) {
    std::vector<int> myBlockIDs;
    std::vector<int> blockIDs;

    CHK_ERR( getConnectivityBlockIDs(myBlockIDs) );
    CHK_ERR( matrixGraph.getConnectivityBlockIDs(blockIDs) );

    for(int i=0; i<numBlocks; ++i) {
      if (myBlockIDs[i] != blockIDs[i]) {
        CHK_ERR( fei::GlobalMax(comm_, localCode, globalCode) );
        equivalent = globalCode > 0 ? false : true;
        return(0);
      }

      const fei::ConnectivityBlock* mycblock = getConnectivityBlock(myBlockIDs[i]);
      const fei::ConnectivityBlock* cblock = matrixGraph.getConnectivityBlock(blockIDs[i]);

      int myNumLists = mycblock->getConnectivityIDs().size();
      int numLists = cblock->getConnectivityIDs().size();

      if (myNumLists != numLists ||
          mycblock->isSymmetric() != cblock->isSymmetric()) {
        CHK_ERR( fei::GlobalMax(comm_, localCode, globalCode) );
        equivalent = globalCode > 0 ? false : true;
        return(0);
      }

      int myNumIDsPerList = mycblock->getRowPattern()->getNumIDs();
      int numIDsPerList = cblock->getRowPattern()->getNumIDs();

      if (myNumIDsPerList != numIDsPerList) {
        CHK_ERR( fei::GlobalMax(comm_, localCode, globalCode) );
        equivalent = globalCode > 0 ? false : true;
        return(0);
      }
    }
  }

  int numMyLagrangeConstraints = getLocalNumLagrangeConstraints();
  int numMySlaveConstraints = getGlobalNumSlaveConstraints();
  int numLagrangeConstraints = matrixGraph.getLocalNumLagrangeConstraints();
  int numSlaveConstraints = matrixGraph.getGlobalNumSlaveConstraints();

  if (numMyLagrangeConstraints != numLagrangeConstraints ||
      numMySlaveConstraints != numSlaveConstraints) {
    CHK_ERR( fei::GlobalMax(comm_, localCode, globalCode) );
    equivalent = globalCode > 0 ? false : true;
    return(0);
  }

  localCode = 0;
  CHK_ERR( fei::GlobalMax(comm_, localCode, globalCode) );
  equivalent = globalCode > 0 ? false : true;
  return(0);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::getNumConnectivityBlocks() const
{
  return(connectivityBlocks_.size());
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::getConnectivityBlockIDs(std::vector<int>& blockIDs) const
{
  blockIDs.resize(connectivityBlocks_.size());

  std::map<int,fei::ConnectivityBlock*>::const_iterator
    cdb_iter = connectivityBlocks_.begin();

  for(size_t i=0; i<blockIDs.size(); ++i, ++cdb_iter) {
    int blockID = (*cdb_iter).first;
    blockIDs[i] = blockID;
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::getNumIDsPerConnectivityList(int blockID) const
{
  const fei::ConnectivityBlock* cblock = getConnectivityBlock(blockID);
  if (cblock == NULL) return(-1);

  const fei::Pattern* pattern = cblock->getRowPattern();
  return(pattern->getNumIDs());
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::getConnectivityNumIndices(int blockID) const
{
  const fei::ConnectivityBlock* cblock = getConnectivityBlock(blockID);
  if (cblock == NULL) return(-1);

  const fei::Pattern* pattern = cblock->getRowPattern();
  return( blockEntryGraph_ ?
    pattern->getNumIDs() : pattern->getNumIndices());
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::getConnectivityNumIndices(int blockID,
                                                  int& numRowIndices,
                                                    int& numColIndices)
{
  fei::ConnectivityBlock* cblock = getConnectivityBlock(blockID);
  if (cblock == NULL) return(-1);

  fei::Pattern* pattern = cblock->getRowPattern();
  numRowIndices = pattern->getNumIndices();
  fei::Pattern* colpattern = cblock->isSymmetric() ?
    pattern : cblock->getColPattern();
  numColIndices = colpattern->getNumIndices();

  return(0);
}

//----------------------------------------------------------------------------
const fei::ConnectivityBlock* fei::MatrixGraph_Impl2::getConnectivityBlock(int blockID) const
{
  std::map<int,fei::ConnectivityBlock*>::const_iterator
    c_iter = connectivityBlocks_.find(blockID);
  if (c_iter == connectivityBlocks_.end()) return(0);

  return( (*c_iter).second );
}

//----------------------------------------------------------------------------
fei::ConnectivityBlock* fei::MatrixGraph_Impl2::getConnectivityBlock(int blockID)
{
  std::map<int,fei::ConnectivityBlock*>::iterator
    c_iter = connectivityBlocks_.find(blockID);
  if (c_iter == connectivityBlocks_.end()) return(0);

  return( (*c_iter).second );
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::getConnectivityIndices(int blockID,
                                             int connectivityID,
                                             int indicesAllocLen,
                                             int* indices,
                                             int& numIndices)
{
  fei::ConnectivityBlock* cblock = getConnectivityBlock(blockID);
  if (cblock == NULL) return(-1);

  std::vector<int>& eqnNums = rowSpace_->getEqnNumbers();
  vspcEqnPtr_ = eqnNums.size() > 0 ? &eqnNums[0] : NULL;

  fei::Pattern* pattern = cblock->getRowPattern();
  numIndices = pattern->getNumIndices();

  int len = numIndices > indicesAllocLen ? indicesAllocLen : numIndices;

  int* records = cblock->getRowConnectivity(connectivityID);
  if (records == NULL) {
    ERReturn(-1);
  }

  fei::Pattern::PatternType pType = pattern->getPatternType();

  if (pType == fei::Pattern::GENERAL || pType == fei::Pattern::SINGLE_IDTYPE) {
    const int* numFieldsPerID = pattern->getNumFieldsPerID();
    const int* fieldIDs = pattern->getFieldIDs();
    int totalNumFields = pattern->getTotalNumFields();

    std::vector<int> fieldSizes(totalNumFields);

    for(int j=0; j<totalNumFields; ++j) {
      fieldSizes[j] = snl_fei::getFieldSize(fieldIDs[j], rowSpace_.get(),
                                            colSpace_.get());
    }

    CHK_ERR( getConnectivityIndices_multiField(pattern->getRecordCollections(),
                                            records, pattern->getNumIDs(),
                                               numFieldsPerID,
                                               fieldIDs, &fieldSizes[0],
                                               len, indices, numIndices) );
  }
  else if (pType == fei::Pattern::SIMPLE) {
    const int* fieldIDs = pattern->getFieldIDs();

    int fieldID = fieldIDs[0];
    unsigned fieldSize = snl_fei::getFieldSize(fieldID,
                                               rowSpace_.get(),
                                               colSpace_.get());

    CHK_ERR( getConnectivityIndices_singleField(pattern->getRecordCollections(),
                                            records, pattern->getNumIDs(),
                                                fieldID, fieldSize,
                                                len, indices, numIndices) );
  }
  else if (pType == fei::Pattern::NO_FIELD) {
    CHK_ERR( getConnectivityIndices_noField(pattern->getRecordCollections(),
                                            records, pattern->getNumIDs(),
                                            len, indices, numIndices) );
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::getConnectivityIndices(int blockID,
                                             int connectivityID,
                                             int rowIndicesAllocLen,
                                             int* rowIndices,
                                             int& numRowIndices,
                                             int colIndicesAllocLen,
                                             int* colIndices,
                                             int& numColIndices)
{
  fei::ConnectivityBlock* cblock = getConnectivityBlock(blockID);
  if (cblock == NULL) return(-1);

  std::vector<int>& eqnNums = rowSpace_->getEqnNumbers();
  vspcEqnPtr_ = eqnNums.size() > 0 ? &eqnNums[0] : NULL;

  fei::Pattern* pattern = cblock->getRowPattern();
  numRowIndices = pattern->getNumIndices();

  int len = numRowIndices > rowIndicesAllocLen ?
    rowIndicesAllocLen : numRowIndices;

  int* records = cblock->getRowConnectivity(connectivityID);
  if (records == NULL) {
    ERReturn(-1);
  }

  fei::Pattern::PatternType pType = pattern->getPatternType();

  if (pType == fei::Pattern::GENERAL || pType == fei::Pattern::SINGLE_IDTYPE) {
    const int* numFieldsPerID = pattern->getNumFieldsPerID();
    const int* fieldIDs = pattern->getFieldIDs();
    int totalNumFields = pattern->getTotalNumFields();

    std::vector<int> fieldSizes(totalNumFields);

    for(int j=0; j<totalNumFields; ++j) {
      fieldSizes[j] = snl_fei::getFieldSize(fieldIDs[j], rowSpace_.get(),
                                            colSpace_.get());
    }

    CHK_ERR( getConnectivityIndices_multiField(pattern->getRecordCollections(),
                                               records, pattern->getNumIDs(),
                                               numFieldsPerID,
                                               fieldIDs, &fieldSizes[0],
                                               len, rowIndices, numRowIndices) );
  }
  else if (pType == fei::Pattern::SIMPLE) {
    const int* fieldIDs = pattern->getFieldIDs();

    int fieldID = fieldIDs[0];
    unsigned fieldSize = snl_fei::getFieldSize(fieldID, rowSpace_.get(),
                                            colSpace_.get());

    CHK_ERR( getConnectivityIndices_singleField(pattern->getRecordCollections(),
                                                records, pattern->getNumIDs(),
                                                fieldID, fieldSize,
                                                len, rowIndices, numRowIndices) );
  }
  else if (pType == fei::Pattern::NO_FIELD) {
    CHK_ERR( getConnectivityIndices_noField(pattern->getRecordCollections(),
                                            records, pattern->getNumIDs(),
                                            len, rowIndices, numRowIndices) );
  }

  fei::Pattern* colpattern = cblock->isSymmetric() ? 
    pattern : cblock->getColPattern();
  pType = colpattern->getPatternType();
  numColIndices = colpattern->getNumIndices();
  len = numColIndices > colIndicesAllocLen ?
    colIndicesAllocLen : numColIndices;

  if (!cblock->isSymmetric()) {
    records = cblock->getColConnectivity(connectivityID);
  }
  if (records == NULL) {
    return(-1);
  }

  if (pType == fei::Pattern::GENERAL || pType == fei::Pattern::SINGLE_IDTYPE) {
    const int* numFieldsPerID = colpattern->getNumFieldsPerID();
    const int* fieldIDs = colpattern->getFieldIDs();
    int totalNumFields = colpattern->getTotalNumFields();

    std::vector<int> fieldSizes(totalNumFields);

    for(int j=0; j<totalNumFields; ++j) {
      fieldSizes[j] = snl_fei::getFieldSize(fieldIDs[j], rowSpace_.get(),
                                            colSpace_.get());
    }

    CHK_ERR( getConnectivityIndices_multiField(colpattern->getRecordCollections(),
                                               records, colpattern->getNumIDs(),
                                               numFieldsPerID,
                                               fieldIDs, &fieldSizes[0],
                                               len, colIndices, numColIndices) );
  }
  else if (pType == fei::Pattern::SIMPLE) {
    const int* fieldIDs = colpattern->getFieldIDs();

    int fieldID = fieldIDs[0];
    unsigned fieldSize = snl_fei::getFieldSize(fieldID, rowSpace_.get(),
                                               colSpace_.get());

    CHK_ERR( getConnectivityIndices_singleField(colpattern->getRecordCollections(),
                                                records, colpattern->getNumIDs(),
                                                fieldID, fieldSize,
                                                len, colIndices, numColIndices) );
  }
  else if (pType == fei::Pattern::NO_FIELD) {
    CHK_ERR( getConnectivityIndices_noField(colpattern->getRecordCollections(),
                                            records, colpattern->getNumIDs(),
                                            len, colIndices, numColIndices) );
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::getLocalNumLagrangeConstraints() const
{
  return( lagrangeConstraints_.size() );
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::addBlockToGraph_multiField_symmetric(fei::Graph* graph,
                                                        fei::ConnectivityBlock* cblock)
{
  fei::Pattern* pattern = cblock->getRowPattern();
  int j;
  int numIDs = pattern->getNumIDs();
  int numIndices = blockEntryGraph_ ?
    pattern->getNumIDs() : pattern->getNumIndices();
  int checkNumIndices = numIndices;
  std::vector<int> indices(numIndices);
  int* indicesPtr = &indices[0];

  const int* numFieldsPerID = pattern->getNumFieldsPerID();
  const int* fieldIDs = pattern->getFieldIDs();
  int totalNumFields = pattern->getTotalNumFields();

  std::vector<int> fieldSizes(totalNumFields);

  for(j=0; j<totalNumFields; ++j) {
    fieldSizes[j] = snl_fei::getFieldSize(fieldIDs[j], rowSpace_.get(),
                                          colSpace_.get());
  }

  std::map<int,int>& connIDs = cblock->getConnectivityIDs();
  std::vector<int>& values = cblock->getRowConnectivities();

  std::map<int,int>::iterator
    c_iter = connIDs.begin(),
    c_end  = connIDs.end();

  for(; c_iter != c_end; ++c_iter) {
    int offset = c_iter->second;
    int* records = &values[offset*numIDs];

    CHK_ERR( getConnectivityIndices_multiField(pattern->getRecordCollections(),
                                                records, numIDs,
                                                 numFieldsPerID,
                                                 fieldIDs,
                                                 &fieldSizes[0],
                                                 numIndices,
                                                 indicesPtr,
                                                 checkNumIndices) );

    if (checkNumIndices != numIndices) {
      ERReturn(-1);
    }

    if (output_level_ > fei::BRIEF_LOGS && output_stream_ != NULL) {
      FEI_OSTREAM& os = *output_stream_;

      const snl_fei::RecordCollection*const* recordColls = pattern->getRecordCollections();
      unsigned thisoffset = 0;
      for(int ii=0; ii<numIDs; ++ii) {
        const fei::Record<int>* record = recordColls[ii]->getRecordWithLocalID(records[ii]);
        int ID = record->getID();
        os << dbgprefix_<<"scatterIndices: ID=" <<ID<<": ";
        int num = pattern->getNumIndicesPerID()[ii];
        for(int jj=0; jj<num; ++jj) {
          os << indicesPtr[thisoffset++] << " ";
        }
        os << FEI_ENDL;
      }

      for(int ii=0; ii<numIndices; ++ii) {
        if (isLogEqn(indicesPtr[ii])) {
          os << "adding Symm inds: ";
          for(int jj=0; jj<numIndices; ++jj) {
            os << indicesPtr[jj] << " ";
          }
          os << FEI_ENDL;
          break;
        }
      }
    }

    //now we have the indices array, so we're ready to push it into
    //the graph container
    if (numIndices == numIDs || !cblock->isDiagonal()) {
      CHK_ERR( graph->addSymmetricIndices(numIndices, indicesPtr,
                                          cblock->isDiagonal()) );
    }
    else {
      int ioffset = 0;
      int* fieldSizesPtr = &fieldSizes[0];
      for(int i=0; i<numIDs; ++i) {
        CHK_ERR( graph->addSymmetricIndices(fieldSizesPtr[i], &(indicesPtr[ioffset])));
        ioffset += fieldSizesPtr[i];
      }
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::addBlockToGraph_multiField_nonsymmetric(fei::Graph* graph,
                                                        fei::ConnectivityBlock* cblock)
{
  fei::Pattern* pattern = cblock->getRowPattern();
  fei::Pattern* colpattern = cblock->getColPattern();
  int j;
  int numIDs = pattern->getNumIDs();
  int numIndices = blockEntryGraph_ ? numIDs : pattern->getNumIndices();
  int checkNumIndices = numIndices;
  std::vector<int> indices(numIndices);
  int* indicesPtr = &indices[0];

  int numColIDs = colpattern->getNumIDs();
  int numColIndices = blockEntryGraph_ ? numColIDs : colpattern->getNumIndices();
  int checkNumColIndices = numColIndices;
  std::vector<int> colindices(numColIndices);
  int* colindicesPtr = &colindices[0];

  const int* numFieldsPerID = pattern->getNumFieldsPerID();
  const int* fieldIDs = pattern->getFieldIDs();
  int totalNumFields = pattern->getTotalNumFields();

  const int* numFieldsPerColID = colpattern->getNumFieldsPerID();
  const int* colfieldIDs = colpattern->getFieldIDs();
  int totalNumColFields = colpattern->getTotalNumFields();

  std::vector<int> fieldSizes(totalNumFields);
  std::vector<int> colfieldSizes(totalNumColFields);

  for(j=0; j<totalNumFields; ++j) {
    fieldSizes[j] = snl_fei::getFieldSize(fieldIDs[j], rowSpace_.get(),
                                          colSpace_.get());
  }
  for(j=0; j<totalNumColFields; ++j) {
    colfieldSizes[j] = snl_fei::getFieldSize(colfieldIDs[j], rowSpace_.get(),
                                           colSpace_.get());
  }

  std::map<int,int>& connIDs = cblock->getConnectivityIDs();
  std::vector<int>& rowrecords = cblock->getRowConnectivities();
  std::vector<int>& colrecords = cblock->getColConnectivities();

  std::map<int,int>::iterator
   c_iter = connIDs.begin(),
   c_end  = connIDs.end();

  for(; c_iter != c_end; ++c_iter) {
    int offset = c_iter->second;
    int* records = &rowrecords[offset*numIDs];

    int* colRecords = &colrecords[offset*numColIDs];

    CHK_ERR( getConnectivityIndices_multiField(pattern->getRecordCollections(),
                                                records, numIDs,
                                                 numFieldsPerID,
                                                 fieldIDs,
                                                 &fieldSizes[0],
                                                 numIndices,
                                                 indicesPtr,
                                                 checkNumIndices) );

    if (checkNumIndices != numIndices) {
      ERReturn(-1);
    }

      
    CHK_ERR( getConnectivityIndices_multiField(colpattern->getRecordCollections(),
                                                colRecords, numColIDs,
                                                 numFieldsPerColID,
                                                 colfieldIDs,
                                                 &colfieldSizes[0],
                                                 numColIndices,
                                                 colindicesPtr,
                                                 checkNumColIndices) );

    if (checkNumColIndices != numColIndices) {
      ERReturn(-1);
    }

    //now we have the indices arrays, so we're ready to push them into
    //the graph container
    for(int r=0; r<numIndices; ++r) {
      CHK_ERR( graph->addIndices(indicesPtr[r], numColIndices, colindicesPtr));
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::getConnectivityIndices_multiField(const snl_fei::RecordCollection*const* recordCollections,
                                                          int* records,
                                                          int numRecords,
                                                          const int* numFieldsPerID,
                                                          const int* fieldIDs,
                                                          const int* fieldSizes,
                                                          int indicesAllocLen,
                                                          int* indices,
                                                          int& numIndices)
{
  numIndices = 0;
  int fld_offset = 0;

  for(int i=0; i<numRecords; ++i) {
    const fei::Record<int>* record = recordCollections[i]->getRecordWithLocalID(records[i]);
    if (record==NULL) continue;

    if (blockEntryGraph_) {
      indices[numIndices++] = record->getNumber();
      continue;
    }

    const fei::FieldMask* fieldMask = record->getFieldMask();
    int* eqnNumbers = vspcEqnPtr_ + record->getOffsetIntoEqnNumbers();

    for(int nf=0; nf<numFieldsPerID[i]; ++nf) {
      int eqnOffset = 0;
      if (!simpleProblem_) {
        fieldMask->getFieldEqnOffset(fieldIDs[fld_offset], eqnOffset);
      }

      for(int fs=0; fs<fieldSizes[fld_offset]; ++fs) {
        indices[numIndices++] = eqnNumbers[eqnOffset+fs];
      }

      ++fld_offset;
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::addBlockToGraph_singleField_symmetric(fei::Graph* graph,
                                                        fei::ConnectivityBlock* cblock)
{
  fei::Pattern* pattern = cblock->getRowPattern();
  int numIDs = pattern->getNumIDs();
  int numIndices = blockEntryGraph_ ? numIDs : pattern->getNumIndices();
  int checkNumIndices = numIndices;
  std::vector<int> indices(numIndices);
  int* indicesPtr = &indices[0];

  const int* fieldIDs = pattern->getFieldIDs();

  int fieldID = fieldIDs[0];
  unsigned fieldSize = snl_fei::getFieldSize(fieldID, rowSpace_.get(),
                                             colSpace_.get());

  std::map<int,int>& connIDs = cblock->getConnectivityIDs();
  std::vector<int>& rowrecords = cblock->getRowConnectivities();

  std::map<int,int>::iterator
    c_iter = connIDs.begin(),
    c_end  = connIDs.end();

  for(; c_iter != c_end; ++c_iter) {
    int offset = c_iter->second;
    int* records = &rowrecords[offset*numIDs];

    CHK_ERR( getConnectivityIndices_singleField(pattern->getRecordCollections(),
                                                records, numIDs,
                                                  fieldID, fieldSize,
                                                  checkNumIndices,
                                                  indicesPtr,
                                                  numIndices) );

    if (checkNumIndices != numIndices) {
      ERReturn(-1);
    }

    if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
        for(int ii=0; ii<numIndices; ++ii) {
          if (isLogEqn(indicesPtr[ii])) {
            FEI_OSTREAM& os = *output_stream_;
            os << "adding Symm inds: ";
            for(int jj=0; jj<numIndices; ++jj) {
              os << indicesPtr[jj] << " ";
            }
            os << FEI_ENDL;
            break;
          }
        }
    }

    //now we have the indices array, so we're ready to push it into
    //the graph container
    if (numIndices == numIDs || !cblock->isDiagonal()) {
      CHK_ERR( graph->addSymmetricIndices(numIndices, indicesPtr,
                                            cblock->isDiagonal()));
    }
    else {
      int ioffset = 0;
      for(int i=0; i<numIDs; ++i) {
        CHK_ERR( graph->addSymmetricIndices(fieldSize, &(indicesPtr[ioffset])));
        ioffset += fieldSize;
      }
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::addBlockToGraph_singleField_nonsymmetric(fei::Graph* graph,
                                                        fei::ConnectivityBlock* cblock)
{
  fei::Pattern* pattern = cblock->getRowPattern();
  fei::Pattern* colpattern = cblock->getColPattern();
  int numIDs = pattern->getNumIDs();
  int numIndices = blockEntryGraph_ ? numIDs : pattern->getNumIndices();
  int checkNumIndices = numIndices;
  std::vector<int> indices(numIndices);
  int* indicesPtr = &indices[0];

  int numColIDs = colpattern->getNumIDs();
  int numColIndices = blockEntryGraph_ ?
    numColIDs : colpattern->getNumIndices();
  int checkNumColIndices = numColIndices;
  std::vector<int> colindices(numColIndices);
  int* colindicesPtr = &colindices[0];

  const int* fieldIDs = pattern->getFieldIDs();

  int rowFieldID = fieldIDs[0];
  int rowFieldSize = snl_fei::getFieldSize(rowFieldID, rowSpace_.get(),
                                           colSpace_.get());

  std::map<int,int>& connIDs = cblock->getConnectivityIDs();
  std::vector<int>& rowrecords = cblock->getRowConnectivities();
  std::vector<int>& colrecords = cblock->getColConnectivities();

  int colFieldID = colpattern->getFieldIDs()[0];
  int colFieldSize = snl_fei::getFieldSize(colFieldID, rowSpace_.get(),
                                           colSpace_.get());

  std::map<int,int>::iterator
    c_iter = connIDs.begin(),
    c_end  = connIDs.end();

  for(; c_iter != c_end; ++c_iter) {
    int offset = c_iter->second;
    int* records = &rowrecords[offset*numIDs];

    int* colRecords = &colrecords[offset*numColIDs];

    if (blockEntryGraph_) {
      rowSpace_->getGlobalBlkIndicesL(numIDs, pattern->getRecordCollections(),
                                      records, checkNumIndices,
                                            indicesPtr, numIndices);

      colSpace_->getGlobalBlkIndicesL(numColIDs, colpattern->getRecordCollections(),
                                      colRecords, checkNumColIndices,
                                            colindicesPtr, numColIndices);
    }
    else {
      rowSpace_->getGlobalIndicesL(numIDs, pattern->getRecordCollections(),
                                      records, rowFieldID,
                                         rowFieldSize, checkNumIndices,
                                         indicesPtr, numIndices);

      colSpace_->getGlobalIndicesL(numColIDs, colpattern->getRecordCollections(),
                                     colRecords, colFieldID,
                                         colFieldSize, checkNumColIndices,
                                         colindicesPtr, numColIndices);
    }

    if (checkNumIndices != numIndices || checkNumColIndices != numColIndices) {
      ERReturn(-1);
    }

    for(int r=0; r<numIndices; ++r) {
      CHK_ERR( graph->addIndices(indicesPtr[r], numColIndices, colindicesPtr));
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::getConnectivityIndices_singleField(const snl_fei::RecordCollection*const* recordCollections,
                                                             int* records,
                                                             int numRecords,
                                                             int fieldID,
                                                             int fieldSize,
                                                             int indicesAllocLen,
                                                             int* indices,
                                                             int& numIndices)
{
  numIndices = 0;

  for(int i=0; i<numRecords; ++i) {
    if (numIndices == indicesAllocLen) break;

    const fei::Record<int>* record = recordCollections[i]->getRecordWithLocalID(records[i]);

    if (blockEntryGraph_) {
      indices[numIndices++] = record->getNumber();
      continue;
    }

    int* eqnNumbers = vspcEqnPtr_+record->getOffsetIntoEqnNumbers();

    int eqnOffset = 0;
    if (!simpleProblem_) {
      const fei::FieldMask* fieldMask = record->getFieldMask();
      fieldMask->getFieldEqnOffset(fieldID, eqnOffset);
    }

    indices[numIndices++] = eqnNumbers[eqnOffset];
    if (fieldSize > 1) {
      for(int fs=1; fs<fieldSize; ++fs) {
        indices[numIndices++] = eqnNumbers[eqnOffset+fs];
      }
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::getConnectivityIndices_noField(const snl_fei::RecordCollection*const* recordCollections,
                                                         int* records,
                                                         int numRecords,
                                                         int indicesAllocLen,
                                                         int* indices,
                                                         int& numIndices)
{
  numIndices = 0;

  for(int i=0; i<numRecords; ++i) {

    const fei::Record<int>* record = recordCollections[i]->getRecordWithLocalID(records[i]);
    int* eqnNumbers = vspcEqnPtr_+record->getOffsetIntoEqnNumbers();

    if (blockEntryGraph_) {
      indices[numIndices++] = record->getNumber();
    }
    else {
      indices[numIndices++] = eqnNumbers[0];
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::addBlockToGraph_noField_symmetric(fei::Graph* graph,
                                                        fei::ConnectivityBlock* cblock)
{
  fei::Pattern* pattern = cblock->getRowPattern();
  int numIDs = pattern->getNumIDs();
  int numIndices = pattern->getNumIndices();
  std::vector<int> indices(numIndices);
  int* indicesPtr = &indices[0];

  std::map<int,int>& connIDs = cblock->getConnectivityIDs();
  std::vector<int>& rowrecords = cblock->getRowConnectivities();

  std::map<int,int>::iterator
    c_iter = connIDs.begin(),
    c_end  = connIDs.end();

  for(; c_iter != c_end; ++c_iter) {
    int offset = c_iter->second;
    int* records = &rowrecords[offset*numIDs];

    int checkNumIndices;
    CHK_ERR( getConnectivityIndices_noField(pattern->getRecordCollections(), records, numIDs, numIndices,
                                            indicesPtr, checkNumIndices) );

    if (checkNumIndices != numIndices) {
      ERReturn(-1);
    }

    if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
      for(int ii=0; ii<numIndices; ++ii) {
        if (isLogEqn(indicesPtr[ii])) {
          FEI_OSTREAM& os = *output_stream_;
          os << "adding Symm inds: ";
          for(int jj=0; jj<numIndices; ++jj) {
            os << indicesPtr[jj] << " ";
          }
          os << FEI_ENDL;
          break;
        }
      }
    }

    //now we have the indices array, so we're ready to push it into
    //the graph container
    CHK_ERR( graph->addSymmetricIndices(numIndices, indicesPtr) );
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::MatrixGraph_Impl2::addBlockToGraph_sparse(fei::Graph* graph,
                                             fei::ConnectivityBlock* cblock)
{
  std::vector<int> row_indices;
  std::vector<int> indices;

  fei::Pattern* pattern = cblock->getRowPattern();
  const snl_fei::RecordCollection*const* recordCollections = pattern->getRecordCollections();

  std::map<int,int>& connIDs = cblock->getConnectivityIDs();
  std::vector<int>& connOffsets = cblock->getConnectivityOffsets();
  int* rowrecords = &(cblock->getRowConnectivities()[0]);
  int* colrecords = &(cblock->getColConnectivities()[0]);
  bool haveField = cblock->haveFieldID();
  int fieldID = cblock->fieldID();
  int fieldSize = 1;

  if (haveField) {
    fieldSize = snl_fei::getFieldSize(fieldID, rowSpace_.get(),
                                      colSpace_.get());
  }

  std::map<int,int>::iterator
    c_iter = connIDs.begin(),
    c_end  = connIDs.end();

  for(; c_iter != c_end; ++c_iter) {
    int offset = c_iter->second;
    int rowlen = connOffsets[offset+1] - offset;

    int* records = &(rowrecords[offset]);

    int checkNumIndices;
    row_indices.resize(fieldSize);

    if (haveField) {
      CHK_ERR( getConnectivityIndices_singleField(recordCollections, records, 1,
                                                  fieldID, fieldSize,
                                                  fieldSize,
                                                  &row_indices[0],
                                                  checkNumIndices) );
    }
    else {
      CHK_ERR( getConnectivityIndices_noField(recordCollections, records, 1, fieldSize,
                                              &row_indices[0],
                                              checkNumIndices) );
    }

    indices.resize(fieldSize*rowlen);
    int* indicesPtr = &indices[0];
    int* crecords = &(colrecords[offset]);

    if (haveField) {
      CHK_ERR( getConnectivityIndices_singleField(recordCollections, crecords, rowlen,
                                                  fieldID, fieldSize,
                                                  fieldSize*rowlen,
                                                  indicesPtr,
                                                  checkNumIndices) );
    }
    else {
      CHK_ERR( getConnectivityIndices_noField(recordCollections, crecords, rowlen, rowlen,
                                              indicesPtr, checkNumIndices) );
    }
    if (checkNumIndices != rowlen) {
      ERReturn(-1);
    }

    //now we have the indices, so we're ready to push them into
    //the graph container
    int* row_ind_ptr = &row_indices[0];
    for(int i=0; i<fieldSize; ++i) {
      CHK_ERR( graph->addIndices(row_ind_ptr[i], fieldSize*rowlen,
                                 indicesPtr) );
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
void fei::MatrixGraph_Impl2::setName(const char* name)
{
  if (name == NULL) return;

  if (name_ == name) return;

  name_ = name;
  dbgprefix_ = "MatGraph_"+name_+": ";
}

//----------------------------------------------------------------------------
void fei::MatrixGraph_Impl2::setIndicesMode(int mode)
{
  if (mode == BLOCK_ENTRY_GRAPH) {
    blockEntryGraph_ = true;
    return;
  }

  if (mode == POINT_ENTRY_GRAPH) {
    blockEntryGraph_ = false;
    return;
  }

  voidERReturn;
}

//----------------------------------------------------------------------------
fei::SharedPtr<fei::FillableMat>
fei::MatrixGraph_Impl2::getSlaveDependencyMatrix()
{
  return( D_ );
}


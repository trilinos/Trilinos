/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#include <fei_trilinos_macros.hpp>
#include <fei_sstream.hpp>

#ifdef HAVE_FEI_AZTECOO

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdexcept>

#include "fei_defs.h"
#include "fei_Data.hpp"
#include "fei_Lookup.hpp"
#include "fei_LinearSystemCore.hpp"
#include "snl_fei_Utils.hpp"
#include "fei_ArrayUtils.hpp"

#undef fei_file
#define fei_file "fei_Aztec_LinSysCore.cpp"
#include "fei_ErrMacros.hpp"

#include "fei_mpi.h"

#ifndef FEI_SER

#define AZTEC_MPI AZTEC_MPI
#define AZ_MPI AZ_MPI
#ifndef MPI
#define MPI MPI
#endif

#endif

#include <az_aztec.h>

#include "fei_Aztec_Map.hpp"
#include "fei_Aztec_BlockMap.hpp"
#include "fei_Aztec_LSVector.hpp"
#include "fei_AztecDMSR_Matrix.hpp"
#include "fei_AztecDVBR_Matrix.hpp"

#include "fei_Aztec_LinSysCore.hpp"

namespace fei_trilinos {

static int azlsc_solveCounter_ = 0;
static int azlsc_debugFileCounter_ = 0;

static std::map<std::string, unsigned> fei_aztec_named_solve_counter;

//=========CONSTRUCTOR==========================================================
Aztec_LinSysCore::Aztec_LinSysCore(MPI_Comm comm)
 : comm_(comm),
   lookup_(NULL),
   haveLookup_(false),
   update_(NULL),
   map_(),
   A_(NULL),
   A_ptr_(NULL),
   x_(NULL),
   b_(NULL),
   bc_(NULL),
   essBCindices_(NULL),
   numEssBCs_(0),
   bcsLoaded_(false),
   explicitDirichletBCs_(false),
   BCenforcement_no_column_mod_(false),
   b_ptr_(NULL),
   matrixAllocated_(false),
   vectorsAllocated_(false),
   blkMatrixAllocated_(false),
   matrixLoaded_(false),
   rhsLoaded_(false),
   needNewPreconditioner_(false),
   tooLateToChooseBlock_(false),
   blockMatrix_(false),
   blkMap_(),
   blkA_(NULL),
   blkA_ptr_(NULL),
   blkUpdate_(NULL),
   azA_(NULL),
   azP_(NULL),
   precondCreated_(false),
   azS_(NULL),
   scalingCreated_(false),
   aztec_options_(NULL),
   aztec_params_(NULL),
   aztec_status_(NULL),
   tmp_x_(NULL),
   tmp_x_touched_(false),
   tmp_b_(NULL),
   tmp_bc_(NULL),
   tmp_b_allocated_(false),
   ML_Vanek_(false),
   numLevels_(9), //arbitrary default...
   rhsIDs_(NULL),
   numRHSs_(0),
   currentRHS_(-1),
   numGlobalEqns_(0),
   localOffset_(0),
   numLocalEqns_(0),
   numGlobalEqnBlks_(0),
   localBlkOffset_(0),
   numLocalEqnBlks_(0),
   localBlockSizes_(NULL),
   outputLevel_(0),
   numParams_(0),
   paramStrings_(NULL),
   name_("dbg"),
   debugOutput_(0),
   debugFileCounter_(0),
   debugPath_(NULL),
   debugFileName_(NULL),
   debugFile_(NULL),
   named_solve_counter_(fei_aztec_named_solve_counter)
{
   masterProc_ = 0;
   numProcs_ = 1;
   thisProc_ = 0;
#ifndef FEI_SER
   MPI_Comm_size(comm_, &numProcs_);
   MPI_Comm_rank(comm_, &thisProc_);
#endif

   aztec_options_ = new int[AZ_OPTIONS_SIZE];
   aztec_params_ = new double[AZ_PARAMS_SIZE];

   AZ_defaults(aztec_options_, aztec_params_);

   aztec_status_ = new double[AZ_STATUS_SIZE];
   for(int i=0; i<AZ_STATUS_SIZE; i++) aztec_status_[i] = 0.0;

   numRHSs_ = 1;
   rhsIDs_ = new int[numRHSs_];
   rhsIDs_[0] = 0;

   std::map<std::string,unsigned>::iterator
     iter = named_solve_counter_.find(name_);
   if (iter == named_solve_counter_.end()) {
     named_solve_counter_.insert(std::make_pair(name_, 0));
   }
}

struct free_any_remaining_az_memory {
 free_any_remaining_az_memory(){}
 ~free_any_remaining_az_memory()
 {
   AZ_manage_memory(0, AZ_CLEAR_ALL, 0, NULL, NULL);
 }
};

//========DESTRUCTOR============================================================
Aztec_LinSysCore::~Aztec_LinSysCore() {
  static free_any_remaining_az_memory when_program_exits;

  if (blkMatrixAllocated_) {
    delete blkA_;
    delete [] blkUpdate_;
    delete [] localBlockSizes_;
    blkMatrixAllocated_ = false;
  }
  if (matrixAllocated_) {
    delete A_;
    matrixAllocated_ = false;
  }

  if (precondCreated_) {
    AZ_precond_destroy(&azP_);
    precondCreated_ = false;
  }

  if (scalingCreated_) {
    AZ_scaling_destroy(&azS_);
    scalingCreated_ = false;
  }

  if (vectorsAllocated_) delete x_;
  delete [] tmp_x_;
  delete [] tmp_bc_;

  for(int i=0; i<numRHSs_; i++) {
    if (vectorsAllocated_) delete b_[i];
    if (tmp_b_allocated_) delete [] tmp_b_[i];
  }
  if (vectorsAllocated_) delete [] b_;
  if (tmp_b_allocated_) delete [] tmp_b_;
  tmp_b_allocated_ = false;
  delete [] update_;

  delete [] aztec_options_;
  delete [] aztec_params_;
  delete [] aztec_status_;

  delete [] rhsIDs_;
  numRHSs_ = 0;

  for(int j=0; j<numParams_; j++) {
    delete [] paramStrings_[j];
  }
  delete [] paramStrings_;
  numParams_ = 0;

  if (debugOutput_) {
    debugOutput_ = 0;
    fclose(debugFile_);
    delete [] debugPath_;
    delete [] debugFileName_;
  }

  delete [] essBCindices_;
  essBCindices_ = NULL;
  numEssBCs_ = 0;
  delete bc_;
}

//==============================================================================
LinearSystemCore* Aztec_LinSysCore::clone() {
   return(new Aztec_LinSysCore(comm_));
}

//==============================================================================
int Aztec_LinSysCore::parameters(int numParams, const char*const * params) {
//
// this function takes parameters for setting internal things like solver
// and preconditioner choice, etc.
//
   debugOutput("parameters");

   if (numParams == 0 || params == NULL) {
      debugOutput("--- no parameters");
      return(0);
   }

   const char* param = NULL;

   snl_fei::mergeStringLists(paramStrings_, numParams_,
                                    params, numParams);

   param = snl_fei::getParamValue("outputLevel",numParams,params);
   if (param != NULL){
      sscanf(param,"%d", &outputLevel_);
   }

   param = snl_fei::getParamValue("AZ_matrix_type", numParams, params);
   if (param != NULL){
      setMatrixType(param);
   }

   param = snl_fei::getParam("EXPLICIT_BC_ENFORCEMENT", numParams, params);
   if (param != NULL){
     explicitDirichletBCs_ = true;
   }
   else {
     explicitDirichletBCs_ = false;
   }

  param = snl_fei::getParam("BC_ENFORCEMENT_NO_COLUMN_MOD",numParams,params);
  if (param != NULL){
    BCenforcement_no_column_mod_ = true;
  }
  else {
    BCenforcement_no_column_mod_ = false;
  }

   param = snl_fei::getParamValue("numLevels", numParams, params);
   if (param != NULL){
      sscanf(param,"%d", &numLevels_);
   }

   bool dbgOutputParam = false;
   char* dbgFileName = NULL;
   char* dbgPath = NULL;

   bool output_level_on = false;
   param = snl_fei::getParamValue("FEI_OUTPUT_LEVEL",numParams,params);
   if (param != NULL) {
     std::string str(param);
     if (str == "ALL" || str == "MATRIX_FILES" || str == "FULL_LOGS") {
       output_level_on = true;
     }
   }

   param = snl_fei::getParamValue("debugOutput",numParams,params);
   if (param != NULL || output_level_on){
      dbgOutputParam = true;
      dbgFileName = new char[128];
      dbgPath = param==NULL ? new char[3] : new char[strlen(param)+1];
      if (param == NULL) sprintf(dbgPath, "./");
      else strcpy(dbgPath, param);

      sprintf(dbgFileName, "AZLSC.%d.%d.%d",
              numProcs_, thisProc_, azlsc_debugFileCounter_);
   }

   param = snl_fei::getParamValue("name",numParams, params);
   if(param != NULL){
      name_ = param;

      std::map<std::string,unsigned>::iterator
        iter = named_solve_counter_.find(name_);
      if (iter == named_solve_counter_.end()) {
        named_solve_counter_.insert(std::make_pair(name_, 0));
      }

      if (dbgOutputParam) {
         if (dbgFileName != NULL) delete [] dbgFileName;
         dbgFileName = new char[256];
         sprintf(dbgFileName, "AZLSC_%s.%d.%d.%d", name_.c_str(),
                 numProcs_, thisProc_, azlsc_debugFileCounter_);
      }
   }

   if (dbgOutputParam) {
     if (azlsc_debugFileCounter_ == azlsc_solveCounter_) {
       setDebugOutput(dbgPath,dbgFileName);
       ++azlsc_debugFileCounter_;
     }
     delete [] dbgFileName;
     delete [] dbgPath;
   }

   if (debugOutput_) {
      fprintf(debugFile_,"--- numParams %d\n",numParams);
      for(int i=0; i<numParams; i++){
         fprintf(debugFile_,"------ paramStrings[%d]: %s\n",i,
                 paramStrings_[i]);
      }
   }

   debugOutput("leaving parameters function");
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::setLookup(Lookup& lookup)
{
   lookup_ = &lookup;
   haveLookup_ = true;
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::setGlobalOffsets(int len, int* nodeOffsets,
                                        int* eqnOffsets, int* blkEqnOffsets)
{
  localOffset_ = eqnOffsets[thisProc_];
  numLocalEqns_ = eqnOffsets[thisProc_+1]-localOffset_;
  numGlobalEqns_ = eqnOffsets[numProcs_];

  localBlkOffset_ = blkEqnOffsets[thisProc_];
  numLocalEqnBlks_ = blkEqnOffsets[thisProc_+1]-localBlkOffset_;
  numGlobalEqnBlks_ = blkEqnOffsets[numProcs_];

  int err = createMiscStuff();
  if (err) return(err);

  if (debugOutput_) {
    fprintf(debugFile_, "setGlobalNodeOffsets, len: %d\n", len);
    for(int i=0; i<len; i++) {
      fprintf(debugFile_, "   nodeOffsets[%d]: %d, eqnOffsets[%d]: %d, blkEqnOffsets[%d], %d\n", i, nodeOffsets[i], i, eqnOffsets[i], i, blkEqnOffsets[i]);
    }
    fflush(debugFile_);
  }
  return(0);
}

//==============================================================================
int Aztec_LinSysCore::setConnectivities(GlobalID elemBlock,
                                  int numElements,
                                  int numNodesPerElem,
                                  const GlobalID* elemIDs,
                                  const int* const* connNodes)
{
   (void) elemBlock;
   (void) numElements;
   (void) numNodesPerElem;
   (void) elemIDs;
   (void) connNodes;
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::setMatrixType(const char* name) {

   if (!strcmp(name, "AZ_VBR_MATRIX")) {
      if (!tooLateToChooseBlock_) {
         blockMatrix_ = true;
         debugOutput("setMatrixType: chose block matrix");
      }
      else {
         fei::console_out() << "Aztec_LinSysCore::setMatrixType: WARNING: Too late to choose"
           << " the DVBR matrix; make this choice before calling "
           << "setMatrixStructure. DMSR will be used." << FEI_ENDL;
      }
   }
   else if (!strcmp(name, "AZ_MSR_MATRIX")) {
      //do nothing, this is the default case.
   }
   else {
      if (thisProc_ == 0) {
         FEI_COUT << "Aztec_LinSysCore: Warning: requested matrix-type <"<<name <<"> not recognized." << FEI_ENDL;
         FEI_COUT << "Aztec_LinSysCore: valid matrix-types are: 'AZ_MSR_MATRIX', 'AZ_VBR_MATRIX'" << FEI_ENDL;
      }
   }
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::setMatrixStructure(int** ptColIndices,
                                          int* ptRowLengths,
                                          int** blkColIndices,
                                          int* blkRowLengths,
                                          int* ptRowsPerBlkRow)
{
  if (debugOutput_) {
    fprintf(debugFile_, "setMatrixStructure\n");
    for(int i=0; i<numLocalEqnBlks_; i++) {
      int blkEqn = localBlkOffset_+i;
      fprintf(debugFile_, "   blkRow %d, ptRowsPerBlkRow %d\n",
              blkEqn, ptRowsPerBlkRow[i]);
    }
    fflush(debugFile_);
  }

  int err = allocateMatrix(ptColIndices, ptRowLengths, blkColIndices,
                           blkRowLengths, ptRowsPerBlkRow);
  return(err);
}

//==============================================================================
int Aztec_LinSysCore::createMiscStuff()
{
//
//This function is where we establish the structures/objects associated
//with the linear algebra library. i.e., do initial allocations, etc.
//
   if (debugOutput_) {
      fprintf(debugFile_,
              "createMiscStuff: numRHSs_: %d\n", numRHSs_);
      fflush(debugFile_);
   }

   if (numLocalEqns_ > numGlobalEqns_)
      messageAbort("createMiscStuff: numLocalEqns > numGlobalEqns");

   if (0 > localOffset_)
      messageAbort("createMiscStuff: localOffset_ out of range");

   update_ = numLocalEqns_ > 0 ? new int[numLocalEqns_] : NULL;
   if (numLocalEqns_ > 0 && update_ == NULL) {
     return(-1);
   }

   int i, j;
   for(j=0; j<numLocalEqns_; j++) update_[j] = localOffset_+j;

   azS_ = AZ_scaling_create();
   if (azS_ != NULL) scalingCreated_ = true;
   else {
     fei::console_out() << "Aztec_LinSysCore::createMiscStuff ERROR: failed to create scaling"
          << FEI_ENDL;
     return(-1);
   }

   if (numRHSs_ <= 0)
      messageAbort("numRHSs_==0. Out of scope or destructor already called?");

    tmp_x_ = numLocalEqns_ > 0 ? new double[numLocalEqns_] : NULL;
    if (numLocalEqns_ > 0 && tmp_x_ == NULL) return(-1);

    tmp_bc_ = numLocalEqns_ > 0 ? new double[numLocalEqns_] : NULL;
    if (numLocalEqns_ > 0 && tmp_bc_ == NULL) return(-1);

    for(i=0; i<numLocalEqns_; i++){
        tmp_x_[i] = 0.0;
        tmp_bc_[i] = 0.0;
    }

    if (!tmp_b_allocated_) {
       tmp_b_ = new double*[numRHSs_];

       for(j=0; j<numRHSs_; j++) {
           tmp_b_[j] = new double[numLocalEqns_];
           for(i=0; i<numLocalEqns_; i++) {
               tmp_b_[j][i] = 0.0;
           }
       }

       tmp_b_allocated_ = true;
    }

   if (currentRHS_ < 0) currentRHS_ = 0;
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::allocateMatrix(int** ptColIndices,
                                      int* ptRowLengths,
                                      int** blkColIndices,
                                      int* blkRowLengths,
                                      int* ptRowsPerBlkRow)
{
  int err;
  if (blockMatrix_) {
    err = createBlockMatrix(blkColIndices, blkRowLengths, ptRowsPerBlkRow);
    return(err);
  }

  tooLateToChooseBlock_ = true;

  map_.reset( new Aztec_Map(numGlobalEqns_, numLocalEqns_, update_, localOffset_, comm_));

  A_ = new AztecDMSR_Matrix(map_);

  if (A_ptr_ == NULL) A_ptr_ = A_;

  //
  //we're going to have to look through the colIndices and see if any
  //of the rows do NOT have an entry on the diagonal. For each row that
  //does have an entry on the diagonal, we subtract 1 from the row-length
  //since the AztecDMSR_Matrix::allocate function wants the length of each
  //row NOT including any entry on the diagonal.
  //

  int* row_lengths = numLocalEqns_ > 0 ? new int[numLocalEqns_] : NULL;

  for (int i = 0; i < numLocalEqns_; i++) {
    if (fei::searchList(localOffset_+i,
                            ptColIndices[i], ptRowLengths[i]) >= 0) {
      row_lengths[i] = ptRowLengths[i] - 1;
    }
    else {
      row_lengths[i] = ptRowLengths[i];
    }
  }

  // so now we know all the row lengths, and can configure our matrix.

  A_ptr_->allocate(row_lengths, ptColIndices);

  delete [] row_lengths;

  matrixAllocated_ = true;
  return(0);
}

//==============================================================================
int Aztec_LinSysCore::createBlockMatrix(int** blkColIndices,
                                         int* blkRowLengths,
                                         int* ptRowsPerBlkRow)
{
  int i;

  blkUpdate_ = new int[numLocalEqnBlks_];

  for(i=0; i<numLocalEqnBlks_; i++) {
    blkUpdate_[i] = localBlkOffset_ + i;
  }

  blkMap_.reset(new Aztec_BlockMap(numGlobalEqns_, numLocalEqns_,
                               update_, localOffset_, comm_,
                               numGlobalEqnBlks_, numLocalEqnBlks_,
                               blkUpdate_,
                               localBlkOffset_, ptRowsPerBlkRow));

  blkA_ = new AztecDVBR_Matrix(blkMap_);

  if (blkA_ptr_ == NULL) blkA_ptr_ = blkA_;

  localBlockSizes_ = new int[numLocalEqnBlks_];

  //now we need to count how many non-zero blocks there are.
  numNonzeroBlocks_ = 0;
  for(i=0; i<numLocalEqnBlks_; i++) {

    numNonzeroBlocks_ += blkRowLengths[i];

    localBlockSizes_[i] = ptRowsPerBlkRow[i];
  }

  //and now we need to flatten the list of block-column-indices into a 1-D
  //list to give to the AztecDVBR_Matrix::allocate function.
  int* blk_col_inds = new int[numNonzeroBlocks_];

  int offset = 0;
  for(i=0; i<numLocalEqnBlks_; i++) {
    for(int j=0; j<blkRowLengths[i]; j++) {
      blk_col_inds[offset++] = blkColIndices[i][j];
    }
  }

  //finally we're ready to allocate our matrix.
  blkA_->allocate(blkRowLengths, blk_col_inds);

  delete [] blk_col_inds;

  blkMatrixAllocated_ = true;
  return(0);
}

//==============================================================================
int Aztec_LinSysCore::resetMatrixAndVector(double s)
{
   if (blkMatrixAllocated_) blkA_ptr_->put(s);
   if (matrixAllocated_) A_ptr_->put(s);

   if (rhsLoaded_) {
      for(int i=0; i<numRHSs_; i++){
         b_[i]->put(s);
      }
   }

   if (b_ptr_ != NULL) b_ptr_->put(s);

   if (bc_ != NULL) bc_->put(s);

   if (!tmp_b_allocated_) {
      for(int j=0; j<numRHSs_; j++) {
         tmp_b_[j] = new double[numLocalEqns_];
      }
   }
   matrixLoaded_ = false;
   rhsLoaded_ = false;
   bcsLoaded_ = false;

   delete [] tmp_bc_;
   tmp_bc_ = new double[numLocalEqns_];

   for(int ii=0; ii<numLocalEqns_; ii++) tmp_bc_[ii] = s;

   delete [] essBCindices_;
   essBCindices_ = NULL;
   numEssBCs_ = 0;

   for(int j=0; j<numRHSs_; j++) {
      for(int i=0; i<numLocalEqns_; i++) {
         tmp_b_[j][i] = s;
      }
   }
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::resetMatrix(double s) {
   if (blkMatrixAllocated_) blkA_ptr_->put(s);
   if (matrixAllocated_) A_ptr_->put(s);

   matrixLoaded_ = false;
   bcsLoaded_ = false;

   if (bc_ != NULL) bc_->put(s);

   delete [] tmp_bc_;
   tmp_bc_ = new double[numLocalEqns_];

   for(int ii=0; ii<numLocalEqns_; ii++) tmp_bc_[ii] = s;

   delete [] essBCindices_;
   essBCindices_ = NULL;
   numEssBCs_ = 0;

   return(0);
}

//==============================================================================
int Aztec_LinSysCore::resetRHSVector(double s) {

   if (rhsLoaded_) {
      for(int i=0; i<numRHSs_; i++){
         b_[i]->put(s);
      }
   }

   if (b_ptr_ != NULL) b_ptr_->put(s);

   if (!tmp_b_allocated_) {
      for(int j=0; j<numRHSs_; j++) {
         tmp_b_[j] = new double[numLocalEqns_];
      }
   }
   rhsLoaded_ = false;
   bcsLoaded_ = false;

   for(int j=0; j<numRHSs_; j++) {
     double* cur_b = tmp_b_[j];
     for(int i=0; i<numLocalEqns_; i++) {
       cur_b[i] = s;
     }
   }
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::sumIntoSystemMatrix(int numPtRows, const int* ptRows,
                            int numPtCols, const int* ptCols,
                            const double* const* values)
{
  matrixLoaded_ = false;

  return(sumIntoPointRow(numPtRows, ptRows, numPtCols, ptCols, values, false));
}

//==============================================================================
int Aztec_LinSysCore::sumIntoSystemMatrix(int numPtRows, const int* ptRows,
                                           int numPtCols, const int* ptCols,
                                           int numBlkRows, const int* blkRows,
                                           int numBlkCols, const int* blkCols,
                                           const double* const* values)
{
   matrixLoaded_ = false;

   if ((A_ptr_ == NULL) && (blkA_ptr_ == NULL))
      messageAbort("sumIntoSystemMatrix: matrix is NULL.");

   if (blockMatrix_) {
      return(sumIntoBlockRow(numBlkRows, blkRows, numBlkCols, blkCols,
                      values, numPtCols, false));
   }
   else {
      return( sumIntoPointRow(numPtRows, ptRows, numPtCols, ptCols, values,
                              false));
   }
}

//==============================================================================
int Aztec_LinSysCore::putIntoSystemMatrix(int numPtRows, const int* ptRows,
                                           int numPtCols, const int* ptCols,
                                           const double* const* values)
{
  matrixLoaded_ = false;

  return(sumIntoPointRow(numPtRows, ptRows, numPtCols, ptCols, values, true));
}

//==============================================================================
int Aztec_LinSysCore::getMatrixRowLength(int row, int& length)
{
  if (!blockMatrix_) {
    if (row >= localOffset_ && row < (localOffset_+numLocalEqns_)) {
      length = A_ptr_->rowLength(row);
      return(0);
    }
  }
  else {
    if (!haveLookup_) return(-1);
    int blkRow = lookup_->ptEqnToBlkEqn(row);
    if (blkRow < 0) return(-1);
    return( blkA_ptr_->getNumNonzerosPerRow(blkRow, length) );
  }

  return(-1);
}

//==============================================================================
int Aztec_LinSysCore::getMatrixRow(int row, double* coefs,
                                   int* indices,
                                   int len, int& rowLength)
{
  int err = getMatrixRowLength(row, rowLength);
  if (err != 0) return(-1);

  int* tmpIndices = indices;
  double* tmpCoefs = coefs;

  if (len < rowLength) {
    tmpIndices = new int[rowLength];
    tmpCoefs = new double[rowLength];
  }

  if (!blockMatrix_) {
    A_ptr_->getRow(row, rowLength, tmpCoefs, tmpIndices);
  }
  else {
    fei::console_out() << "Aztec_LinSysCore::getMatrixRow: not implemented." << FEI_ENDL;
    return(-1);
  }

  if (len < rowLength) {
    for(int i=0; i<len; i++) {
      indices[i] = tmpIndices[i];
      coefs[i] = tmpCoefs[i];
    }
    delete [] tmpIndices;
    delete [] tmpCoefs;
  }

  return(0);
}

//==============================================================================
int Aztec_LinSysCore::sumIntoBlockRow(int numBlkRows, const int* blkRows,
                                       int numBlkCols, const int* blkCols,
                                       const double* const* values,
                                       int numPtCols,
                                      bool overwriteInsteadOfAccumulate)
{
  int i;
  //first, let's figure out which of the incoming blkRows is the biggest --
  //i.e., which contains the most point-rows.
  int maxBlkSize = 0;

  for(i=0; i<numBlkRows; i++) {
    int thisSize = lookup_->getBlkEqnSize(blkRows[i]);

    if (maxBlkSize < thisSize) maxBlkSize = thisSize;
  }

  //now we can allocate an array to hold the values for passing each block
  //row into the block-matrix.
  double* coefs = new double[numPtCols*maxBlkSize];

  int rowOffset = 0;
  for(i=0; i<numBlkRows; i++) {
    //now, copy the values for this block row into the coefs array.
    copyBlockRow(i, blkRows, numBlkCols, blkCols, &(values[rowOffset]), coefs);

    //now shove it into the matrix.
    if (overwriteInsteadOfAccumulate) {
      int err = blkA_ptr_->putBlockRow(blkRows[i], coefs,
                                       (int*)blkCols, numBlkCols);
      if (err != 0) {
        fei::console_out() << thisProc_ << " DVBR::putBlockRow failed." << FEI_ENDL;
        return(err);
      }
    }
    else {
      int err = blkA_ptr_->sumIntoBlockRow(blkRows[i], coefs,
                                           (int*)blkCols, numBlkCols);
      if (err != 0) {
        fei::console_out() << thisProc_ << " DVBR::sumIntoBlockRow failed." << FEI_ENDL;
        return(err);
      }
    }

    rowOffset += lookup_->getBlkEqnSize(blkRows[i]);
  }

  delete [] coefs;
  return(0);
}

//==============================================================================
int Aztec_LinSysCore::copyBlockRow(int i, const int* blkRows,
                                    int numBlkCols, const int* blkCols,
                                    const double* const* values,
                                    double* coefs){

   int rowSize = 0, colSize = 0;
   //coefs needs to contain the entries for each block together, and those
   //entries need to be in column-major order
   int colOffset = 0;
   int coefOffset = 0;
   for(int b=0; b<numBlkCols; b++) {

      int err = blkA_ptr_->getBlockSize(blkRows[i], blkCols[b],
                                            rowSize, colSize);
      if (err) {
         fei::console_out() << "Aztec_LSC:copyBlockRow: ERROR in getBlockSize" << FEI_ENDL;
         return(-1);
      }

      for(int j=colOffset; j<colOffset+colSize; j++) {
         for(int r=0; r<rowSize; r++) {
            coefs[coefOffset++] = values[r][j];
         }
      }
      colOffset += colSize;
   }
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::getBlockSize(int blkInd) {

   int localBlkRow = blkInd - localBlkOffset_;
   if (localBlkRow >= 0 && localBlkRow < numLocalEqnBlks_) {
      //if this is a local blkInd, we have its size in an array.
      return(localBlockSizes_[localBlkRow]);
   }
   else {
      //else it ain't local, so we'll get its size from the matrix because
      //we know that the matrix obtained it when it allocated itself.

      int numRemoteBlocks = blkA_ptr_->getNumRemoteBlocks();
      int* remoteInds = blkA_ptr_->getRemoteBlockIndices();
      int* remoteSizes = blkA_ptr_->getRemoteBlockSizes();

      int ins;
      int index = fei::binarySearch(blkInd, remoteInds, numRemoteBlocks, ins);
      if (index < 0)
         messageAbort("getBlockSize: can't find blkInd.");

      return(remoteSizes[index]);
   }
}

//==============================================================================
int Aztec_LinSysCore::sumIntoPointRow(int numPtRows, const int* ptRows,
                                       int numPtCols, const int* ptColIndices,
                                       const double* const* values,
                                      bool overwriteInsteadOfAccumulate)
{
  int i, j;
  if (debugOutput_) {
    fprintf(debugFile_,"sumIntoPointRow, %d rows\n", numPtRows);
    for(i=0; i<numPtRows; ++i) {
      for(j=0; j<numPtCols; ++j) {
        fprintf(debugFile_,"  sipr row %d, col %d, value: %e\n", ptRows[i],
                ptColIndices[j], values[i][j]);
      }
    }
    fflush(debugFile_);
  }

  if (!blockMatrix_) {
    if (overwriteInsteadOfAccumulate) {
      for(i=0; i<numPtRows; ++i) {
        CHK_ERR( A_ptr_->putRow(ptRows[i], numPtCols, values[i], ptColIndices) );
      }
    }
    else {
      int err = A_ptr_->sumIntoRow(numPtRows, ptRows, numPtCols, ptColIndices,
                                   values);
      if (err != 0) {
        FEI_OSTRINGSTREAM osstr;
        osstr << "Aztec_LinSysCore::sumIntoPointRow ERROR calling A_ptr->sumIntoRow";
        throw std::runtime_error(osstr.str());
      }
    }

    return(0);
  }

  if (!haveLookup_) {
    messageAbort("sumIntoPointRow: need lookup object, don't have it.");
  }

  int* blkIntData = new int[numPtRows*2 + numPtCols*2];
  int* blkRows = blkIntData;
  int* blkRowOffsets = blkIntData+numPtRows;
  int* blkCols = blkIntData+2*numPtRows;
  int* blkColOffsets = blkIntData+2*numPtRows+numPtCols;;

  //now fill the blkRows and blkCols lists.

  for(i=0; i<numPtRows; i++) {
    blkRows[i] = lookup_->ptEqnToBlkEqn(ptRows[i]);
    blkRowOffsets[i] = lookup_->getOffsetIntoBlkEqn(blkRows[i], ptRows[i]);
  }

  for(i=0; i<numPtCols; i++) {
    blkCols[i] = lookup_->ptEqnToBlkEqn(ptColIndices[i]);
    if (blkCols[i] < 0) {
      fei::console_out() << thisProc_ << " lookup ptEqnToBlkEqn("<<ptColIndices[i]<<"): "
           << blkCols[i] << FEI_ENDL;
      messageAbort("negative blk-col");
    }

    blkColOffsets[i] = lookup_->getOffsetIntoBlkEqn(blkCols[i], ptColIndices[i]);
  }

  for(i=0; i<numPtRows; i++) {
    for(j=0; j<numPtCols; j++) {
      sumPointIntoBlockRow(blkRows[i], blkRowOffsets[i],
                           blkCols[j], blkColOffsets[j], values[i][j]);
    }
  }

  delete [] blkIntData;
  return(0);
}

//==============================================================================
int Aztec_LinSysCore::sumPointIntoBlockRow(int blkRow, int rowOffset,
                                            int blkCol, int colOffset,
                                            double value)
{
   int rowSize = lookup_->getBlkEqnSize(blkRow);
   int colSize = lookup_->getBlkEqnSize(blkCol);

   int len = rowSize*colSize;
   if (len <= 0) {
      fei::console_out() << thisProc_ << ", ALSC::spibr: blkRow: " << blkRow << ", blkCol: " << blkCol << ", rowSize: " << rowSize << ", colSize: " << colSize << FEI_ENDL;
      messageAbort("sumPointIntoBlockRow: len <= 0");
   }

   double* val = new double[rowSize*colSize];

   for(int i=0; i<len; i++) val[i] = 0.0;

   val[colOffset*rowSize + rowOffset] = value;

   int err = blkA_ptr_->sumIntoBlockRow(blkRow, val, &blkCol, 1);
   if (err != 0) {
     fei::console_out() << thisProc_ << " DVBR::sumIntoBlockRow failed" << FEI_ENDL;
     return(err);
   }

   delete [] val;
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::sumIntoRHSVector(int num,
                                       const double* values,
                                       const int* indices)
{
  //
  //This function scatters (accumulates) values into the linear-system's
  //currently selected RHS vector.
  //
  // num is how many values are being passed,
  // indices holds the global 'row-numbers' into which the values go,
  // and values holds the actual coefficients to be scattered.
  //
  rhsLoaded_ = false;
  double* cur_b = tmp_b_[currentRHS_];

  if (debugOutput_) {
    for(int i=0; i<num; ++i) {
      fprintf(debugFile_,"sumIntoRHS %d, %e\n", indices[i], values[i]);
      fflush(debugFile_);
    }
  }

  for(int i=0; i<num; i++){
    int localRow = indices[i] - localOffset_;
    if (localRow < 0 || localRow > numLocalEqns_) continue;

    cur_b[localRow] += values[i];
  }
  return(0);
}

//==============================================================================
int Aztec_LinSysCore::putIntoRHSVector(int num, const double* values,
                                       const int* indices)
{
//
//This function scatters (puts) values into the linear-system's
//currently selected RHS vector.
//
// num is how many values are being passed,
// indices holds the global 'row-numbers' into which the values go,
// and values holds the actual coefficients to be scattered.
//
   rhsLoaded_ = false;

   for(int i=0; i<num; i++){
     int localRow = indices[i] - localOffset_;
     if (localRow < 0 || localRow > numLocalEqns_) continue;

     if (debugOutput_) {
       fprintf(debugFile_,"putIntoRHS %d, %e\n", indices[i], values[i]);
       fflush(debugFile_);
     }

     tmp_b_[currentRHS_][localRow] = values[i];
   }
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::getFromRHSVector(int num, double* values,
                                       const int* indices)
{
  //
  //This function retrieves values from the linear-system's
  //currently selected RHS vector.
  //
  // num is how many values are being requested,
  // indices holds the global 'row-numbers' for which values are requested,
  // and values holds the actual coefficients to be returned.
  //
  //if non-local indices are supplied, the corresponding positions in the values
  //array are not referenced.
  //

  for(int i=0; i<num; i++){
    int localRow = indices[i] - localOffset_;
    if (localRow < 0 || localRow > numLocalEqns_) continue;

    values[i] = tmp_b_[currentRHS_][localRow];
  }
  return(0);
}

//==============================================================================
int Aztec_LinSysCore::matrixLoadComplete() {

   if (debugOutput_) {
      fprintf(debugFile_,"matrixLoadComplete\n");
      fflush(debugFile_);
   }

   if (matrixLoaded_ && rhsLoaded_) return(0);

   fei::SharedPtr<Aztec_Map> tmpMap;
   int* data_org = NULL;

   if (blockMatrix_) {
      tmpMap = blkMap_;
   }
   else {
      tmpMap = map_;
   }

   if (!matrixLoaded_) {
     if (blockMatrix_) {
       if (!blkA_ptr_->isLoaded()) blkA_ptr_->loadComplete();
     }
     else {
       if (!A_ptr_->isFilled()) {
         A_ptr_->fillComplete();
       }
     }

     matrixLoaded_ = true;

     needNewPreconditioner_ = true;
   }

   if (blockMatrix_) data_org = blkA_ptr_->getData_org();
   else data_org = A_ptr_->getAZ_MATRIX_PTR()->data_org;

   Aztec_LSVector* tmp = NULL;

   //if x_ is not null, then it probably has a previous solution in it that
   //we might as well keep for the next initial guess unless the user has
   //specifically set an initial guess.
   if (x_ != NULL) {
      tmp = new Aztec_LSVector(*x_);
      *tmp = *x_;
   }

   //if x_ hasn't been allocated yet, we better do that now.
   if (x_ == NULL) x_ = new Aztec_LSVector(tmpMap, data_org);
   if (bc_ == NULL) bc_ = new Aztec_LSVector(tmpMap, data_org);

   //if we did save off a copy of x_ above, let's put it back now.
   if (tmp != NULL) {
      *x_ = *tmp;
      delete tmp;
   }

   //if we've put boundary-condition data in tmp_bc_ then we better move it
   //into bc_ now.
   if (tmp_bc_ != NULL) {
     for(int j=0; j<numEssBCs_; j++) {
       int index = essBCindices_[j];
       (*bc_)[index] = tmp_bc_[index-localOffset_];
     }
   }

   if (tmp_x_touched_) {
      //and now we can fill x_ with the stuff we've been holding in
      //tmp_x_, if tmp_x_ has been touched... i.e., if the user loaded an
      //initial guess into it.
      for(int i=0; i<numLocalEqns_; i++){
         (*x_)[i+localOffset_] = tmp_x_[i];
      }
   }

   //now we're going to get the AZ_MATRIX ptr out of the AztecDMSR matrix
   //wrapper class.

   if (blockMatrix_) azA_ = blkA_ptr_->getAZ_MATRIX_Ptr();
   else azA_ = A_ptr_->getAZ_MATRIX_PTR();

   if (rhsLoaded_) return(0);

   if (b_ != NULL) {
      for(int i=0; i<numRHSs_; i++) delete b_[i];
      delete [] b_;
   }

   b_ = new Aztec_LSVector*[numRHSs_];
   for(int j=0; j<numRHSs_; j++) {
      b_[j] = new Aztec_LSVector(tmpMap, data_org);
      if (b_[j] == NULL) return(-1);

      //now fill b_[j] with the stuff we've been holding in tmp_b_[j].
      for(int i=0; i<numLocalEqns_; i++){
         (*(b_[j]))[i+localOffset_] = tmp_b_[j][i];
      }
   }

   b_ptr_ = b_[currentRHS_];
   vectorsAllocated_ = true;

   if (!bcsLoaded_) modifyRHSforBCs();
   bcsLoaded_ = false;

   rhsLoaded_ = true;

   if (debugOutput_) {
      fprintf(debugFile_,"leaving matrixLoadComplete\n");
      fflush(debugFile_);
   }
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::getBlkEqnsAndOffsets(int* ptEqns, int* blkEqns,
                                            int* blkOffsets, int numEqns)
{
   if (!haveLookup_) messageAbort("getBlkEqnsAndOffsets: don't have Lookup");

   for(int i=0; i<numEqns; i++) {
      blkEqns[i] = lookup_->ptEqnToBlkEqn(ptEqns[i]);
      if (blkEqns[i] < 0) {
        fei::console_out() << thisProc_ << "ptEqn: " << ptEqns[i] << ", blkEqn: " << blkEqns[i]
             << FEI_ENDL;
        messageAbort("getBlkEqnsAndOffsets: blk-eqn lookup failed.");
      }

      blkOffsets[i] = lookup_->getOffsetIntoBlkEqn(blkEqns[i], ptEqns[i]);
      if (blkOffsets[i] < 0) {
         messageAbort("getBlkEqnsAndOffsets: blk-offset lookup failed.");
      }
   }
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::enforceEssentialBC(int* globalEqn,
                                         double* alpha,
                                         double* gamma, int len)
{
  //
  //This function must enforce an essential boundary condition on each local
  //equation in 'globalEqn'. This means that the following modification
  //should be made to A and b, for each globalEqn:
  //
  //for(each locally-owned equation i in the globalEqn array){
  //   zero row i and put 1.0 on the diagonal
  //   for(each row j in column i) {
  //      if (i==j) b[i] = gamma/alpha;
  //      else b[j] -= (gamma/alpha)*A[j,i];
  //      A[j,i] = 0.0;
  //   }
  //}
  //

  if (len == 0) return(0);

  std::vector<int> bcEqns; bcEqns.reserve(len);
  std::vector<int> indirect; indirect.reserve(len);
  for(int i=0; i<len; ++i) {
    bcEqns.push_back(globalEqn[i]);
    indirect.push_back(i);
  }

  fei::insertion_sort_with_companions<int>(len, &bcEqns[0], &indirect[0]);

  if (debugOutput_) {
    fprintf(debugFile_,"numEssBCs: %d\n", len);
    for(int ii=0; ii<len; ++ii) {
      fprintf(debugFile_, "   EssBC eqn %d, alpha %e gamma %e\n",
              bcEqns[ii], alpha[indirect[ii]], gamma[indirect[ii]]);
    }
    fflush(debugFile_);
  }

  bcsLoaded_ = true;

  int localEnd = localOffset_ + numLocalEqns_ - 1;

  if (debugOutput_) {
    fprintf(debugFile_,"localOffset_: %d, localEnd: %d\n", localOffset_, localEnd);
    fflush(debugFile_);
  }

  {
    int* newBCindices = new int[numEssBCs_+len];
    int ii;
    for(ii=0; ii<numEssBCs_; ii++) newBCindices[ii] = essBCindices_[ii];
    int offset = 0;
    for(ii=0; ii<len; ii++) {
      if ((localOffset_ <= globalEqn[ii]) && (globalEqn[ii] <= localEnd)){
        newBCindices[numEssBCs_+offset++] = globalEqn[ii];
      }
    }

    if (offset > 0) {
      delete [] essBCindices_;
      essBCindices_ = newBCindices;
      numEssBCs_ += offset;
    }
    else {
      delete [] newBCindices;
    }
  }

  if (blockMatrix_) {
    int* blkIntData = new int[len*2];
    int* blkEqns = blkIntData;
    int* blkOffsets = blkIntData+len;

    getBlkEqnsAndOffsets(globalEqn, blkEqns, blkOffsets, len);

    enforceBlkEssentialBC(blkEqns, blkOffsets, alpha, gamma, len);

    delete [] blkIntData;

    return(0);
  }

  for(int i=0; i<len; i++) {

    int globalEqn_i = bcEqns[i];

    //if globalEqn[i] isn't local, then the processor that owns it
    //should be running this code too. Otherwise there's trouble...

    if ((localOffset_ > globalEqn_i) || (globalEqn_i > localEnd)) continue;

    //zero this row, except for the diagonal coefficient.
    A_ptr_->setDiagEntry(globalEqn_i, 1.0);
    int* offDiagIndices = NULL;
    double* offDiagCoefs = NULL;
    int offDiagLength = 0;
    A_ptr_->getOffDiagRowPointers(globalEqn_i, offDiagIndices,
                                  offDiagCoefs, offDiagLength);

    for(int jjj=0; jjj<offDiagLength; jjj++) offDiagCoefs[jjj] = 0.0;

    //also, make the rhs modification here.
    double bc_term = gamma[indirect[i]]/alpha[indirect[i]];
    double rhs_term = bc_term;
    if (explicitDirichletBCs_) rhs_term = 0.0;

    if (rhsLoaded_) {
      (*b_ptr_)[globalEqn_i] = rhs_term;
      (*bc_)[globalEqn_i] = bc_term;
    }
    else {
      tmp_b_[currentRHS_][globalEqn_i -localOffset_] = rhs_term;
      tmp_bc_[globalEqn_i -localOffset_] = bc_term;
    }
  }

  if (BCenforcement_no_column_mod_ == true) {
    return(0);
  }

  for(int row=localOffset_; row<=localEnd; row++) {

    int insertPoint = -1;
    int index = fei::binarySearch(row, &bcEqns[0], len, insertPoint);
    if (index >= 0) continue;

    int* offDiagIndices2 = NULL;
    double* offDiagCoefs2 = NULL;
    int offDiagLength2 = 0;
    A_ptr_->getOffDiagRowPointers(row, offDiagIndices2,
                                  offDiagCoefs2, offDiagLength2);

    //look through this row to find the non-zeros in positions that
    //correspond to eqns in bcEqns and make the appropriate modification.

    for(int j=0; j<offDiagLength2; j++) {

      int col_index = A_ptr_->getAztec_Map()->getTransformedEqn(offDiagIndices2[j]);

      int idx = fei::binarySearch(col_index, &bcEqns[0], len, insertPoint);
      if (idx < 0) continue;

      double bc_term = gamma[indirect[idx]]/alpha[indirect[idx]];

      double value = offDiagCoefs2[j]*bc_term;

      if (rhsLoaded_) {
        (*b_ptr_)[row] -= value;
        (*bc_)[row] -= value;
      }
      else {
        tmp_b_[currentRHS_][row-localOffset_] -= value;
        tmp_bc_[row-localOffset_] -= value;
      }

      if (debugOutput_) {
        fprintf(debugFile_,"BC mod, rhs %d  -= %e\n", row, value);
        fprintf(debugFile_,"BC, set A(%d,%d)==%e, to 0.0\n",
                row, bcEqns[idx], offDiagCoefs2[j]);
      }

      offDiagCoefs2[j] = 0.0;

    }//for offDiagLength2

  }//for localEnd

  return(0);
}

//==============================================================================
int Aztec_LinSysCore::enforceBlkEssentialBC(int* blkEqn,
                                            int* blkOffset,
                                            double* alpha,
                                            double* gamma,
                                            int len)
{
//
//This function must enforce an essential boundary condition on each local
//equation specified by the pair 'blkEqn' and 'blkOffset'. A point-equation
//resides within a block-equation. The blkOffset gives the distance from the
//beginning of the block-equation to the point-equation.

//The following modification should be made to A and b, for each incoming
//point-equation:
//
//for(each local equation i){
//   for(each column j in row i) {
//      if (i==j) b[i] = gamma/alpha;
//      else b[j] -= (gamma/alpha)*A[j,i];
//   }
//}
//
//all of row and column 'i' in A should be zeroed,
//except for 1.0 on the diagonal.
//
   int val_length = 0;
   double* val = NULL;
   int colInd_length = 0;
   int* blkCols = NULL;
   int rowNNZs = 0, numCols = 0, err = 0;

   double* val2 = NULL;
   int val2Len = val_length;
   int* cols2 = NULL;
   int cols2Len = colInd_length;

   int localEnd = localBlkOffset_ + numLocalEqnBlks_ - 1;

   for(int i=0; i<len; i++) {

      //if blkEqn[i] isn't local, then the processor that owns it
      //should be running this code too. Otherwise there's trouble...

      if ((localBlkOffset_ > blkEqn[i]) || (blkEqn[i] > localEnd)){
         continue;
      }

      err = getBlockRow(blkEqn[i], val, val_length, blkCols, colInd_length,
                        numCols, rowNNZs);
      if (err) {
         fei::console_out() << "Aztec_LSC: ERROR in getBlockRow" << FEI_ENDL;
         return(-1);
      }

      //now let's do the BC modification to this block-row.

      err = blkRowEssBCMod(blkEqn[i], blkOffset[i], val, blkCols, numCols,
                             rowNNZs, alpha[i], gamma[i]);

      blkA_ptr_->putBlockRow(blkEqn[i], val, blkCols, numCols);

      //now let's take advantage of the symmetry of element-contributions to
      //do the column-wise part of the BC modification. Since the structure of
      //the matrix arising from element contributions is symmetric, we know that
      //the column indices in row 'blkEqn' correspond to rows which have entries
      //in column 'blkEqn'. So we can modify the column in those rows rather
      //than searching all rows in the matrix looking for that column.

      //so let's loop over the block-columns and do the column-wise essBC mod
      //to those rows.
      
      for(int j=0; j<numCols; j++) {

         int col_row = blkCols[j];

         //if this column-index doesn't correspond to a local row, skip it.
         if ((localOffset_ > col_row) || (col_row > localEnd)) continue;

         //if this column is the diagonal column, skip it.
         if (col_row == blkEqn[i]) continue;

         int thisNNZ = 0;
         int thisNumBlks = 0;
         err = getBlockRow(col_row, val2, val2Len, cols2, cols2Len,
                           thisNumBlks, thisNNZ);

         err = blkColEssBCMod(col_row, blkEqn[i], blkOffset[i], val2, cols2,
                              thisNumBlks, thisNNZ, alpha[i], gamma[i]);

         blkA_ptr_->putBlockRow(col_row, val2, cols2, thisNumBlks);

      }// end for(j<rowLength) loop
   }

   delete [] val;
   delete [] blkCols;
   delete [] val2;
   delete [] cols2;
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::blkRowEssBCMod(int blkEqn, int blkOffset, double* val,
                                       int* blkCols, int numCols, int numPtNNZ,
                                       double alpha, double gamma)
{
//
//This function performs an essential boundary-condition modification for a
//single block-row.
//
   //we need to know which point-row this block-row corresponds to, so
   //we can do stuff to the rhs vector.
         
   int pointRow = blockRowToPointRow(blkEqn);

   int offset = 0;

   for(int j=0; j<numCols; j++) {

      int err, ptRows = 0, ptCols = 0;
      err = blkA_ptr_->getBlockSize(blkEqn, blkCols[j], ptRows, ptCols);

      if (err) {
         fei::console_out() << "Aztec_LSC::blkRowEssBCMod: error in getBlockSize" << FEI_ENDL;
         return(1);
      }

      if (blkCols[j] == blkEqn) {
         //this is the diagonal block, so we need to diagonalize the
         //'blkOffset'th point-row and point-column, leaving a 1.0 on the
         //diagonal.
         double bc_term = gamma/alpha;
         double rhs_term = bc_term;
         if (explicitDirichletBCs_) rhs_term = 0.0;

         int thisOffset = offset;

         for(int jj=0; jj<ptCols; jj++) {
            if (jj==blkOffset) {
               //if this is the point-column of interest, we move the 
               //contents over to the rhs in the appropriate way, except
               //for the entry on the row-of-interest, which we just
               //set to 1.0.

               for(int row=0; row<ptRows; row++) {
                  if (row==blkOffset) {
                     if (rhsLoaded_) {
                        (*b_ptr_)[pointRow+row] = rhs_term;
                        (*bc_)[pointRow+row] = bc_term;
                     }
                     else {
                        tmp_b_[currentRHS_][pointRow+row-localOffset_] = rhs_term;
                        tmp_bc_[pointRow+row-localOffset_] = bc_term;
                     }
                     val[thisOffset+row] = 1.0;
                  }
                  else {
                     if (rhsLoaded_) {
                        (*b_ptr_)[pointRow+row] -= val[thisOffset+row]
                                                * bc_term;
                        (*bc_)[pointRow+row] -= val[thisOffset+row]
                                                * bc_term;
                     }
                     else {
                        tmp_b_[currentRHS_][pointRow+row-localOffset_] -= val[thisOffset+row]*bc_term;
                        tmp_bc_[pointRow+row-localOffset_] -= val[thisOffset+row]*bc_term;
                     }
                     val[thisOffset+row] = 0.0;
                  }
               }
            }
            else {
               val[thisOffset+blkOffset] = 0.0;
            }

            thisOffset += ptRows;
         }
      }
      else {
         //this isn't the diagonal block, so we just want to zero the
         //whole 'blkOffset'th point-row in this block.
         int thisOffset = offset + blkOffset;
         for(int ii=0; ii<ptCols; ii++) {
            val[thisOffset] = 0.0;
            thisOffset += ptRows;
         }
      }

      offset += ptRows*ptCols;
   }

   return(0);
}

//==============================================================================
int Aztec_LinSysCore::blkColEssBCMod(int blkRow, int blkEqn, int blkOffset,
                                     double* val, int* blkCols, int numCols,
                                     int numPtNNZ, double alpha, double gamma)
{
//
//This function does the column-wise modification for an Essential BC, to 
//block column 'blkEqn', to the point-equation 'blkOffset' equations into the
//block, for the block-row contained in val,blkCols.
//
//NOTE: blkEqn is a 0-based equation number.

   int thisPtRow = blockRowToPointRow(blkRow);
   int thisRowSize = 0, thisColSize = 0;

   //look through the row to find the non-zero blk in position
   //blkEqn and make the appropriate modification.
   int err, offset = 0;
   for(int j=0; j<numCols; j++) {
      err = blkA_ptr_->getBlockSize(blkRow, blkCols[j],
                                    thisRowSize, thisColSize);

      if (err) {
         fei::console_out() << "Aztec_LinSysCore::blkColEssBCMod: ERROR in getBlockSize" << FEI_ENDL;
         return(1);
      }

      if (blkCols[j] == blkEqn) {
         double bc_term = gamma/alpha;

         int thisOffset = offset + blkOffset*thisRowSize;

         for(int row=0; row<thisRowSize; row++) {
            if (rhsLoaded_) {
               (*b_ptr_)[thisPtRow+row] -= val[thisOffset+row] * bc_term;
               (*bc_)[thisPtRow+row] -= val[thisOffset+row] * bc_term;
            }
            else {
               tmp_b_[currentRHS_][thisPtRow+row-localOffset_] -=
                                         val[thisOffset+row]*bc_term;
               tmp_bc_[thisPtRow+row-localOffset_] -=
                                         val[thisOffset+row]*bc_term;
            }
            val[thisOffset+row] = 0.0;
         }

         break;
      }

      offset += thisRowSize*thisColSize;
   }// end for(j<numCols) loop

   return(0);
}

//==============================================================================
int Aztec_LinSysCore::blockRowToPointRow(int blkRow) {
//
//This function returns a (global 0-based) point-equation which corresponds to
//the first point-equation in block-row 'blkRow'.
//
   int localBlkRow = blkRow - localBlkOffset_;

   if (localBlkRow < 0 || localBlkRow > numLocalEqnBlks_) return(-1);
 
   int localPtRow = 0;
   for(int i=0; i<localBlkRow; i++) {
      localPtRow += localBlockSizes_[i];
   }

   int pointRow = localPtRow + localOffset_;
   return(pointRow);
}

//==============================================================================
int Aztec_LinSysCore::getBlockRow(int blkRow, double*& val, int& valLen,
                                  int*& blkColInds, int& blkColIndLen,
                                  int& numNzBlks, int& numNNZ) {
//
//This function gets a block row from the VBR matrix. The val array and the
//blkColInds array are of lengths valLen and blkColIndLen, respectively.
//On output, val and blkColInds are lengthened if necessary, with the new
//lengths updated in valLen and blkColIndLen. The actual number of nonzero
//blocks and nonzero point-entries are returned in numNzBlks and numNNZ.
//
   numNNZ = 0;
   int err = blkA_ptr_->getNumNonzerosPerRow(blkRow, numNNZ);
   if (err) {
      fei::console_out() << "Aztec_LSC::getBlockRow: ERROR in getNumNonzerosPerRow" << FEI_ENDL;
      return(1);
   }

   numNzBlks = 0;
   err = blkA_ptr_->getNumBlocksPerRow(blkRow, numNzBlks);
   if (err) {
      fei::console_out() << "Aztec_LSC::getBlockRow: ERROR in getNumBlocksPerRow" << FEI_ENDL;
      return(1);
   }

   if (numNNZ > valLen) {
      double* newVals = new double[numNNZ];
      delete [] val;
      val = newVals;
      valLen = numNNZ;
   }

   if (numNzBlks > blkColIndLen) {
      int* newCols = new int[numNzBlks];
      delete [] blkColInds;
      blkColInds = newCols;
      blkColIndLen = numNzBlks;
   }

   err = blkA_ptr_->getBlockRow(blkRow, val, blkColInds, numNzBlks);

   return(err);
}

//==============================================================================
int Aztec_LinSysCore::enforceRemoteEssBCs(int numEqns, int* globalEqns,
                                          int** colIndices, int* colIndLen,
                                          double** coefs)
{
  bcsLoaded_ = true;

  //writeA("preRemBC");

  if (debugOutput_) {
    for(int i=0; i<numEqns; ++i) {
      fprintf(debugFile_,"remBC row %d, (cols,coefs): ", globalEqns[i]);
      for(int j=0; j<colIndLen[i]; ++j) {
        fprintf(debugFile_, "(%d,%e) ",colIndices[i][j], coefs[i][j]);
      }
      fprintf(debugFile_,"\n");
    }
    fflush(debugFile_);
  }

   if (blockMatrix_) {
      int* blkEqns = new int[numEqns];
      int* blkOffsets = new int[numEqns];

      getBlkEqnsAndOffsets(globalEqns, blkEqns, blkOffsets, numEqns);

      int** blkCols = new int*[numEqns];
      int** blkColOffsets = new int*[numEqns];

      for(int i=0; i<numEqns; i++) {
         blkCols[i] = new int[colIndLen[i]];
         blkColOffsets[i] = new int[colIndLen[i]];

         getBlkEqnsAndOffsets(colIndices[i], blkCols[i], blkColOffsets[i],
                              colIndLen[i]);
      }

      enforceBlkRemoteEssBCs(numEqns, blkEqns, blkCols, blkColOffsets,
                             colIndLen, coefs);

      delete [] blkEqns;
      delete [] blkOffsets;

      for(int j=0; j<numEqns; j++) {
         delete [] blkCols[j];
         delete [] blkColOffsets[j];
      }
      delete [] blkCols;
      delete [] blkColOffsets;

      return(0);
   }

   int localEnd = localOffset_ + numLocalEqns_ - 1;

   for(int i=0; i<numEqns; i++) {
     int globalEqn_i = globalEqns[i];

     if ((localOffset_ > globalEqn_i) || (globalEqn_i > localEnd)){
       continue;
     }

     int rowLen = 0;
     int* AcolInds = NULL;
     double* Acoefs = NULL;

     A_ptr_->getOffDiagRowPointers(globalEqn_i, AcolInds, Acoefs, rowLen);

     for(int j=0; j<colIndLen[i]; j++) {
       for(int k=0; k<rowLen; k++) {
         if (A_ptr_->getAztec_Map()->getTransformedEqn(AcolInds[k]) == colIndices[i][j]) {
           double value = Acoefs[k]*coefs[i][j];

           double old_rhs_val = 0.0;
           if (rhsLoaded_) {
             old_rhs_val = (*b_ptr_)[globalEqn_i];
             (*b_ptr_)[globalEqn_i] -= value;
             (*bc_)[globalEqn_i] -= value;
           }
           else {
             old_rhs_val = tmp_b_[currentRHS_][globalEqn_i -localOffset_];
             tmp_b_[currentRHS_][globalEqn_i -localOffset_] -= value;
             tmp_bc_[globalEqn_i -localOffset_] -= value;
           }

           if (debugOutput_) {
             fprintf(debugFile_,"remBC mod, rhs %d (%e) -= %e\n",
                     globalEqn_i, old_rhs_val, value);
             fprintf(debugFile_,"remBC, set A(%d,%d)==%e, to 0.0\n",
                     globalEqn_i, AcolInds[k], Acoefs[k]);
           }

           Acoefs[k] = 0.0;
         }
       }
     }

   }

   return(0);
}

//==============================================================================
int Aztec_LinSysCore::enforceBlkRemoteEssBCs(int numEqns, int* blkEqns,
                               int** blkColInds, int** blkColOffsets,
                               int* blkColLens, double** remEssBCCoefs) {
   int val_length = 0;
   double* val = NULL;
   int colInd_length = 0;
   int* blkCols = NULL;
   int rowNNZs = 0, numCols = 0, err = 0;

   int localEnd = localBlkOffset_ + numLocalEqnBlks_ - 1;

   for(int i=0; i<numEqns; i++) {
      if ((localBlkOffset_ > blkEqns[i]) || (blkEqns[i] > localEnd)){
         continue;
      }

      err = getBlockRow(blkEqns[i], val, val_length, blkCols, colInd_length,
                        numCols, rowNNZs);
      if (err) {
         fei::console_out() << "Aztec_LSC:enforceBlkRemoteEssBC ERROR in getBlockRow" << FEI_ENDL;
         return(-1);
      }

      //we need to know which point-row this block-row corresponds to, so
      //we can do stuff to the rhs vector.

      int pointRow = blockRowToPointRow(blkEqns[i]);

      int offset = 0;

      for(int j=0; j<numCols; j++) {
         int ptRows = 0, ptCols = 0;
         err = blkA_ptr_->getBlockSize(blkEqns[i], blkCols[j],
                                       ptRows, ptCols);
         if (err) {
            fei::console_out() << "Aztec_LSC::enforceBlkRemoteEssBCs: error in getBlockSize"
                << FEI_ENDL;
            return(-1);
         }

         for(int k=0; k<blkColLens[i]; k++) {
            if (blkColInds[i][k] == blkCols[j]) {
               int thisOffset = offset + blkColOffsets[i][k] * ptRows;
               double rhsTerm = remEssBCCoefs[i][k];

               double* bvec = &(tmp_b_[currentRHS_][pointRow-localOffset_]);
               double* bcvec = &(tmp_bc_[pointRow-localOffset_]);
               for(int row=0; row<ptRows; row++) {
                  double& coef = val[thisOffset+row];
                  bvec[row] -= coef*rhsTerm;
                  bcvec[row] -= coef*rhsTerm;
                  coef = 0.0;
               }

               blkA_ptr_->putBlockRow(blkEqns[i], &val[offset],
                                      &blkCols[j], 1);
            }
         }

         offset += ptRows*ptCols;
      }
   }

   delete [] val;
   delete [] blkCols;
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::getMatrixPtr(Data& data)
{
   if (!matrixLoaded_) matrixLoadComplete();

   if (blockMatrix_) {
      data.setTypeName("AztecDVBR_Matrix");
      data.setDataPtr((void*)blkA_ptr_);
   }
   else {
      data.setTypeName("AztecDMSR_Matrix");
      data.setDataPtr((void*)A_ptr_);
   }
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::copyInMatrix(double scalar, const Data& data) {
//
//Overwrites the current internal matrix with a scaled copy of the
//input argument.
//
   if (blockMatrix_) {
      if (strcmp("AztecDVBR_Matrix", data.getTypeName()))
         messageAbort("copyInMatrix: data type string not 'AztecDVBR_Matrix'.");

      AztecDVBR_Matrix* source = (AztecDVBR_Matrix*)data.getDataPtr();

      if (blkA_ != NULL) delete blkA_;
      blkA_ = new AztecDVBR_Matrix(*source);
      blkA_ptr_ = blkA_;

      VBRmatPlusScaledMat(blkA_ptr_, scalar, source);

      azA_ = blkA_ptr_->getAZ_MATRIX_Ptr();
   }
   else {
      if (strcmp("AztecDMSR_Matrix", data.getTypeName()))
         messageAbort("copyInMatrix: data type string not 'AztecDMSR_Matrix'.");

      AztecDMSR_Matrix* source = (AztecDMSR_Matrix*)data.getDataPtr();
      A_ptr_->copyStructure(*source);

      MSRmatPlusScaledMat(A_ptr_, scalar, source);

      azA_ = A_ptr_->getAZ_MATRIX_PTR();
   }

   needNewPreconditioner_ = true;
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::copyOutMatrix(double scalar, Data& data) {
//
//Passes out a scaled copy of the current internal matrix.
//

   if (!matrixLoaded_) matrixLoadComplete();

   if (blockMatrix_) {
      AztecDVBR_Matrix* outmat = new AztecDVBR_Matrix(*blkA_ptr_);

      //the copy-constructor above took all structural info from blkA_ptr_,
      //but not the data coefficients.

      VBRmatPlusScaledMat(outmat, scalar, blkA_ptr_);

      data.setTypeName("AztecDVBR_Matrix");
      data.setDataPtr((void*)outmat);
   }
   else {
      AztecDMSR_Matrix* outmat = new AztecDMSR_Matrix(*A_ptr_);

      outmat->scale(scalar);

      data.setTypeName("AztecDMSR_Matrix");
      data.setDataPtr((void*)outmat);
   }
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::sumInMatrix(double scalar, const Data& data) {

   if (!matrixLoaded_) matrixLoadComplete();

   if (blockMatrix_) {
     if (strcmp("AztecDVBR_Matrix", data.getTypeName())) {
       fei::console_out() << "Aztec_LinSysCore::sumInMatrix ERROR, incoming type-string: "
            << data.getTypeName() << ", expected AztecDVBR_Matrix." << FEI_ENDL;
       messageAbort("Aborting.");
     }
     AztecDVBR_Matrix* source = (AztecDVBR_Matrix*)data.getDataPtr();

     VBRmatPlusScaledMat(blkA_ptr_, scalar, source);
   }
   else {
      if (strcmp("AztecDMSR_Matrix", data.getTypeName()))
         messageAbort("sumInMatrix: data type string not 'AztecDMSR_Matrix'.");

      AztecDMSR_Matrix* source = (AztecDMSR_Matrix*)data.getDataPtr();

      MSRmatPlusScaledMat(A_ptr_, scalar, source);
   }

   needNewPreconditioner_ = true;
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::getRHSVectorPtr(Data& data) {

   if (!matrixLoaded_) matrixLoadComplete();

   data.setTypeName("Aztec_LSVector");
   data.setDataPtr((void*)b_ptr_);
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::copyInRHSVector(double scalar, const Data& data) {

   if (!rhsLoaded_) matrixLoadComplete();

   if (strcmp("Aztec_LSVector", data.getTypeName()))
      messageAbort("copyInRHSVector: data's type string not 'Aztec_LSVector'.");

   Aztec_LSVector* sourcevec = (Aztec_LSVector*)data.getDataPtr();

   *b_ptr_ = *sourcevec;

   if (scalar != 1.0) b_ptr_->scale(scalar);
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::copyOutRHSVector(double scalar, Data& data) {

   if (!rhsLoaded_) matrixLoadComplete();

   Aztec_LSVector* outvec = new Aztec_LSVector(*b_ptr_);

   outvec->put(0.0);

   outvec->addVec(scalar, *b_ptr_);

   data.setTypeName("Aztec_LSVector");
   data.setDataPtr((void*)outvec);
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::sumInRHSVector(double scalar, const Data& data) {

   if (!rhsLoaded_) matrixLoadComplete();

   if (strcmp("Aztec_LSVector", data.getTypeName()))
      messageAbort("sumInRHSVector: data's type string not 'Aztec_LSVector'.");

   Aztec_LSVector* source = (Aztec_LSVector*)data.getDataPtr();

   b_ptr_->addVec(scalar, *source);
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::destroyMatrixData(Data& data) {

   if (blockMatrix_) {
      if (strcmp("AztecDVBR_Matrix", data.getTypeName()))
         messageAbort("destroyMatrixData: data doesn't contain a AztecDVBR_Matrix.");

      AztecDVBR_Matrix* mat = (AztecDVBR_Matrix*)data.getDataPtr();
      delete mat;
   }
   else {
      if (strcmp("AztecDMSR_Matrix", data.getTypeName()))
         messageAbort("destroyMatrixData: data doesn't contain a AztecDMSR_Matrix.");

      AztecDMSR_Matrix* mat = (AztecDMSR_Matrix*)data.getDataPtr();
      delete mat;
   }
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::destroyVectorData(Data& data) {

   if (strcmp("Aztec_LSVector", data.getTypeName()))
      messageAbort("destroyVectorData: data doesn't contain a Aztec_LSVector.");

   Aztec_LSVector* vec = (Aztec_LSVector*)data.getDataPtr();
   delete vec;
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::setNumRHSVectors(int numRHSs, const int* rhsIDs) {

   if (numRHSs < 0)
      messageAbort("setNumRHSVectors: numRHSs < 0.");

   if (numRHSs == 0) return(0);

   delete [] rhsIDs_;
   numRHSs_ = numRHSs;
   rhsIDs_ = new int[numRHSs_];

   for(int i=0; i<numRHSs_; i++) rhsIDs_[i] = rhsIDs[i];
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::setRHSID(int rhsID) {

   for(int i=0; i<numRHSs_; i++){
      if (rhsIDs_[i] == rhsID){
         currentRHS_ = i;
         if (rhsLoaded_) b_ptr_ = b_[currentRHS_];
         return(0);
      }
   }

   messageAbort("setRHSID: rhsID not found.");
   return(-1);
}

//==============================================================================
int Aztec_LinSysCore::putInitialGuess(const int* eqnNumbers,
                                      const double* values,
                                      int len) {
//
//This function scatters (puts) values into the linear-system's soln vector.
//
// num is how many values are being passed,
// indices holds the global 'row-numbers' into which the values go,
// and values holds the actual coefficients to be scattered.
//

   int localEnd = localOffset_ + numLocalEqns_ -1;
   if (matrixLoaded_) {
      for(int i=0; i<len; i++){
         if ((localOffset_ > eqnNumbers[i]) || (eqnNumbers[i] > localEnd))
            continue;

         (*x_)[eqnNumbers[i]] = values[i];
      }
   }
   else {
      tmp_x_touched_ = true;
      for(int i=0; i<len; i++){
         if ((localOffset_ > eqnNumbers[i]) || (eqnNumbers[i] > localEnd))
            continue;
         tmp_x_[eqnNumbers[i]-localOffset_] = values[i];
      }
   }
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::getSolution(double* answers, int len) {
//
//The caller must allocate the memory for 'answers',
//and len must be set to the right value -- i.e., len should equal
//numLocalEqns_.
//
   if (len != numLocalEqns_)
      messageAbort("getSolution: len != numLocalEqns_.");

   for(int i=0; i<numLocalEqns_; i++) {
      answers[i] = (*x_)[localOffset_ + i];
   }
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::getSolnEntry(int eqnNumber, double& answer) {
//
//This function returns a single solution entry, the coefficient for
//equation number eqnNumber.
//
   int localEnd = localOffset_ + numLocalEqns_ -1;
   if ((localOffset_ > eqnNumber) || (eqnNumber > localEnd))
      messageAbort("getSolnEntry: eqnNumber out of range.");

   answer = (*x_)[eqnNumber];
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::formResidual(double* values, int len)
{
   if (len != numLocalEqns_) {
      messageAbort("formResidual: len != numLocalEqns_.");
   }

   if (!matrixLoaded_ || !rhsLoaded_) {
      matrixLoadComplete();
   }

   //after the solve x_ and b_ptr_ were 'inv_order'd back into user-ordering.
   //Before we can do this residual calculation, we have to put them into Aztec
   //ordering again.

   int* update_index = NULL;
   if (blockMatrix_) {
      update_index = blkA_ptr_->getUpdate_index();
   }
   else {
      update_index = A_ptr_->getAztec_Map()->update_index;
   }

   AZ_reorder_vec((double*)(x_->startPointer()), azA_->data_org, update_index,
                  azA_->rpntr);

   AZ_reorder_vec((double*)(b_ptr_->startPointer()), azA_->data_org,
                  update_index, azA_->rpntr);

   Aztec_LSVector* r = new Aztec_LSVector(*x_);

   if (blockMatrix_) blkA_ptr_->matvec(*x_, *r); //form r = A*x
   else A_ptr_->matvec(*x_, *r);

   r->addVec(-1.0, *b_ptr_); //form r = A*x - b

   r->scale(-1.0); //form r = b - A*x (not really necessary, but...)

   //now let's get the residual r into user ordering...

   Aztec_LSVector* rtmp = new Aztec_LSVector(*x_);

   AZ_invorder_vec((double*)(r->startPointer()), azA_->data_org, update_index,
                   azA_->rpntr, (double*)rtmp->startPointer());

   *r = *rtmp;

   const double* rptr = r->startPointer();

   for(int i=0; i<numLocalEqns_; i++) {
      values[i] = rptr[i];
   }

   //finally, let's put x_ and b_ptr_ back into user ordering before we exit...
   AZ_invorder_vec((double*)(x_->startPointer()), azA_->data_org, update_index,
                   azA_->rpntr, (double*)rtmp->startPointer());

   *x_ = *rtmp;

   AZ_invorder_vec((double*)(b_ptr_->startPointer()), azA_->data_org,
                   update_index, azA_->rpntr, (double*)rtmp->startPointer());

   *b_ptr_ = *rtmp;

   delete rtmp;
   delete r;
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::selectSolver(const char* name) {
//
// This function takes a string naming the solver and sets the solver choice
// in the options array accordingly.
//
   char msg[64];

   if (!strcmp(name, "AZ_gmres")) {
      aztec_options_[AZ_solver] = AZ_gmres;
      sprintf(msg, "AZ_gmres solver.");
   }
   else if (!strcmp(name, "AZ_cg")) {
      aztec_options_[AZ_solver] = AZ_cg;
      sprintf(msg, "AZ_cg solver.");
   }
   else if (!strcmp(name, "AZ_bicgstab")) {
      aztec_options_[AZ_solver] = AZ_bicgstab;
      sprintf(msg, "AZ_bicgstab solver.");
   }
   else if (!strcmp(name, "AZ_cgs")) {
      aztec_options_[AZ_solver] = AZ_cgs;
      sprintf(msg, "AZ_cgs solver.");
   }
   else if (!strcmp(name, "AZ_tfqmr")) {
      aztec_options_[AZ_solver] = AZ_tfqmr;
      sprintf(msg, "AZ_tfqmr solver.");
   }
   else if (!strcmp(name, "AZ_lu")) {
      aztec_options_[AZ_solver] = AZ_lu;
      sprintf(msg, "AZ_lu solver.");
   }
   else {
      aztec_options_[AZ_solver] = AZ_gmres;
      sprintf(msg, "AZ_gmres default solver.");
      if (thisProc_ == 0) {
         FEI_COUT << "Aztec_LinSysCore: Warning: requested solver <" << name << "> not recognized, defaulting to AZ_gmres." << FEI_ENDL;
      }
   }

   debugOutput(msg);

   return(0);
}

//==============================================================================
int Aztec_LinSysCore::selectPreconditioner(const char* name) {

   char msg[128];
   sprintf(msg, "selectPreconditioner(%s)", name);
   debugOutput(msg);

   if (!strcmp(name, "AZ_none")) {
      aztec_options_[AZ_precond] = AZ_none;
      sprintf(msg, " -- selected: AZ_none.");
   }
   else if (!strcmp(name, "AZ_Jacobi")) {
      aztec_options_[AZ_precond] = AZ_Jacobi;
      sprintf(msg, " -- selected: AZ_Jacobi.");
   }
   else if (!strcmp(name, "AZ_Neumann")) {
      aztec_options_[AZ_precond] = AZ_Neumann;
      sprintf(msg, " -- selected: AZ_Neumann.");
   }
   else if (!strcmp(name, "AZ_ls")) {
      aztec_options_[AZ_precond] = AZ_ls;
      sprintf(msg, " -- selected: AZ_ls.");
   }
   else if (!strcmp(name, "AZ_sym_GS")) {
      aztec_options_[AZ_precond] = AZ_sym_GS;
      sprintf(msg, " -- selected: AZ_sym_GS.");
   }
   else if (!strcmp(name, "AZ_dom_decomp")) {
      aztec_options_[AZ_precond] = AZ_dom_decomp;
      sprintf(msg, " -- selected: AZ_dom_decomp.");
   }
//#ifndef NOT_USING_ML
#ifdef USING_ML
   else if (!strcmp(name, "ML_Vanek")) {
      ML_Vanek_ = true;
      aztec_options_[AZ_precond] = AZ_user_precond;
      sprintf(msg, " -- selected: AZ_user_precond.");
   }
#endif
   else {
      aztec_options_[AZ_precond] = AZ_none;
      sprintf(msg," -- selected: Default, AZ_none.\n");
      if (thisProc_ == 0) {
         FEI_COUT << "Aztec_LinSysCore: Warning: requested preconditioner <" << name << "> not recognized, defaulting to AZ_none." << FEI_ENDL;
      }
   }

   debugOutput(msg);

   return(0);
}

//==============================================================================
void Aztec_LinSysCore::setSubdomainSolve(const char* name) {

   char msg[128];
   sprintf(msg, "setSubdomainSolve(%s)", name);
   debugOutput(msg);

   if (!strcmp(name, "AZ_lu")) {
      aztec_options_[AZ_subdomain_solve] = AZ_lu;
      sprintf(msg, " -- selected AZ_lu");
   }
   else if (!strcmp(name, "AZ_ilu")) {
      aztec_options_[AZ_subdomain_solve] = AZ_ilu;
      sprintf(msg, " -- selected AZ_ilu");
   }
   else if (!strcmp(name, "AZ_ilut")) {
      aztec_options_[AZ_subdomain_solve] = AZ_ilut;
      sprintf(msg, " -- selected AZ_ilut");
   }
   else if (!strcmp(name, "AZ_rilu")) {
      aztec_options_[AZ_subdomain_solve] = AZ_rilu;
      sprintf(msg, " -- selected AZ_rilu");
   }
   else if (!strcmp(name, "AZ_bilu")) {
      aztec_options_[AZ_subdomain_solve] = AZ_bilu;
      sprintf(msg, " -- selected AZ_bilu");
   }
   else if (!strcmp(name, "AZ_icc")) {
      aztec_options_[AZ_subdomain_solve] = AZ_icc;
      sprintf(msg, " -- selected AZ_icc");
   }
   else {
      if (thisProc_ == 0) {
         FEI_COUT << "Aztec_LinSysCore: Warning: requested subdomain-solve <" << name << "> not recognized." << FEI_ENDL;
      }
   }

   debugOutput(msg);
}

//==============================================================================
void Aztec_LinSysCore::setTypeOverlap(const char* name) {

   char msg[128];
   sprintf(msg, "setTypeOverlap(%s)", name);
   debugOutput(msg);

   if (!strcmp(name, "AZ_standard")) {
      aztec_options_[AZ_type_overlap] = AZ_standard;
      sprintf(msg, " -- selected AZ_standard");
   }
   else if (!strcmp(name, "AZ_ilu")) {
      aztec_options_[AZ_type_overlap] = AZ_symmetric;
      sprintf(msg, " -- selected AZ_symmetric");
   }
   else {
      if (thisProc_ == 0) {
         FEI_COUT << "Aztec_LinSysCore: Warning: requested type-overlap <" << name << "> not recognized." << FEI_ENDL;
      }
   }

   debugOutput(msg);
}

//==============================================================================
int Aztec_LinSysCore::writeSystem(const char* name)
{
   writeA(name);
   return(0);
}

//==============================================================================
void Aztec_LinSysCore::recordUserParams()
{
   checkForParam("AZ_tol", numParams_, paramStrings_,
                  aztec_params_[AZ_tol]);

   checkForParam("AZ_drop", numParams_, paramStrings_,
                  aztec_params_[AZ_drop]);

   checkForParam("AZ_ilut_fill", numParams_, paramStrings_,
                 aztec_params_[AZ_ilut_fill]);

   checkForParam("AZ_omega", numParams_, paramStrings_,
                 aztec_params_[AZ_omega]);

   checkForParam("AZ_weights", numParams_, paramStrings_,
                 aztec_params_[AZ_weights]);
}

//==============================================================================
void Aztec_LinSysCore::recordUserOptions()
{
//
//Private function, to be called from launchSolver.
//
   const char* param = NULL;

   param = snl_fei::getParamValue("AZ_solver", numParams_, paramStrings_);
   if (param != NULL){
      selectSolver(param);
   }
 
   param = snl_fei::getParamValue("AZ_precond", numParams_, paramStrings_);
   if (param != NULL){
      selectPreconditioner(param);
   }

   param = snl_fei::getParamValue("AZ_subdomain_solve",
                                         numParams_, paramStrings_);
   if (param != NULL){
      setSubdomainSolve(param);
   }

   param = snl_fei::getParamValue("AZ_scaling", numParams_, paramStrings_);
   if (param != NULL){
      setScalingOption(param);
   }

   param = snl_fei::getParamValue("AZ_conv", numParams_, paramStrings_);
   if (param != NULL){
      setConvTest(param);
   }

   param = snl_fei::getParamValue("AZ_pre_calc",
                                         numParams_, paramStrings_);
   if (param != NULL){
      setPreCalc(param);
   }

   param = snl_fei::getParamValue("AZ_overlap", numParams_, paramStrings_);
   if (param != NULL){
      setOverlap(param);
   }

   param = snl_fei::getParamValue("AZ_type_overlap",
                                         numParams_, paramStrings_);
   if (param != NULL){
      setTypeOverlap(param);
   }

   param = snl_fei::getParamValue("AZ_orthog", numParams_, paramStrings_);
   if (param != NULL){
      setOrthog(param);
   }

   param = snl_fei::getParamValue("AZ_aux_vec", numParams_, paramStrings_);
   if (param != NULL){
      setAuxVec(param);
   }

   param = snl_fei::getParamValue("AZ_output", numParams_, paramStrings_);
   if (param != NULL){
      setAZ_output(param);
   }
   else aztec_options_[AZ_output] = outputLevel_;

   checkForOption("AZ_poly_ord", numParams_, paramStrings_,
                  aztec_options_[AZ_poly_ord]);

   checkForOption("AZ_kspace", numParams_, paramStrings_,
                  aztec_options_[AZ_kspace]);

   checkForOption("AZ_max_iter", numParams_, paramStrings_,
                  aztec_options_[AZ_max_iter]);

   checkForOption("AZ_reorder", numParams_, paramStrings_,
                  aztec_options_[AZ_reorder]);

   checkForOption("AZ_graph_fill", numParams_, paramStrings_,
                  aztec_options_[AZ_graph_fill]);

   checkForOption("AZ_keep_info", numParams_, paramStrings_,
                  aztec_options_[AZ_keep_info]);
}

//==============================================================================
int Aztec_LinSysCore::writeA(const char* name)
{
  if (name == NULL) {
    return(-1);
  }

  FEI_OSTRINGSTREAM osstr;

  if (debugPath_ != NULL) {
    osstr << debugPath_;
  }
  else osstr << ".";

  osstr << "/A_" << name << ".mtx";

  if (blockMatrix_) osstr << ".vbr";

  std::string str = osstr.str();
  const char* matname = str.c_str();

  if (blockMatrix_) {
    blkA_ptr_->writeToFile(matname);
  }
  else {
    A_ptr_->writeToFile(matname);
  }

  return(0);
}

//==============================================================================
int Aztec_LinSysCore::writeVec(Aztec_LSVector* v, const char* name)
{
  if (name == NULL || v == NULL) {
    return(-1);
  }

  FEI_OSTRINGSTREAM osstr;

  if (debugPath_ != NULL) {
    osstr << debugPath_;
  }
  else osstr << ".";

  osstr << "/" << name << ".vec";

  std::string str = osstr.str();

  const char* vecname = str.c_str();

  v->writeToFile(vecname);

  return(0);
}

//==============================================================================
int Aztec_LinSysCore::modifyRHSforBCs()
{
  for(int i=0; i<numLocalEqns_; i++) {
    (*b_ptr_)[i+localOffset_] += tmp_bc_[i];
  }

  if (explicitDirichletBCs_) {
    for(int j=0; j<numEssBCs_; j++) {
      int index = essBCindices_[j];
      (*b_ptr_)[index] = 0.0;
    }
  }
  else {
    for(int j=0; j<numEssBCs_; j++) {
      int index = essBCindices_[j];
      (*b_ptr_)[index] = tmp_bc_[index-localOffset_];
    }
  }

  return(0);
}

//==============================================================================
int Aztec_LinSysCore::explicitlySetDirichletBCs()
{
  for(int j=0; j<numEssBCs_; j++) {
    int index = essBCindices_[j];
    if (rhsLoaded_) {
      (*b_ptr_)[index] = (*bc_)[index];
      (*x_)[index] = (*bc_)[index];
    }
    else {
      (*b_ptr_)[index] = tmp_bc_[index-localOffset_];
      (*x_)[index] = tmp_bc_[index-localOffset_];
    }
  }
  return(0);
}

//==============================================================================
int Aztec_LinSysCore::launchSolver(int& solveStatus, int& iterations) {
//
//This function does any last-second setup required for the
//linear solver, then goes ahead and launches the solver to get
//the solution vector.
//Also, if possible, the number of iterations that were performed
//is stored in the iterations_ variable.
//

   unsigned counter = 0;
   std::map<std::string,unsigned>::iterator 
     iter = named_solve_counter_.find(name_);
   if (iter == named_solve_counter_.end()) {
     fei::console_out() << "fei: Aztec_LinSysCore::LaunchSolver: internal error."
      << FEI_ENDL;
   }
   else {
     counter = iter->second++;
   }

   if (debugOutput_ && outputLevel_ > 1) {
      FEI_OSTRINGSTREAM osstr;
      osstr << name_ << "_Aztec.np"<<numProcs_<<".slv"<<counter;
      std::string str = osstr.str();

      writeA(str.c_str());

      FEI_OSTRINGSTREAM x0_osstr;
      x0_osstr << "x0_" << str;
      std::string x0_str = x0_osstr.str();

      writeVec(x_, x0_str.c_str());

      FEI_OSTRINGSTREAM b_osstr;
      b_osstr << "b_" << str;
      std::string b_str = b_osstr.str();

      writeVec(b_ptr_, b_str.c_str());
   }

   if (needNewPreconditioner_) {
     if (precondCreated_) {
       AZ_precond_destroy(&azP_);
     }

//#ifndef NOT_USING_ML
#ifdef USING_ML
     ML* ml = NULL;

     if (ML_Vanek_) {
       int numFineSweeps = 2;
       int numCoarseSweeps = 2;
       double omega = 0.67;

       initialize_ML(azA_, &azP_, numLevels_,
                     numFineSweeps, numCoarseSweeps, omega,
                     map_->getProcConfig(), &ml);
     }
     else {
       //set the preconditioning matrix Pmat to point to azA_ (Amat).
       azA_->data_org[AZ_name] = 0;
       azP_ = AZ_precond_create(azA_, AZ_precondition, NULL);
     }
#else
     azA_->data_org[AZ_name] = 0;
     azP_ = AZ_precond_create(azA_, AZ_precondition, NULL);
#endif
     needNewPreconditioner_ = false;
     aztec_options_[AZ_pre_calc] = AZ_calc;
     aztec_options_[AZ_conv] = AZ_rhs;
     aztec_options_[AZ_orthog] = AZ_modified;

     recordUserParams();
     recordUserOptions();

     aztec_options_[AZ_keep_info] = 1;
   }
   else {
     aztec_options_[AZ_pre_calc] = AZ_reuse;
   }

   precondCreated_ = true;

   int* proc_config = NULL;
   int* update_index = NULL;
   if (blockMatrix_) {
      proc_config = blkMap_->getProcConfig();
      update_index = blkA_ptr_->getUpdate_index();
   }
   else {
      proc_config = map_->getProcConfig();
      update_index = A_ptr_->getAztec_Map()->update_index;
   }

   AZ_reorder_vec((double*)(x_->startPointer()), azA_->data_org, update_index,
                  azA_->rpntr);

   AZ_reorder_vec((double*)(b_ptr_->startPointer()), azA_->data_org,
                  update_index, azA_->rpntr);

   AZ_iterate((double*)(x_->startPointer()),
              (double*)(b_ptr_->startPointer()),
              aztec_options_, aztec_params_, aztec_status_,
              proc_config, azA_, azP_, azS_);

   iterations = (int)aztec_status_[AZ_its];

   solveStatus = (int)aztec_status_[AZ_why];

   azlsc_solveCounter_++;

   Aztec_LSVector* xtmp = new Aztec_LSVector(*x_);

   //now we need to put x_ back into user-ordering for when we're asked to
   //hand out solution entries.

   AZ_invorder_vec((double*)(x_->startPointer()), azA_->data_org, update_index,
                   azA_->rpntr, (double*)xtmp->startPointer());

   *x_ = *xtmp;

   //let's put b back into user-ordering too...
   AZ_invorder_vec((double*)(b_ptr_->startPointer()), azA_->data_org,
                   update_index, azA_->rpntr, (double*)xtmp->startPointer());

   *b_ptr_ = *xtmp;

   delete xtmp;

   if (explicitDirichletBCs_) explicitlySetDirichletBCs();

   if (debugOutput_) {
      FEI_OSTRINGSTREAM osstr;
      osstr << name_ << "_Aztec.np"<<numProcs_<<".slv"<<counter;
      std::string str = osstr.str();

      FEI_OSTRINGSTREAM x_osstr;
      x_osstr << "x_" << str;
      std::string x_str = x_osstr.str();

      writeVec(x_, x_str.c_str());
   }
   return(0);
}

//==============================================================================
void Aztec_LinSysCore::setScalingOption(const char* param) {

   char* msg = new char[128];
   for(int i=0; i<128; i++) msg[i] = '\0';

   if (!strcmp(param, "AZ_none")) {
      aztec_options_[AZ_scaling] = AZ_none;
      sprintf(msg, "No scaling");
   }
   else if (!strcmp(param, "AZ_Jacobi")) {
      aztec_options_[AZ_scaling] = AZ_Jacobi;
      sprintf(msg, "AZ_Jacobi scaling");
   }
   else if (!strcmp(param, "AZ_BJacobi")) {
      aztec_options_[AZ_scaling] = AZ_BJacobi;
      sprintf(msg, "AZ_BJacobi scaling");
   }
   else if (!strcmp(param, "AZ_row_sum")) {
      aztec_options_[AZ_scaling] = AZ_row_sum;
      sprintf(msg, "AZ_row_sum scaling");
   }
   else if (!strcmp(param, "AZ_sym_diag")) {
      aztec_options_[AZ_scaling] = AZ_sym_diag;
      sprintf(msg, "AZ_sym_diag scaling");
   }
   else if (!strcmp(param, "AZ_sym_row_sum")) {
      aztec_options_[AZ_scaling] = AZ_sym_row_sum;
      sprintf(msg, "AZ_sym_row_sum scaling");
   }
   else {
      if (thisProc_ == 0) {
         FEI_COUT << "Aztec_LinSysCore: Warning: requested scaling <" << param << "> not recognized." << FEI_ENDL;
      }
   }

   debugOutput(msg);

   delete [] msg;

   return;
}

//==============================================================================
void Aztec_LinSysCore::setConvTest(const char* param) {

   char msg[64];

   if (!strcmp(param, "AZ_r0")) {
      aztec_options_[AZ_conv] = AZ_r0;
      sprintf(msg, "AZ_conv AZ_r0");
   }
   else if (!strcmp(param, "AZ_rhs")) {
      aztec_options_[AZ_conv] = AZ_rhs;
      sprintf(msg, "AZ_conv AZ_rhs");
   }
   else if (!strcmp(param, "AZ_Anorm")) {
      aztec_options_[AZ_conv] = AZ_Anorm;
      sprintf(msg, "AZ_conv AZ_Anorm");
   }
   else if (!strcmp(param, "AZ_sol")) {
      aztec_options_[AZ_conv] = AZ_sol;
      sprintf(msg, "AZ_conv AZ_sol");
   }
   else if (!strcmp(param, "AZ_weighted")) {
      aztec_options_[AZ_conv] = AZ_weighted;
      sprintf(msg, "AZ_conv AZ_weighted");
   }
   else if (!strcmp(param, "AZ_noscaled")) {
      aztec_options_[AZ_conv] = AZ_noscaled;
      sprintf(msg, "AZ_conv AZ_noscaled");
   }
   else {
      if (thisProc_ == 0) {
         FEI_COUT << "Aztec_LinSysCore: Warning: requested convergence test <" << param << "> not recognized." << FEI_ENDL;
      }
   }

   debugOutput(msg);

    return;
}

//==============================================================================
void Aztec_LinSysCore::setPreCalc(const char* param)
{
   char msg[64];

   if (!strcmp(param, "AZ_calc")) {
      aztec_options_[AZ_pre_calc] = AZ_calc;
      sprintf(msg, "AZ_pre_calc AZ_calc");
   }
   else if (!strcmp(param, "AZ_recalc")) {
      aztec_options_[AZ_pre_calc] = AZ_recalc;
      sprintf(msg, "AZ_pre_calc AZ_recalc");
   }
   else if (!strcmp(param, "AZ_reuse")) {
      aztec_options_[AZ_pre_calc] = AZ_reuse;
      sprintf(msg, "AZ_pre_calc AZ_reuse");
   }
   else {
      if (thisProc_ == 0) {
         FEI_COUT << "Aztec_LinSysCore: Warning: requested pre_calc <" << param << "> not recognized." << FEI_ENDL;
      }
   }

   debugOutput(msg);

   return;
}

//==============================================================================
void Aztec_LinSysCore::setOverlap(const char* param)
{
   char msg[64];

   if (!strcmp(param, "AZ_none")) {
      aztec_options_[AZ_overlap] = AZ_none;
      sprintf(msg, "AZ_overlap AZ_none");
   }
   else if (!strcmp(param, "AZ_diag")) {
      aztec_options_[AZ_overlap] = AZ_diag;
      sprintf(msg, "AZ_overlap AZ_diag");
   }
   else if (!strcmp(param, "AZ_full")) {
      aztec_options_[AZ_overlap] = AZ_full;
      sprintf(msg, "AZ_overlap AZ_full");
   }
   else {
      checkForOption("AZ_overlap", numParams_, paramStrings_,
                     aztec_options_[AZ_overlap]);
   }

   debugOutput(msg);

   return;
}

//==============================================================================
void Aztec_LinSysCore::setOrthog(const char* param)
{
   char msg[64];

   if (!strcmp(param, "AZ_classic")) {
      aztec_options_[AZ_orthog] = AZ_classic;
      sprintf(msg, "AZ_orthog AZ_classic");
   }
   else if (!strcmp(param, "AZ_modified")) {
      aztec_options_[AZ_orthog] = AZ_modified;
      sprintf(msg, "AZ_orthog AZ_modified");
   }
   else {
      if (thisProc_ == 0) {
         FEI_COUT << "Aztec_LinSysCore: Warning: requested orthog. <" << param << "> not recognized." << FEI_ENDL;
      }
   }

   debugOutput(msg);

   return;
}

//==============================================================================
void Aztec_LinSysCore::setAuxVec(const char* param)
{
   char msg[64];

   if (!strcmp(param, "AZ_resid")) {
      aztec_options_[AZ_aux_vec] = AZ_resid;
      sprintf(msg, "AZ_aux_vec AZ_resid");
   }
   else if (!strcmp(param, "AZ_rand")) {
      aztec_options_[AZ_aux_vec] = AZ_rand;
      sprintf(msg, "AZ_aux_vec AZ_rand");
   }
   else {
      if (thisProc_ == 0) {
         FEI_COUT << "Aztec_LinSysCore: Warning: requested aux_vec <" << param << "> not recognized." << FEI_ENDL;
      }
   }

   debugOutput(msg);

   return;
}

//==============================================================================
void Aztec_LinSysCore::setAZ_output(const char* param)
{
   char msg[64];
   int out = -1;
   int num = sscanf(param, "%d", &out);
   if (num == 1 && out > -1) {
     sprintf(msg, "AZ_output %d", out);
     aztec_options_[AZ_output] = out;
   }
   else if (!strcmp(param, "AZ_all")) {
      aztec_options_[AZ_output] = AZ_all;
      sprintf(msg, "AZ_output AZ_all");
   }
   else if (!strcmp(param, "AZ_none")) {
      aztec_options_[AZ_output] = AZ_none;
      sprintf(msg, "AZ_output AZ_none");
   }
   else if (!strcmp(param, "AZ_warnings")) {
      aztec_options_[AZ_output] = AZ_warnings;
      sprintf(msg, "AZ_output AZ_warnings");
   }
   else if (!strcmp(param, "AZ_last")) {
      aztec_options_[AZ_output] = AZ_last;
      sprintf(msg, "AZ_output AZ_last");
   }
   else {
      if (thisProc_ == 0) {
         FEI_COUT << "Aztec_LinSysCore: Warning: requested AZ_output <" << param << "> not recognized." << FEI_ENDL;
      }
   }

   debugOutput(msg);

   return;
}

//==============================================================================
void Aztec_LinSysCore::checkForParam(const char* paramName,
                                     int numParams, char** paramStrings,
                                     double& param) {
   const char* parameter =
     snl_fei::getParamValue(paramName, numParams, paramStrings);
   if (parameter != NULL) {
      sscanf(parameter, "%le", &param);
   }
}

//==============================================================================
void Aztec_LinSysCore::checkForOption(const char* paramName,
                                      int numParams, char** paramStrings,
                                      int& param) {
  const char* parameter =
    snl_fei::getParamValue(paramName, numParams, paramStrings);
   if (parameter != NULL) {
      sscanf(parameter, "%d", &param);
   }
}

//==============================================================================
void Aztec_LinSysCore::setDebugOutput(const char* path, const char* name){
//
//This function turns on debug output, and opens a file to put it in.
//
   if (debugOutput_) {
      fprintf(debugFile_,"setDebugOutput closing this file.");
      fflush(debugFile_);
      fclose(debugFile_);
      debugFile_ = NULL;
   }

   int pathLength = strlen(path);
   if (path != debugPath_) {
      delete [] debugPath_;
      debugPath_ = new char[pathLength + 1];
      sprintf(debugPath_, "%s", path);
   }

   int nameLength = strlen(name);
   if (name != debugFileName_) {
      delete [] debugFileName_;
      debugFileName_ = new char[nameLength + 1];
      sprintf(debugFileName_,"%s",name);
   }

   char* dbFileName = new char[pathLength + nameLength + 3];

   sprintf(dbFileName, "%s/%s", path, name);

   debugOutput_ = 1;
   debugFile_ = fopen(dbFileName, "w");

   if (!debugFile_){
      fei::console_out() << "couldn't open debug output file: " << dbFileName << FEI_ENDL;
      debugOutput_ = 0;
      delete [] debugPath_;
      debugPath_ = NULL;
      delete [] debugFileName_;
      debugFileName_ = NULL;
   }

   delete [] dbFileName;
}

//==============================================================================
int Aztec_LinSysCore::VBRmatPlusScaledMat(AztecDVBR_Matrix* A,
                                           double scalar,
                                           AztecDVBR_Matrix* source)
{
   int* nnz = new int[numLocalEqnBlks_];
   int* nblk = new int[numLocalEqnBlks_];
   int* src_nnz = new int[numLocalEqnBlks_];
   int* src_nblk = new int[numLocalEqnBlks_];

   if (nnz == NULL || nblk == NULL || src_nnz==NULL || src_nblk==NULL) {
      messageAbort("VBRMatPlusScaledMat: allocation failed");
   }

   A->getNumNonzerosPerRow(nnz);
   A->getNumBlocksPerRow(nblk);
   source->getNumNonzerosPerRow(src_nnz);
   source->getNumBlocksPerRow(src_nblk);

   int i, max_nnz = 0, max_nblk = 0;
   for(i=0; i<numLocalEqnBlks_; i++) {
      if (nnz[i] != src_nnz[i] || nblk[i] != src_nblk[i]) {
         messageAbort("VBRmatPlusScaledMat: matrix sizes don't match.");
      }
      if (max_nnz < nnz[i]) max_nnz = nnz[i];
      if (max_nblk < nblk[i]) max_nblk = nblk[i];
   }
   
   delete [] nnz;
   delete [] nblk;
   delete [] src_nnz;
   delete [] src_nblk;

   double* val = new double[max_nnz];
   int* colInds = new int[max_nblk];
   if (val==NULL || colInds==NULL) {
      messageAbort("VBRmatPlusScaledMat: allocation failed");
   }
   int len, nnzBlks;

   for(i=0; i<numLocalEqnBlks_; i++) {
      int row = localBlkOffset_+i;
      int err = source->getNumBlocksPerRow(row, nnzBlks);
      err += source->getNumNonzerosPerRow(row, len);
      err += source->getBlockRow(row, val, colInds, nnzBlks);

      if (err) messageAbort("VBRmatPlusScaledMat: error getting src row");

      for(int j=0; j<len; j++) val[j] *= scalar;

      err = A->sumIntoBlockRow(row, val, colInds, nnzBlks);
      if (err) messageAbort("VBRmatPlusScaledMat: error summing in row");
   }

   delete [] val;
   delete [] colInds;
   return(0);
}

//==============================================================================
int Aztec_LinSysCore::MSRmatPlusScaledMat(AztecDMSR_Matrix* A,
                                           double scalar,
                                           AztecDMSR_Matrix* source)
{
  return(A->addScaledMatrix(scalar, *source));
}

//==============================================================================
void Aztec_LinSysCore::debugOutput(const char* msg) const {
   if (debugOutput_) {
      fprintf(debugFile_, "%s\n", msg);
      fflush(debugFile_);
   }
}

//==============================================================================
int Aztec_LinSysCore::messageAbort(const char* msg) const {
   fei::console_out() << "Aztec_LinSysCore: " << msg << " Aborting." << FEI_ENDL;
#ifndef FEI_SER
   MPI_Abort(comm_, -1);
#else
   abort();
#endif
   return(-1);
}

}//namespace fei_trilinos

#endif
//HAVE_FEI_AZTECOO

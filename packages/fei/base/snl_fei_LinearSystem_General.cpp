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


#include <fei_macros.hpp>
#include <fei_utils.hpp>

#include <cmath>

#include <snl_fei_LinearSystem_General.hpp>
#include <fei_MatrixReducer.hpp>
#include <fei_Matrix_Impl.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_MatrixGraph.hpp>
#include <fei_SparseRowGraph.hpp>
#include <snl_fei_Constraint.hpp>
#include <fei_Record.hpp>
#include <fei_utils.hpp>
#include <fei_impl_utils.hpp>
#include <fei_LogManager.hpp>

#include <fei_DirichletBCRecord.hpp>
#include <fei_DirichletBCManager.hpp>
#include <fei_EqnBuffer.hpp>
#include <fei_LinSysCoreFilter.hpp>

#undef fei_file
#define fei_file "snl_fei_LinearSystem_General.cpp"
#include <fei_ErrMacros.hpp>

//----------------------------------------------------------------------------
snl_fei::LinearSystem_General::LinearSystem_General(fei::SharedPtr<fei::MatrixGraph>& matrixGraph)
  : fei::LinearSystem(matrixGraph),
    comm_(matrixGraph->getRowSpace()->getCommunicator()),
    essBCvalues_(NULL),
    resolveConflictRequested_(false),
    bcs_trump_slaves_(false),
    explicitBCenforcement_(false),
    BCenforcement_no_column_mod_(false),
    localProc_(0),
    numProcs_(1),
    name_(),
    named_loadcomplete_counter_(),
    iwork_(),
    dwork_(),
    dbgprefix_("LinSysG: ")
{
  localProc_ = fei::localProc(comm_);
  numProcs_  = fei::numProcs(comm_);

  fei::SharedPtr<fei::VectorSpace> vecSpace = matrixGraph->getRowSpace();

  std::vector<int> offsets;
  vecSpace->getGlobalIndexOffsets(offsets);

  firstLocalOffset_ = offsets[localProc_];
  lastLocalOffset_ = offsets[localProc_+1]-1;

  setName("dbg");
}

//----------------------------------------------------------------------------
snl_fei::LinearSystem_General::~LinearSystem_General()
{
  delete essBCvalues_;
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_General::parameters(int numParams,
				   const char* const* paramStrings)
{
  if (numParams == 0 || paramStrings == NULL) return(0);

  const char* param = snl_fei::getParam("name", numParams, paramStrings);
  if (param != NULL) {
    if (strlen(param) < 6) ERReturn(-1);

    setName(&(param[5]));
  }

  param = snl_fei::getParam("resolveConflict",numParams,paramStrings);
  if (param != NULL){
    resolveConflictRequested_ = true;
  }

  param = snl_fei::getParam("BCS_TRUMP_SLAVE_CONSTRAINTS",
                            numParams,paramStrings);
  if (param != NULL) {
    bcs_trump_slaves_ = true;
  }

  param = snl_fei::getParam("EXPLICIT_BC_ENFORCEMENT",numParams,paramStrings);
  if (param != NULL){
    explicitBCenforcement_ = true;
  }

  param = snl_fei::getParam("BC_ENFORCEMENT_NO_COLUMN_MOD",numParams,paramStrings);
  if (param != NULL){
    BCenforcement_no_column_mod_ = true;
  }

  param = snl_fei::getParamValue("FEI_OUTPUT_LEVEL",numParams,paramStrings);
  if (param != NULL) {
    setOutputLevel(fei::utils::string_to_output_level(param));
  }

  if (matrix_.get() != NULL) {
    fei::Matrix* matptr = matrix_.get();
    fei::MatrixReducer* matred = dynamic_cast<fei::MatrixReducer*>(matptr);
    if (matred != NULL) {
      matptr = matred->getTargetMatrix().get();
    }
    fei::Matrix_Impl<LinearSystemCore>* lscmatrix =
      dynamic_cast<fei::Matrix_Impl<LinearSystemCore>*>(matptr);
    if (lscmatrix != NULL) {
      lscmatrix->getMatrix()->parameters(numParams, (char**)paramStrings);
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_General::parameters(const fei::ParameterSet& params)
{
  int numParams = 0;
  const char** paramStrings = NULL;
  std::vector<std::string> stdstrings;
  fei::utils::convert_ParameterSet_to_strings(&params, stdstrings);
  fei::utils::strings_to_char_ptrs(stdstrings, numParams, paramStrings);

  int err = parameters(numParams, paramStrings);

  delete [] paramStrings;

  return(err);
}

//----------------------------------------------------------------------------
void snl_fei::LinearSystem_General::setName(const char* name)
{
  if (name == NULL) return;

  if (name_ == name) return;

  name_ = name;

  std::map<std::string,unsigned>::iterator
    iter = named_loadcomplete_counter_.find(name_);
  if (iter == named_loadcomplete_counter_.end()) {
    static int counter = 0;
    named_loadcomplete_counter_.insert(std::make_pair(name_, counter));
    ++counter;
  }
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_General::loadEssentialBCs(int numIDs,
                         const int* IDs,
                         int idType,
                         int fieldID,
                         int offsetIntoField,
                         const double* prescribedValues)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != 0) {
    FEI_OSTREAM& os = *output_stream_;
    os << "loadEssentialBCs, numIDs: "<<numIDs<<", idType: " <<idType
    <<", fieldID: "<<fieldID<<", offsetIntoField: "<<offsetIntoField<<FEI_ENDL;
  }

  return fei::LinearSystem::loadEssentialBCs(numIDs, IDs, idType, fieldID,
                                           offsetIntoField, prescribedValues);
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_General::loadEssentialBCs(int numIDs,
                         const int* IDs,
                         int idType,
                         int fieldID,
                         const int* offsetsIntoField,
                         const double* prescribedValues)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != 0) {
    FEI_OSTREAM& os = *output_stream_;
    for(int i=0; i<numIDs; ++i) {
      os << "loadEssentialBCs idType: " <<idType
        <<", fieldID: "<<fieldID<<", ID: " << IDs[i]<<", offsetIntoField: "<<offsetsIntoField[i]<<", val: " << prescribedValues[i] << FEI_ENDL;
    }
  }

  return fei::LinearSystem::loadEssentialBCs(numIDs, IDs, idType, fieldID,
                                           offsetsIntoField, prescribedValues);
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_General::loadComplete(bool applyBCs,
                                                bool globalAssemble)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != 0) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"loadComplete"<<FEI_ENDL;
  }

  if (dbcManager_ == NULL) {
    dbcManager_ = new fei::DirichletBCManager(matrixGraph_->getRowSpace());
  }

  if (globalAssemble) {

    if (matrix_.get() != NULL) {
      CHK_ERR( matrix_->gatherFromOverlap() );
    }

    if (rhs_.get() != NULL) {
      CHK_ERR( rhs_->gatherFromOverlap() );
    }

  }

  unsigned counter = 0;

  std::map<std::string,unsigned>::iterator
    iter = named_loadcomplete_counter_.find(name_);
  if (iter == named_loadcomplete_counter_.end()) {
    FEI_COUT << "fei::LinearSystem::loadComplete internal error, name "
      << name_ << " not found." << FEI_ENDL;
  }
  else {
    counter = iter->second++;
  }

  if (output_level_ >= fei::FULL_LOGS) {
    std::string opath = fei::LogManager::getLogManager().getOutputPath();
    if (opath == "") opath = ".";

    FEI_OSTRINGSTREAM Aname;
    FEI_OSTRINGSTREAM bname;
    FEI_OSTRINGSTREAM xname;
    Aname << opath << "/";
    bname << opath << "/";
    xname << opath << "/";

    Aname << "A_"<<name_<<".preBC.np"<<numProcs_<<".slv"<<counter<< ".mtx";

    bname << "b_"<<name_<<".preBC.np"<<numProcs_<<".slv"<<counter<< ".vec";

    std::string Aname_str = Aname.str();
    const char* Aname_c_str = Aname_str.c_str();
    CHK_ERR( matrix_->writeToFile(Aname_c_str) );

    std::string bname_str = bname.str();
    const char* bname_c_str = bname_str.c_str();
    CHK_ERR( rhs_->writeToFile(bname_c_str) );
  }

  CHK_ERR( implementBCs(applyBCs) );

  if (globalAssemble) {
    CHK_ERR( matrix_->globalAssemble() );
  }

  if (output_level_ == fei::STATS || output_level_ == fei::ALL) {
    int globalNumSlaveCRs = matrixGraph_->getGlobalNumSlaveConstraints();
    if (localProc_ == 0) {
      FEI_COUT << "Global Neqns: " << matrix_->getGlobalNumRows();
      if (globalNumSlaveCRs > 0) {
	FEI_COUT << ", Global NslaveCRs: " << globalNumSlaveCRs;
      }
      FEI_COUT << FEI_ENDL;
    }
  }

  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"Neqns=" << matrix_->getGlobalNumRows();
    int globalNumSlaveCRs = matrixGraph_->getGlobalNumSlaveConstraints();
    if (globalNumSlaveCRs > 0) {
      os << ", Global NslaveCRs=" << globalNumSlaveCRs;
    }
    os << FEI_ENDL;
  }

  if (output_level_ >= fei::MATRIX_FILES) {
    std::string opath = fei::LogManager::getLogManager().getOutputPath();
    if (opath == "") opath = ".";

    FEI_OSTRINGSTREAM Aname;
    FEI_OSTRINGSTREAM bname;
    FEI_OSTRINGSTREAM xname;
    Aname << opath << "/";
    bname << opath << "/";
    xname << opath << "/";

    Aname << "A_" <<name_<<".np"<<numProcs_<< ".slv" << counter << ".mtx";

    bname << "b_" <<name_<<".np"<<numProcs_<< ".slv" << counter << ".vec";

    xname << "x0_" <<name_<<".np"<<numProcs_<< ".slv" << counter << ".vec";

    std::string Aname_str = Aname.str();
    const char* Aname_c_str = Aname_str.c_str();
    CHK_ERR( matrix_->writeToFile(Aname_c_str) );

    std::string bname_str = bname.str();
    const char* bname_c_str = bname_str.c_str();
    CHK_ERR( rhs_->writeToFile(bname_c_str) );

    std::string xname_str = xname.str();
    const char* xname_c_str = xname_str.c_str();
    CHK_ERR( soln_->writeToFile(xname_c_str) );
  }

  return(0);
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_General::setBCValuesOnVector(fei::Vector* vector)
{
  if (essBCvalues_ == NULL) {
    return(0);
  }

  if (essBCvalues_->size() == 0) {
    return(0);
  }

  CHK_ERR( vector->copyIn(essBCvalues_->size(),
			  &(essBCvalues_->indices())[0],
			  &(essBCvalues_->coefs())[0]) );

  return(0);
}

//----------------------------------------------------------------------------
bool snl_fei::LinearSystem_General::eqnIsEssentialBC(int globalEqnIndex) const
{
  if (essBCvalues_ == NULL) return(false);

  std::vector<int>& indices = essBCvalues_->indices();
  int offset = fei::binarySearch(globalEqnIndex, indices);
  return( offset < 0 ? false : true);
}

//----------------------------------------------------------------------------
void snl_fei::LinearSystem_General::getEssentialBCs(std::vector<int>& bcEqns,
                                             std::vector<double>& bcVals) const
{
  bcEqns.clear();
  bcVals.clear();
  if (essBCvalues_ == NULL) return;

  int num = essBCvalues_->size();
  bcEqns.resize(num);
  bcVals.resize(num);
  int* essbcs = &(essBCvalues_->indices())[0];
  double* vals = &(essBCvalues_->coefs())[0];
  for(int i=0; i<num; ++i) {
    bcEqns[i] = essbcs[i];
    bcVals[i] = vals[i];
  }
}

//----------------------------------------------------------------------------
void snl_fei::LinearSystem_General::getConstrainedEqns(std::vector<int>& crEqns) const
{
  matrixGraph_->getConstrainedIndices(crEqns);
}

//----------------------------------------------------------------------------
int extractDirichletBCs(fei::DirichletBCManager* bcManager,
                fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                fei::CSVec* essBCvalues,
                bool resolveConflictRequested,
                bool bcs_trump_slaves)
{
//  int numLocalBCs = bcManager->getNumBCRecords();
//  int globalNumBCs = 0;
//  MPI_Comm comm = matrixGraph->getRowSpace()->getCommunicator();
//  fei::GlobalSum(comm, numLocalBCs, globalNumBCs);
//  if (globalNumBCs == 0) {
//    return(0);
//  }

  fei::SharedPtr<fei::FillableMat> localBCeqns(new fei::FillableMat);
  fei::SharedPtr<fei::Matrix_Impl<fei::FillableMat> > bcEqns;
//  matrixGraph->getRowSpace()->initComplete();
  int numSlaves = matrixGraph->getGlobalNumSlaveConstraints();
  fei::SharedPtr<fei::Reducer> reducer = matrixGraph->getReducer();

  int numIndices = numSlaves>0 ?
    reducer->getLocalReducedEqns().size() :
    matrixGraph->getRowSpace()->getNumIndices_Owned();

  bool zeroSharedRows = false;
  bcEqns.reset(new fei::Matrix_Impl<fei::FillableMat>(localBCeqns, matrixGraph, numIndices, zeroSharedRows));
  fei::SharedPtr<fei::Matrix> bcEqns_reducer;
  if (numSlaves > 0) {
    bcEqns_reducer.reset(new fei::MatrixReducer(reducer, bcEqns));
  }

  fei::Matrix& bcEqns_mat = bcEqns_reducer.get()==NULL ?
      *bcEqns : *bcEqns_reducer;

  CHK_ERR( bcManager->finalizeBCEqns(bcEqns_mat, bcs_trump_slaves) );

  if (resolveConflictRequested) {
    fei::SharedPtr<fei::FillableMat> mat = bcEqns->getMatrix();
    std::vector<int> bcEqnNumbers;
    fei::get_row_numbers(*mat, bcEqnNumbers);
    CHK_ERR( snl_fei::resolveConflictingCRs(*matrixGraph, bcEqns_mat,
                                            bcEqnNumbers) );
  }

  std::vector<int> essEqns;
  std::vector<double> values;

  std::map<int,fei::FillableMat*>& remotes = bcEqns->getRemotelyOwnedMatrices();
  std::map<int,fei::FillableMat*>::iterator
    it = remotes.begin(),
    it_end = remotes.end();
  for(; it!=it_end; ++it) {
    fei::impl_utils::separate_BC_eqns( *(it->second), essEqns, values);
  }

//  CHK_ERR( bcEqns->gatherFromOverlap(false) );

  fei::impl_utils::separate_BC_eqns( *(bcEqns->getMatrix()), essEqns, values);

  if (essEqns.size() > 0) {
    int* essEqnsPtr = &essEqns[0];
    double* valuesPtr = &values[0];

    for(unsigned i=0; i<essEqns.size(); ++i) {
      int eqn = essEqnsPtr[i];
      double value = valuesPtr[i];
      fei::put_entry(*essBCvalues, eqn, value);
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_General::implementBCs(bool applyBCs)
{
  if (essBCvalues_ != NULL) {
    delete essBCvalues_;
  }

  essBCvalues_ = new fei::CSVec;

  CHK_ERR( extractDirichletBCs(dbcManager_, matrixGraph_,
                       essBCvalues_,  resolveConflictRequested_,
                      bcs_trump_slaves_) );

  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    std::vector<int>& indices = essBCvalues_->indices();
    std::vector<double>& coefs= essBCvalues_->coefs();
    for(size_t i=0; i<essBCvalues_->size(); ++i) {
      os << "essBCeqns["<<i<<"]: "<<indices[i]<<", "<<coefs[i]<<FEI_ENDL;
    }
  }

  //If the underlying matrix is a LinearSystemCore instance, then this
  //function will return 0, and we're done. A non-zero return-code means
  //we should continue and enforce the BCs assuming a general matrix.

  int returncode = enforceEssentialBC_LinSysCore();
  if (returncode == 0) {
    return(0);
  }

  fei::CSVec allEssBCs;
  if (!BCenforcement_no_column_mod_) {
    fei::impl_utils::global_union(comm_, *essBCvalues_, allEssBCs);

    if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
      FEI_OSTREAM& os = *output_stream_;
      os << "  implementBCs, essBCvalues_.size(): "<<essBCvalues_->size()
         << ", allEssBCs.size(): " << allEssBCs.size()<<FEI_ENDL;
    }
  }

  if (essBCvalues_->size() > 0) {
    enforceEssentialBC_step_1(*essBCvalues_);
  }

  if (!BCenforcement_no_column_mod_ && allEssBCs.size() > 0) {
    enforceEssentialBC_step_2(allEssBCs);
  }

  return(0);
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_General::enforceEssentialBC_LinSysCore()
{
  fei::Matrix* matptr = matrix_.get();
  fei::MatrixReducer* matred = dynamic_cast<fei::MatrixReducer*>(matptr);
  if (matred != NULL) {
    matptr = matred->getTargetMatrix().get();
  }

  fei::Matrix_Impl<LinearSystemCore>* lscmatrix =
    dynamic_cast<fei::Matrix_Impl<LinearSystemCore>*>(matptr);
  if (lscmatrix == 0) {
    return(-1);
  }

  int localsize = matrixGraph_->getRowSpace()->getNumIndices_Owned();
  fei::SharedPtr<fei::Reducer> reducer = matrixGraph_->getReducer();
  if (matrixGraph_->getGlobalNumSlaveConstraints() > 0) {
    localsize = reducer->getLocalReducedEqns().size();
  }

  fei::SharedPtr<fei::FillableMat> inner(new fei::FillableMat);
  bool zeroSharedRows = false;
  fei::SharedPtr<fei::Matrix_Impl<fei::FillableMat> > matrix;
  matrix.reset(new fei::Matrix_Impl<fei::FillableMat>(inner, matrixGraph_, localsize, zeroSharedRows));

  fei::SharedPtr<fei::SparseRowGraph> remoteGraph =
    matrixGraph_->getRemotelyOwnedGraphRows();

  if (!BCenforcement_no_column_mod_) {
    CHK_ERR( snl_fei::gatherRemoteEssBCs(*essBCvalues_, remoteGraph.get(), *matrix) );
  }

  unsigned numBCRows = inner->getNumRows();

  if (output_stream_ != NULL && output_level_ >= fei::BRIEF_LOGS) {
    FEI_OSTREAM& os = *output_stream_;
    os << "#enforceEssentialBC_LinSysCore RemEssBCs to enforce: "
       << numBCRows << FEI_ENDL;
  }

  if (numBCRows > 0 && !BCenforcement_no_column_mod_) {
    std::vector<int*> colIndices(numBCRows);
    std::vector<double*> coefs(numBCRows);
    std::vector<int> colIndLengths(numBCRows);

    fei::CSRMat csrmat(*inner);
    fei::SparseRowGraph& srg = csrmat.getGraph();

    int numEqns = csrmat.getNumRows();
    int* eqns = &(srg.rowNumbers[0]);
    int* rowOffsets = &(srg.rowOffsets[0]);

    for(int i=0; i<numEqns; ++i) {
      colIndices[i] = &(srg.packedColumnIndices[rowOffsets[i]]);
      coefs[i] = &(csrmat.getPackedCoefs()[rowOffsets[i]]);
      colIndLengths[i] = rowOffsets[i+1] - rowOffsets[i];
    }

    int** colInds = &colIndices[0];
    int* colIndLens = &colIndLengths[0];
    double** BCcoefs = &coefs[0];

    if (output_stream_ != NULL && output_level_ > fei::BRIEF_LOGS) {
      FEI_OSTREAM& os = *output_stream_;
      for(int i=0; i<numEqns; ++i) {
        os << "remBCeqn: " << eqns[i] << ", inds/coefs: ";
        for(int j=0; j<colIndLens[i]; ++j) {
          os << "("<<colInds[i][j]<<","<<BCcoefs[i][j]<<") ";
        }
        os << FEI_ENDL;
      }
    }

    int errcode = lscmatrix->getMatrix()->enforceRemoteEssBCs(numEqns,
							      eqns,
							      colInds,
							      colIndLens,
							      BCcoefs);
    if (errcode != 0) {
      return(errcode);
    }
  }

  int numEqns = essBCvalues_->size();
  if (numEqns > 0) {
    int* eqns = &(essBCvalues_->indices())[0];
    double* bccoefs = &(essBCvalues_->coefs())[0];
    std::vector<double> ones(numEqns, 1.0);

    return(lscmatrix->getMatrix()->enforceEssentialBC(eqns, &ones[0],
						    bccoefs, numEqns));
  }

  return(0);
}

//----------------------------------------------------------------------------
void snl_fei::LinearSystem_General::enforceEssentialBC_step_1(fei::CSVec& essBCs)
{
  //to enforce essential boundary conditions, we do the following:
  //
  //  1.  for each eqn (== essBCs.indices()[n]), {
  //        put zeros in row A[eqn], but leave 1.0 on the diagonal
  //        set b[eqn] = essBCs.coefs()[n]
  //      }
  //
  //  2.  for i in 1..numRows (i.e., all rows) {
  //        if (i in bcEqns) continue;
  //        b[i] -= A[i,eqn] * essBCs.coefs()[n]
  //        A[i,eqn] = 0.0;
  //      }
  //
  //It is important to note that for step 1, essBCs need only contain
  //local eqns, but for step 2 it should contain *ALL* bc eqns.
  //
  //This function performs step 1.

  int numEqns = essBCs.size();
  int* eqns = &(essBCs.indices())[0];
  double* bcCoefs = &(essBCs.coefs())[0];

  std::vector<double> coefs;
  std::vector<int> indices;

  fei::SharedPtr<fei::Reducer> reducer = matrixGraph_->getReducer();
  bool haveSlaves = reducer.get()!=NULL;

  try {
  for(int i=0; i<numEqns; i++) {
    int eqn = eqns[i];

    //if slave-constraints are present, the incoming bc-eqns are in
    //the reduced equation space. so we actually have to translate them back
    //to the unreduced space before passing them into the fei::Matrix object,
    //because the fei::Matrix object has machinery to translate unreduced eqns
    //to the reduced space.
    //Also, our firstLocalOffset_ and lastLocalOffset_ attributes are in the
    //unreduced space.
    if (haveSlaves) {
      eqn = reducer->translateFromReducedEqn(eqn);
    }

    if (eqn < firstLocalOffset_ || eqn > lastLocalOffset_) continue;

    //put gamma/alpha on the rhs for this ess-BC equation.
    double bcValue = bcCoefs[i];
    int err = rhs_->copyIn(1, &eqn, &bcValue);
    if (err != 0) {
      FEI_OSTRINGSTREAM osstr;
      osstr <<"snl_fei::LinearSystem_General::enforceEssentialBC_step_1 ERROR: "
	    << "err="<<err<<" returned from rhs_->copyIn row="<<eqn;
      throw std::runtime_error(osstr.str());
    }

    err = getMatrixRow(matrix_.get(), eqn, coefs, indices);
    if (err != 0 || indices.size() < 1) {
      continue;
    }

    int rowLen = indices.size();
    int* indPtr = &indices[0];

    //first, put zeros in the row and 1.0 on the diagonal...
    for(int j=0; j<rowLen; j++) {
      if (indPtr[j] == eqn) coefs[j] = 1.0;
      else coefs[j] = 0.0;
    }

    double* coefPtr = &coefs[0];

    err = matrix_->copyIn(1, &eqn, rowLen, indPtr, &coefPtr);
    if (err != 0) {
      FEI_OSTRINGSTREAM osstr;
      osstr <<"snl_fei::LinearSystem_General::enforceEssentialBC_step_1 ERROR: "
	    << "err="<<err<<" returned from matrix_->copyIn row="<<eqn;
      throw std::runtime_error(osstr.str());
    }
  }//for i
  }
  catch(std::runtime_error& exc) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei::LinearSystem::enforceEssentialBC: ERROR, caught exception: "
        << exc.what();
    throw std::runtime_error(osstr.str());
  }
}

//----------------------------------------------------------------------------
void snl_fei::LinearSystem_General::enforceEssentialBC_step_2(fei::CSVec& essBCs)
{
  //to enforce essential boundary conditions, we do the following:
  //
  //  1.  for each eqn (== essBCs.indices()[n]), {
  //        put zeros in row A[eqn], but leave 1.0 on the diagonal
  //        set b[eqn] = essBCs.coefs()[n]
  //      }
  //
  //  2.  for i in 1..numRows (i.e., all rows) {
  //        if (i in bcEqns) continue;
  //        b[i] -= A[i,eqn] * essBCs.coefs()[n]
  //        A[i,eqn] = 0.0;
  //      }
  //
  //It is important to note that for step 1, essBCs need only contain
  //local eqns, but for step 2 it should contain *ALL* bc eqns.
  //
  //This function performs step 2.

  int numBCeqns = essBCs.size();
  if (numBCeqns < 1) {
    return;
  }

  int* bcEqns = &(essBCs.indices())[0];
  double* bcCoefs = &(essBCs.coefs())[0];

  fei::SharedPtr<fei::Reducer> reducer = matrixGraph_->getReducer();
  bool haveSlaves = reducer.get()!=NULL;
  if (haveSlaves) {
    for(int i=0; i<numBCeqns; ++i) {
      bcEqns[i] = reducer->translateFromReducedEqn(bcEqns[i]);
    }
  }

  int firstBCeqn = bcEqns[0];
  int lastBCeqn = bcEqns[numBCeqns-1];

  std::vector<double> coefs;
  std::vector<int> indices;

  int insertPoint;

  int nextBCeqnOffset = 0;
  int nextBCeqn = bcEqns[nextBCeqnOffset];

  for(int i=firstLocalOffset_; i<=lastLocalOffset_; ++i) {
    if (haveSlaves) {
      if (reducer->isSlaveEqn(i)) continue;
    }

    bool should_continue = false;
    if (i >= nextBCeqn) {
      if (i == nextBCeqn) {
	++nextBCeqnOffset;
	if (nextBCeqnOffset < numBCeqns) {
	  nextBCeqn = bcEqns[nextBCeqnOffset];
	}
	else {
	  nextBCeqn = lastLocalOffset_+1;
	}

	should_continue = true;
      }
      else {
	while(nextBCeqn <= i) {
	  if (nextBCeqn == i) should_continue = true;
	  ++nextBCeqnOffset;
	  if (nextBCeqnOffset < numBCeqns) {
	    nextBCeqn = bcEqns[nextBCeqnOffset];
	  }
	  else {
	    nextBCeqn = lastLocalOffset_+1;
	  }
	}
      }
    }

    if (should_continue) continue;

    int err = getMatrixRow(matrix_.get(), i, coefs, indices);
    if (err != 0 || indices.size() < 1) {
      continue;
    }

    int numIndices = indices.size();
    int* indicesPtr = &indices[0];
    double* coefsPtr = &coefs[0];
    bool modifiedCoef = false;

    fei::insertion_sort_with_companions(numIndices, indicesPtr, coefsPtr);

    if (indicesPtr[0] > lastBCeqn || indicesPtr[numIndices-1] < firstBCeqn) {
      continue;
    }

    double value = 0.0;
    int offset = 0;

    for(int j=0; j<numIndices; ++j) {
      int idx = indicesPtr[j];
      offset = fei::binarySearch(idx, bcEqns, numBCeqns,
				     insertPoint);
      if (offset > -1) {
	value -= bcCoefs[offset]*coefsPtr[j];

	coefsPtr[j] = 0.0;
	modifiedCoef = true;
      }
    }

    if (modifiedCoef) {
      err = matrix_->copyIn(1, &i, numIndices, indicesPtr, &coefsPtr);
      if (err != 0) {
	FEI_OSTRINGSTREAM osstr;
	osstr <<"snl_fei::LinearSystem_General::enforceEssentialBC_step_2 ERROR: "
	      << "err="<<err<<" returned from matrix_->copyIn, row="<<i;
	throw std::runtime_error(osstr.str());
      }
    }

    const double fei_eps = 1.e-49;
    if (std::abs(value) > fei_eps) {
      rhs_->sumIn(1, &i, &value);

      if (output_level_ >= fei::FULL_LOGS && output_stream_ != 0) {
	FEI_OSTREAM& os = *output_stream_;
	os << "enfEssBC_step2: rhs["<<i<<"] += "<<value<<FEI_ENDL;
      }
    }
  }
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_General::getMatrixRow(fei::Matrix* matrix, int row,
						std::vector<double>& coefs,
						std::vector<int>& indices)
{
  int len = 0;
  int err = matrix->getRowLength(row, len);
  if (err != 0) {
    coefs.resize(0);
    indices.resize(0);
    return(err);
  }

  if ((int)coefs.size() != len) {
    coefs.resize(len);
  }
  if ((int)indices.size() != len) {
    indices.resize(len);
  }

  CHK_ERR( matrix->copyOutRow(row, len, &coefs[0], &indices[0]));

  return(0);
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_General::loadLagrangeConstraint(int constraintID,
							  const double *weights,
							  double rhsValue)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << "loadLagrangeConstraint crID: "<<constraintID<<FEI_ENDL;
  }

  Constraint<fei::Record<int>*>* cr =
    matrixGraph_->getLagrangeConstraint(constraintID);
  if (cr == NULL) {
    return(-1);
  }

  CHK_ERR( matrixGraph_->getConstraintConnectivityIndices(cr, iwork_) );

  //Let's attach the weights to the constraint-record now.
  std::vector<double>& cr_weights = cr->getMasterWeights();
  cr_weights.resize(iwork_.size());
  for(unsigned i=0; i<iwork_.size(); ++i) {
    cr_weights.push_back(weights[i]);
  }

  fei::SharedPtr<fei::VectorSpace> vecSpace = matrixGraph_->getRowSpace();

  int crEqn = -1;
  CHK_ERR( vecSpace->getGlobalIndex(cr->getIDType(),
				     cr->getConstraintID(),
				     crEqn) );

  //now add the row contribution to the matrix and rhs
  int numIndices = iwork_.size();
  int* indicesPtr = &(iwork_[0]);

  CHK_ERR( matrix_->sumIn(1, &crEqn, numIndices, indicesPtr, &weights) );

  CHK_ERR( rhs_->sumIn(1, &crEqn, &rhsValue) );

  //now add the column contributions to the matrix
  for(int k=0; k<numIndices; ++k) {
    double* thisWeight = (double*)(&(weights[k]));
    CHK_ERR( matrix_->sumIn(1, &(indicesPtr[k]), 1, &crEqn, &thisWeight) );
  }

  return(0);
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_General::loadPenaltyConstraint(int constraintID,
							 const double *weights,
							 double penaltyValue,
							 double rhsValue)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << "loadPenaltyConstraint crID: "<<constraintID<<FEI_ENDL;
  }

  Constraint<fei::Record<int>*>* cr =
    matrixGraph_->getPenaltyConstraint(constraintID);
  if (cr == NULL) {
    return(-1);
  }

  CHK_ERR( matrixGraph_->getConstraintConnectivityIndices(cr, iwork_) );

  fei::SharedPtr<fei::VectorSpace> vecSpace = matrixGraph_->getRowSpace();

  int numIndices = iwork_.size();
  int* indicesPtr = &(iwork_[0]);

  //now add the contributions to the matrix and rhs
  std::vector<double> coefs(numIndices);
  double* coefPtr = &coefs[0];
  for(int i=0; i<numIndices; ++i) {
    for(int j=0; j<numIndices; ++j) {
      coefPtr[j] = weights[i]*weights[j]*penaltyValue;
    }
    CHK_ERR( matrix_->sumIn(1, &(indicesPtr[i]), numIndices, indicesPtr,
			    &coefPtr) );

    double rhsCoef = weights[i]*penaltyValue*rhsValue;
    CHK_ERR( rhs_->sumIn(1, &(indicesPtr[i]), &rhsCoef) );
  }

  return(0);
}


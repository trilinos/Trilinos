/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_trilinos_macros.hpp"
#include <fei_Factory_Trilinos.hpp>

#include <fei_VectorReducer.hpp>
#include <fei_MatrixReducer.hpp>
#include <fei_Matrix_Local.hpp>
#include <fei_Vector_Local.hpp>

#include <fei_Solver_AztecOO.hpp>
#include <fei_Solver_Amesos.hpp>

Factory_Trilinos::Factory_Trilinos(MPI_Comm comm)
  : fei::Factory(comm),
    comm_(comm),
    reducer_(),
    lpm_epetrabasic_(),
    use_lpm_epetrabasic_(false),
    useAmesos_(false),
    use_feiMatrixLocal_(false),
    blockEntryMatrix_(false),
    orderRowsWithLocalColsFirst_(false),
    outputLevel_(0)
{
}

Factory_Trilinos::~Factory_Trilinos()
{
}

int Factory_Trilinos::parameters(int numParams,
                                  const char* const* paramStrings)
{
  std::vector<std::string> stdstrings;
  fei::utils::char_ptrs_to_strings(numParams, paramStrings, stdstrings);

  fei::ParameterSet paramset;
  fei::utils::parse_strings(stdstrings, " ", paramset);

  parameters(paramset);
  return(0);
}

void Factory_Trilinos::parameters(const fei::ParameterSet& parameterset)
{
  fei::Factory::parameters(parameterset);

  parameterset.getIntParamValue("outputLevel", outputLevel_);

  bool blkGraph = false;
  bool blkMatrix = false;

  parameterset.getBoolParamValue("BLOCK_GRAPH", blkGraph);
  parameterset.getBoolParamValue("BLOCK_MATRIX", blkMatrix);

  blockEntryMatrix_ = (blkGraph || blkMatrix);

  const fei::Param* param = parameterset.get("Trilinos_Solver");
  if (param != 0) {
    if (param->getType() == fei::Param::STRING) {
      std::string strval = param->getStringValue();
      std::string::size_type ii = strval.find("Amesos");
      if (ii != std::string::npos) {
        useAmesos_ = true;
      }
    }
  }

  parameterset.getBoolParamValue("LPM_EpetraBasic", use_lpm_epetrabasic_);
  if (use_lpm_epetrabasic_) {
    create_LinProbMgr();
  }

  parameterset.getBoolParamValue("USE_FEI_MATRIX_LOCAL", use_feiMatrixLocal_);

  parameterset.getBoolParamValue("ORDER_ROWS_WITH_LOCAL_COLS_FIRST",
                                 orderRowsWithLocalColsFirst_);
}

fei::SharedPtr<fei::MatrixGraph>
Factory_Trilinos::createMatrixGraph(fei::SharedPtr<fei::VectorSpace> rowSpace,
                                   fei::SharedPtr<fei::VectorSpace> colSpace,
                                   const char* name)
{
  static fei::MatrixGraph_Impl2::Factory factory2;
  if (!use_lpm_epetrabasic_) {
    return(factory2.createMatrixGraph(rowSpace, colSpace, name));
  }

  return(factory2.createMatrixGraph(rowSpace, colSpace, name));
}

fei::SharedPtr<fei::Vector>
Factory_Trilinos::wrapVector(fei::SharedPtr<fei::VectorSpace> vecSpace,
                             fei::SharedPtr<Epetra_MultiVector> multiVec)
{
  fei::SharedPtr<fei::Vector> vecptr;
  if (vecSpace->getNumIndices_Owned() != multiVec->Map().NumMyPoints()) {
    //vecSpace isn't compatible with multiVec's map, so return empty vecptr
    return(vecptr);
  }

  vecptr.reset(new fei::Vector_Impl<Epetra_MultiVector>(vecSpace,
                                                       multiVec.get(),
                                     multiVec->Map().NumMyPoints(), false));
  return(vecptr);
}

fei::SharedPtr<fei::Vector>
Factory_Trilinos::wrapVector(fei::SharedPtr<fei::MatrixGraph> matGraph,
                             fei::SharedPtr<Epetra_MultiVector> multiVec)
{
  fei::SharedPtr<fei::Vector> vecptr;
  if (matGraph->getRowSpace()->getNumIndices_Owned() !=
      multiVec->Map().NumMyPoints()) {
    //vector-space isn't compatible with multiVec's map, so return empty vecptr
    return(vecptr);
  }

  vecptr.reset(new fei::Vector_Impl<Epetra_MultiVector>(matGraph->getRowSpace(),
                                                       multiVec.get(),
                                        multiVec->Map().NumMyPoints(), false));
  return(vecptr);
}

fei::SharedPtr<fei::Vector>
Factory_Trilinos::createVector(fei::SharedPtr<fei::VectorSpace> vecSpace,
                               bool isSolutionVector,
                               int numVectors)
{
  if (use_feiMatrixLocal_) {
    fei::SharedPtr<fei::Vector> feivec(new fei::Vector_Local(vecSpace));
    return(feivec);
  }

//  if (blockEntryMatrix_) {
//    throw std::runtime_error("Factory_Trilinos: fei ERROR, block-entry matrices/vectors not currently supported.");
//  }

  std::vector<int> indices;
  int err = 0, localSize;
  if (reducer_.get() != NULL) {
    indices = reducer_->getLocalReducedEqns();
    localSize = indices.size();
  }
  else {
    if (blockEntryMatrix_) {
      localSize = vecSpace->getNumBlkIndices_Owned();
      indices.resize(localSize*2);
      err = vecSpace->getBlkIndices_Owned(localSize, &indices[0], &indices[localSize], localSize);
    }
    else {
      localSize = vecSpace->getNumIndices_Owned();
      indices.resize(localSize);
      err = vecSpace->getIndices_Owned(localSize, &indices[0], localSize);
    }
  }
  if (err != 0) {
    throw std::runtime_error("error in vecSpace->getIndices_Owned");
  }

  fei::SharedPtr<fei::Vector> feivec, tmpvec;
  if (!use_lpm_epetrabasic_) {
    try {
      Epetra_BlockMap emap = blockEntryMatrix_ ? 
        Trilinos_Helpers::create_Epetra_BlockMap(vecSpace) :
        Trilinos_Helpers::create_Epetra_Map(comm_, indices);

      Epetra_MultiVector* emvec = new Epetra_MultiVector(emap, numVectors);

      tmpvec.reset(new fei::Vector_Impl<Epetra_MultiVector>(vecSpace,
                                              emvec, localSize,
                                        isSolutionVector, true));
    }
    catch(std::runtime_error& exc) {
      FEI_CERR << "Factory_Trilinos::createVector: caught exception '"
               << exc.what() << "', re-throwing..." << FEI_ENDL;
      throw exc;
    }
  }
  else {
    tmpvec.reset(new fei::Vector_Impl<fei::LinearProblemManager>(vecSpace,
                                    lpm_epetrabasic_.get(), localSize,
                                        isSolutionVector));
  }

  if (reducer_.get() != NULL) {
    feivec.reset(new fei::VectorReducer(reducer_,
                                        tmpvec, isSolutionVector));
  }
  else {
    feivec = tmpvec;
  }

  return(feivec);
}

fei::SharedPtr<fei::Vector>
Factory_Trilinos::createVector(fei::SharedPtr<fei::VectorSpace> vecSpace,
                               int numVectors)
{
  bool isSolnVector = false;
  return(createVector(vecSpace, isSolnVector, numVectors));
}

fei::SharedPtr<fei::Vector>
Factory_Trilinos::createVector(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                              int numVectors)
{
  bool isSolnVector = false;
  return(createVector(matrixGraph, isSolnVector, numVectors));
}

fei::SharedPtr<fei::Vector>
Factory_Trilinos::createVector(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                               bool isSolutionVector,
                               int numVectors)
{
  if (use_feiMatrixLocal_) {
    fei::SharedPtr<fei::Vector> feivec(new fei::Vector_Local(matrixGraph->getRowSpace()));
    return(feivec);
  }

//  if (blockEntryMatrix_) {
//    throw std::runtime_error("Factory_Trilinos: fei ERROR, block-entry matrices/vectors not currently supported.");
//  }

  int globalNumSlaves = matrixGraph->getGlobalNumSlaveConstraints();

  if (globalNumSlaves > 0 && reducer_.get()==NULL) {
    reducer_ = matrixGraph->getReducer();
  }

  fei::SharedPtr<fei::Vector> feivec, tmpvec;

  std::vector<int> indices;
  int err = 0, localSize;
  fei::SharedPtr<fei::VectorSpace> vecSpace = matrixGraph->getRowSpace();
  if (reducer_.get() != NULL) {
    indices = reducer_->getLocalReducedEqns();
    localSize = indices.size();
  }
  else {
    if (blockEntryMatrix_) {
      localSize = vecSpace->getNumBlkIndices_Owned();
      indices.resize(localSize*2);
      err = vecSpace->getBlkIndices_Owned(localSize, &indices[0], &indices[localSize], localSize);
    }
    else {
      localSize = vecSpace->getNumIndices_Owned();
      indices.resize(localSize);
      err = vecSpace->getIndices_Owned(localSize, &indices[0], localSize);
    }
  }
  if (err != 0) {
    throw std::runtime_error("error in vecSpace->getIndices_Owned");
  }

  if (!use_lpm_epetrabasic_) {
    try {
      Epetra_BlockMap emap = blockEntryMatrix_ ?
        Trilinos_Helpers::create_Epetra_BlockMap(vecSpace) :
        Trilinos_Helpers::create_Epetra_Map(comm_, indices);

      Epetra_MultiVector* emvec = new Epetra_MultiVector(emap, numVectors);

      tmpvec.reset(new fei::Vector_Impl<Epetra_MultiVector>(matrixGraph->getRowSpace(), emvec,
                                                  localSize,
                                                       isSolutionVector, true));
    }
    catch(std::runtime_error& exc) {
      FEI_CERR << "Factory_Trilinos::createVector: caught exception '"
               << exc.what() << "', re-throwing..." << FEI_ENDL;
      throw exc;
    }
  }
  else {
    vecSpace = matrixGraph->getRowSpace();

    lpm_epetrabasic_->setRowDistribution(indices);
    tmpvec.reset(new fei::Vector_Impl<fei::LinearProblemManager>(vecSpace,
                                   lpm_epetrabasic_.get(),
                                   localSize, isSolutionVector));
  }

  if (reducer_.get() != NULL) {
    feivec.reset(new fei::VectorReducer(reducer_, tmpvec, isSolutionVector));
  }
  else {
    feivec = tmpvec;
  }

  return(feivec);
}

fei::SharedPtr<fei::Matrix>
Factory_Trilinos::createMatrix(fei::SharedPtr<fei::MatrixGraph> matrixGraph)
{
  int globalNumSlaves = matrixGraph->getGlobalNumSlaveConstraints();

  if (globalNumSlaves > 0 && reducer_.get()==NULL) {
    reducer_ = matrixGraph->getReducer();
  }

  if (use_lpm_epetrabasic_) {
    return(
      Trilinos_Helpers::create_from_LPM_EpetraBasic(matrixGraph,
                                                    blockEntryMatrix_,
                                                    reducer_,
                                                    lpm_epetrabasic_)
    );
  }

  if (use_feiMatrixLocal_) {
    return(fei::Matrix_Local::create_Matrix_Local(matrixGraph, blockEntryMatrix_));
  }

  return(
      Trilinos_Helpers::create_from_Epetra_Matrix(matrixGraph,
                                                  blockEntryMatrix_,
                                                  reducer_,
                                                  orderRowsWithLocalColsFirst_)
  );
}

fei::SharedPtr<fei::Solver>
Factory_Trilinos::createSolver(const char* name)
{
  fei::SharedPtr<fei::Solver> solver;

  if (useAmesos_ || name != 0) {
    if (name != 0) {
      std::string sname(name);
      std::string::size_type ii = sname.find("Amesos");
      if (ii == std::string::npos) {
        solver.reset(new Solver_AztecOO);
        return(solver);
      }
    }
#ifdef HAVE_AMESOS
    solver.reset(new Solver_Amesos);
#else
    FEI_CERR << "Factory_Trilinos::createSolver: ERROR, Amesos requested,"
      << " but HAVE_AMESOS is not defined so Amesos is not available."
      <<FEI_ENDL;
    return(solver);
#endif
  }
  else {
    solver.reset(new Solver_AztecOO);
  }

  return(solver);
}

void Factory_Trilinos::create_LinProbMgr(bool replace_if_already_created)
{
  if (!use_lpm_epetrabasic_) {
    return;
  }

  bool need_to_create_lpm = false;

  if (lpm_epetrabasic_.get() == NULL) {
    need_to_create_lpm = true;
  }

  if (replace_if_already_created) {
    need_to_create_lpm = true;
  }

  if (need_to_create_lpm) {
    fei::SharedPtr<fei::LinearProblemManager>
      newlpm(new LinProbMgr_EpetraBasic(comm_));

    lpm_epetrabasic_ = newlpm;
  }
}


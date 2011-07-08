/*--------------------------------------------------------------------*/
/*    Copyright 2006 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <fei_macros.hpp>

#include <fei_Include_Trilinos.hpp>

#ifdef HAVE_FEI_EPETRA
#include <fei_VectorTraits_Epetra.hpp>
#include <fei_MatrixTraits_Epetra.hpp>
#include <fei_LinProbMgr_EpetraBasic.hpp>
#endif

#include <fei_Trilinos_Helpers.hpp>
#include <fei_ParameterSet.hpp>
#include <fei_Vector_Impl.hpp>
#include <fei_VectorReducer.hpp>
#include <fei_Matrix_Impl.hpp>
#include <fei_MatrixReducer.hpp>
//#include <EpetraExt_BlockMapOut.h>

namespace Trilinos_Helpers {

#ifdef HAVE_FEI_EPETRA

Epetra_Map
create_Epetra_Map(MPI_Comm comm,
                  const std::vector<int>& local_eqns)
{
#ifndef FEI_SER
  Epetra_MpiComm EComm(comm);
#else
  Epetra_SerialComm EComm;
#endif

  int localSize = local_eqns.size();
  int globalSize = 0;
  EComm.SumAll(&localSize, &globalSize, 1);

  if (localSize < 0 || globalSize < 0) {
    throw std::runtime_error("Trilinos_Helpers::create_Epetra_Map: negative local or global size.");
  }

  Epetra_Map EMap(globalSize, localSize, 0, EComm);
  return(EMap);
}

Epetra_BlockMap
create_Epetra_BlockMap(const fei::SharedPtr<fei::VectorSpace>& vecspace)
{
  if (vecspace.get() == 0) {
    throw std::runtime_error("create_Epetra_Map needs non-null fei::VectorSpace");
  }

#ifndef FEI_SER
  MPI_Comm comm = vecspace->getCommunicator();
  Epetra_MpiComm EComm(comm);
#else
  Epetra_SerialComm EComm;
#endif

  int localSizeBlk = vecspace->getNumBlkIndices_Owned();
  int globalSizeBlk = vecspace->getGlobalNumBlkIndices();

  if (localSizeBlk < 0 || globalSizeBlk < 0) {
    throw std::runtime_error("Trilinos_Helpers::create_Epetra_BlockMap: fei::VectorSpace has negative local or global size.");
  }

  std::vector<int> blkEqns(localSizeBlk*2);
  int* blkEqnsPtr = &(blkEqns[0]);

  int chkNum = 0;
  int errcode = vecspace->getBlkIndices_Owned(localSizeBlk,
					      blkEqnsPtr, blkEqnsPtr+localSizeBlk,
					      chkNum);
  if (errcode != 0) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "create_Epetra_BlockMap ERROR, nonzero errcode="<<errcode
	  << " returned by vecspace->getBlkIndices_Owned.";
    throw std::runtime_error(osstr.str());
  }

  Epetra_BlockMap EBMap(globalSizeBlk, localSizeBlk,
			blkEqnsPtr, blkEqnsPtr+localSizeBlk, 0, EComm);

  return(EBMap);
}

Epetra_CrsGraph
create_Epetra_CrsGraph(const fei::SharedPtr<fei::MatrixGraph>& matgraph,
		       bool blockEntries, bool orderRowsWithLocalColsFirst)
{
  fei::SharedPtr<fei::SparseRowGraph> localSRGraph =
    matgraph->createGraph(blockEntries);
  if (localSRGraph.get() == NULL) {
    throw std::runtime_error("create_Epetra_CrsGraph ERROR in fei::MatrixGraph::createGraph");
  }

  int numLocallyOwnedRows = localSRGraph->rowNumbers.size();
  int* rowNumbers = &(localSRGraph->rowNumbers[0]);
  int* rowOffsets = &(localSRGraph->rowOffsets[0]);
  int* packedColumnIndices = &(localSRGraph->packedColumnIndices[0]);

  fei::SharedPtr<fei::VectorSpace> vecspace = matgraph->getRowSpace();
  MPI_Comm comm = vecspace->getCommunicator();
  std::vector<int>& local_eqns = localSRGraph->rowNumbers;

  Epetra_BlockMap emap = blockEntries ?
    create_Epetra_BlockMap(vecspace) : create_Epetra_Map(comm, local_eqns);

  if (orderRowsWithLocalColsFirst == true &&
      emap.Comm().NumProc() > 2 && !blockEntries) {
    bool* used_row = new bool[local_eqns.size()];
    for(unsigned ii=0; ii<local_eqns.size(); ++ii) used_row[ii] = false;

    int offset = 0;
    std::vector<int> ordered_local_eqns(local_eqns.size());
    for(unsigned ii=0; ii<local_eqns.size(); ++ii) {
      bool row_has_off_proc_cols = false;
      for(int j=rowOffsets[ii]; j<rowOffsets[ii+1]; ++j) {
        if (emap.MyGID(packedColumnIndices[j]) == false) {
          row_has_off_proc_cols = true;
          break;
        }
      }

      if (row_has_off_proc_cols == false) {
        ordered_local_eqns[offset++] = rowNumbers[ii];
        used_row[ii] = true;
      }
    }

    for(unsigned ii=0; ii<local_eqns.size(); ++ii) {
      if (used_row[ii] == true) continue;
      ordered_local_eqns[offset++] = rowNumbers[ii];
    }

    emap = create_Epetra_Map(comm, ordered_local_eqns);
    delete [] used_row;
  }

//  EpetraExt::BlockMapToMatrixMarketFile("EBMap.np12.mm",emap,"AriaTest");

  std::vector<int> rowLengths; rowLengths.reserve(numLocallyOwnedRows);
  for(int ii=0; ii<numLocallyOwnedRows; ++ii) {
    rowLengths.push_back(rowOffsets[ii+1]-rowOffsets[ii]);
  }

  bool staticProfile = true;
  Epetra_CrsGraph egraph(Copy, emap, &(rowLengths[0]), staticProfile);

  const Epetra_Comm& ecomm = emap.Comm();
  int localProc = ecomm.MyPID();

  int firstLocalEqn = numLocallyOwnedRows > 0 ? rowNumbers[0] : -1;

  int offset = 0;
  for(int i=0; i<numLocallyOwnedRows; ++i) {
    int err = egraph.InsertGlobalIndices(firstLocalEqn+i,
					 rowLengths[i],
					 &(packedColumnIndices[offset]));
    if (err != 0) {
      fei::console_out() << "proc " << localProc << " err-return " << err
               << " inserting row " << firstLocalEqn+i<<", cols ";
      for(int ii=0; ii<rowLengths[i]; ++ii) {
	fei::console_out() << packedColumnIndices[offset+ii]<<",";
      }
      fei::console_out() << FEI_ENDL;
      throw std::runtime_error("... occurred in create_Epetra_CrsGraph");
    }

    offset += rowLengths[i];
  }

  //Epetra_BlockMap* domainmap = const_cast<Epetra_BlockMap*>(&(epetraGraph_->DomainMap()));
  //Epetra_BlockMap* rangemap = const_cast<Epetra_BlockMap*>(&(epetraGraph_->RangeMap()));
  egraph.FillComplete();

  return(egraph);
}

fei::SharedPtr<fei::Matrix>
create_from_Epetra_Matrix(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                          bool blockEntryMatrix,
                          fei::SharedPtr<fei::Reducer> reducer,
                          bool orderRowsWithLocalColsFirst)
{
  fei::SharedPtr<fei::VectorSpace> vecSpace = matrixGraph->getRowSpace();
  //localSize is in point-equations, because we will only use it for constructing
  //the fei::Matrix_Impl, and those always deal in point-equations.
  int localSize = vecSpace->getNumIndices_Owned();

  fei::SharedPtr<fei::Matrix> tmpmat;
  fei::SharedPtr<fei::Matrix> feimat;
  try {
    Epetra_CrsGraph egraph =
      Trilinos_Helpers::create_Epetra_CrsGraph(matrixGraph, blockEntryMatrix,
                                               orderRowsWithLocalColsFirst);

    if (blockEntryMatrix) {
      fei::SharedPtr<Epetra_VbrMatrix>
        epetraMatrix(new Epetra_VbrMatrix(Copy, egraph));

      tmpmat.reset(new fei::Matrix_Impl<Epetra_VbrMatrix>(epetraMatrix,
                                                   matrixGraph, localSize));
      zero_Epetra_VbrMatrix(epetraMatrix.get());
    }
    else {
      fei::SharedPtr<Epetra_CrsMatrix>
        epetraMatrix(new Epetra_CrsMatrix(Copy, egraph));

      tmpmat.reset(new fei::Matrix_Impl<Epetra_CrsMatrix>(epetraMatrix,
                                          matrixGraph, localSize));
    }
  }
  catch(std::runtime_error& exc) {
    fei::console_out() << "Trilinos_Helpers::create_from_Epetra_Matrix ERROR, "
           << "caught exception: '" << exc.what() << "', rethrowing..."
           << FEI_ENDL;
    throw exc;
  }

  if (reducer.get() != NULL) {
    feimat.reset(new fei::MatrixReducer(reducer, tmpmat));
  }
  else {
    feimat = tmpmat;
  }

  return(feimat);
}

fei::SharedPtr<fei::Matrix>
create_from_LPM_EpetraBasic(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                            bool blockEntryMatrix,
                            fei::SharedPtr<fei::Reducer> reducer,
                            fei::SharedPtr<fei::LinearProblemManager>
                              lpm_epetrabasic)
{
  fei::SharedPtr<fei::SparseRowGraph> srgraph =
    matrixGraph->createGraph(blockEntryMatrix);
  if (srgraph.get() == NULL) {
    throw std::runtime_error("create_LPM_EpetraBasic ERROR in fei::MatrixGraph::createGraph");
  }

  fei::SharedPtr<fei::SparseRowGraph> sharedsrg(srgraph);
  int localSize;
  if (reducer.get() != NULL) {
    std::vector<int>& reduced_eqns = reducer->getLocalReducedEqns();
    lpm_epetrabasic->setRowDistribution(reduced_eqns);
    lpm_epetrabasic->setMatrixGraph(sharedsrg);
    localSize = reduced_eqns.size();
  }
  else {
    fei::SharedPtr<fei::VectorSpace> vecSpace = matrixGraph->getRowSpace();
    int err = 0,chkNum;
    std::vector<int> indices;
    if (blockEntryMatrix) {
      localSize = vecSpace->getNumBlkIndices_Owned();
      indices.resize(localSize*2);
      err = vecSpace->getBlkIndices_Owned(localSize, &indices[0],
                                  &indices[localSize], chkNum);
    }
    else {
      localSize = vecSpace->getNumIndices_Owned();
      indices.resize(localSize);
      err = vecSpace->getIndices_Owned(indices);
    }
    if (err != 0) {
      throw std::runtime_error("Factory_Trilinos: createMatrix: error in vecSpace->getIndices_Owned");
    }

    lpm_epetrabasic->setRowDistribution(indices);
    lpm_epetrabasic->setMatrixGraph(sharedsrg);
  }

  fei::SharedPtr<fei::Matrix> tmpmat(new fei::Matrix_Impl<fei::LinearProblemManager>(lpm_epetrabasic, matrixGraph, localSize));

  fei::SharedPtr<fei::Matrix> feimat;

  if (reducer.get() != NULL) {
    feimat.reset(new fei::MatrixReducer(reducer, tmpmat));
  }
  else {
    feimat = tmpmat;
  }

  return(feimat);
}
#endif //HAVE_FEI_EPETRA

void copy_parameterset(const fei::ParameterSet& paramset,
                       Teuchos::ParameterList& paramlist)
{
  fei::ParameterSet::const_iterator
    iter = paramset.begin(),
    iter_end = paramset.end();

  for(; iter != iter_end; ++iter) {
    const fei::Param& param = *iter;
    fei::Param::ParamType ptype = param.getType();
    switch(ptype) {
    case fei::Param::STRING:
      paramlist.set(param.getName(), param.getStringValue());
      break;
    case fei::Param::DOUBLE:
      paramlist.set(param.getName(), param.getDoubleValue());
      break;
    case fei::Param::INT:
      paramlist.set(param.getName(), param.getIntValue());
      break;
    case fei::Param::BOOL:
      paramlist.set(param.getName(), param.getBoolValue());
      break;
    default:
      break;
    }
  }
}

void copy_parameterlist(const Teuchos::ParameterList& paramlist,
                        fei::ParameterSet& paramset)
{
  Teuchos::ParameterList::ConstIterator
    iter = paramlist.begin(),
    iter_end = paramlist.end();

  for(; iter != iter_end; ++iter) {
    const Teuchos::ParameterEntry& param = paramlist.entry(iter);
    if (param.isType<std::string>()) {
      paramset.add(fei::Param(paramlist.name(iter).c_str(), Teuchos::getValue<std::string>(param).c_str()));
    }
    else if (param.isType<double>()) {
      paramset.add(fei::Param(paramlist.name(iter).c_str(), Teuchos::getValue<double>(param)));
    }
    else if (param.isType<int>()) {
      paramset.add(fei::Param(paramlist.name(iter).c_str(), Teuchos::getValue<int>(param)));
    }
    else if (param.isType<bool>()) {
      paramset.add(fei::Param(paramlist.name(iter).c_str(), Teuchos::getValue<bool>(param)));
    }
  }
}

#ifdef HAVE_FEI_EPETRA

Epetra_MultiVector*
get_Epetra_MultiVector(fei::Vector* feivec, bool soln_vec)
{
  fei::Vector* vecptr = feivec;
  fei::VectorReducer* feireducer = dynamic_cast<fei::VectorReducer*>(feivec);
  if (feireducer != NULL) vecptr = feireducer->getTargetVector().get();

  fei::Vector_Impl<Epetra_MultiVector>* fei_epetra_vec =
    dynamic_cast<fei::Vector_Impl<Epetra_MultiVector>*>(vecptr);
  fei::Vector_Impl<fei::LinearProblemManager>* fei_lpm =
    dynamic_cast<fei::Vector_Impl<fei::LinearProblemManager>*>(vecptr);

  if (fei_epetra_vec == NULL && fei_lpm == NULL) {
    throw std::runtime_error("failed to obtain Epetra_MultiVector from fei::Vector.");
  }

  if (fei_epetra_vec != NULL) {
    return( fei_epetra_vec->getUnderlyingVector());
  }

  LinProbMgr_EpetraBasic* lpm_epetrabasic =
    dynamic_cast<LinProbMgr_EpetraBasic*>(fei_lpm->getUnderlyingVector());
  if (lpm_epetrabasic == 0) {
    throw std::runtime_error("fei Trilinos_Helpers: ERROR getting LinProbMgr_EpetraBasic");
  }

  if (soln_vec) {
    return(lpm_epetrabasic->get_solution_vector().get());
  }

  return(lpm_epetrabasic->get_rhs_vector().get());
}

Epetra_VbrMatrix* get_Epetra_VbrMatrix(fei::Matrix* feimat)
{
  fei::Matrix_Impl<Epetra_VbrMatrix>* fei_epetra_vbr =
    dynamic_cast<fei::Matrix_Impl<Epetra_VbrMatrix>*>(feimat);
  fei::MatrixReducer* feireducer =
    dynamic_cast<fei::MatrixReducer*>(feimat);

  if (feireducer != NULL) {
    fei::SharedPtr<fei::Matrix> feimat2 = feireducer->getTargetMatrix();
    fei_epetra_vbr =
      dynamic_cast<fei::Matrix_Impl<Epetra_VbrMatrix>*>(feimat2.get());
  }

  if (fei_epetra_vbr == NULL) {
    throw std::runtime_error("failed to obtain Epetra_VbrMatrix from fei::Matrix.");
  }

  return(fei_epetra_vbr->getMatrix().get());
}

Epetra_CrsMatrix* get_Epetra_CrsMatrix(fei::Matrix* feimat)
{
  fei::Matrix_Impl<Epetra_CrsMatrix>* fei_epetra_crs =
    dynamic_cast<fei::Matrix_Impl<Epetra_CrsMatrix>*>(feimat);
  fei::Matrix_Impl<fei::LinearProblemManager>* fei_lpm =
    dynamic_cast<fei::Matrix_Impl<fei::LinearProblemManager>*>(feimat);
  fei::MatrixReducer* feireducer =
    dynamic_cast<fei::MatrixReducer*>(feimat);

  if (feireducer != NULL) {
    fei::SharedPtr<fei::Matrix> feimat2 = feireducer->getTargetMatrix();
    fei_epetra_crs =
      dynamic_cast<fei::Matrix_Impl<Epetra_CrsMatrix>*>(feimat2.get());
    fei_lpm =
      dynamic_cast<fei::Matrix_Impl<fei::LinearProblemManager>*>(feimat2.get());
  }

  if (fei_epetra_crs == NULL && fei_lpm == NULL) {
    throw std::runtime_error("failed to obtain Epetra_CrsMatrix from fei::Matrix.");
  }

  if (fei_epetra_crs != NULL) {
    return(fei_epetra_crs->getMatrix().get());
  }

  LinProbMgr_EpetraBasic* lpm_epetrabasic =
    dynamic_cast<LinProbMgr_EpetraBasic*>(fei_lpm->getMatrix().get());
  if (lpm_epetrabasic == 0) {
    throw std::runtime_error("fei Trilinos_Helpers ERROR getting LinProbMgr_EpetraBasic");
  }

  return(lpm_epetrabasic->get_A_matrix().get());
}


void get_Epetra_pointers(fei::SharedPtr<fei::Matrix> feiA,
                         fei::SharedPtr<fei::Vector> feix,
                         fei::SharedPtr<fei::Vector> feib,
                         Epetra_CrsMatrix*& crsA,
                         Epetra_Operator*& opA,
                         Epetra_MultiVector*& x,
                         Epetra_MultiVector*& b)
{
  x = get_Epetra_MultiVector(feix.get(), true);
  b = get_Epetra_MultiVector(feib.get(), false);

  const char* matname = feiA->typeName();
  if (!strcmp(matname, "Epetra_VbrMatrix")) {
    Epetra_VbrMatrix* A = get_Epetra_VbrMatrix(feiA.get());
    opA = A;
  }
  else {
    crsA = get_Epetra_CrsMatrix(feiA.get());
    opA = crsA;
  }
}

int zero_Epetra_VbrMatrix(Epetra_VbrMatrix* mat)
{
  const Epetra_CrsGraph& crsgraph = mat->Graph();
  const Epetra_BlockMap& rowmap = crsgraph.RowMap();
  const Epetra_BlockMap& colmap = crsgraph.ColMap();
  mat->RowMatrixRowMap();//generate point objects
  int maxBlkRowSize = mat->GlobalMaxRowDim();
  int maxBlkColSize = mat->GlobalMaxColDim();
  std::vector<double> zeros(maxBlkRowSize*maxBlkColSize, 0);
  int numMyRows = rowmap.NumMyElements();
  int* myRows = rowmap.MyGlobalElements();
  for(int i=0; i<numMyRows; ++i) {
    int row = myRows[i];
    int rowlength = 0;
    int* colindicesView = NULL;
    int localrow = rowmap.LID(row);
    int err = crsgraph.ExtractMyRowView(localrow, rowlength, colindicesView);
    if (err != 0) {
      return err;
    }
    err = mat->BeginReplaceMyValues(localrow, rowlength, colindicesView);
    if (err != 0) {
      return err;
    }
    int blkRowSize = rowmap.ElementSize(localrow);
    for(int j=0; j<rowlength; ++j) {
      int blkColSize = colmap.ElementSize(colindicesView[j]);
      err = mat->SubmitBlockEntry(&zeros[0], maxBlkRowSize,
                            blkRowSize, blkColSize);
      if (err != 0) {
        return err;
      }
    }
    err = mat->EndSubmitEntries();
    if (err != 0) {
      return err;
    }
  }

  return 0;
}

#endif //HAVE_FEI_EPETRA

}//namespace Trilinos_Helpers


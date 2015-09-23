/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <cmath>
#include <fei_mpi.h>
#include <test_utils/fei_test_utils.hpp>
#include <test_utils/LibraryFactory.hpp>
#include <test_utils/test_Matrix.hpp>
#include <test_utils/test_VectorSpace.hpp>
#include <test_utils/test_MatrixGraph.hpp>
#include <fei_MatrixGraph_Impl2.hpp>
#include <fei_Factory.hpp>
#include <fei_defs.h>
#include <snl_fei_Factory.hpp>
#include <fei_Vector_Impl.hpp>
#include <fei_Matrix_Impl.hpp>

#ifdef HAVE_FEI_AZTECOO
#include <fei_Aztec_LinSysCore.hpp>
#endif
#include <fei_Factory_Trilinos.hpp>

#ifdef HAVE_FEI_FETI
#include <FETI_DP_FiniteElementData.h>
#endif

#undef fei_file
#define fei_file "test_Matrix.cpp"
#include <fei_ErrMacros.hpp>

int test_matrix_unit1()
{
  std::vector<double> data1D;
  std::vector<const double*>data2D;

  int numRows = 2;
  int numCols = 3;
  double** values = new double*[numRows];
  int i, j;
  for(i=0; i<numRows; ++i) {
    values[i] = new double[numCols];
    for(j=0; j<numCols; ++j) {
      values[i][j] = i*1.0;
    }
  }

  fei::Matrix_core::copyTransposeToWorkArrays(numRows, numCols,
					   values, data1D, data2D);

  for(i=0; i<numRows; ++i) {
    for(j=0; j<numCols; ++j) {
      if (std::abs(values[i][j] - (data2D[j])[i]) > 1.e-49) {
	ERReturn(-1);
      }
    }
    delete [] values[i];
  }

  delete [] values;

  return(0);
}

void test_Matrix_unit2(MPI_Comm comm, int numProcs, int localProc)
{
  if (numProcs > 1) {
    return;
  }

  FEI_COUT << "testing fei::Matrix_Impl...";

  fei::SharedPtr<fei::VectorSpace> rowspace(new fei::VectorSpace(comm));
  fei::SharedPtr<fei::VectorSpace> colspace;

  int rowfield = 0, rowfieldsize = 1;
  int idType = 0;
  rowspace->defineFields(1, &rowfield, &rowfieldsize);
  rowspace->defineIDTypes(1, &idType);

  fei::SharedPtr<fei::MatrixGraph> mgraph(new fei::MatrixGraph_Impl2(rowspace, colspace));

  int patternID1 = mgraph->definePattern(2, idType, rowfield);

  fei::Pattern* rowpattern = mgraph->getPattern(patternID1);

  mgraph->initConnectivityBlock(0, 1, patternID1);

  std::vector<int> ids(2);
  ids[0] = 0; ids[1] = 1;

  int err = mgraph->initConnectivity(0, 0, &ids[0]);
  if (err) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "test_Matrix_unit2, initConnectivity returned err="<<err;
    throw std::runtime_error(osstr.str());
  }

  err = mgraph->initComplete();
  if (err) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "test_Matrix_unit2, initComplete returned err="<<err;
    throw std::runtime_error(osstr.str());
  }

  bool factory_created = false;
  fei::SharedPtr<fei::Factory> factory;
  try {
    factory = fei::create_fei_Factory(comm, "Trilinos");
    factory_created = true;
  }
  catch(...) {}

  if (!factory_created) {
    FEI_COUT << "failed to create Trilinos factory."<<FEI_ENDL;
    return;
  }

  fei::SharedPtr<fei::Matrix> feimat = factory->createMatrix(mgraph);

  int numrowindices = rowpattern->getNumIndices();

  std::vector<double> coefs(numrowindices*numrowindices, 1.0);
  std::vector<double*> coefs_2D(numrowindices);
  for(int i=0; i<numrowindices; ++i) {
    coefs_2D[i] = &(coefs[i*numrowindices]);
  }

  err = feimat->sumIn(0, 0, &coefs_2D[0]);
  if (err) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "test_Matrix_unit2, feimat->sumIn returned err="<<err;
    throw std::runtime_error(osstr.str());
  }

  err = feimat->globalAssemble();
  if (err) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "test_Matrix_unit2, feimat->globalAssemble returned err="<<err;
    throw std::runtime_error(osstr.str());
  }

  err = feimat->writeToFile("feimat2.mtx", false);
  if (err) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "test_Matrix_unit2, feimat->writeToFile returned err="<<err;
    throw std::runtime_error(osstr.str());
  }

  fei::FillableMat feimat_ss;
  err = fei_test_utils::copy_feiMatrix_to_FillableMat(*feimat, feimat_ss);
  if (err) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "test_Matrix_unit2, copy_feiMatrix_to_FillableMat returned err="<<err;
    throw std::runtime_error(osstr.str());
  }

  fei_test_utils::writeMatrix("feimat_ss2.mtx", feimat_ss);

  FEI_COUT << "ok"<<FEI_ENDL;
}

void test_Matrix_unit4(MPI_Comm comm, int numProcs, int localProc)
{
  if (numProcs > 1) {
    return;
  }

  FEI_COUT << "testing fei::Matrix_Impl with FEI_BLOCK_DIAGONAL_ROW...";

  fei::SharedPtr<fei::VectorSpace> rowspace(new fei::VectorSpace(comm));
  fei::SharedPtr<fei::VectorSpace> colspace;

  int rowfield = 0, rowfieldsize = 2;
  int idType = 0;
  rowspace->defineFields(1, &rowfield, &rowfieldsize);
  rowspace->defineIDTypes(1, &idType);

  fei::SharedPtr<fei::MatrixGraph> mgraph(new fei::MatrixGraph_Impl2(rowspace, colspace));

  int patternID1 = mgraph->definePattern(2, idType, rowfield);

  fei::Pattern* rowpattern = mgraph->getPattern(patternID1);

  bool diagonal = true;
  mgraph->initConnectivityBlock(0, 1, patternID1, diagonal);

  std::vector<int> ids(2);
  ids[0] = 0; ids[1] = 1;

  int err = mgraph->initConnectivity(0, 0, &ids[0]);
  if (err) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "test_Matrix_unit4, initConnectivity returned err="<<err;
    throw std::runtime_error(osstr.str());
  }

  err = mgraph->initComplete();
  if (err) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "test_Matrix_unit4, initComplete returned err="<<err;
    throw std::runtime_error(osstr.str());
  }

  fei::SharedPtr<fei::Factory> factory;
  try {
    factory = fei::create_fei_Factory(comm, "Trilinos");
  }
  catch(...) {
    FEI_COUT << "Trilinos not available."<<FEI_ENDL;
    return;
  }

  fei::Param blktrue("BLOCK_MATRIX", true);
  fei::Param blkfalse("BLOCK_MATRIX", false);

  fei::ParameterSet paramset;
  paramset.add(blktrue);
  factory->parameters(paramset);

  fei::SharedPtr<fei::Matrix> feiblkmat = factory->createMatrix(mgraph);

  paramset.add(blkfalse);
  factory->parameters(paramset);

  fei::SharedPtr<fei::Matrix> feimat = factory->createMatrix(mgraph);

  int numrowindices = rowpattern->getNumIndices();

  std::vector<double> coefs(numrowindices*rowfieldsize*rowfieldsize, 1.0);
  std::vector<double*> coefs_2D(numrowindices*rowfieldsize);
  int offset = 0;
  for(int i=0; i<numrowindices*rowfieldsize; ++i) {
    coefs_2D[i] = &(coefs[offset]);
    offset += rowfieldsize;
  }

  err = feimat->sumIn(0, 0, &coefs_2D[0], FEI_BLOCK_DIAGONAL_ROW);
  if (err) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "test_Matrix_unit4, feimat->sumIn returned err="<<err;
    throw std::runtime_error(osstr.str());
  }

  err = feiblkmat->sumIn(0, 0, &coefs_2D[0], FEI_BLOCK_DIAGONAL_ROW);
  if (err) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "test_Matrix_unit4, feiblkmat->sumIn returned err="<<err;
    throw std::runtime_error(osstr.str());
  }

  err = feimat->globalAssemble();
  if (err) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "test_Matrix_unit4, feimat->globalAssemble returned err="<<err;
    throw std::runtime_error(osstr.str());
  }

  err = feiblkmat->globalAssemble();
  if (err) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "test_Matrix_unit4, feimat->globalAssemble returned err="<<err;
    throw std::runtime_error(osstr.str());
  }

  feimat->writeToFile("feimat_blkdiag.mtx");
  feiblkmat->writeToFile("feiblkmat_blkdiag.mtx");

  FEI_COUT << "ok"<<FEI_ENDL;
}

test_Matrix::test_Matrix(MPI_Comm comm)
  : tester(comm)
{
}

test_Matrix::~test_Matrix()
{
}

int test_Matrix::runtests()
{

  //-------------------------------
  // Test a Factory_Trilinos matrix.
  if (localProc_==0) FEI_COUT << "Testing Factory_Trilinos fei::Matrix..." << FEI_ENDL;
  fei::SharedPtr<fei::Factory> factory_trilinos(new Factory_Trilinos(comm_));

  fei::SharedPtr<fei::Matrix> mat = create_matrix(factory_trilinos);

  matrix_test1(mat);

  if (localProc_==0) FEI_COUT << FEI_ENDL;

#ifdef HAVE_FEI_AZTECOO
  //-------------------------------
  // Test a matrix from a "snl_fei::Factory", which was created by the
  // test-util function create_fei_Factory for library "Aztec", which causes
  // an "old" LinearSystemCore object to be used as the underlying assembly layer.
  // (It still ends up using Aztec if you create a solver and solve the linear-
  // system, but it assembles data into non-Trilinos matrix/vector data structures.)
  if (localProc_==0) FEI_COUT << "Testing 'old' factory_aztec fei::Matrix..." << FEI_ENDL;
  fei::SharedPtr<fei::Factory> factory_aztec =
    fei::create_fei_Factory(comm_, "Aztec");

  fei::SharedPtr<fei::Matrix> mat1 = create_matrix(factory_aztec);

  matrix_test1(mat1);
#endif

  return(0);
}

fei::SharedPtr<fei::Matrix>
test_Matrix::create_matrix(fei::SharedPtr<fei::Factory> factory)
{
  testData test_data(localProc_, numProcs_);

  fei::SharedPtr<fei::VectorSpace> vspace =
    test_VectorSpace::create_VectorSpace(comm_, &test_data, localProc_, numProcs_,
					 false, false, (const char*)0, factory);
  int err = vspace->initComplete();
  if (err != 0) {
    FEI_COUT << "ERROR, failed to create valid fei::VectorSpace." << FEI_ENDL;
    throw std::runtime_error("test_Vector::vector_test1: ERROR, failed to create valid fei::VectorSpace.");
  }

  fei::SharedPtr<fei::MatrixGraph> mgraph =
    factory->createMatrixGraph(vspace, vspace, NULL);

  std::vector<int>& fieldIDs = test_data.fieldIDs;
  std::vector<int>& idTypes = test_data.idTypes;
  std::vector<int>& ids = test_data.ids;

  int numIDs = ids.size();
  int fieldID = fieldIDs[0];
  int idType = idTypes[0];

  int patternID = mgraph->definePattern(numIDs, idType, fieldID);

  mgraph->initConnectivityBlock(0, 1, patternID);

  mgraph->initConnectivity(0, 0, &ids[0]);

  mgraph->initComplete();

  fei::SharedPtr<fei::Matrix> matrix = factory->createMatrix(mgraph);
  return(matrix);
}

void test_Matrix::matrix_test1(fei::SharedPtr<fei::Matrix> mat)
{
  if (localProc_==0)
    FEI_COUT << "  matrix_test1: testing fei::Matrix with type '"
	   << mat->typeName() << "':"<<FEI_ENDL;

  fei::SharedPtr<fei::MatrixGraph> mgraph = mat->getMatrixGraph();

  fei::SharedPtr<fei::VectorSpace> rspace = mgraph->getRowSpace();

  if (localProc_==0)
    FEI_COUT << "   testing get{Global/Local}NumRows,getRowLength...";

  int mglobalrows = mat->getGlobalNumRows();
  int vglobaleqns = rspace->getGlobalNumIndices();

  if (mglobalrows != vglobaleqns) {
    throw std::runtime_error("mat reports different num rows than vector-space eqns");
  }

  std::vector<int> global_offsets;
  rspace->getGlobalIndexOffsets(global_offsets);

  int my_num_rows = mat->getLocalNumRows();
  if (my_num_rows != global_offsets[localProc_+1]-global_offsets[localProc_]) {
    throw std::runtime_error("num-local-rows mis-match between mat and vector-space");
  }

  int i, my_first_row = global_offsets[localProc_];
  std::vector<int> row_lengths(my_num_rows);

  for(i=0; i<my_num_rows; ++i) {
    int errcode = mat->getRowLength(i+my_first_row, row_lengths[i]);
    if (errcode != 0) {
      throw std::runtime_error("nonzero errcode from mat->getRowLength");
    }
  }

  if (localProc_==0) FEI_COUT << "ok" << FEI_ENDL;
}

int test_Matrix::serialtest1()
{
  fei::SharedPtr<testData> testdata(new testData(localProc_, numProcs_));
  std::vector<int>& fieldIDs = testdata->fieldIDs;
  std::vector<int>& fieldSizes = testdata->fieldSizes;
  std::vector<int>& idTypes = testdata->idTypes;
  std::vector<int>& ids = testdata->ids;

  fei::SharedPtr<fei::VectorSpace> vspc(new fei::VectorSpace(comm_, "sU_Mat"));

  vspc->defineFields(fieldIDs.size(),
		     &fieldIDs[0],
		     &fieldSizes[0]);

  vspc->defineIDTypes(idTypes.size(), &idTypes[0]);

  fei::SharedPtr<fei::MatrixGraph>
    matgraph(new fei::MatrixGraph_Impl2(vspc, vspc, "sU_Mat"));

  int numIDs = ids.size();
  int fieldID = fieldIDs[0];
  int idType = idTypes[0];

  int patternID = matgraph->definePattern(numIDs, idType, fieldID);

  CHK_ERR( matgraph->initConnectivityBlock(0, 1, patternID) );

  CHK_ERR( matgraph->initConnectivity(0, 0, &ids[0]) );

  CHK_ERR( matgraph->initComplete() );

  fei::SharedPtr<fei::FillableMat> ssmat(new fei::FillableMat), ssmatT(new fei::FillableMat);
  int localsize = matgraph->getRowSpace()->getNumIndices_Owned();
  fei::SharedPtr<fei::Matrix> matrix(new fei::Matrix_Impl<fei::FillableMat>(ssmat, matgraph, localsize));

  fei::SharedPtr<fei::Matrix> matrixT(new fei::Matrix_Impl<fei::FillableMat>(ssmatT, matgraph, localsize));

  std::vector<int> indices(numIDs);
  CHK_ERR( matgraph->getConnectivityIndices(0, 0, numIDs, &indices[0], numIDs) );

  std::vector<double> data1(numIDs*numIDs);
  std::vector<double*> data2d(numIDs);

  int i;
  for(i=0; i<numIDs; ++i) {
    data2d[i] = &(data1[i*numIDs]);
  }

  for(i=0; i<numIDs*numIDs; ++i) {
    data1[i] = 1.0*i;
  }

  CHK_ERR( matrix->sumIn(numIDs, &indices[0], numIDs, &indices[0],
			 &data2d[0], 0) );

  CHK_ERR( matrix->sumIn(0, 0, &data2d[0], 0) );

  CHK_ERR( matrixT->sumIn(numIDs, &indices[0],
			 numIDs, &indices[0], &data2d[0], 3) );

  CHK_ERR( matrixT->sumIn(0, 0, &data2d[0], 3) );

  if (*ssmat != *ssmatT) {
    ERReturn(-1);
  }

  return(0);
}

int test_Matrix::serialtest2()
{
  testData* testdata = new testData(localProc_, numProcs_);
  std::vector<int>& idTypes = testdata->idTypes;
  std::vector<int>& ids = testdata->ids;

  fei::SharedPtr<fei::VectorSpace> vspc(new fei::VectorSpace(comm_, "sU_Mat"));

  vspc->defineIDTypes(idTypes.size(), &idTypes[0]);

  fei::SharedPtr<fei::MatrixGraph> matgraph(new fei::MatrixGraph_Impl2(vspc, vspc, "sU_Mat"));

  int numIDs = ids.size();
  int idType = idTypes[0];

  int patternID = matgraph->definePattern(numIDs, idType);

  CHK_ERR( matgraph->initConnectivityBlock(0, 1, patternID) );

  CHK_ERR( matgraph->initConnectivity(0, 0, &ids[0]) );

  CHK_ERR( matgraph->initComplete() );

  fei::SharedPtr<fei::FillableMat> ssmat(new fei::FillableMat), ssmatT(new fei::FillableMat);
  int localsize = matgraph->getRowSpace()->getNumIndices_Owned();
  fei::Matrix* matrix = new fei::Matrix_Impl<fei::FillableMat>(ssmat, matgraph, localsize);

  fei::Matrix* matrixT = new fei::Matrix_Impl<fei::FillableMat>(ssmatT, matgraph, localsize);

  std::vector<int> indices(numIDs);
  CHK_ERR( matgraph->getConnectivityIndices(0, 0, numIDs, &indices[0], numIDs) );

  std::vector<double> data1(numIDs*numIDs);
  std::vector<double*> data2d(numIDs);

  int i;
  for(i=0; i<numIDs; ++i) {
    data2d[i] = &(data1[i*numIDs]);
  }

  for(i=0; i<numIDs*numIDs; ++i) {
    data1[i] = 1.0*i;
  }

  CHK_ERR( matrix->sumIn(numIDs, &indices[0],
			 numIDs, &indices[0], &data2d[0], 0) );

  CHK_ERR( matrixT->sumIn(numIDs, &indices[0],
			 numIDs, &indices[0], &data2d[0], 3) );

  if (*ssmat != *ssmatT) {
    ERReturn(-1);
  }

  delete matrix;
  delete matrixT;
  delete testdata;

  return(0);
}

int test_Matrix::serialtest3()
{
  testData* testdata = new testData(localProc_, numProcs_);
  std::vector<int>& fieldIDs = testdata->fieldIDs;
  std::vector<int>& fieldSizes = testdata->fieldSizes;
  std::vector<int>& idTypes = testdata->idTypes;
  std::vector<int>& ids = testdata->ids;

  fei::SharedPtr<fei::VectorSpace> vspc(new fei::VectorSpace(comm_, "sU_Mat3"));

  vspc->defineFields(fieldIDs.size(), &fieldIDs[0], &fieldSizes[0]);

  vspc->defineIDTypes(idTypes.size(), &idTypes[0]);

  fei::SharedPtr<fei::MatrixGraph>
    matgraph(new fei::MatrixGraph_Impl2(vspc, vspc, "sU_Mat3"));

  int numIDs = ids.size();
  int fieldID = fieldIDs[0];
  int idType = idTypes[0];

  int patternID = matgraph->definePattern(numIDs, idType, fieldID);

  CHK_ERR( matgraph->initConnectivityBlock(0, 1, patternID) );

  CHK_ERR( matgraph->initConnectivity(0, 0, &ids[0]) );

  //set up a slave constraint that defines id 2, field 0 to be equal to
  //id 1, field 0.
  int offsetOfSlave = 1;
  int offsetIntoSlaveField = 0;
  std::vector<double> weights(2);
  weights[0] = 1.0;
  weights[1] = -1.0;
  double rhsValue = 0.0;
  std::vector<int> cr_idtypes(2, idTypes[0]);
  std::vector<int> cr_fieldIDs(2, fieldIDs[0]);

  CHK_ERR( matgraph->initSlaveConstraint(2, //numIDs
					&cr_idtypes[0],
					&ids[1],
					&cr_fieldIDs[0],
					offsetOfSlave,
					offsetIntoSlaveField,
					&weights[0],
					rhsValue) );

  CHK_ERR( matgraph->initComplete() );

  fei::SharedPtr<fei::FillableMat> ssmat(new fei::FillableMat);
  int localsize = matgraph->getRowSpace()->getNumIndices_Owned();
  localsize -= 1;//subtract the slave
  fei::Matrix* matrix = new fei::Matrix_Impl<fei::FillableMat>(ssmat, matgraph, localsize);

  if (matrix == NULL) {
    ERReturn(-1);
  }

  std::vector<int> indices(numIDs);
  CHK_ERR( matgraph->getConnectivityIndices(0, 0, numIDs,
					   &indices[0], numIDs) );

  std::vector<double> data1(numIDs*numIDs);
  std::vector<double*> data2d(numIDs);

  int i;
  for(i=0; i<numIDs; ++i) {
    data2d[i] = &(data1[i*numIDs]);
  }

  for(i=0; i<numIDs*numIDs; ++i) {
    data1[i] = 1.0*i;
  }

  CHK_ERR( matrix->sumIn(numIDs, &indices[0],
			 numIDs, &indices[0], &data2d[0], 0) );

  CHK_ERR( matrix->sumIn(0, 0, &data2d[0], 0) );

  delete matrix;
  delete testdata;

  return(0);
}

int test_Matrix::test1()
{
  return(0);
}

int test_Matrix::test2()
{
 return(0);
}

int test_Matrix::test3()
{
#ifdef HAVE_FEI_FETI
  testData* testdata = new testData(localProc_, numProcs_);
  std::vector<int>& idTypes = testdata->idTypes;
  std::vector<int>& ids = testdata->ids;

  fei::SharedPtr<FiniteElementData> fedata(new FETI_DP_FiniteElementData(comm_));

  std::string paramstr("debugOutput .");
  char* param = const_cast<char*>(paramstr.c_str());

  CHK_ERR( fedata->parameters(1, &param) );

  fei::SharedPtr<fei::Factory> factory(new snl_fei::Factory(fedata, idTypes[0]));

  fei::SharedPtr<fei::VectorSpace> vectorSpacePtr =
    test_VectorSpace::create_VectorSpace(comm_,
					 testdata, localProc_, numProcs_,
					 false, false, "U_FEMat", factory);

  fei::SharedPtr<fei::MatrixGraph> matrixGraphPtr =
    test_MatrixGraph::create_MatrixGraph(testdata, localProc_, numProcs_,
					 false, false, "U_FEMat", vectorSpacePtr,
					 factory);

  CHK_ERR( matrixGraphPtr->initComplete() );

  fei::SharedPtr<fei::Vector> vec_fed = factory->createVector(vectorSpacePtr);

  fei::SharedPtr<fei::Matrix> mat_fed = factory->createMatrix(matrixGraphPtr);

  fei::Matrix_Impl<FiniteElementData>* smat2 = 
    dynamic_cast<fei::Matrix_Impl<FiniteElementData>*>(mat_fed.get());
  if (smat2 == NULL) {
    ERReturn(-1);
  }

  int blockID=0;
  int numIndices = matrixGraphPtr->getConnectivityNumIndices(blockID);

  std::vector<int> indicesArray(numIndices);
  int* indicesPtr = &indicesArray[0];

  int checkNumIndices = 0;
  CHK_ERR( matrixGraphPtr->getConnectivityIndices(blockID, 0,
					     numIndices, indicesPtr,
					     checkNumIndices) );

  std::vector<double> data(ids.size(), 1.0);
  double* dptr = &data[0];
  std::vector<double*> coefPtrs(ids.size(), dptr);

  CHK_ERR( mat_fed->sumIn(blockID, 0, &coefPtrs[0]) );

  CHK_ERR( vec_fed->sumIn(blockID, 0, &data[0]) );

  CHK_ERR( mat_fed->gatherFromOverlap() );

  CHK_ERR( fedata->loadComplete() );


  delete testdata;

  MPI_Barrier(comm_);

#endif  //HAVE_FEI_FETI

  return(0);
}

int test_Matrix::test4()
{
  return(0);
}


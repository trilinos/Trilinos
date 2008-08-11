/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/


#include <fei_macros.hpp>

#include <test_utils/test_LinearSystem.hpp>
#include <test_utils/test_VectorSpace.hpp>
#include <test_utils/test_MatrixGraph.hpp>
#include <snl_fei_Factory.hpp>
#include <fei_Vector.hpp>
#include <fei_Matrix.hpp>
#include <snl_fei_LinearSystem_General.hpp>

#include <test_utils/LibraryFactory.hpp>

#ifdef FEI_HAVE_TRILINOS
#include <support-Trilinos/fei-aztec.hpp>
#endif

#undef fei_file
#define fei_file "test_LinearSystem.cpp"
#include <fei_ErrMacros.hpp>

test_LinearSystem::test_LinearSystem(MPI_Comm comm)
 : tester(comm)
{
}

test_LinearSystem::~test_LinearSystem()
{
}

int test_LS_test6(MPI_Comm comm, int numProcs, int localProc,
                  const std::string& path)
{
#ifdef FEI_HAVE_TRILINOS
  FEI_COUT << "testing LinearSystem::loadEssentialBCs with idType != 0...";

  fei::SharedPtr<testData> testdata(new testData(localProc, numProcs));
  std::vector<int>& fieldIDs = testdata->fieldIDs;
  std::vector<int>& fieldSizes=testdata->fieldSizes;
  std::vector<int>& idTypes = testdata->idTypes;
  std::vector<int>& ids = testdata->ids;

  idTypes[0] = 1;

  fei::SharedPtr<LinearSystemCore> az_lsc(new Aztec_LinSysCore(comm));

  FEI_OSTRINGSTREAM osstr;
  osstr << "debugOutput ";
  if (path_.empty()) osstr << ".";
  else osstr << path_;

  std::string str = osstr.str();
  const char* charptr = str.c_str();

  CHK_ERR( az_lsc->parameters(1, &charptr) );

  fei::SharedPtr<fei::Factory> factory(new snl_fei::Factory(comm, az_lsc));

  fei::SharedPtr<fei::VectorSpace> vectorSpacePtr =
    test_VectorSpace::create_VectorSpace(comm,
					 testdata.get(), localProc, numProcs,
					 true, true, "U_LS6", factory);

  fei::SharedPtr<fei::MatrixGraph> matrixGraphPtr =
    test_MatrixGraph::create_MatrixGraph(testdata.get(), localProc, numProcs,
				      true, false, "U_LS6", vectorSpacePtr,
				      factory, path);

  CHK_ERR( matrixGraphPtr->initComplete() );

  fei::SharedPtr<fei::Vector> vec_lsc = factory->createVector(vectorSpacePtr);

  fei::SharedPtr<fei::Vector> vec_lsc2 = factory->createVector(vectorSpacePtr, true);

  fei::SharedPtr<fei::Matrix> mat_lsc = factory->createMatrix(matrixGraphPtr);

  fei::SharedPtr<fei::LinearSystem> linsys = factory->createLinearSystem(matrixGraphPtr);
  linsys->setMatrix(mat_lsc);
  linsys->setSolutionVector(vec_lsc2);
  linsys->setRHS(vec_lsc);

  int blockID=0;
  int numIndices = matrixGraphPtr->getConnectivityNumIndices(blockID);

  feiArray<int> indicesArray(numIndices);
  int* indicesPtr = indicesArray.dataPtr();

  int checkNumIndices = 0;
  CHK_ERR( matrixGraphPtr->getConnectivityIndices(blockID, 0,
					     numIndices, indicesPtr,
					     checkNumIndices) );

  feiArray<int> bcids(2);
  feiArray<int> bcidTypes(2);
  if (localProc == 0) {
    bcids = ids[0];
  }
  else {
    bcids = ids[ids.size()-1];
  }

  bcidTypes = idTypes[1];

  feiArray<double> data1(numIndices);
  feiArray<double> data2(numIndices);
  data1 = 0.0;
  data2 = 0.0;

  data1[0] = 1.0; data1[1] = 0.0; data1[2] = 0.0;
  data2[0] = 0.0; data2[1] = 1.0; data2[2] = 0.0;

  feiArray<double*> coefPtrs(numIndices);
  feiArray<double*> bcdata(2);
  bcdata[0] = data1.dataPtr();
  bcdata[1] = data2.dataPtr();

  for(int ii=0; ii<numIndices; ++ii) coefPtrs[ii] = data1.dataPtr();

  CHK_ERR( mat_lsc->sumIn(numIndices, indicesPtr, numIndices, indicesPtr,
			  coefPtrs.dataPtr()) );

  CHK_ERR( vec_lsc->sumInFieldData(fieldIDs[0], idTypes[0],
				    ids.size(), &ids[0],
				    data1.dataPtr()) );

  CHK_ERR( linsys->loadEssentialBCs(1, bcids.dataPtr(), idTypes[0],
				   fieldIDs[1], fieldSizes[1],
				   bcdata.dataPtr(), bcdata.dataPtr()) );

  CHK_ERR( mat_lsc->gatherFromOverlap() );

  CHK_ERR( az_lsc->matrixLoadComplete() );

  CHK_ERR( linsys->loadComplete() );

  std::vector<int> bceqns;
  std::vector<double> bcvals;
  linsys->getEssentialBCs(bceqns, bcvals);
  if (bceqns.size() != 1) {
    ERReturn(-7);
  }

  CHK_ERR( linsys->setBCValuesOnVector(vec_lsc2.get()) );

  MPI_Barrier(comm);

  FEI_COUT << "ok"<<FEI_ENDL;
#endif
  return(0);
}

int test_LinearSystem::runtests()
{
  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  CHK_ERR( test5() );

  CHK_ERR( test_LS_test6(comm_, numProcs_, localProc_, path_) );

  return(0);
}

int test_LinearSystem::test1()
{
#ifdef FEI_HAVE_TRILINOS
  fei::SharedPtr<testData> testdata(new testData(localProc_, numProcs_));
  std::vector<int>& fieldIDs = testdata->fieldIDs;
  std::vector<int>& fieldSizes=testdata->fieldSizes;
  std::vector<int>& idTypes = testdata->idTypes;
  std::vector<int>& ids = testdata->ids;

  fei::SharedPtr<LinearSystemCore> az_lsc(new Aztec_LinSysCore(comm_));

  FEI_OSTRINGSTREAM osstr;
  osstr << "debugOutput ";
  if (path_.empty()) osstr << ".";
  else osstr << path_;

  std::string str = osstr.str();
  const char* charptr = str.c_str();

  CHK_ERR( az_lsc->parameters(1, &charptr) );

  fei::SharedPtr<fei::Factory> factory(new snl_fei::Factory(comm_, az_lsc));

  fei::SharedPtr<fei::VectorSpace> vectorSpacePtr =
    test_VectorSpace::create_VectorSpace(comm_,
					 testdata.get(), localProc_, numProcs_,
					 true, true, "U_LS", factory);

  fei::SharedPtr<fei::MatrixGraph> matrixGraphPtr =
    test_MatrixGraph::create_MatrixGraph(testdata.get(), localProc_, numProcs_,
				      true, false, "U_LS", vectorSpacePtr,
				      factory, path_);

  CHK_ERR( matrixGraphPtr->initComplete() );

  fei::SharedPtr<fei::Vector> vec_lsc = factory->createVector(vectorSpacePtr);

  fei::SharedPtr<fei::Vector> vec_lsc2 = factory->createVector(vectorSpacePtr, true);

  fei::SharedPtr<fei::Matrix> mat_lsc = factory->createMatrix(matrixGraphPtr);

  fei::SharedPtr<fei::LinearSystem> linsys = factory->createLinearSystem(matrixGraphPtr);
  linsys->setMatrix(mat_lsc);
  linsys->setSolutionVector(vec_lsc2);
  linsys->setRHS(vec_lsc);

  int blockID=0;
  int numIndices = matrixGraphPtr->getConnectivityNumIndices(blockID);

  feiArray<int> indicesArray(numIndices);
  int* indicesPtr = indicesArray.dataPtr();

  int checkNumIndices = 0;
  CHK_ERR( matrixGraphPtr->getConnectivityIndices(blockID, 0,
					     numIndices, indicesPtr,
					     checkNumIndices) );

  feiArray<int> bcids(2);
  feiArray<int> bcidTypes(2);
  bcids = ids[1];
  bcidTypes = idTypes[1];

  feiArray<double> data1(numIndices);
  feiArray<double> data2(numIndices);
  data1 = 0.0;
  data2 = 0.0;

  data1[0] = 1.0; data1[1] = 0.0; data1[2] = 0.0;
  data2[0] = 0.0; data2[1] = 1.0; data2[2] = 0.0;

  feiArray<double*> coefPtrs(numIndices);
  feiArray<double*> bcdata(2);
  bcdata[0] = data1.dataPtr();
  bcdata[1] = data2.dataPtr();

  for(int ii=0; ii<numIndices; ++ii) coefPtrs[ii] = data1.dataPtr();

  CHK_ERR( mat_lsc->sumIn(numIndices, indicesPtr, numIndices, indicesPtr,
			  coefPtrs.dataPtr()) );

  CHK_ERR( vec_lsc->sumInFieldData(fieldIDs[0], idTypes[0],
				    ids.size(), &ids[0],
				    data1.dataPtr()) );

  CHK_ERR( linsys->loadEssentialBCs(2, bcids.dataPtr(), idTypes[0],
				   fieldIDs[1], fieldSizes[1],
				   bcdata.dataPtr(), bcdata.dataPtr()) );

  CHK_ERR( mat_lsc->gatherFromOverlap() );

  CHK_ERR( az_lsc->matrixLoadComplete() );

  CHK_ERR( linsys->loadComplete() );

  CHK_ERR( linsys->setBCValuesOnVector(vec_lsc2.get()) );

  CHK_ERR( az_lsc->writeSystem("U_LS") );

  linsys->putAttribute("solnvec", vec_lsc2.get());
  fei::Vector* solnvec = NULL;
  CHK_ERR( linsys->getAttribute("solnvec", (void*&)solnvec) );

  CHK_ERR( solnvec->writeToFile("tls_solnvec") );

  MPI_Barrier(comm_);
#endif

  return(0);
}

int test_LinearSystem::test2()
{
#ifdef FEI_HAVE_TRILINOS
  fei::SharedPtr<testData> testdata(new testData(localProc_, numProcs_));
  std::vector<int>& fieldIDs = testdata->fieldIDs;
  std::vector<int>& idTypes = testdata->idTypes;
  std::vector<int>& ids = testdata->ids;

  fei::SharedPtr<LinearSystemCore> az_lsc(new Aztec_LinSysCore(comm_));

  char* param = new char[64];
  sprintf(param,"debugOutput .");

  CHK_ERR( az_lsc->parameters(1, &param) );
  delete [] param;

  fei::SharedPtr<fei::Factory> factory(new snl_fei::Factory(comm_, az_lsc));

  fei::SharedPtr<fei::VectorSpace> vectorSpacePtr =
    test_VectorSpace::create_VectorSpace(comm_,
					 testdata.get(), localProc_, numProcs_,
					 false, false, "U_LS2", factory);

  fei::SharedPtr<fei::MatrixGraph> matrixGraphPtr =
    test_MatrixGraph::create_MatrixGraph(testdata.get(), localProc_, numProcs_,
					 false, false, "U_LS2", vectorSpacePtr,
					 factory, path_);

  feiArray<int> crIDTypes(2);
  feiArray<int> crFieldIDs(2);
  crIDTypes[0] = idTypes[0]; crIDTypes[1] = idTypes[0];
  crFieldIDs[0] = fieldIDs[0]; crFieldIDs[1] = fieldIDs[0];

  CHK_ERR( matrixGraphPtr->initLagrangeConstraint(0, idTypes[1],
						  2, //numIDs
						  crIDTypes.dataPtr(),
						  &(ids[1]),
						  crFieldIDs.dataPtr()) );

  CHK_ERR( matrixGraphPtr->initComplete() );

  fei::SharedPtr<fei::Vector> vec_lsc = factory->createVector(vectorSpacePtr);

  fei::SharedPtr<fei::Vector> vec_lsc2 = factory->createVector(vectorSpacePtr, true);

  fei::SharedPtr<fei::Matrix> mat_lsc = factory->createMatrix(matrixGraphPtr);

  fei::SharedPtr<fei::LinearSystem> linsys = factory->createLinearSystem(matrixGraphPtr);
  linsys->setMatrix(mat_lsc);
  linsys->setSolutionVector(vec_lsc2);
  linsys->setRHS(vec_lsc);

  int blockID=0;
  int numIndices = matrixGraphPtr->getConnectivityNumIndices(blockID);

  feiArray<int> indicesArray(numIndices);
  int* indicesPtr = indicesArray.dataPtr();

  int checkNumIndices = 0;
  CHK_ERR( matrixGraphPtr->getConnectivityIndices(blockID, 0,
					     numIndices, indicesPtr,
					     checkNumIndices) );

  feiArray<double> data(ids.size());
  data = 1.0;
  double* dptr = data.dataPtr();
  feiArray<double*> coefPtrs(ids.size());
  feiArray<double> crdata(2);
  crdata[0] = 1.0;
  crdata[1] = -1.0;

  for(unsigned ii=0; ii<ids.size(); ++ii) coefPtrs[ii] = dptr;

  CHK_ERR( mat_lsc->sumIn(numIndices, indicesPtr, numIndices, indicesPtr,
			  coefPtrs.dataPtr()) );

  CHK_ERR( vec_lsc->sumInFieldData(fieldIDs[0], idTypes[0],
				    ids.size(), &ids[0],
				    data.dataPtr()) );

  CHK_ERR( linsys->loadLagrangeConstraint(0, crdata.dataPtr(), 0.0) );

  CHK_ERR( mat_lsc->gatherFromOverlap() );

  CHK_ERR( az_lsc->matrixLoadComplete() );

   CHK_ERR( linsys->loadComplete() );

  std::vector<int> crindices;
  linsys->getConstrainedEqns(crindices);
  if (crindices.size() != 2) {
    ERReturn(-7);
  }

  CHK_ERR( az_lsc->writeSystem("U_LS2") );

  MPI_Barrier(comm_);
#endif

  return(0);
}

int test_LinearSystem::test3()
{
#ifdef FEI_HAVE_TRILINOS
  fei::SharedPtr<testData> testdata(new testData(localProc_, numProcs_));
  std::vector<int>& fieldIDs = testdata->fieldIDs;
  std::vector<int>& idTypes = testdata->idTypes;
  std::vector<int>& ids = testdata->ids;

  fei::SharedPtr<LinearSystemCore> az_lsc(new Aztec_LinSysCore(comm_));

  char* param = new char[64];
  sprintf(param,"debugOutput .");

  CHK_ERR( az_lsc->parameters(1, &param) );

  fei::SharedPtr<fei::Factory> factory(new snl_fei::Factory(comm_, az_lsc));

  fei::SharedPtr<fei::VectorSpace> vectorSpacePtr =
    test_VectorSpace::create_VectorSpace(comm_,
					 testdata.get(), localProc_, numProcs_,
					 false, false, "U_LS3", factory);

  fei::SharedPtr<fei::MatrixGraph> matrixGraphPtr =
    test_MatrixGraph::create_MatrixGraph(testdata.get(), localProc_, numProcs_,
					 false, false, "U_LS3", vectorSpacePtr,
					 factory, path_);

  feiArray<int> crIDTypes(2);
  feiArray<int> crFieldIDs(2);
  crIDTypes[0] = idTypes[0]; crIDTypes[1] = idTypes[0];
  crFieldIDs[0] = fieldIDs[0]; crFieldIDs[1] = fieldIDs[0];

  CHK_ERR( matrixGraphPtr->initPenaltyConstraint(0, idTypes[1],
						  2, //numIDs
						  crIDTypes.dataPtr(),
						  &(ids[1]),
						  crFieldIDs.dataPtr()) );

  CHK_ERR( matrixGraphPtr->initComplete() );

  fei::SharedPtr<fei::Vector> vec_lsc = factory->createVector(vectorSpacePtr);

  fei::SharedPtr<fei::Vector> vec_lsc2 = factory->createVector(vectorSpacePtr, true);

  fei::SharedPtr<fei::Matrix> mat_lsc = factory->createMatrix(matrixGraphPtr);

  fei::SharedPtr<fei::LinearSystem> linsys = factory->createLinearSystem(matrixGraphPtr);
  CHK_ERR( linsys->parameters(1, &param) );
  delete [] param;

  linsys->setMatrix(mat_lsc);
  linsys->setSolutionVector(vec_lsc2);
  linsys->setRHS(vec_lsc);

  int blockID=0;
  int numIndices = matrixGraphPtr->getConnectivityNumIndices(blockID);

  feiArray<int> indicesArray(numIndices);
  int* indicesPtr = indicesArray.dataPtr();

  int checkNumIndices = 0;
  CHK_ERR( matrixGraphPtr->getConnectivityIndices(blockID, 0,
					     numIndices, indicesPtr,
					     checkNumIndices) );

  feiArray<double> data(ids.size());
  data = 1.0;
  double* dptr = data.dataPtr();
  feiArray<double*> coefPtrs(ids.size());
  feiArray<double> crdata(2);
  crdata[0] = 1.0;
  crdata[1] = -1.0;

  for(unsigned ii=0; ii<ids.size(); ++ii) coefPtrs[ii] = dptr;

  CHK_ERR( mat_lsc->sumIn(numIndices, indicesPtr, numIndices, indicesPtr,
			  coefPtrs.dataPtr()) );

  CHK_ERR( vec_lsc->sumInFieldData(fieldIDs[0], idTypes[0],
				    ids.size(), &ids[0],
				    data.dataPtr()) );

  CHK_ERR( linsys->loadPenaltyConstraint(0, crdata.dataPtr(), 100.0, 0.0) );

  CHK_ERR( mat_lsc->gatherFromOverlap() );

  CHK_ERR( az_lsc->matrixLoadComplete() );

   CHK_ERR( linsys->loadComplete() );

  CHK_ERR( az_lsc->writeSystem("U_LS3") );

  MPI_Barrier(comm_);
#endif

  return(0);
}

int test_LinearSystem::test4()
{
#ifdef FEI_HAVE_TRILINOS
  fei::SharedPtr<testData> testdata(new testData(localProc_, numProcs_));
  std::vector<int>& fieldIDs = testdata->fieldIDs;
  std::vector<int>& idTypes = testdata->idTypes;
  std::vector<int>& ids = testdata->ids;

  fei::SharedPtr<LinearSystemCore> az_lsc(new Aztec_LinSysCore(comm_));

  char* param = new char[64];
  sprintf(param,"debugOutput .");

  CHK_ERR( az_lsc->parameters(1, &param) );
  delete [] param;

  fei::SharedPtr<fei::Factory> factory(new snl_fei::Factory(comm_, az_lsc));

  fei::SharedPtr<fei::VectorSpace> vectorSpacePtr =
    test_VectorSpace::create_VectorSpace(comm_,
					 testdata.get(), localProc_, numProcs_,
					 false, false, "U_LS4", factory);

  fei::SharedPtr<fei::MatrixGraph> matrixGraphPtr =
    test_MatrixGraph::create_MatrixGraph(testdata.get(), localProc_, numProcs_,
					 false, false, "U_LS4", vectorSpacePtr,
					 factory, path_);

  feiArray<int> crIDTypes(2);
  feiArray<int> crFieldIDs(2);
  crIDTypes[0] = idTypes[0]; crIDTypes[1] = idTypes[0];
  crFieldIDs[0] = fieldIDs[0]; crFieldIDs[1] = fieldIDs[0];

  CHK_ERR( matrixGraphPtr->initLagrangeConstraint(0, idTypes[1],
						  2, //numIDs
						  crIDTypes.dataPtr(),
						  &(ids[1]),
						  crFieldIDs.dataPtr()) );

  CHK_ERR( matrixGraphPtr->initComplete() );

  fei::SharedPtr<fei::Vector> vec_lsc = factory->createVector(vectorSpacePtr);

  fei::SharedPtr<fei::Vector> vec_lsc2 = factory->createVector(vectorSpacePtr, true);

  fei::SharedPtr<fei::Matrix> mat_lsc = factory->createMatrix(matrixGraphPtr);

  fei::SharedPtr<fei::LinearSystem> linsys = factory->createLinearSystem(matrixGraphPtr);
  linsys->setMatrix(mat_lsc);
  linsys->setSolutionVector(vec_lsc2);
  linsys->setRHS(vec_lsc);

  int blockID=0;
  int numIndices = matrixGraphPtr->getConnectivityNumIndices(blockID);

  feiArray<int> indicesArray(numIndices);
  int* indicesPtr = indicesArray.dataPtr();

  int checkNumIndices = 0;
  CHK_ERR( matrixGraphPtr->getConnectivityIndices(blockID, 0,
					     numIndices, indicesPtr,
					     checkNumIndices) );

  feiArray<double> data(ids.size());
  data = 1.0;
  double* dptr = data.dataPtr();
  feiArray<double*> coefPtrs(ids.size());
  feiArray<double> crdata(2);
  crdata[0] = 1.0;
  crdata[1] = -1.0;

  for(unsigned ii=0; ii<ids.size(); ++ii) coefPtrs[ii] = dptr;

  CHK_ERR( mat_lsc->sumIn(numIndices, indicesPtr, numIndices, indicesPtr,
			  coefPtrs.dataPtr()) );

  CHK_ERR( vec_lsc->sumInFieldData(fieldIDs[0], idTypes[0],
				   ids.size(), &ids[0],
				   data.dataPtr()) );

  CHK_ERR( linsys->loadLagrangeConstraint(0, crdata.dataPtr(), 0.0) );

  CHK_ERR( mat_lsc->gatherFromOverlap() );

  CHK_ERR( az_lsc->matrixLoadComplete() );

  CHK_ERR( linsys->loadComplete() );

  CHK_ERR( az_lsc->writeSystem("U_LS4") );

  MPI_Barrier(comm_);
#endif

  return(0);
}

int test_LinearSystem::test5()
{
#ifdef FEI_HAVE_TRILINOS
  testData* testdata = new testData(localProc_, numProcs_);
  std::vector<int>& fieldIDs = testdata->fieldIDs;
  std::vector<int>& fieldSizes=testdata->fieldSizes;
  std::vector<int>& idTypes = testdata->idTypes;
  std::vector<int>& ids = testdata->ids;

  fei::SharedPtr<LinearSystemCore> az_lsc(new Aztec_LinSysCore(comm_));

  char* param = new char[64];
  sprintf(param,"debugOutput .");

  CHK_ERR( az_lsc->parameters(1, &param) );
  delete [] param;

  fei::SharedPtr<fei::Factory> factory(new snl_fei::Factory(comm_, az_lsc));

  fei::SharedPtr<fei::VectorSpace> vectorSpacePtr =
    test_VectorSpace::create_VectorSpace(comm_,
					 testdata, localProc_, numProcs_,
					 true, true, "U_LS5", factory);

  fei::SharedPtr<fei::MatrixGraph> matrixGraphPtr =
    test_MatrixGraph::create_MatrixGraph(testdata, localProc_, numProcs_,
					 true, false, "U_LS5", vectorSpacePtr,
					 factory, path_);

  CHK_ERR( matrixGraphPtr->initComplete() );

  fei::SharedPtr<fei::Vector> vec_lsc = factory->createVector(vectorSpacePtr);

  fei::SharedPtr<fei::Vector> vec_lsc2 = factory->createVector(vectorSpacePtr, true);

  fei::SharedPtr<fei::Matrix> mat_lsc = factory->createMatrix(matrixGraphPtr);

  fei::SharedPtr<fei::LinearSystem> linsys = factory->createLinearSystem(matrixGraphPtr);

  linsys->setMatrix(mat_lsc);
  linsys->setSolutionVector(vec_lsc2);
  linsys->setRHS(vec_lsc);

  int blockID=0;
  int numIndices = matrixGraphPtr->getConnectivityNumIndices(blockID);

  feiArray<int> indicesArray(numIndices);
  int* indicesPtr = indicesArray.dataPtr();

  int checkNumIndices = 0;
  CHK_ERR( matrixGraphPtr->getConnectivityIndices(blockID, 0,
					     numIndices, indicesPtr,
					     checkNumIndices) );

  feiArray<int> bcids(2);
  feiArray<int> bcidTypes(2);
  bcids = ids[1];
  bcidTypes = idTypes[1];

  feiArray<double> data1(numIndices);
  feiArray<double> data2(numIndices);
  data1 = 0.0;
  data2 = 0.0;

  data1[0] = 1.0; data1[1] = 0.0; data1[2] = 0.0;
  data2[0] = 0.0; data2[1] = 1.0; data2[2] = 0.0;

  feiArray<double*> coefPtrs(numIndices);
  feiArray<double*> bcdata(2);
  bcdata[0] = data1.dataPtr();
  bcdata[1] = data2.dataPtr();

  for(int ii=0; ii<numIndices; ++ii) coefPtrs[ii] = data1.dataPtr();

  CHK_ERR( mat_lsc->sumIn(numIndices, indicesPtr, numIndices, indicesPtr,
			  coefPtrs.dataPtr()) );

  CHK_ERR( vec_lsc->sumInFieldData(fieldIDs[0], idTypes[0],
				    ids.size(), &ids[0],
				    data1.dataPtr()) );

  CHK_ERR( linsys->loadEssentialBCs(2, bcids.dataPtr(), idTypes[0],
				   fieldIDs[1], fieldSizes[1],
				   bcdata.dataPtr(), bcdata.dataPtr()) );

  CHK_ERR( mat_lsc->gatherFromOverlap() );

  CHK_ERR( az_lsc->matrixLoadComplete() );

  CHK_ERR( linsys->loadComplete() );

  CHK_ERR( linsys->setBCValuesOnVector(vec_lsc2.get()) );

  CHK_ERR( az_lsc->writeSystem("U_LS5") );

  delete testdata;

  MPI_Barrier(comm_);
#endif

  return(0);
}


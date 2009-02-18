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

#ifdef HAVE_FEI_AZTECOO
#include <fei_Aztec_LinSysCore.hpp>
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

int test_LinearSystem::runtests()
{
  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  CHK_ERR( test5() );

  return(0);
}

int test_LinearSystem::test1()
{
  return(0);
}

int test_LinearSystem::test2()
{
#ifdef HAVE_FEI_AZTECOO
  fei::SharedPtr<testData> testdata(new testData(localProc_, numProcs_));
  std::vector<int>& fieldIDs = testdata->fieldIDs;
  std::vector<int>& idTypes = testdata->idTypes;
  std::vector<int>& ids = testdata->ids;

  fei::SharedPtr<LinearSystemCore> az_lsc(new fei_trilinos::Aztec_LinSysCore(comm_));

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

  std::vector<int> crIDTypes(2);
  std::vector<int> crFieldIDs(2);
  crIDTypes[0] = idTypes[0]; crIDTypes[1] = idTypes[0];
  crFieldIDs[0] = fieldIDs[0]; crFieldIDs[1] = fieldIDs[0];

  CHK_ERR( matrixGraphPtr->initLagrangeConstraint(0, idTypes[1],
						  2, //numIDs
						  &crIDTypes[0],
						  &(ids[1]),
						  &crFieldIDs[0]) );

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

  std::vector<int> indicesArray(numIndices);
  int* indicesPtr = &indicesArray[0];

  int checkNumIndices = 0;
  CHK_ERR( matrixGraphPtr->getConnectivityIndices(blockID, 0,
					     numIndices, indicesPtr,
					     checkNumIndices) );

  std::vector<double> data(ids.size(), 1.0);
  double* dptr = &data[0];
  std::vector<double*> coefPtrs(ids.size());
  std::vector<double> crdata(2);
  crdata[0] = 1.0;
  crdata[1] = -1.0;

  for(unsigned ii=0; ii<ids.size(); ++ii) coefPtrs[ii] = dptr;

  CHK_ERR( mat_lsc->sumIn(numIndices, indicesPtr, numIndices, indicesPtr,
			  &coefPtrs[0]) );

  CHK_ERR( vec_lsc->sumInFieldData(fieldIDs[0], idTypes[0],
				    ids.size(), &ids[0],
				    &data[0]) );

  CHK_ERR( linsys->loadLagrangeConstraint(0, &crdata[0], 0.0) );

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
#ifdef HAVE_FEI_AZTECOO
  fei::SharedPtr<testData> testdata(new testData(localProc_, numProcs_));
  std::vector<int>& fieldIDs = testdata->fieldIDs;
  std::vector<int>& idTypes = testdata->idTypes;
  std::vector<int>& ids = testdata->ids;

  fei::SharedPtr<LinearSystemCore> az_lsc(new fei_trilinos::Aztec_LinSysCore(comm_));

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

  std::vector<int> crIDTypes(2);
  std::vector<int> crFieldIDs(2);
  crIDTypes[0] = idTypes[0]; crIDTypes[1] = idTypes[0];
  crFieldIDs[0] = fieldIDs[0]; crFieldIDs[1] = fieldIDs[0];

  CHK_ERR( matrixGraphPtr->initPenaltyConstraint(0, idTypes[1],
						  2, //numIDs
						  &crIDTypes[0],
						  &(ids[1]),
						  &crFieldIDs[0]) );

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

  std::vector<int> indicesArray(numIndices);
  int* indicesPtr = &indicesArray[0];

  int checkNumIndices = 0;
  CHK_ERR( matrixGraphPtr->getConnectivityIndices(blockID, 0,
					     numIndices, indicesPtr,
					     checkNumIndices) );

  std::vector<double> data(ids.size(), 1.0);
  double* dptr = &data[0];
  std::vector<double*> coefPtrs(ids.size());
  std::vector<double> crdata(2);
  crdata[0] = 1.0;
  crdata[1] = -1.0;

  for(unsigned ii=0; ii<ids.size(); ++ii) coefPtrs[ii] = dptr;

  CHK_ERR( mat_lsc->sumIn(numIndices, indicesPtr, numIndices, indicesPtr,
			  &coefPtrs[0]) );

  CHK_ERR( vec_lsc->sumInFieldData(fieldIDs[0], idTypes[0],
				    ids.size(), &ids[0],
				    &data[0]) );

  CHK_ERR( linsys->loadPenaltyConstraint(0, &crdata[0], 100.0, 0.0) );

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
#ifdef HAVE_FEI_AZTECOO
  fei::SharedPtr<testData> testdata(new testData(localProc_, numProcs_));
  std::vector<int>& fieldIDs = testdata->fieldIDs;
  std::vector<int>& idTypes = testdata->idTypes;
  std::vector<int>& ids = testdata->ids;

  fei::SharedPtr<LinearSystemCore> az_lsc(new fei_trilinos::Aztec_LinSysCore(comm_));

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

  std::vector<int> crIDTypes(2);
  std::vector<int> crFieldIDs(2);
  crIDTypes[0] = idTypes[0]; crIDTypes[1] = idTypes[0];
  crFieldIDs[0] = fieldIDs[0]; crFieldIDs[1] = fieldIDs[0];

  CHK_ERR( matrixGraphPtr->initLagrangeConstraint(0, idTypes[1],
						  2, //numIDs
						  &crIDTypes[0],
						  &(ids[1]),
						  &crFieldIDs[0]) );

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

  std::vector<int> indicesArray(numIndices);
  int* indicesPtr = &indicesArray[0];

  int checkNumIndices = 0;
  CHK_ERR( matrixGraphPtr->getConnectivityIndices(blockID, 0,
					     numIndices, indicesPtr,
					     checkNumIndices) );

  std::vector<double> data(ids.size(), 1.0);
  double* dptr = &data[0];
  std::vector<double*> coefPtrs(ids.size());
  std::vector<double> crdata(2);
  crdata[0] = 1.0;
  crdata[1] = -1.0;

  for(unsigned ii=0; ii<ids.size(); ++ii) coefPtrs[ii] = dptr;

  CHK_ERR( mat_lsc->sumIn(numIndices, indicesPtr, numIndices, indicesPtr,
			  &coefPtrs[0]) );

  CHK_ERR( vec_lsc->sumInFieldData(fieldIDs[0], idTypes[0],
				   ids.size(), &ids[0], &data[0]) );

  CHK_ERR( linsys->loadLagrangeConstraint(0, &crdata[0], 0.0) );

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
  return(0);
}


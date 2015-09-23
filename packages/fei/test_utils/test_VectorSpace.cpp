/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/


#include <fei_macros.hpp>

#include <test_utils/test_VectorSpace.hpp>

#include <snl_fei_Factory.hpp>
#include <fei_ParameterSet.hpp>
#include <fei_LibraryWrapper.hpp>
#include <assert.h>

#undef fei_file
#define fei_file "test_VectorSpace.cpp"
#include <fei_ErrMacros.hpp>

test_VectorSpace::test_VectorSpace(MPI_Comm comm)
 : tester(comm)
{
}

test_VectorSpace::~test_VectorSpace()
{
}

int test_VectorSpace::runtests()
{
  CHK_ERR( test0() );
  CHK_ERR( test1() );
  CHK_ERR( test2() );
  CHK_ERR( test3() );
  CHK_ERR( test4() );
  return(0);
}

int test_VectorSpace::test0()
{
  int i, ID = 0;
  int idType = 0;
  std::vector<int> fieldIDs(numProcs_);
  int fieldSize = 2;
  std::vector<int> fieldSizes(numProcs_, fieldSize);

  for(i=0; i<numProcs_; ++i) fieldIDs[i] = i;

  fei::VectorSpace vspace(comm_);

  vspace.defineIDTypes(1, &idType);
  vspace.defineFields(fieldIDs.size(), &fieldIDs[0], &fieldSizes[0]);

  CHK_ERR( vspace.addDOFs(idType, 1, &ID) );

  for(i=0; i<numProcs_; ++i) {
    if (i == localProc_) continue;

    int numSharingProcsPerID = 1;
    int sharingProc = i;

    CHK_ERR( vspace.initSharedIDs(1, idType, &ID,
				  &numSharingProcsPerID, &sharingProc) );
  }

  CHK_ERR( vspace.addDOFs(fieldIDs[localProc_], idType,
				      1, &ID) );

  CHK_ERR( vspace.initComplete() );

  std::vector<int> globalIndexOffsets;

  vspace.getGlobalIndexOffsets(globalIndexOffsets);

  return(0);
}

int test_VectorSpace::test1()
{
  //A general test of fei::VectorSpace, including defining fields and
  //idTypes, initializing solution entries, shared ids if numProcs_ > 1, then
  //testing getGlobalIndexOffsets and getGlobalIndex several times.

  testData* testdata = new testData(localProc_, numProcs_);
  std::vector<int>& fieldIDs = testdata->fieldIDs;
  std::vector<int>& idTypes = testdata->idTypes;
  std::vector<int>& ids = testdata->ids;

  int numDOFsPerID = 4; //1 scalar field and 1 3-vector field

  fei::SharedPtr<LibraryWrapper> wrapper;
  fei::SharedPtr<fei::Factory> factory(new snl_fei::Factory(comm_, wrapper));

  fei::SharedPtr<fei::VectorSpace> vectorSpacePtr =
    test_VectorSpace::create_VectorSpace(comm_,
                                         testdata, localProc_, numProcs_,
					 true, //defineBothFields
					 true, //initSolnBothFields
					 "U", //name
					 factory);

  CHK_ERR( vectorSpacePtr->initComplete() );

  fei::VectorSpace& vectorSpace = *vectorSpacePtr;

  std::vector<int> globalOffsets;

  vectorSpace.getGlobalIndexOffsets(globalOffsets);

  if (localProc_ > 0) {
    if (globalOffsets[localProc_] != (localProc_+1)*2*numDOFsPerID) {
      ERReturn(-1);
    }
  }

  if (globalOffsets[localProc_+1] != (localProc_+2)*2*numDOFsPerID) {
    ERReturn(-1);
  }

  std::vector<int> globalBlkOffsets;
  vectorSpace.getGlobalBlkIndexOffsets(globalBlkOffsets);
  if (localProc_ > 0) {
    if (globalBlkOffsets[localProc_] != (localProc_+1)*2) {
      ERReturn(-1);
    }
  }

  if (globalBlkOffsets[localProc_+1] != (localProc_+2)*2) {
    ERReturn(-1);
  }

  int len = ids.size();
  int globalIndex = 0;
  CHK_ERR( vectorSpace.getGlobalIndex(idTypes[0], ids[0], fieldIDs[0],
				    0, 0, globalIndex) );

  int correctIndex = globalOffsets[localProc_];
  int correctBlkIndex = globalBlkOffsets[localProc_];
  if (localProc_ != 0) correctIndex -= 2*numDOFsPerID;
  if (localProc_ != 0) correctBlkIndex -= 2;

  if (globalIndex != correctIndex) {
    ERReturn(-1);
  }

  int globalBlkIndex = 0;
  CHK_ERR( vectorSpace.getGlobalBlkIndex(idTypes[0], ids[0], globalBlkIndex) );
  if (globalBlkIndex != correctBlkIndex) {
    fei::console_out() << "localProc_: " << localProc_ << ", globalBlkIndex: "
	 << globalBlkIndex << ", correctBlkIndex: " << correctBlkIndex << FEI_ENDL;
    ERReturn(-1);
  }

  CHK_ERR( vectorSpace.getGlobalIndex(idTypes[0],  ids[1], fieldIDs[0],
				    0, 0, globalIndex) );

  correctIndex = globalOffsets[localProc_] + 4;
  if (localProc_ != 0) correctIndex -= 2*numDOFsPerID;

  if (globalIndex != correctIndex) {
    ERReturn(-1);
  }

  CHK_ERR( vectorSpace.getGlobalBlkIndex(idTypes[0], ids[1], globalBlkIndex) );
  correctBlkIndex = globalBlkOffsets[localProc_]+1;
  if (localProc_ != 0) correctBlkIndex -= 2;

  if (globalBlkIndex != correctBlkIndex) {
   fei::console_out() << "2localProc_: " << localProc_ << ", globalBlkIndex: "
	 << globalBlkIndex << ", correctBlkIndex: " << correctBlkIndex << FEI_ENDL;
    ERReturn(-1); 
  }

  CHK_ERR( vectorSpace.getGlobalIndex(idTypes[0], ids[len-1], fieldIDs[0],
				    0, 0, globalIndex) );
  correctIndex = globalOffsets[localProc_] + 12;
  if (localProc_ != 0) correctIndex -= 2*numDOFsPerID;

  if (globalIndex != correctIndex) {
    ERReturn(-1);
  }

  CHK_ERR( vectorSpace.getGlobalIndex(idTypes[0], ids[0], fieldIDs[1],
				    0, 0, globalIndex) );
  correctIndex = globalOffsets[localProc_] + 1;
  if (localProc_ != 0) correctIndex -= 2*numDOFsPerID;

  if (globalIndex != correctIndex) {
    ERReturn(-1);
  }

  CHK_ERR( vectorSpace.getGlobalIndex(idTypes[0], ids[1], fieldIDs[1],
				    0, 0, globalIndex) );
  correctIndex = globalOffsets[localProc_] + 5;
  if (localProc_ != 0) correctIndex -= 2*numDOFsPerID;

  if (globalIndex != correctIndex) {
    ERReturn(-1);
  }

  CHK_ERR( vectorSpace.getGlobalIndex(idTypes[0], ids[len-1], fieldIDs[1],
				    0, 0, globalIndex) );
  correctIndex = globalOffsets[localProc_] + 13;
  if (localProc_ != 0) correctIndex -= 2*numDOFsPerID;

  if (globalIndex != correctIndex) {
    ERReturn(-1);
  }

  std::vector<int> globalIndices(ids.size()*numDOFsPerID);

  CHK_ERR( vectorSpace.getGlobalIndices(ids.size(),
					&(ids[0]),
					idTypes[0], fieldIDs[0],
					&globalIndices[0] ));

  std::vector<int> idFieldIDs(ids.size(), fieldIDs[1]);
  std::vector<int> idIDTypes(ids.size(), idTypes[0]);

  CHK_ERR( vectorSpace.getGlobalIndices(ids.size(),
					&ids[0],
					idTypes[0], fieldIDs[0],
					&globalIndices[0] ));

  CHK_ERR( vectorSpace.getGlobalBlkIndices(ids.size(),
					   &ids[0],
					   idTypes[0],
					   &globalIndices[0] ));

  CHK_ERR( vectorSpace.getGlobalIndices(ids.size(),
					&ids[0],
					&idIDTypes[0],
					&idFieldIDs[0],
					&globalIndices[0]) );

  unsigned numFields = vectorSpace.getNumFields();
  if (numFields != fieldIDs.size()) {
    ERReturn(-1);
  }

  std::vector<int> testgetfields;
  vectorSpace.getFields(testgetfields);
  if (testgetfields.size() != numFields) {
    ERReturn(-1);
  }

  if (fieldIDs != testgetfields) {
    ERReturn(-1);
  }

  delete testdata;

  return(0);
}

int test_VectorSpace::test2()
{
  //only run this test if asserts are enabled (i.e., if NDEBUG is not defined)
  bool run_this_test = false;
  assert( (run_this_test = true) == true);
  if (!run_this_test) return(0);

  //Initializes a fei::VectorSpace object using the same data as used
  //above, but then adds a globally inconsistent shared-id
  //before calling initComplete. The goal is to make sure that the shared-id
  //safety-check in fei::VectorSpace will catch the incorrect shared-id
  //data and return an error-code.

  testData* testdata = new testData(localProc_, numProcs_);
  std::vector<int>& idTypes = testdata->idTypes;
  std::vector<int>& ids = testdata->ids;

  fei::SharedPtr<LibraryWrapper> wrapper;
  fei::SharedPtr<fei::Factory> factory(new snl_fei::Factory(comm_, wrapper));

  fei::SharedPtr<fei::VectorSpace> vectorSpacePtr =
    test_VectorSpace::create_VectorSpace(comm_,
                                         testdata, localProc_,
					 numProcs_,
					 true, //defineBothFields
					 true, //initSolnBothFields
					 "U2", factory);

  if (localProc_ < numProcs_-1) {
    int numSharingProcsPerID = 1;
    int sharingProc = localProc_+1;

    CHK_ERR( vectorSpacePtr->initSharedIDs(1, idTypes[0], &ids[0],
				     &numSharingProcsPerID,
				     &sharingProc) );
  }

  fei::ParameterSet paramset;
  paramset.add(fei::Param("FEI_CHECK_SHARED_IDS", true));

  vectorSpacePtr->setParameters(paramset);

  //we just provided globally inconsistent shared-id info, so we expect the
  //initComplete method to return(-1) after the shared-id safety-check fails,
  //if we're running on more than one processor.
  int err = vectorSpacePtr->initComplete();
  if (numProcs_ > 1) if (err == 0) return(-1);
  delete testdata;

  return(0);
}

int test_VectorSpace::test3()
{
  testData* testdata = new testData(localProc_, numProcs_);

  fei::SharedPtr<LibraryWrapper> wrapper;
  fei::SharedPtr<fei::Factory> factory(new snl_fei::Factory(comm_, wrapper));

  fei::SharedPtr<fei::VectorSpace> vectorSpacePtr =
    test_VectorSpace::create_VectorSpace(comm_,
                                         testdata, localProc_, numProcs_,
					 true, //defineBothFields
					 true, //initSolnBothFields
					 "U3", factory);

  fei::SharedPtr<fei::VectorSpace> copy =
    factory->createVectorSpace(comm_, "U3copy");

  fei::Param param("debugOutput", ".");
  fei::ParameterSet paramset;
  paramset.add(param);
  copy->setParameters(paramset);


  CHK_ERR( copy->addVectorSpace(vectorSpacePtr.get()) );

  CHK_ERR( vectorSpacePtr->initComplete() );

  CHK_ERR( copy->initComplete() );

  std::vector<int> globalOffsets;
  std::vector<int> globalOffsetsCopy;

  vectorSpacePtr->getGlobalIndexOffsets(globalOffsets);

  copy->getGlobalIndexOffsets(globalOffsetsCopy);

  for(size_t i=0; i<globalOffsets.size(); ++i) {
    if (globalOffsets[i] != globalOffsetsCopy[i]) {
      ERReturn(-1);
    }
  }

  CHK_ERR( copy->initComplete() );

  copy->getGlobalIndexOffsets(globalOffsetsCopy);

  if (globalOffsetsCopy[numProcs_] != globalOffsets[numProcs_]) {
    ERReturn(-1);
  }

  delete testdata;

  return(0);
}

int test_VectorSpace::test4()
{
  return(0);
}

fei::SharedPtr<fei::VectorSpace>
test_VectorSpace::create_VectorSpace(MPI_Comm comm)
{
  int localProc = 0, numProcs = 1;
#ifndef FEI_SER
  MPI_Comm_rank(comm, &localProc);
  MPI_Comm_size(comm, &numProcs);
#endif

  testData test_data(localProc, numProcs);
  fei::SharedPtr<fei::Factory> factory;

  fei::SharedPtr<fei::VectorSpace> vspace =
    test_VectorSpace::create_VectorSpace(comm, &test_data, localProc, numProcs,
					 false, false, (const char*)0, factory);
  int err = vspace->initComplete();
  if (err != 0) {
    FEI_COUT << "ERROR, failed to create valid fei::VectorSpace." << FEI_ENDL;
    throw std::runtime_error("test_Vector::vector_test1: ERROR, failed to create valid fei::VectorSpace.");
  }

  return(vspace);
}

fei::SharedPtr<fei::VectorSpace>
test_VectorSpace::create_VectorSpace(MPI_Comm comm,
				     testData* testdata,
				     int localProc,
				     int numProcs,
				     bool defineBothFields,
				     bool initSolnBothFields,
				     const char* name,
				     fei::SharedPtr<fei::Factory> factory,
				     bool turnOnDebugOutput)
{
  //
  //This function creates a VectorSpace object, then initializes it as follows:
  //
  //defineFields testdata->fieldIDs, testdata->fieldSizes
  //defineIDTypes testdata->idTypes
  //addDOFs testdata->ids with associated testdata->fieldIDs[0]
  //addDOFs testdata->ids with testdata->fieldIDs[1] if bothFields
  //initSharedIDs  testdata->[shared-id-data] with testdata->idTypes[0]
  //
  fei::SharedPtr<fei::VectorSpace> vsptr;

  if (factory.get() == NULL) {
    vsptr.reset(new fei::VectorSpace(comm, name));
  }
  else {
    vsptr = factory->createVectorSpace(comm, name);
  }

  fei::VectorSpace& vectorSpace = *vsptr;

  std::vector<int>& fieldIDs = testdata->fieldIDs;
  std::vector<int>& fieldSizes = testdata->fieldSizes;
  std::vector<int>& idTypes = testdata->idTypes;
  std::vector<int>& ids = testdata->ids;
  std::vector<int>& sharedIDs = testdata->sharedIDs;
  std::vector<int>& numSharingProcsPerID = testdata->numSharingProcsPerID;
  std::vector<int>& sharingProcs = testdata->sharingProcs;

  fei::ParameterSet paramset;
  fei::Param param1("name", name);
  paramset.add(param1);
  if (turnOnDebugOutput) {
    fei::Param param2("debugOutput", ".");
    paramset.add(param2);
  }

  vectorSpace.setParameters(paramset);

  int numFields = defineBothFields ? 2 : 1;
  vectorSpace.defineFields(numFields,
			   &fieldIDs[0],
			   &fieldSizes[0]);

  vectorSpace.defineIDTypes(idTypes.size(),
			    &idTypes[0]);

  vectorSpace.addDOFs(fieldIDs[0],
				  idTypes[0],
				  ids.size(),
				  &ids[0]);

  if (initSolnBothFields) {
    vectorSpace.addDOFs(fieldIDs[1],
				    idTypes[0],
				    ids.size(),
				    &ids[0]);
  }

  vectorSpace.initSharedIDs(sharedIDs.size(),
    idTypes[0],
    sharedIDs.size() ? &sharedIDs[0] : 0,
    numSharingProcsPerID.size() ? &numSharingProcsPerID[0] : 0,
    sharingProcs.size() ? &sharingProcs[0] : 0
    );

  return(vsptr);
}

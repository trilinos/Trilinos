
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <fei_mpi.h>
#include <fei_iostream.hpp>
#include <fei_ReverseMapper.hpp>
#include <fei_VectorSpace.hpp>

#include <vector>
#include <cmath>

namespace {

TEUCHOS_UNIT_TEST(ReverseMapper, test1)
{
  MPI_Comm comm = MPI_COMM_WORLD;

  int numProcs = 1;
#ifndef FEI_SER
  MPI_Comm_size(comm, &numProcs);
#endif

  if (numProcs > 1) return;

  fei::VectorSpace vspace(comm);

  int fieldIDs[2] = {1, 2};
  int fieldSizes[2] = {1, 2};

  vspace.defineFields(2, &fieldIDs[0], &fieldSizes[0]);

  int idTypes[2] = {1, 2};
  vspace.defineIDTypes(2, &idTypes[0]);

  std::vector<int> ids(10);
  for(size_t i=0; i<ids.size(); ++i) ids[i] = i;

  vspace.addDOFs(fieldIDs[0], idTypes[0], ids.size(), &ids[0]);
  vspace.addDOFs(fieldIDs[1], idTypes[0], ids.size(), &ids[0]);
  vspace.addDOFs(fieldIDs[0], idTypes[1], ids.size(), &ids[0]);

  vspace.initComplete();

  fei::ReverseMapper rm(vspace);

  //fei::VectorSpace produces a set of equation-numbers grouped by
  //field, then by id, then by id-type. In other words, all components of
  //a field are contiguous in the equation-space, then all equations for
  //fields at an id (such as a node) are contiguous, etc.
  //The VectorSpace initialized above was given both fields for ids with
  //type idTypes[0], and one field for ids with type idTypes[1].
  //There should be 50 equations total, numbered 0 .. 49. ids[0] with
  //idTypes[0] should have 3 equations. The 5th equation, eqn==4, should
  //be the first scalar component of the 2nd field on ids[1].

  fei::EqnRecord er1 = rm.getEqnRecord(4);

  TEUCHOS_TEST_EQUALITY(er1.IDType, idTypes[0], out, success);
  TEUCHOS_TEST_EQUALITY(er1.ID, ids[1], out, success);

  TEUCHOS_TEST_EQUALITY(er1.fieldID, fieldIDs[1], out, success);

  TEUCHOS_TEST_EQUALITY(er1.offset, 0, out, success);
}

}//namespace <anonymous>


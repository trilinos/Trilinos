// @HEADER
// ***********************************************************************
//                Copyright message goes here.   
// ***********************************************************************
// @HEADER

/*! \file Util.cpp
 *  \brief Tests methods in Zoltan2_Util.hpp
 *
 *   \todo Some of the tests require that you look at the output
 *          to know if they did the right thing.  Enhance this so
 *          the test itself determines correctness.
 */

#include <Zoltan2_Util.hpp>   
#include <Zoltan2_TestHelpers.hpp>   
#include <Zoltan2_PartitioningSolution.hpp>   
#include <Zoltan2_BasicIdentifierInput.hpp>   

using Teuchos::RCP;
using Teuchos::Comm;
using Zoltan2::Environment;
using Zoltan2::PartitioningSolution;
using namespace std;

void convertSolutionToImportListTest(RCP<const Comm<int> > &comm)
{
  int numProcs = comm->getSize();
  int rank = comm->getRank();

  typedef Zoltan2::BasicUserTypes<scalar_t, gno_t, lno_t, gno_t> userTypes_t;
  typedef Zoltan2::IdentifierMap<userTypes_t> idMap_t;

  // Set num_local_parts to 3.

  Teuchos::ParameterList myParams("testParameterList");
  Teuchos::ParameterList &parParams = myParams.sublist("partitioning");

  parParams.set("num_local_parts", 3);

  zoltan2_partId_t myMinPart = rank * 3;
  zoltan2_partId_t myMaxPart = (rank+1) * 3 - 1;

  // Create an environment

  RCP<const Environment> env = rcp(new Environment(myParams, comm));

  // Assume a contiguous distribution of gids to procs.
  
  lno_t localNumObjects = 100;
  gno_t myBaseId = localNumObjects * rank;

  gno_t *mygids  = new gno_t [localNumObjects];
  for (lno_t i=0, base=myBaseId; i < localNumObjects; i++, base++)
    mygids[i] = base;

  ArrayRCP<const gno_t> gidArray(mygids, 0, localNumObjects, true);

  RCP<const idMap_t> idMap = rcp(new idMap_t(env, comm, gidArray, false));

  // An identifier adapter type for Solution

  typedef Zoltan2::BasicIdentifierInput<userTypes_t> idAdapter_t;

  // Create a partitioning solution.  Update with a solution
  // where gids are assigned to parts in a round robin fashion.

  RCP<PartitioningSolution<idAdapter_t> > solution =
    rcp(new PartitioningSolution<idAdapter_t>(env, comm, idMap, 1));

  size_t numGlobalParts = solution->getGlobalNumberOfParts();

  TEST_FAIL_AND_EXIT(*comm, (numGlobalParts == 3 * numProcs),
    "numGlobalParts", 1);

  zoltan2_partId_t *partList = new zoltan2_partId_t [localNumObjects];

  for (lno_t i=0; i < localNumObjects; i++){
    partList[i] = mygids[i] % numGlobalParts;
  }
    
  ArrayRCP<zoltan2_partId_t> partArray(partList, 0, localNumObjects, true);

  ArrayRCP<Zoltan2::MetricValues<scalar_t> > metrics;

  solution->setParts(gidArray,
                     partArray,
                     metrics);

  // Check if gids assigned to me are correct.

  ArrayRCP<int> dummyInfoOut;
  ArrayRCP<int> dummyInfoIn;
  ArrayRCP<gno_t> myImports;

  size_t numNewObjects = 
    Zoltan2::convertSolutionToImportList<idAdapter_t, int>(
      *solution, dummyInfoOut, myImports, dummyInfoIn);

  int fail = 0;

  for (lno_t i=0; !fail && i < numNewObjects; i++){
    gno_t mapsToPart = myImports[i] % numGlobalParts;
    if (mapsToPart < myMinPart || mapsToPart > myMaxPart)
      fail = 1;
  }
  
  TEST_FAIL_AND_EXIT(*comm, fail == 0, "myImports", 1);
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  convertSolutionToImportListTest(comm);

  if (comm->getRank() == 0)
    std::cout << "PASS" << std::endl;
}

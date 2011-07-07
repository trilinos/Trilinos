
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <fei_iostream.hpp>
#include <fei_Pattern.hpp>
#include <fei_MatrixGraph_Impl2.hpp>

#include <vector>
#include <cmath>

TEUCHOS_UNIT_TEST(MatGraph, MatGraph_test1)
{
  int numprocs = fei::numProcs(MPI_COMM_WORLD);
  if (numprocs > 1) return;

  //two id-types: 0, 1:
  int idT[] = {0, 1};
  snl_fei::RecordCollection* recColls[] = {NULL,NULL};

  //set up a pattern for an element that has 6 ids: 3 nodes and 3 edges.
  const int numIDs = 6;

  //assume 0 is the node-type and 1 is the edge-type:
  int idTypes[] = {0, 0, 0, 1, 1, 1};

  //for the first pattern, only the edge ids will have a field attached:
  int fieldsPerID[] = {0, 0, 0, 1, 1, 1};

  int fieldID = 0;
  int fieldSize = 1;

  int fieldIDs[] = {fieldID, fieldID, fieldID};
  int fieldSizes[] = {fieldSize, fieldSize, fieldSize};

  fei::Pattern pattern1(numIDs, &idTypes[0], &recColls[0], &fieldsPerID[0], &fieldIDs[0], &fieldSizes[0]);

  //declare a vector-space, do some rudimentary initializations:

  fei::SharedPtr<fei::VectorSpace> rowspace(new fei::VectorSpace(MPI_COMM_WORLD));
  rowspace->defineIDTypes(2, &idT[0]);
  rowspace->defineFields(1, &fieldID, &fieldSize);

  fei::SharedPtr<fei::VectorSpace> colspace;

  //declare a matrix-graph:
  fei::MatrixGraph_Impl2 mgraph(rowspace, colspace);

  int patternID1 = mgraph.definePattern(numIDs, &idTypes[0], &fieldsPerID[0], &fieldIDs[0]);

  //unit-test: make sure the matrix-graph's pattern is the same as our
  //explicitly-declared pattern:
  fei::Pattern* pttn1 = mgraph.getPattern(patternID1);

  TEUCHOS_TEST_EQUALITY(pattern1 == *pttn1, true, out, success);

  //now declare a second pattern which is the same except now fields are
  //attached to nodes instead of edges:
  fieldsPerID[0] = 1; fieldsPerID[1] = 1; fieldsPerID[2] = 1;
  fieldsPerID[3] = 0; fieldsPerID[4] = 0; fieldsPerID[5] = 0;

  fei::Pattern pattern2(numIDs, &idTypes[0], &recColls[0], &fieldsPerID[0], &fieldIDs[0], &fieldSizes[0]);

  int patternID2 = mgraph.definePattern(numIDs, &idTypes[0], &fieldsPerID[0], &fieldIDs[0]);
  fei::Pattern* pttn2 = mgraph.getPattern(patternID2);
  TEUCHOS_TEST_EQUALITY(pattern2 == *pttn2, true, out, success);

  //declare two element-blocks, one for each pattern. each element block will have
  //just one element:
  mgraph.initConnectivityBlock(0, 1, patternID1);
  mgraph.initConnectivityBlock(1, 1, patternID2);

//Two-element mesh, each element has 3 vertex-nodes, and 3 edges:
/*
   0      1
   o---4--o
    \     | \
     \    |  \
      5   6   7
       \  |    \
        \ |     \
          o---8--o
          2      3
*/

  //the first element has nodes 0, 1, 2 and edges 4, 5, 6:
  int ids0[] = {0, 1, 2, 4, 5, 6};

  mgraph.initConnectivity(0, 0, &ids0[0]);

  //the second element has nodes 1, 2, 3 and edges 6, 7, 8:
  int ids1[] = {1, 2, 3, 6, 7, 8};

  mgraph.initConnectivity(1, 0, &ids1[0]);

  mgraph.initComplete();

  fei::SharedPtr<fei::SparseRowGraph> srg = mgraph.createGraph(false);

//The way we set things up, the graph should have 6 rows,
//with 3 nonzeros per row:

  TEUCHOS_TEST_EQUALITY(srg->rowNumbers.size(), 6, out, success);
  TEUCHOS_TEST_EQUALITY(srg->packedColumnIndices.size(), 18, out, success);
}


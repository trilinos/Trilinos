
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <fei_iostream.hpp>
#include <fei_Pattern.hpp>

#include <vector>
#include <cmath>

TEUCHOS_UNIT_TEST(Pattern, Pattern_test1)
{
  int numIDs = 6;
  std::vector<int> idTypes(numIDs);
  std::vector<snl_fei::RecordCollection*> recColls(numIDs,(snl_fei::RecordCollection*)NULL);
  std::vector<int> fieldsPerID(numIDs);
  std::vector<int> fieldIDs(3);
  std::vector<int> fieldSizes(3, 1);

  idTypes[0] = 0;
  idTypes[1] = 0;
  idTypes[2] = 0;
  idTypes[3] = 1;
  idTypes[4] = 1;
  idTypes[5] = 1;

  fieldsPerID[0] = 0;
  fieldsPerID[1] = 0;
  fieldsPerID[2] = 0;
  fieldsPerID[3] = 1;
  fieldsPerID[4] = 1;
  fieldsPerID[5] = 1;

  fieldIDs[0] = 0;
  fieldIDs[1] = 0;
  fieldIDs[2] = 0;

  fei::Pattern pattern1(numIDs, 0, recColls[0], &fieldsPerID[0], &fieldIDs[0], &fieldSizes[0]);

  TEUCHOS_TEST_EQUALITY(pattern1.getTotalNumFields(), 3, out, success);
  TEUCHOS_TEST_EQUALITY(pattern1.getNumIndices(), 3, out, success);

  fei::Pattern pattern2(numIDs, &idTypes[0], &recColls[0], &fieldsPerID[0], &fieldIDs[0], &fieldSizes[0]);

  TEUCHOS_TEST_EQUALITY(pattern2.getTotalNumFields(), 3, out, success);
  TEUCHOS_TEST_EQUALITY(pattern2.getNumIndices(), 3, out, success);
}


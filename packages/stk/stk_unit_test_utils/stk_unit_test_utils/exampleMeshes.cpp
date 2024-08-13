#include <stk_unit_test_utils/exampleMeshes.h>

namespace unitTestUtils
{
namespace exampleMeshes
{

namespace simple_fields {

void fillDataForUnitCube(std::vector<double> &coordinates) {
  unitTestUtils::exampleMeshes::fillDataForUnitCube(coordinates);
}

void fillDataForRectangloid(std::vector<double> &coordinates) {
  unitTestUtils::exampleMeshes::fillDataForRectangloid(coordinates);
}

Iogn::ExodusData createExodusDataForDisconnectedHex8s(int numberOfHexes) {
  return unitTestUtils::exampleMeshes::createExodusDataForDisconnectedHex8s(numberOfHexes);
}

} // namespace simple_fields

}
}

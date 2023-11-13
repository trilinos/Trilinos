#ifndef KRINO_KRINO_KRINO_LIB_AKRI_VOLUMEPRESERVINGSNAPPINGLIMITER_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_VOLUMEPRESERVINGSNAPPINGLIMITER_HPP_
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <Akri_FieldRef.hpp>
#include <stk_math/StkVector.hpp>

namespace krino {

class AuxMetaData;

class VolumePreservingSnappingLimiter
{
public:
  typedef std::function<stk::mesh::Part*(const stk::mesh::BulkData &, const stk::mesh::Entity)> ElementToBlockConverter;

  VolumePreservingSnappingLimiter(
    const stk::mesh::BulkData & mesh,
    const FieldRef coordsField,
    const ElementToBlockConverter & elementToBlockConverter,
    const double volumeConservationTol);
  bool is_snap_allowed(const stk::mesh::Entity node, const stk::math::Vector3d & snapLocation) const;
private:
  std::set<stk::mesh::Part*> get_blocks_to_consider(const stk::mesh::Entity node) const;
  const stk::mesh::BulkData & myMesh;
  const AuxMetaData & myAuxMeta;
  ElementToBlockConverter myElementToBlockConverter;
  FieldRef myCoordsField;
  double myVolumeConservationTol;
};

}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_VOLUMEPRESERVINGSNAPPINGLIMITER_HPP_ */

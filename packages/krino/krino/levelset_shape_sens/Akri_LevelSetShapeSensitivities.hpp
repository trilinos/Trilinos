#ifndef KRINO_KRINO_KRINO_LIB_AKRI_LEVELSETSHAPESENSITIVITIES_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_LEVELSETSHAPESENSITIVITIES_HPP_

#include <array>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <Akri_LS_Field.hpp>
#include <Akri_LevelSetShapeSensitivityTypes.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Types.hpp>

namespace krino {

std::string output_sensitivity(const LevelSetShapeSensitivity & sens);

// Decomposes a mesh to conform to level sets and computes the shape sensitivities of the
// interface.  Construct it before the MetaData is committed: it registers the fields for the
// decomposition and for the sensitivities of either method.  With snapping, the decomposition
// snaps a copy of each level set, so the sensitivities can read the original, and uses predictable
// snap selection, so that the snapped sensitivities are exact (see Snap_Selection).  Without
// snapping, the CDFEM edge tolerance is set to its minimum.  The mesh must outlive this object.
class ConformingDecompositionWithSensitivities
{
public:
  // If interfaceRestriction is given, each level set's interface is restricted to it (for example
  // to the sides and nodes of the material blocks, so that it is the material boundary); otherwise
  // it is everywhere the level set changes sign.
  ConformingDecompositionWithSensitivities(stk::mesh::MetaData & meta,
      const std::vector<LS_Field> & lsFields,
      const bool doSnapping,
      const stk::mesh::Selector * interfaceRestriction = nullptr);

  // Decomposes the mesh for the current level set values; can be repeated as they change.
  void decompose();

  // Removes the decomposition: deletes the child elements and undoes the snapping of nodes.  Call
  // it before setting the level sets of a new design, so that they are set and cut on the original
  // mesh rather than the previous decomposition.
  void reset_to_undecomposed_mesh();

  const std::vector<LS_Field> & levelset_fields() const { return myLSFields; }
  bool does_snapping() const { return myDoSnapping; }

  std::vector<LevelSetShapeSensitivity>
  sensitivities(const unsigned lsIndex, const ShapeSensitivityMethod method) const;
  std::vector<std::vector<LevelSetShapeSensitivity>>
  all_sensitivities(const ShapeSensitivityMethod method) const;
  template <size_t NNODES>
  InterfaceFacetsWithSensitivities<NNODES>
  facets_and_sensitivities(const unsigned lsIndex, const ShapeSensitivityMethod method) const;

private:
  const LS_Field & levelset_field(const unsigned lsIndex) const;
  stk::mesh::Selector interface_selector(const LS_Field & lsField) const;

  stk::mesh::MetaData & myMeta;
  std::vector<LS_Field> myLSFields;
  bool myDoSnapping;
  std::optional<stk::mesh::Selector> myInterfaceRestriction;
};

// The position (0 to 1) along the unsnapped edge from unsnappedEdgeLoc0 to unsnappedEdgeLoc1 of a
// node snapped by snapDisp from unsnappedEdgeLoc0, and the square of its distance from that
// position.  Exposed for testing.
std::pair<double, double> snapped_node_edge_position_and_square_error(const unsigned dim,
    const stk::math::Vector3d & snapDisp,
    const stk::math::Vector3d & unsnappedEdgeLoc0,
    const stk::math::Vector3d & unsnappedEdgeLoc1);
}



#endif /* KRINO_KRINO_KRINO_LIB_AKRI_LEVELSETSHAPESENSITIVITIES_HPP_ */

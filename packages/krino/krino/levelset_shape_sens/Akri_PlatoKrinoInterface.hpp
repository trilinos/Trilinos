#ifndef KRINO_KRINO_LEVELSET_SHAPE_SENS_AKRI_PLATOKRINOINTERFACE_HPP_
#define KRINO_KRINO_LEVELSET_SHAPE_SENS_AKRI_PLATOKRINOINTERFACE_HPP_

#include <array>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <mpi.h>

#include <Akri_LevelSetShapeSensitivityTypes.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Types.hpp> // stk::mesh::EntityId

//------------------------------------------------------------------------------
// PlatoKrinoInterface: a facade over krino's level set decomposition and shape
// sensitivities for optimization drivers such as plato.
//
// Design:
//   * This is the only krino header a consumer needs; it includes only krino's
//     sensitivity value types (Akri_LevelSetShapeSensitivityTypes.hpp) and STK.
//   * The background mesh can be read from a file or generated as a bounding box
//     by BackgroundMesh, or built by the consumer.
//   * Each design gets a new ConformingMeshWithSensitivities, which decomposes a
//     copy of the background mesh; the background mesh is never modified.
//
// Conventions:
//   * The level set is positive in the material and negative in the void.  With
//     several level sets, the material is where every level set is positive, and
//     each level set's interface, sensitivities and facets are on the material
//     boundary.
//     Geometry describes the material: the union of the insides of its spheres,
//     of the half-spaces behind its planes, whose normals point out of the
//     material, and of the material sides of its facet file surfaces.
//   * The sensitivities' parent node ids are ids of the background mesh's nodes,
//     which the copy keeps.  The copy has all the background mesh's fields; on
//     the nodes the decomposition adds, the level sets are interpolated and the
//     other fields have their initial values.
//   * Parallel: a shared interface node appears on every processor that has it;
//     each interface facet appears on exactly one processor.  Everything is
//     collective except the accessors (bulk_data, num_level_sets,
//     level_set_field and the selectors); invalid arguments on any processor
//     throw on every processor.  Level set indices must be the same on every
//     processor.
//------------------------------------------------------------------------------

namespace PlatoKrinoInterface {

// Implicit-geometry primitives of the material region (see the conventions above).
struct SphereSpec
{
  std::array<double,3> center{{0.,0.,0.}};
  double radius = 0.;
};

// The half-space normal.x + offset < 0; normal need not be a unit vector.
struct PlaneSpec
{
  std::array<double,3> normal{{1.,0.,0.}};
  double offset = 0.;
};

// A closed surface read from a facet file: STL, FAC, PLY or EXO, from the file's
// extension.  The file's coordinates are multiplied by scale.  The material is
// inside the surface, or outside it if !materialInside.  3D only.
struct FacetFileSpec
{
  std::string fileName;
  std::array<double, 3> scale{{1., 1., 1.}};
  bool materialInside = true;
};

struct Geometry
{
  std::vector<SphereSpec> spheres{};
  std::vector<PlaneSpec> planes{};
  std::vector<FacetFileSpec> facetFiles{};
};

// A caller-supplied level-set value as a function of a node's coordinates.
using LevelSetValueFunction = std::function<double(const std::array<double, 3> & coords)>;

// The level set fields of a BackgroundMesh.
struct BackgroundMeshOptions
{
  // Nodal level set fields "LS", or "LS0", "LS1", ... for several, on every node.
  unsigned numLevelSets = 1;
};

// A background mesh read from a file or generated as a bounding box, with the
// level set fields of its options.  Owns the mesh.  Move-only.
class BackgroundMesh
{
public:
  // decompMethod is the parallel decomposition of the file, or "" if it is
  // already split per processor.
  static BackgroundMesh read(const std::string & fileName,
      MPI_Comm comm,
      const std::string & decompMethod = "rib",
      const BackgroundMeshOptions & options = {});
  // A triangle (2D) or tetrahedron (3D) mesh of the box with the given corners, with flat walls
  // whose sides are in the sidesets Xminus, Xplus, Yminus, Yplus (and Zminus, Zplus in 3D).
  static BackgroundMesh bounding_box(unsigned dim,
      const std::array<double, 3> & minCorner,
      const std::array<double, 3> & maxCorner,
      double meshSize,
      MPI_Comm comm,
      const BackgroundMeshOptions & options = {});

  ~BackgroundMesh();
  BackgroundMesh(BackgroundMesh &&) noexcept;
  BackgroundMesh & operator=(BackgroundMesh &&) noexcept;
  BackgroundMesh(const BackgroundMesh &) = delete;
  BackgroundMesh & operator=(const BackgroundMesh &) = delete;

  stk::mesh::BulkData & bulk_data();
  const stk::mesh::BulkData & bulk_data() const;
  unsigned num_level_sets() const;
  stk::mesh::FieldBase & level_set_field(unsigned lsIndex = 0) const;

private:
  class Impl;
  explicit BackgroundMesh(std::unique_ptr<Impl> impl);

  std::unique_ptr<Impl> myImpl;
};

// Set a level set field, a nodal double field of mesh on every node: to the
// signed distance to the geometry, to valueFn at each node's coordinates, or to
// values at the given nodes of this processor.  Only owned nodes' values are
// kept, and they are communicated to the other processors.
void initialize_field_from_geometry(const stk::mesh::BulkData & mesh,
    stk::mesh::FieldBase & levelSetField,
    const Geometry & geometry);
void initialize_field_from_function(const stk::mesh::BulkData & mesh,
    stk::mesh::FieldBase & levelSetField,
    const LevelSetValueFunction & valueFn);
void set_field_values(const stk::mesh::BulkData & mesh,
    stk::mesh::FieldBase & levelSetField,
    const std::vector<stk::mesh::EntityId> & nodeIds,
    const std::vector<double> & values);

// How a ConformingMeshWithSensitivities decomposes its mesh and computes its
// sensitivities.
struct DecompositionOptions
{
  // Snap nodes near the interface onto it.
  bool doSnapping = false;
  // Described in Akri_LevelSetShapeSensitivityTypes.hpp.
  krino::ShapeSensitivityMethod sensitivityMethod =
      krino::ShapeSensitivityMethod::EDGE_DIRECTION_USING_EDGE_GRADIENT;
};

// The interface-conforming decomposition of a copy of backgroundMesh for the
// level sets in levelSetFields (nodal double fields of backgroundMesh on all its
// nodes), and its shape sensitivities.  Move-only.
class ConformingMeshWithSensitivities
{
public:
  ConformingMeshWithSensitivities(const stk::mesh::BulkData & backgroundMesh,
      const std::vector<const stk::mesh::FieldBase *> & levelSetFields,
      const DecompositionOptions & options = {});

  ~ConformingMeshWithSensitivities();
  ConformingMeshWithSensitivities(ConformingMeshWithSensitivities &&) noexcept;
  ConformingMeshWithSensitivities & operator=(ConformingMeshWithSensitivities &&) noexcept;
  ConformingMeshWithSensitivities(const ConformingMeshWithSensitivities &) = delete;
  ConformingMeshWithSensitivities & operator=(const ConformingMeshWithSensitivities &) = delete;

  // The conforming mesh.  Use the selectors below for its material, void and
  // interfaces; changing it invalidates the sensitivities.
  stk::mesh::BulkData & bulk_data();
  const stk::mesh::BulkData & bulk_data() const;
  unsigned num_level_sets() const;

  // The active elements of the material, with their sides and nodes.
  stk::mesh::Selector material_selector() const;
  // The active elements of the void, with their sides and nodes.
  stk::mesh::Selector void_selector() const;
  // The material boundary sides of level set lsIndex, where it changes sign
  // between the material and the void; the same sides as interface_facets.
  stk::mesh::Selector interface_selector(unsigned lsIndex = 0) const;

  // Shape sensitivities of level set lsIndex, one per owned or shared interface
  // node.  all_shape_sensitivities()[i] is shape_sensitivities(i).
  std::vector<krino::LevelSetShapeSensitivity> shape_sensitivities(unsigned lsIndex = 0) const;
  std::vector<std::vector<krino::LevelSetShapeSensitivity>> all_shape_sensitivities() const;

  // The interface facets of level set lsIndex that this processor owns, with
  // their nodes' sensitivities; NNODES is the spatial dimension (2 or 3).
  template <size_t NNODES>
  krino::InterfaceFacetsWithSensitivities<NNODES> interface_facets(unsigned lsIndex = 0) const;

  // Outputs the conforming mesh: its material, and its void too if includeVoidPhase.
  void write_decomposed_mesh(
      const std::string & fileName, bool includeVoidPhase, int step = 1, double time = 0.0);

private:
  class Impl;
  std::unique_ptr<Impl> myImpl;
};

} // namespace PlatoKrinoInterface

#endif /* KRINO_KRINO_LEVELSET_SHAPE_SENS_AKRI_PLATOKRINOINTERFACE_HPP_ */

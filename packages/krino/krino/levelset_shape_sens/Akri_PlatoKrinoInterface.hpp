#ifndef KRINO_KRINO_LEVELSET_SHAPE_SENS_AKRI_PLATOKRINOINTERFACE_HPP_
#define KRINO_KRINO_LEVELSET_SHAPE_SENS_AKRI_PLATOKRINOINTERFACE_HPP_

#include <array>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <mpi.h>

#include <stk_mesh/base/Types.hpp> // stk::mesh::EntityId

//------------------------------------------------------------------------------
// PlatoKrinoInterface: a narrow, stable facade over the krino level-set / CDFEM
// decomposition + shape-sensitivity workflow, intended as the single supported
// entry point for external optimization drivers (e.g. plato).
//
// Design goals:
//   * This is the ONLY krino header a consumer needs to include. No krino type
//     (LS_Field, MeshInterface, FieldRef, Surface, LevelSetShapeSensitivity, ...)
//     appears below; the krino mesh state is hidden behind a pimpl and geometry
//     / sensitivities are expressed with std types plus stk::mesh::EntityId.
//     Consumers still use their own MPI, STK, and standard library, but no
//     <Akri_*.hpp>. (STK is not krino, and every intended consumer is already an
//     STK application, so EntityId is the natural node-id currency; coordinate
//     deltas stay std::array to avoid depending on <stk_math/StkVector.hpp>.)
//   * It packages workflows that consumers currently re-implement by hand: the
//     multi-step "read a mesh and register the fields needed to decompose it"
//     recipe, and the entire cut-surface shape-sensitivity computation. Both
//     already exist inside krino; this facade just exposes intent-level verbs.
//   * doSnapping / doComputeClosestPointSensitivities must be threaded
//     identically through field setup, decomposition, and sensitivity
//     extraction. They are captured once in LevelSetMesh so callers cannot make
//     them inconsistent.
//   * Multiple level sets are supported throughout; single-level-set callers use
//     the lsIndex == 0 defaults.
//------------------------------------------------------------------------------

namespace PlatoKrinoInterface {

// Shape sensitivity of one cut-surface (interface) node's coordinates with
// respect to the level-set values at its two background parent nodes. This is
// the facade's krino-free stand-in for krino::LevelSetShapeSensitivity.
//   dCoordsdParentLevelSets[i] = d(interface node coords) / d(level set at parentNodeIds[i])
struct ShapeSensitivity
{
  stk::mesh::EntityId interfaceNodeId = 0;
  std::vector<stk::mesh::EntityId> parentNodeIds;
  std::vector<std::array<double,3>> dCoordsdParentLevelSets;
};

// Implicit-geometry primitives used to initialize a level-set field. The set of
// spheres and planes handed to initialize_levelset_field_from_geometry is
// unioned (CSG) to form the initial zero level set.
struct SphereSpec
{
  std::array<double,3> center{{0.,0.,0.}};
  double radius = 0.;
};

struct PlaneSpec
{
  std::array<double,3> normal{{1.,0.,0.}};
  double offset = 0.;
};

struct Geometry
{
  std::vector<SphereSpec> spheres;
  std::vector<PlaneSpec> planes;
};

// Opaque owning handle for a krino mesh prepared for level-set decomposition.
// All krino state lives in the pimpl, so this header stays krino-free. Move-only.
class LevelSetMesh
{
public:
  LevelSetMesh();
  ~LevelSetMesh();
  LevelSetMesh(LevelSetMesh &&) noexcept;
  LevelSetMesh & operator=(LevelSetMesh &&) noexcept;
  LevelSetMesh(const LevelSetMesh &) = delete;
  LevelSetMesh & operator=(const LevelSetMesh &) = delete;

  unsigned spatial_dimension() const;
  unsigned num_level_sets() const;

  class Impl;                 // defined in the .cpp; holds the krino mesh + fields
  Impl & impl();
  const Impl & impl() const;

private:
  std::unique_ptr<Impl> myImpl;
};

// A caller-supplied level-set value as a function of a node's coordinates. Used
// by initialize_levelset_field_from_function so an optimizer can stamp its own
// design field without touching any krino field API.
using LevelSetValueFunction = std::function<double(const std::array<double,3> & coords)>;

//--- Workflow A: read + setup ------------------------------------------------

// Read a mesh from file and prepare it for level-set decomposition: register
// numLevelSets level sets (per-interface policy, void phase for any negative
// level set), register every CDFEM field required for the requested snapping /
// closest-point-sensitivity modes, then populate and activate the mesh.
// decompMethod is the parallel decomposition ("rib" by default).
LevelSetMesh read_and_setup_mesh_for_levelset_decomposition(
    const std::string & fileName,
    MPI_Comm comm,
    unsigned numLevelSets,
    bool doSnapping,
    bool doComputeClosestPointSensitivities,
    const std::string & decompMethod = "rib");

//--- Workflow B: initialize level-set fields ---------------------------------

// Initialize level set lsIndex from a union of implicit primitives.
void initialize_levelset_field_from_geometry(
    LevelSetMesh & mesh,
    const Geometry & geometry,
    unsigned lsIndex = 0);

// Initialize level set lsIndex by evaluating valueFn at each node's coordinates.
void initialize_levelset_field_from_function(
    LevelSetMesh & mesh,
    const LevelSetValueFunction & valueFn,
    unsigned lsIndex = 0);

//--- Workflow B: cut ---------------------------------------------------------

// Decompose the mesh so element boundaries conform to the current level-set
// field values, using the options captured at setup time.
void decompose_mesh(LevelSetMesh & mesh);

//--- Workflow C: shape sensitivities -----------------------------------------

// Shape sensitivities of the cut-surface node coordinates for level set lsIndex.
// First-class replacement for hand-rolled child-node-stencil differentiation.
std::vector<ShapeSensitivity> get_shape_sensitivities(
    const LevelSetMesh & mesh,
    unsigned lsIndex = 0);

// Shape sensitivities for every level set, indexed [lsIndex][interfaceNode].
std::vector<std::vector<ShapeSensitivity>> get_all_shape_sensitivities(
    const LevelSetMesh & mesh);

//--- Workflow B: write -------------------------------------------------------

// Output the decomposed mesh. When includeVoidPhase is false, "*_void" phase
// parts are excluded (the convention consumers currently reproduce with ad hoc
// part-name string matching).
void write_decomposed_mesh(
    LevelSetMesh & mesh,
    const std::string & fileName,
    bool includeVoidPhase,
    int step = 1,
    double time = 0.0);

} // namespace PlatoKrinoInterface

#endif /* KRINO_KRINO_LEVELSET_SHAPE_SENS_AKRI_PLATOKRINOINTERFACE_HPP_ */

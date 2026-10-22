#include <Akri_PlatoKrinoInterface.hpp>

#include <memory>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <Akri_AuxMetaData.hpp>
#include <Akri_AnalyticSurf.hpp>
#include <Akri_CSG_Surface.hpp>
#include <Akri_CSG_Tree.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_LevelSetShapeSensitivities.hpp>
#include <Akri_MeshFromFile.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodalSurfaceDistance.hpp>
#include <Akri_OutputUtils.hpp>
#include <Akri_Phase_Support.hpp>

namespace PlatoKrinoInterface {

//--- pimpl: the actual krino mesh state --------------------------------------

class LevelSetMesh::Impl
{
public:
  std::unique_ptr<krino::MeshInterface> mesh;
  std::vector<krino::LS_Field> lsFields;
  bool doSnapping = false;
  bool doComputeClosestPointSensitivities = false;

  stk::mesh::BulkData & bulk() { return mesh->bulk_data(); }
  const stk::mesh::BulkData & bulk() const { return mesh->bulk_data(); }
  stk::mesh::MetaData & meta() { return mesh->meta_data(); }
  const stk::mesh::MetaData & meta() const { return mesh->meta_data(); }

  void check_index(const unsigned lsIndex) const
  {
    STK_ThrowRequireMsg(lsIndex < lsFields.size(),
        "PlatoKrinoInterface: level set index " << lsIndex << " out of range (" << lsFields.size() << " level sets).");
  }
};

LevelSetMesh::LevelSetMesh() : myImpl(std::make_unique<Impl>()) {}
LevelSetMesh::~LevelSetMesh() = default;
LevelSetMesh::LevelSetMesh(LevelSetMesh &&) noexcept = default;
LevelSetMesh & LevelSetMesh::operator=(LevelSetMesh &&) noexcept = default;

LevelSetMesh::Impl & LevelSetMesh::impl() { return *myImpl; }
const LevelSetMesh::Impl & LevelSetMesh::impl() const { return *myImpl; }

unsigned LevelSetMesh::spatial_dimension() const { return myImpl->meta().spatial_dimension(); }
unsigned LevelSetMesh::num_level_sets() const { return myImpl->lsFields.size(); }

//--- local conversions -------------------------------------------------------

namespace {

stk::math::Vector3d to_vec(const std::array<double,3> & a) { return stk::math::Vector3d(a[0], a[1], a[2]); }
std::array<double,3> to_array(const stk::math::Vector3d & v) { return {{v[0], v[1], v[2]}}; }

ShapeSensitivity convert(const krino::LevelSetShapeSensitivity & sens)
{
  ShapeSensitivity out;
  out.interfaceNodeId = sens.interfaceNodeId;
  out.parentNodeIds = sens.parentNodeIds;
  out.dCoordsdParentLevelSets.reserve(sens.dCoordsdParentLevelSets.size());
  for (auto & d : sens.dCoordsdParentLevelSets)
    out.dCoordsdParentLevelSets.push_back(to_array(d));
  return out;
}

std::vector<ShapeSensitivity> convert(const std::vector<krino::LevelSetShapeSensitivity> & sensitivities)
{
  std::vector<ShapeSensitivity> out;
  out.reserve(sensitivities.size());
  for (auto & sens : sensitivities)
    out.push_back(convert(sens));
  return out;
}

stk::mesh::Selector build_output_selector(const stk::mesh::MetaData & meta, const stk::mesh::Part & activePart, const bool includeVoidPhase)
{
  if (includeVoidPhase)
    return activePart;

  stk::mesh::PartVector outputParts;
  for (auto * part : meta.get_parts())
    if (stk::io::is_part_io_part(*part) && part->name().find("_void") == std::string::npos)
      outputParts.push_back(part);
  return activePart & stk::mesh::selectUnion(outputParts);
}

} // anonymous namespace

//--- Workflow A: read + setup ------------------------------------------------

LevelSetMesh read_and_setup_mesh_for_levelset_decomposition(
    const std::string & fileName,
    MPI_Comm comm,
    const unsigned numLevelSets,
    const bool doSnapping,
    const bool doComputeClosestPointSensitivities,
    const std::string & decompMethod)
{
  LevelSetMesh result;
  LevelSetMesh::Impl & impl = result.impl();
  impl.doSnapping = doSnapping;
  impl.doComputeClosestPointSensitivities = doComputeClosestPointSensitivities;

  impl.mesh = std::make_unique<krino::MeshFromFile>(fileName, comm, decompMethod);

  stk::mesh::MetaData & meta = impl.meta();
  impl.lsFields = krino::LSPerInterfacePolicy::setup_levelsets_on_all_blocks_with_void_phase_for_any_negative_levelset(meta, numLevelSets);

  krino::setup_fields_for_conforming_decomposition(meta, impl.lsFields, doSnapping, doComputeClosestPointSensitivities);

  impl.mesh->populate_mesh();
  krino::activate_all_entities(impl.bulk(), krino::AuxMetaData::get(meta).active_part());

  return result;
}

//--- Workflow B: initialize --------------------------------------------------

void initialize_levelset_field_from_geometry(LevelSetMesh & mesh, const Geometry & geometry, const unsigned lsIndex)
{
  LevelSetMesh::Impl & impl = mesh.impl();
  impl.check_index(lsIndex);

  std::vector<std::shared_ptr<krino::Surface>> surfaces;
  for (const auto & sphere : geometry.spheres)
    surfaces.emplace_back(std::make_shared<krino::Sphere>(to_vec(sphere.center), sphere.radius));
  for (const auto & plane : geometry.planes)
    surfaces.emplace_back(std::make_shared<krino::Plane>(to_vec(plane.normal), plane.offset));

  STK_ThrowRequireMsg(!surfaces.empty(), "initialize_levelset_field_from_geometry: geometry has no spheres or planes.");

  krino::CSG_Surface csgSurface(krino::CSGTree::unite(surfaces));
  krino::compute_nodal_surface_distance(impl.bulk(), impl.meta().coordinate_field(), impl.lsFields[lsIndex].isovar, csgSurface);
}

void initialize_levelset_field_from_function(LevelSetMesh & mesh, const LevelSetValueFunction & valueFn, const unsigned lsIndex)
{
  LevelSetMesh::Impl & impl = mesh.impl();
  impl.check_index(lsIndex);

  const unsigned nDim = impl.meta().spatial_dimension();
  const krino::FieldRef coordsField = impl.meta().coordinate_field();
  const krino::FieldRef distanceField = impl.lsFields[lsIndex].isovar;

  const stk::mesh::BucketVector & buckets = impl.bulk().get_buckets(stk::topology::NODE_RANK, stk::mesh::selectField(distanceField));
  auto distData = distanceField.data<double, stk::mesh::ReadWrite>();
  auto coordData = coordsField.data<double>();

  for (auto && bucketPtr : buckets)
  {
    const stk::mesh::Bucket & b = *bucketPtr;
    auto dist = distData.bucket_values(b);
    for (stk::mesh::EntityIdx n : b.entities())
    {
      const stk::math::Vector3d x(coordData.entity_values(b[n]), nDim);
      dist(n) = valueFn({{x[0], x[1], x[2]}});
    }
  }

  stk::mesh::communicate_field_data(impl.bulk(), {&distanceField.field()});
}

//--- Workflow B: cut ---------------------------------------------------------

void decompose_mesh(LevelSetMesh & mesh)
{
  LevelSetMesh::Impl & impl = mesh.impl();
  krino::decompose_mesh_to_conform_to_levelsets(impl.bulk(), impl.lsFields, impl.doSnapping, impl.doComputeClosestPointSensitivities);
}

//--- Workflow C: shape sensitivities -----------------------------------------

std::vector<ShapeSensitivity> get_shape_sensitivities(const LevelSetMesh & mesh, const unsigned lsIndex)
{
  const LevelSetMesh::Impl & impl = mesh.impl();
  impl.check_index(lsIndex);
  return convert(krino::get_levelset_shape_sensitivities(impl.bulk(), impl.lsFields[lsIndex], impl.doComputeClosestPointSensitivities));
}

std::vector<std::vector<ShapeSensitivity>> get_all_shape_sensitivities(const LevelSetMesh & mesh)
{
  const LevelSetMesh::Impl & impl = mesh.impl();
  const std::vector<std::vector<krino::LevelSetShapeSensitivity>> allSens =
      krino::get_levelset_shape_sensitivities(impl.bulk(), impl.lsFields, impl.doComputeClosestPointSensitivities);

  std::vector<std::vector<ShapeSensitivity>> out;
  out.reserve(allSens.size());
  for (auto & perLevelSet : allSens)
    out.push_back(convert(perLevelSet));
  return out;
}

//--- Workflow B: write -------------------------------------------------------

void write_decomposed_mesh(LevelSetMesh & mesh, const std::string & fileName, const bool includeVoidPhase, const int step, const double time)
{
  LevelSetMesh::Impl & impl = mesh.impl();
  const stk::mesh::Selector outputSelector = build_output_selector(impl.meta(), krino::AuxMetaData::get(impl.meta()).active_part(), includeVoidPhase);
  krino::fix_ownership_and_output_composed_mesh_with_fields(impl.bulk(), outputSelector, fileName, step, time);
}

} // namespace PlatoKrinoInterface

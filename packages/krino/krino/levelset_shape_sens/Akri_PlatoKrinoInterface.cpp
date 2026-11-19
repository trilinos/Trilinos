#include <Akri_PlatoKrinoInterface.hpp>

#include <algorithm>
#include <cctype>
#include <memory>
#include <set>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <Akri_AuxMetaData.hpp>
#include <Akri_AnalyticSurf.hpp>
#include <Akri_BoundingBoxMesh.hpp>
#include <Akri_CSG_Surface.hpp>
#include <Akri_CSG_Tree.hpp>
#include <Akri_FacetFileReader.hpp>
#include <Akri_FacetFileSurface.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_LevelSetPolicy.hpp>
#include <Akri_LevelSetShapeSensitivities.hpp>
#include <Akri_MeshCloneUtils.hpp>
#include <Akri_MeshFromFile.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_NodalSurfaceDistance.hpp>
#include <Akri_OutputUtils.hpp>
#include <Akri_ParallelErrorMessage.hpp>
#include <Akri_Phase_Support.hpp>

namespace PlatoKrinoInterface {

namespace {

stk::math::Vector3d to_vec(const std::array<double,3> & a) { return stk::math::Vector3d(a[0], a[1], a[2]); }

// Logs to err why levelSetField can't be a level set field of mesh: it must be a nodal double field
// of mesh, on every node.
void check_level_set_field(krino::ParallelErrorMessage & err,
    const stk::mesh::FieldBase & levelSetField,
    const stk::mesh::BulkData & mesh)
{
  const std::string fieldDescription = "Level set field " + levelSetField.name();
  if (&levelSetField.get_mesh() != &mesh)
  {
    err << fieldDescription << " is a field of another mesh.\n";
    return;
  }
  if (levelSetField.entity_rank() != stk::topology::NODE_RANK)
  {
    err << fieldDescription << " is not a nodal field.\n";
    return;
  }
  if (!levelSetField.type_is<double>())
  {
    err << fieldDescription << " is not a double field.\n";
    return;
  }
  const stk::mesh::BucketVector & bucketsWithoutField =
      mesh.get_buckets(stk::topology::NODE_RANK, !stk::mesh::selectField(levelSetField));
  if (!bucketsWithoutField.empty())
    err << fieldDescription << " is not on every node, e.g. node "
        << mesh.identifier((*bucketsWithoutField.front())[0]) << " on processor "
        << mesh.parallel_rank() << ".\n";
}

std::vector<stk::mesh::FieldBase *>
declare_level_set_fields(stk::mesh::MetaData & meta, const unsigned numLevelSets)
{
  std::vector<stk::mesh::FieldBase *> fields;
  for (unsigned i = 0; i < numLevelSets; ++i)
  {
    const std::string name = (numLevelSets == 1) ? "LS" : "LS" + std::to_string(i);
    auto & field = meta.declare_field<double>(stk::topology::NODE_RANK, name);
    stk::mesh::put_field_on_mesh(field, meta.universal_part(), nullptr);
    fields.push_back(&field);
  }
  return fields;
}

// The upper-case extension of fileName, which names its facet file format.
std::string facet_file_format(const std::string & fileName)
{
  const size_t dot = fileName.find_last_of('.');
  std::string format = (dot == std::string::npos) ? "" : fileName.substr(dot + 1);
  std::transform(format.begin(),
      format.end(),
      format.begin(),
      [](const unsigned char c) { return static_cast<char>(std::toupper(c)); });
  return format;
}

// The conforming blocks of the material, where every level set is positive.
stk::mesh::Selector material_blocks_selector(
    const stk::mesh::MetaData & meta, const std::vector<krino::LS_Field> & lsFields)
{
  krino::PhaseTag materialPhase;
  for (const auto & lsField : lsFields)
    materialPhase.add(lsField.identifier, +1);
  return krino::Phase_Support::get(meta).get_selector_of_blocks_with_phase(materialPhase);
}

} // anonymous namespace

//--- background mesh ---------------------------------------------------------

class BackgroundMesh::Impl
{
public:
  std::unique_ptr<krino::MeshInterface> mesh;
  std::vector<stk::mesh::FieldBase *> levelSetFields;
};

BackgroundMesh::BackgroundMesh(std::unique_ptr<Impl> impl) : myImpl(std::move(impl)) {}
BackgroundMesh::~BackgroundMesh() = default;
BackgroundMesh::BackgroundMesh(BackgroundMesh &&) noexcept = default;
BackgroundMesh & BackgroundMesh::operator=(BackgroundMesh &&) noexcept = default;

BackgroundMesh BackgroundMesh::read(const std::string & fileName,
    MPI_Comm comm,
    const std::string & decompMethod,
    const BackgroundMeshOptions & options)
{
  auto impl = std::make_unique<Impl>();
  impl->mesh = std::make_unique<krino::MeshFromFile>(fileName, comm, decompMethod);
  impl->levelSetFields = declare_level_set_fields(impl->mesh->meta_data(), options.numLevelSets);
  impl->mesh->populate_mesh();
  return BackgroundMesh(std::move(impl));
}

BackgroundMesh BackgroundMesh::bounding_box(const unsigned dim,
    const std::array<double, 3> & minCorner,
    const std::array<double, 3> & maxCorner,
    const double meshSize,
    MPI_Comm comm,
    const BackgroundMeshOptions & options)
{
  STK_ThrowRequireMsg(
      dim == 2 || dim == 3, "BackgroundMesh::bounding_box: dimension must be 2 or 3.");
  auto impl = std::make_unique<Impl>();
  auto bboxMesh = std::make_unique<krino::BoundingBoxMesh>(
      (dim == 2) ? stk::topology::TRIANGLE_3_2D : stk::topology::TETRAHEDRON_4, comm);
  bboxMesh->set_mesh_structure_type((dim == 2) ?
          krino::FLAT_WALLED_TRIANGULAR_LATTICE_BOUNDING_BOX_MESH :
          krino::FLAT_WALLED_BCC_BOUNDING_BOX_MESH);
  bboxMesh->set_domain(
      krino::BoundingBoxMesh::BoundingBoxType(to_vec(minCorner), to_vec(maxCorner)), meshSize);
  impl->levelSetFields = declare_level_set_fields(bboxMesh->meta_data(), options.numLevelSets);
  bboxMesh->populate_mesh();
  impl->mesh = std::move(bboxMesh);
  return BackgroundMesh(std::move(impl));
}

stk::mesh::BulkData & BackgroundMesh::bulk_data() { return myImpl->mesh->bulk_data(); }
const stk::mesh::BulkData & BackgroundMesh::bulk_data() const { return myImpl->mesh->bulk_data(); }
unsigned BackgroundMesh::num_level_sets() const { return myImpl->levelSetFields.size(); }

stk::mesh::FieldBase & BackgroundMesh::level_set_field(const unsigned lsIndex) const
{
  STK_ThrowRequireMsg(lsIndex < myImpl->levelSetFields.size(),
      "BackgroundMesh: level set index " << lsIndex << " out of range ("
                                         << myImpl->levelSetFields.size() << " level sets).");
  return *myImpl->levelSetFields[lsIndex];
}

//--- level set initializers --------------------------------------------------

void initialize_field_from_geometry(const stk::mesh::BulkData & mesh,
    stk::mesh::FieldBase & levelSetField,
    const Geometry & geometry)
{
  krino::ParallelErrorMessage err(mesh.parallel());
  check_level_set_field(err, levelSetField, mesh);
  if (geometry.spheres.empty() && geometry.planes.empty() && geometry.facetFiles.empty())
    err << "The geometry has no spheres, planes or facet files.\n";
  for (const auto & sphere : geometry.spheres)
    if (!(sphere.radius > 0.)) err << "Sphere radius " << sphere.radius << " is not positive.\n";
  for (const auto & plane : geometry.planes)
    if (to_vec(plane.normal).length() == 0.) err << "Plane normal is zero.\n";
  for (const auto & facetFile : geometry.facetFiles)
    if (!krino::is_facet_file_format(facet_file_format(facetFile.fileName)))
      err << "Facet file " << facetFile.fileName
          << " is not STL, FAC, PLY or EXO (from its extension).\n";
  if (!geometry.facetFiles.empty() && mesh.mesh_meta_data().spatial_dimension() != 3)
    err << "Facet files need a 3D mesh.\n";
  err.require_empty("initialize_field_from_geometry: invalid arguments.");

  std::vector<std::shared_ptr<krino::Surface>> surfaces;
  for (const auto & sphere : geometry.spheres)
    surfaces.emplace_back(std::make_shared<krino::Sphere>(to_vec(sphere.center), sphere.radius));
  // krino::Plane unitizes its normal, so scale the offset with it.
  for (const auto & plane : geometry.planes)
  {
    const stk::math::Vector3d normal = to_vec(plane.normal);
    const double normalLength = normal.length();
    surfaces.emplace_back(
        std::make_shared<krino::Plane>(normal / normalLength, plane.offset / normalLength));
  }
  // Like the other surfaces, a facet file surface of sign 1 is negative inside.
  for (const auto & facetFile : geometry.facetFiles)
    surfaces.emplace_back(std::make_shared<krino::FacetFileSurface>(mesh.parallel(),
        facet_file_format(facetFile.fileName),
        facetFile.fileName,
        facetFile.materialInside ? 1 : -1,
        to_vec(facetFile.scale)));

  // Negate the signed distance to the union, which is negative inside it, so that the material is
  // positive.
  krino::CSG_Surface csgSurface(krino::CSGTree::unite(surfaces));
  krino::compute_nodal_surface_distance(
      mesh, mesh.mesh_meta_data().coordinate_field(), krino::FieldRef(levelSetField), csgSurface);
  stk::mesh::field_scale(-1., levelSetField);
}

void initialize_field_from_function(const stk::mesh::BulkData & mesh,
    stk::mesh::FieldBase & levelSetField,
    const LevelSetValueFunction & valueFn)
{
  krino::ParallelErrorMessage err(mesh.parallel());
  check_level_set_field(err, levelSetField, mesh);
  err.require_empty("initialize_field_from_function: invalid arguments.");

  const unsigned nDim = mesh.mesh_meta_data().spatial_dimension();
  const krino::FieldRef coordsField = mesh.mesh_meta_data().coordinate_field();
  const krino::FieldRef distanceField(levelSetField);

  const stk::mesh::BucketVector & buckets =
      mesh.get_buckets(stk::topology::NODE_RANK, stk::mesh::selectField(distanceField));
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

  stk::mesh::communicate_field_data(mesh, {&levelSetField});
}

void set_field_values(const stk::mesh::BulkData & mesh,
    stk::mesh::FieldBase & levelSetField,
    const std::vector<stk::mesh::EntityId> & nodeIds,
    const std::vector<double> & values)
{
  krino::ParallelErrorMessage err(mesh.parallel());
  check_level_set_field(err, levelSetField, mesh);
  if (nodeIds.size() != values.size())
  {
    err << nodeIds.size() << " node ids but " << values.size() << " values on processor "
        << mesh.parallel_rank() << ".\n";
  }
  else
  {
    constexpr size_t maxNumMissingNodesToReport = 5;
    size_t numMissingNodes = 0;
    for (const stk::mesh::EntityId nodeId : nodeIds)
    {
      if (mesh.is_valid(mesh.get_entity(stk::topology::NODE_RANK, nodeId))) continue;
      if (numMissingNodes++ < maxNumMissingNodesToReport)
        err << "Node " << nodeId << " is not on processor " << mesh.parallel_rank() << ".\n";
    }
    if (numMissingNodes > maxNumMissingNodesToReport)
      err << numMissingNodes - maxNumMissingNodesToReport << " more nodes are not on processor "
          << mesh.parallel_rank() << ".\n";
  }
  err.require_empty("set_field_values: invalid arguments.");

  auto levelSetData = krino::FieldRef(levelSetField).data<double, stk::mesh::ReadWrite>();
  for (size_t i = 0; i < nodeIds.size(); ++i)
    levelSetData.entity_values(mesh.get_entity(stk::topology::NODE_RANK, nodeIds[i]))() = values[i];

  stk::mesh::communicate_field_data(mesh, {&levelSetField});
}

//--- conforming mesh ---------------------------------------------------------

class ConformingMeshWithSensitivities::Impl
{
public:
  std::shared_ptr<stk::mesh::MetaData> meta;
  std::unique_ptr<stk::mesh::BulkData> mesh;
  std::unique_ptr<krino::ConformingDecompositionWithSensitivities> decomposition;
  krino::ShapeSensitivityMethod sensitivityMethod =
      krino::ShapeSensitivityMethod::EDGE_DIRECTION_USING_EDGE_GRADIENT;

  void check_index(const unsigned lsIndex) const
  {
    const size_t numLevelSets = decomposition->levelset_fields().size();
    STK_ThrowRequireMsg(lsIndex < numLevelSets,
        "PlatoKrinoInterface: level set index " << lsIndex << " out of range (" << numLevelSets
                                                << " level sets).");
  }
};

// Copies the level set values of the copied background field into the krino level set field, if
// they are different fields (the background field can have the krino field's name).
static void copy_level_set_values(stk::mesh::BulkData & mesh,
    const stk::mesh::FieldBase & copiedBackgroundField,
    const krino::FieldRef levelSetField)
{
  if (&copiedBackgroundField == &levelSetField.field()) return;

  const krino::FieldRef fromField(copiedBackgroundField);
  auto fromData = fromField.data<double>();
  auto toData = levelSetField.data<double, stk::mesh::ReadWrite>();
  for (auto * bucketPtr : mesh.get_buckets(stk::topology::NODE_RANK,
           stk::mesh::selectField(fromField) & stk::mesh::selectField(levelSetField)))
    for (auto node : *bucketPtr)
      toData.entity_values(node)() = fromData.entity_values(node)();
  stk::mesh::communicate_field_data(mesh, {&levelSetField.field()});
}

ConformingMeshWithSensitivities::ConformingMeshWithSensitivities(
    const stk::mesh::BulkData & backgroundMesh,
    const std::vector<const stk::mesh::FieldBase *> & levelSetFields,
    const DecompositionOptions & options)
    : myImpl(std::make_unique<Impl>())
{
  krino::ParallelErrorMessage err(backgroundMesh.parallel());
  if (levelSetFields.empty()) err << "There are no level set fields.\n";
  for (size_t i = 0; i < levelSetFields.size(); ++i)
  {
    if (levelSetFields[i] == nullptr)
      err << "Level set field " << i << " is null.\n";
    else
      check_level_set_field(err, *levelSetFields[i], backgroundMesh);
  }
  err.require_empty("ConformingMeshWithSensitivities: invalid level set fields.");

  Impl & impl = *myImpl;
  impl.sensitivityMethod = options.sensitivityMethod;

  stk::mesh::MeshBuilder builder(backgroundMesh.parallel());
  builder.set_aura_option(backgroundMesh.is_automatic_aura_on() ?
          stk::mesh::BulkData::AUTO_AURA :
          stk::mesh::BulkData::NO_AUTO_AURA);
  impl.meta = builder.create_meta_data();
  krino::clone_meta_data_parts_and_fields(backgroundMesh.mesh_meta_data(), *impl.meta);
  krino::clone_io_attributes_and_surface_to_block_mapping(
      backgroundMesh.mesh_meta_data(), *impl.meta);
  krino::AuxMetaData & auxMeta = krino::AuxMetaData::create(*impl.meta);

  const std::vector<krino::LS_Field> lsFields = krino::LSPerInterfacePolicy::
      setup_levelsets_on_all_blocks_with_void_phase_for_any_negative_levelset(
          *impl.meta, levelSetFields.size());
  // With several level sets, each level set's interface is restricted to the material boundary.
  const stk::mesh::Selector materialBlocks = material_blocks_selector(*impl.meta, lsFields);
  impl.decomposition = std::make_unique<krino::ConformingDecompositionWithSensitivities>(
      *impl.meta, lsFields, options.doSnapping, &materialBlocks);

  impl.mesh = builder.create(impl.meta);
  impl.meta->commit();
  krino::clone_bulk_data(backgroundMesh, *impl.mesh);
  // krino needs the exposed and interior block boundary sides, which a mesh krino reads has but the
  // background mesh need not.
  stk::mesh::create_exposed_block_boundary_sides(
      *impl.mesh, impl.meta->universal_part(), {&auxMeta.exposed_boundary_part()});
  stk::mesh::create_interior_block_boundary_sides(
      *impl.mesh, impl.meta->universal_part(), {&auxMeta.block_boundary_part()});
  krino::activate_all_entities(*impl.mesh, auxMeta.active_part());

  for (size_t i = 0; i < levelSetFields.size(); ++i)
  {
    const stk::mesh::FieldBase * copiedField =
        impl.meta->get_field(stk::topology::NODE_RANK, levelSetFields[i]->name());
    STK_ThrowRequireMsg(copiedField != nullptr,
        "ConformingMeshWithSensitivities: the copy has no level set field "
            << levelSetFields[i]->name() << ".");
    copy_level_set_values(*impl.mesh, *copiedField, lsFields[i].isovar);
  }

  impl.decomposition->decompose();
}

ConformingMeshWithSensitivities::~ConformingMeshWithSensitivities() = default;
ConformingMeshWithSensitivities::ConformingMeshWithSensitivities(
    ConformingMeshWithSensitivities &&) noexcept = default;
ConformingMeshWithSensitivities &
ConformingMeshWithSensitivities::operator=(ConformingMeshWithSensitivities &&) noexcept = default;

stk::mesh::BulkData & ConformingMeshWithSensitivities::bulk_data() { return *myImpl->mesh; }
const stk::mesh::BulkData & ConformingMeshWithSensitivities::bulk_data() const
{
  return *myImpl->mesh;
}
unsigned ConformingMeshWithSensitivities::num_level_sets() const
{
  return myImpl->decomposition->levelset_fields().size();
}

std::vector<krino::LevelSetShapeSensitivity>
ConformingMeshWithSensitivities::shape_sensitivities(const unsigned lsIndex) const
{
  const Impl & impl = *myImpl;
  impl.check_index(lsIndex);
  return impl.decomposition->sensitivities(lsIndex, impl.sensitivityMethod);
}

std::vector<std::vector<krino::LevelSetShapeSensitivity>>
ConformingMeshWithSensitivities::all_shape_sensitivities() const
{
  const Impl & impl = *myImpl;
  return impl.decomposition->all_sensitivities(impl.sensitivityMethod);
}

stk::mesh::Selector ConformingMeshWithSensitivities::material_selector() const
{
  const Impl & impl = *myImpl;
  return krino::AuxMetaData::get(*impl.meta).active_part() &
      material_blocks_selector(*impl.meta, impl.decomposition->levelset_fields());
}

stk::mesh::Selector ConformingMeshWithSensitivities::void_selector() const
{
  const Impl & impl = *myImpl;
  const krino::Phase_Support & phaseSupport = krino::Phase_Support::get(*impl.meta);
  std::set<unsigned> voidBlockOrdinals;
  for (const auto & lsField : impl.decomposition->levelset_fields())
    for (const unsigned blockOrdinal :
        phaseSupport.get_negative_levelset_block_ordinals(lsField.identifier))
      voidBlockOrdinals.insert(blockOrdinal);
  stk::mesh::PartVector voidBlocks;
  for (const unsigned blockOrdinal : voidBlockOrdinals)
    voidBlocks.push_back(&impl.meta->get_part(blockOrdinal));
  return krino::AuxMetaData::get(*impl.meta).active_part() & stk::mesh::selectUnion(voidBlocks);
}

stk::mesh::Selector
ConformingMeshWithSensitivities::interface_selector(const unsigned lsIndex) const
{
  const Impl & impl = *myImpl;
  impl.check_index(lsIndex);
  return material_selector() &
      krino::Phase_Support::get(*impl.meta)
          .get_negative_levelset_interface_selector(
              impl.decomposition->levelset_fields()[lsIndex].identifier);
}

template <size_t NNODES>
krino::InterfaceFacetsWithSensitivities<NNODES>
ConformingMeshWithSensitivities::interface_facets(const unsigned lsIndex) const
{
  const Impl & impl = *myImpl;
  impl.check_index(lsIndex);
  return impl.decomposition->facets_and_sensitivities<NNODES>(lsIndex, impl.sensitivityMethod);
}

template krino::InterfaceFacetsWithSensitivities<2>
ConformingMeshWithSensitivities::interface_facets<2>(const unsigned lsIndex) const;
template krino::InterfaceFacetsWithSensitivities<3>
ConformingMeshWithSensitivities::interface_facets<3>(const unsigned lsIndex) const;

void ConformingMeshWithSensitivities::write_decomposed_mesh(
    const std::string & fileName, const bool includeVoidPhase, const int step, const double time)
{
  Impl & impl = *myImpl;
  const stk::mesh::Selector outputSelector = includeVoidPhase ?
      stk::mesh::Selector(krino::AuxMetaData::get(*impl.meta).active_part()) :
      material_selector();
  krino::fix_ownership_and_output_composed_mesh_with_fields(
      *impl.mesh, outputSelector, fileName, step, time);
}

} // namespace PlatoKrinoInterface

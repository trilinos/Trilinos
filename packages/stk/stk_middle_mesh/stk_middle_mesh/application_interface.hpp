#ifndef APPLICATION_INTERFACE_H
#define APPLICATION_INTERFACE_H

#include "incremental_mesh_boundary_snapper_opts.hpp"
#include "mesh_quality_improver_opts.hpp"
#include "mesh_snapper_opts.hpp"
#include "nonconformal4_opts.hpp"
#include "bounding_box_search_opts.hpp"

namespace stk {
namespace middle_mesh {

enum class ParallelSearchType
{
  StkSearch = 0
};

struct ParallelSearchOpts
{
    ParallelSearchType type = ParallelSearchType::StkSearch;

    // TODO: are there options?
};

enum class VolumeSnapType
{
  None = 0,
  Standard
};

struct VolumeSnapOpts
{
    VolumeSnapType type = VolumeSnapType::None;

    MeshSnapperOpts standardOpts;
};

enum class BoundarySnapAndQualityImprovementType
{
  None = 0,
  SnapThenQuality,
  IncrementalBoundarySnap,
};

struct BoundarySnapAndQualityImprovementOpts
{
    BoundarySnapAndQualityImprovementType type = BoundarySnapAndQualityImprovementType::IncrementalBoundarySnap;

    MeshQualityImproverOpts snapThenQualityOpts;
    IncrementalBoundarySnapperOpts incrementalMeshBoundarySnapOpts;
};

enum class MiddleGridType
{
  NormalProjection = 0
};

struct MiddleGridOpts
{
    MiddleGridType type = MiddleGridType::NormalProjection;

    NormalProjectionOpts normalProjectionOpts;
    search::BoundingBoxSearchOpts searchOpts;
};

class XiCoordinates
{
  public:
    virtual ~XiCoordinates() = default;

    virtual const std::vector<utils::Point>& get_xi_coords(mesh::MeshEntityType type) = 0;

    virtual std::pair<double, double> get_xi_coord_range(mesh::MeshEntityType type) = 0;
};


// Class for constructing a middle grid when a given set of processes has
// two meshes
class ApplicationInterface
{
  public:
    virtual ~ApplicationInterface() {}

    // creates the middle grid
    virtual void create_middle_grid() = 0;

    // returns a parallel middle grid.  The local part of the middle grid
    // overlaps the local part of mesh1
    virtual std::shared_ptr<mesh::Mesh> get_middle_grid_for_mesh1() = 0;

    // returns a parallel middle grid.  The local part of the middle grid
    // overlaps the local part of mesh2
    virtual std::shared_ptr<mesh::Mesh> get_middle_grid_for_mesh2() = 0;

    // returns a field that maps the elements of mesh_in to the elements of mesh1.
    virtual mesh::FieldPtr<mesh::MeshEntityPtr> get_mesh1_classification() = 0;

    // returns a field that maps the elements of mesh_in to the elements of mesh2.
    virtual mesh::FieldPtr<mesh::MeshEntityPtr> get_mesh2_classification() = 0;

    // returns a field that maps the elements of mesh1 to the elements of the middle grid.
    virtual mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> compute_mesh1_inverse_classification() = 0;

    // returns a field that maps the elements of mesh2 to the elements of the middle grid.
    virtual mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> compute_mesh2_inverse_classification() = 0;

    // returns a field that maps the elements of the mesh returned by get_middle_grid_for_mesh1()
    // to the corresponding elements on the mesh returned by get_middle_grid_for_mesh2()
    virtual mesh::FieldPtr<mesh::RemoteSharedEntity> get_remote_info_mesh_one_to_two() = 0;

    // returns a field that maps the elements of the mesh returned by get_middle_grid_for_mesh2()
    // to the corresponding elements on the mesh returned by get_middle_grid_for_mesh1()
    virtual mesh::FieldPtr<mesh::RemoteSharedEntity> get_remote_info_mesh_two_to_one() = 0;

    // if a XiCoordinates object was passed into the constructor, returns a field
    // containing the xi coordinates in each middle mesh element  projected onto mesh1
    virtual mesh::FieldPtr<utils::Point> get_xi_points_on_mesh1() = 0;

    // similar to get_xi_points_on_mesh1(), but returns the points projected onto mesh2
    virtual mesh::FieldPtr<utils::Point> get_xi_points_on_mesh2() = 0;
};

enum class ApplicationInterfaceType
{
  FakeParallel = 0,
  Parallel
};

std::shared_ptr<ApplicationInterface> application_interface_factory(
    ApplicationInterfaceType type, std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
    MPI_Comm unionComm,
    std::shared_ptr<XiCoordinates> xiPts = nullptr,
    const ParallelSearchOpts& parallelSearchOpts                  = ParallelSearchOpts(),
    const VolumeSnapOpts& volumeSnapOpts                          = VolumeSnapOpts(),
    const BoundarySnapAndQualityImprovementOpts& boundarySnapOpts = BoundarySnapAndQualityImprovementOpts(),
    const MiddleGridOpts& middleGridOpts                          = MiddleGridOpts());

} // namespace middle_mesh
} // namespace stk

#endif

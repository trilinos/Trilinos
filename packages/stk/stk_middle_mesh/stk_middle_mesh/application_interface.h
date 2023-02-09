#ifndef APPLICATION_INTERFACE_H
#define APPLICATION_INTERFACE_H

// TODO: consider moving the options structure to their own files, so
//       we can include only the struct and not the class
#include "incremental_mesh_boundary_snapper_opts.h"
#include "mesh_quality_improver_opts.h"
#include "mesh_snapper_opts.h"
#include "nonconformal4_opts.h"

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
};

// Class for constructing a middle grid when two non-overlapping sets of
// processes have one mesh each.
class ApplicationInterfaceMPMD
{
  public:
    virtual ~ApplicationInterfaceMPMD() {}

    // returns a parallel mesh that contains the middle grid elements
    // classified on either mesh1 or mesh2, which ever was non-null when passed into
    // the constructor
    // If both are null, returns a nullptr
    virtual std::shared_ptr<mesh::Mesh> create_middle_grid() = 0;

    // returns a field that maps the elements of mesh_in to the elements of mesh1 or mesh 2, whichever
    // was not null when passed into the constructor.
    virtual mesh::FieldPtr<mesh::MeshEntityPtr> get_mesh_classification() = 0;

    // returns a field that maps the elements of mesh1 or mesh 2, whichever was not null when passed
    // into the constructor, to the elements of the middle grid.
    virtual mesh::VariableSizeFieldPtr<mesh::MeshEntityPtr> compute_mesh_inverse_classification() = 0;
};

// Class for constructing a middle grid when a given set of processes has
// two meshes
class ApplicationInterfaceSPMD
{
  public:
    virtual ~ApplicationInterfaceSPMD() {}

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
};

enum class ApplicationInterfaceType
{
  FakeParallel = 0
};

std::shared_ptr<ApplicationInterfaceMPMD> application_interface_mpmd_factory(
    ApplicationInterfaceType type, std::shared_ptr<mesh::Mesh> mesh, bool isMesh1, MPI_Comm unionComm,
    const ParallelSearchOpts& parallelSearchOpts                  = ParallelSearchOpts(),
    const VolumeSnapOpts& volumeSnapOpts                          = VolumeSnapOpts(),
    const BoundarySnapAndQualityImprovementOpts& boundarySnapOpts = BoundarySnapAndQualityImprovementOpts(),
    const MiddleGridOpts& middleGridOpts                          = MiddleGridOpts());

std::shared_ptr<ApplicationInterfaceSPMD> application_interface_spmd_factory(
    ApplicationInterfaceType type, std::shared_ptr<mesh::Mesh> mesh1, std::shared_ptr<mesh::Mesh> mesh2,
    const ParallelSearchOpts& parallelSearchOpts                  = ParallelSearchOpts(),
    const VolumeSnapOpts& volumeSnapOpts                          = VolumeSnapOpts(),
    const BoundarySnapAndQualityImprovementOpts& boundarySnapOpts = BoundarySnapAndQualityImprovementOpts(),
    const MiddleGridOpts& middleGridOpts                          = MiddleGridOpts());

} // namespace middle_mesh
} // namespace stk

#endif
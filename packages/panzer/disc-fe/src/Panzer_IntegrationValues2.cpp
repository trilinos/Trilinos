// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_IntegrationValues2.hpp"
#include "Panzer_UtilityAlgs.hpp"

#include "Shards_CellTopology.hpp"

#include "Kokkos_DynRankView.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_RealSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_CubatureControlVolume.hpp"
#include "Intrepid2_CubatureControlVolumeSide.hpp"
#include "Intrepid2_CubatureControlVolumeBoundary.hpp"

#include "Panzer_CommonArrayFactories.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_SubcellConnectivity.hpp"
#include "Panzer_ConvertNormalToRotationMatrix.hpp"

// FIXME: There are some calls in Intrepid2 that require non-const arrays when they should be const - search for PHX::getNonConstDynRankViewFromConstMDField
#include "Phalanx_GetNonConstDynRankViewFromConstMDField.hpp"

namespace panzer {

namespace {

Teuchos::RCP<Intrepid2::Cubature<PHX::Device::execution_space,double,double>>
getIntrepidCubature(const panzer::IntegrationRule & ir)
{
  typedef panzer::IntegrationDescriptor ID;
  Teuchos::RCP<Intrepid2::Cubature<PHX::Device::execution_space,double,double> > ic;

  Intrepid2::DefaultCubatureFactory cubature_factory;

  if(ir.getType() == ID::CV_SIDE){
    ic = Teuchos::rcp(new Intrepid2::CubatureControlVolumeSide<PHX::Device::execution_space,double,double>(*ir.topology));
  } else if(ir.getType() == ID::CV_VOLUME){
    ic = Teuchos::rcp(new Intrepid2::CubatureControlVolume<PHX::Device::execution_space,double,double>(*ir.topology));
  } else if(ir.getType() == ID::CV_BOUNDARY){
    TEUCHOS_ASSERT(ir.isSide());
    ic = Teuchos::rcp(new Intrepid2::CubatureControlVolumeBoundary<PHX::Device::execution_space,double,double>(*ir.topology,ir.getSide()));
  } else if(ir.getType() == ID::VOLUME){
    ic = cubature_factory.create<PHX::Device::execution_space,double,double>(*(ir.topology),ir.getOrder());
  } else if(ir.getType() == ID::SIDE){
    ic = cubature_factory.create<PHX::Device::execution_space,double,double>(*(ir.side_topology),ir.getOrder());
  } else if(ir.getType() == ID::SURFACE){
    // closed surface integrals don't exist in intrepid.
  } else {
    TEUCHOS_ASSERT(false);
  }

  return ic;
}

template<typename Scalar>
void
correctVirtualNormals(PHX::MDField<Scalar,Cell,IP,Dim> normals,
                      const int num_real_cells,
                      const int num_virtual_cells,
                      const shards::CellTopology & cell_topology,
                      const SubcellConnectivity & face_connectivity)
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::correctVirtualNormals()",corr_virt_norms);

  // What we want is for the adjoining face of the virtual cell to have normals that are the negated real cell's normals.
  // we correct the normals here:
  const int space_dim           = cell_topology.getDimension();
  const int faces_per_cell      = cell_topology.getSubcellCount(space_dim-1);
  const int points              = normals.extent_int(1);
  const int points_per_face     = points / faces_per_cell;

  Kokkos::parallel_for(num_virtual_cells, KOKKOS_LAMBDA(const int & virtual_cell_ordinal){

    const int virtual_cell = virtual_cell_ordinal+num_real_cells;

    // Find the face and local face ids for the given virtual cell
    // Note that virtual cells only connect to the owned cells through a single face
    int face, virtual_lidx;
    for (int local_face_id=0; local_face_id<faces_per_cell; local_face_id++){
      // Faces that exist have positive indexes
      face = face_connectivity.subcellForCell(virtual_cell, local_face_id);
      if (face >= 0){
        virtual_lidx = local_face_id;
        break;
      }
    }

    // Indexes for real cell
    const int real_side = (face_connectivity.cellForSubcell(face, 0) == virtual_cell) ? 1 : 0;
    const int real_cell = face_connectivity.cellForSubcell(face,real_side);
    const int real_lidx = face_connectivity.localSubcellForSubcell(face,real_side);

    // Make sure it is a real cell (it should actually be an owned cell)
    KOKKOS_ASSERT(real_cell < num_real_cells);

    for(int point_ordinal=0; point_ordinal<points_per_face; point_ordinal++){

      // Point indexes for virtual and real point on face
      const int virtual_cell_point = points_per_face * virtual_lidx + point_ordinal;
      const int real_cell_point = points_per_face * real_lidx + point_ordinal;

      for (int d=0; d<space_dim; d++)
        normals(virtual_cell,virtual_cell_point,d) = -normals(real_cell,real_cell_point,d);

    }

    // Clear other normals
    for (int local_face_id=0; local_face_id<faces_per_cell; local_face_id++){

      if (local_face_id == virtual_lidx) continue;

      for (int point_ordinal=0; point_ordinal<points_per_face; point_ordinal++){
        const int point = local_face_id * points_per_face + point_ordinal;
        for (int dim=0; dim<space_dim; dim++)
          normals(virtual_cell,point,dim) = 0.0;
      }
    }
  });
  PHX::Device::execution_space().fence();
}


template<typename Scalar>
void
correctVirtualRotationMatrices(PHX::MDField<Scalar,Cell,IP,Dim,Dim> rotation_matrices,
                               const int num_real_cells,
                               const int num_virtual_cells,
                               const shards::CellTopology & cell_topology,
                               const SubcellConnectivity & face_connectivity)
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::correctVirtualRotationMatrices()",corr_virt_rotmat);

  // What we want is for the adjoining face of the virtual cell to have normals that are the negated real cell's normals.
  // we correct the normals here:
  const int space_dim           = cell_topology.getDimension();
  const int faces_per_cell      = cell_topology.getSubcellCount(space_dim-1);
  const int points              = rotation_matrices.extent_int(1);
  const int points_per_face     = points / faces_per_cell;

  Kokkos::parallel_for(num_virtual_cells, KOKKOS_LAMBDA(const int & virtual_cell_ordinal){

    const int virtual_cell = virtual_cell_ordinal+num_real_cells;

    // Find the face and local face ids for the given virtual cell
    // Note that virtual cells only connect to the owned cells through a single face
    int face, virtual_lidx;
    for (int local_face_id=0; local_face_id<faces_per_cell; local_face_id++){
      // Faces that exist have positive indexes
      face = face_connectivity.subcellForCell(virtual_cell, local_face_id);
      if (face >= 0){
        virtual_lidx = local_face_id;
        break;
      }
    }

    // The normals already have the correction applied, so we just need to zero out the rotation matrices on the other faces

    // Clear other rotation matrices
    for (int local_face_id=0; local_face_id<faces_per_cell; local_face_id++){

      if (local_face_id == virtual_lidx) continue;

      for (int point_ordinal=0; point_ordinal<points_per_face; point_ordinal++){
        const int point = local_face_id * points_per_face + point_ordinal;
        for (int dim=0; dim<3; dim++)
          for (int dim2=0; dim2<3; dim2++)
            rotation_matrices(virtual_cell,point,dim,dim2) = 0.0;
      }
    }
  });
  PHX::Device::execution_space().fence();
}

template<typename Scalar>
void
applyBasePermutation(PHX::MDField<Scalar,IP> field,
                     PHX::MDField<const int,Cell,IP> permutations)
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::applyBasePermutation(rank 1)",app_base_perm_r1);
  MDFieldArrayFactory af("",true);

  const int num_ip = field.extent(0);

  auto scratch = af.template buildStaticArray<Scalar,IP>("scratch", num_ip);
  scratch.deep_copy(field);
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int & ){
    for(int ip=0; ip<num_ip; ++ip)
      if (ip != permutations(0,ip))
        field(ip) = scratch(permutations(0,ip));
  });
  PHX::Device::execution_space().fence();
}

template<typename Scalar>
void
applyBasePermutation(PHX::MDField<Scalar,IP,Dim> field,
                     PHX::MDField<const int,Cell,IP> permutations)
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::applyBasePermutation(rank 2)",app_base_perm_r2);
  MDFieldArrayFactory af("",true);

  const int num_ip = field.extent(0);
  const int num_dim = field.extent(1);

  auto scratch = af.template buildStaticArray<Scalar,IP,Dim>("scratch", num_ip,num_dim);
  scratch.deep_copy(field);
  Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int & ){
    for(int ip=0; ip<num_ip; ++ip)
      if (ip != permutations(0,ip))
        for(int dim=0; dim<num_dim; ++dim)
        field(ip,dim) = scratch(permutations(0,ip),dim);
  });
  PHX::Device::execution_space().fence();
}

template<typename Scalar>
void
applyPermutation(PHX::MDField<Scalar,Cell,IP> field,
                 PHX::MDField<const int,Cell,IP> permutations)
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::applyPermutation(rank 2)",app_perm_r2);
  MDFieldArrayFactory af("",true);

  const int num_cells = field.extent(0);
  const int num_ip = field.extent(1);

  auto scratch = af.template buildStaticArray<Scalar,Cell,IP>("scratch", num_cells, num_ip);
  scratch.deep_copy(field);
  Kokkos::parallel_for(num_cells, KOKKOS_LAMBDA(const int & cell){
    for(int ip=0; ip<num_ip; ++ip)
      if (ip != permutations(cell,ip))
        field(cell,ip) = scratch(cell,permutations(cell,ip));
  });
  PHX::Device::execution_space().fence();
}

template<typename Scalar>
void
applyPermutation(PHX::MDField<Scalar,Cell,IP,Dim> field,
                 PHX::MDField<const int,Cell,IP> permutations)
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::applyPermutation(rank 3)",app_perm_r3);
  MDFieldArrayFactory af("",true);

  const int num_cells = field.extent(0);
  const int num_ip = field.extent(1);
  const int num_dim = field.extent(2);

  auto scratch = af.template buildStaticArray<Scalar,Cell,IP,Dim>("scratch", num_cells, num_ip, num_dim);
  scratch.deep_copy(field);
  Kokkos::parallel_for(num_cells, KOKKOS_LAMBDA(const int & cell){
    for(int ip=0; ip<num_ip; ++ip)
      if (ip != permutations(cell,ip))
        for(int dim=0; dim<num_dim; ++dim)
          field(cell,ip,dim) = scratch(cell,permutations(cell,ip),dim);
  });
  PHX::Device::execution_space().fence();
}

template<typename Scalar>
void
applyPermutation(PHX::MDField<Scalar,Cell,IP,Dim,Dim> field,
                 PHX::MDField<const int,Cell,IP> permutations)
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::applyPermutation(rank 4)",app_perm_r4);
  MDFieldArrayFactory af("",true);

  const int num_cells = field.extent(0);
  const int num_ip = field.extent(1);
  const int num_dim = field.extent(2);
  const int num_dim2 = field.extent(3);

  auto scratch = af.template buildStaticArray<Scalar,Cell,IP,Dim,Dim>("scratch", num_cells, num_ip, num_dim, num_dim2);
  scratch.deep_copy(field);
  Kokkos::parallel_for(num_cells, KOKKOS_LAMBDA(const int & cell){
    for(int ip=0; ip<num_ip; ++ip)
      if (ip != permutations(cell,ip))
        for(int dim=0; dim<num_dim; ++dim)
          for(int dim2=0; dim2<num_dim2; ++dim2)
            field(cell,ip,dim,dim2) = scratch(cell,permutations(cell,ip),dim,dim2);
  });
  PHX::Device::execution_space().fence();
}


// Find the permutation that maps the set of points coords to other_coords. To
// avoid possible finite precision issues, == is not used, but rather
// min(norm(.)).
template <typename Scalar>
PHX::MDField<const int,Cell,IP>
generatePermutations(const int num_cells,
                     PHX::MDField<const Scalar,Cell,IP,Dim> coords,
                     PHX::MDField<const Scalar,Cell,IP,Dim> other_coords)
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::generatePermutations()",gen_perms);

  const int num_ip = coords.extent(1);
  const int num_dim = coords.extent(2);

  MDFieldArrayFactory af("",true);

  // This will store the permutations to go from coords to other_coords
  auto permutation = af.template buildStaticArray<int,Cell,IP>("permutation", num_cells, num_ip);
  permutation.deep_copy(0);

  // This is scratch space - there is likely a better way to do this
  auto taken = af.template buildStaticArray<int,Cell,IP>("taken", num_cells, num_ip);
  taken.deep_copy(0);

  Kokkos::parallel_for(num_cells, KOKKOS_LAMBDA(const int & cell){

    for (int ip = 0; ip < num_ip; ++ip) {
      // Find an other point to associate with ip.
      int i_min = 0;
      Scalar d_min = -1;
      for (int other_ip = 0; other_ip < num_ip; ++other_ip) {
        // For speed, skip other points that are already associated.
        if(taken(cell,other_ip)) continue;
        // Compute the distance between the two points.
        Scalar d(0);
        for (int dim = 0; dim < num_dim; ++dim) {
          const Scalar diff = coords(cell, ip, dim) - other_coords(cell, other_ip, dim);
          d += diff*diff;
        }
        if (d_min < 0 || d < d_min) {
          d_min = d;
          i_min = other_ip;
        }
      }
      // Record the match.
      permutation(cell,ip) = i_min;
      // This point is now taken.
      taken(cell,i_min) = 1;
    }
  });
  PHX::Device::execution_space().fence();

  return permutation;

}

template <typename Scalar>
PHX::MDField<const int,Cell,IP>
generateSurfacePermutations(const int num_cells,
                            const SubcellConnectivity face_connectivity,
                            PHX::MDField<const Scalar,Cell,IP,Dim> surface_points,
                            PHX::MDField<const Scalar,Cell,IP,Dim,Dim> surface_rotation_matrices)

{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::generateSurfacePermutations()",gen_surf_perms);

  // The challenge for this call is handling wedge-based periodic boundaries
  // We need to make sure that we can align points along faces that are rotated with respect to one another.
  // Below we will see an algorithm that rotates two boundaries into a shared frame, then figures out
  // how to permute the points on one face to line up with the other.

  const int num_points = surface_points.extent_int(1);
  const int num_faces_per_cell = face_connectivity.numSubcellsOnCellHost(0);
  const int num_points_per_face = num_points / num_faces_per_cell;
  const int cell_dim = surface_points.extent(2);

  MDFieldArrayFactory af("",true);

  // This will store the permutations
  auto permutation = af.template buildStaticArray<int,Cell,IP>("permutation", num_cells, num_points);
  permutation.deep_copy(0);

  // Fill permutations with trivial values (i.e. no permutation - this will get overwritten for some cells)
  Kokkos::parallel_for(num_cells,KOKKOS_LAMBDA (const int & cell) {
    for(int point=0; point<num_points; ++point)
      permutation(cell,point) = point;
  });

  // Now we need to align the cubature points for each face
  // If there is only one point there is no need to align things
  if(num_points_per_face != 1) {
    // To enforce that quadrature points on faces are aligned properly we will iterate through faces,
    // map points to a plane shared by the faces, then re-order quadrature points on the "1" face to
    // line up with the "0" face

    // Utility calls
#define PANZER_DOT(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
#define PANZER_CROSS(a,b,c) {c[0] = a[1]*b[2] - a[2]*b[1]; c[1] = a[2]*b[0] - a[0]*b[2]; c[2] = a[0]*b[1] - a[1]*b[0];}

    // Iterate through faces
    Kokkos::parallel_for("face iteration",face_connectivity.numSubcells(),KOKKOS_LAMBDA (const int face) {
      // Cells for sides 0 and 1
      const int cell_0 = face_connectivity.cellForSubcell(face,0);
      const int cell_1 = face_connectivity.cellForSubcell(face,1);

      // If this face doesn't connect to anything we don't need to worry about alignment
      if(cell_1 < 0)
        return;

      // Local face index for sides 0 and 1
      const int lidx_0 = face_connectivity.localSubcellForSubcell(face,0);
      const int lidx_1 = face_connectivity.localSubcellForSubcell(face,1);

      // If the cell exists, then the local face index needs to exist
      KOKKOS_ASSERT(lidx_1 >= 0);

      // To compare points on planes, it is good to center the planes around the origin
      // Find the face center for the left and right cell (may not be the same point - e.g. periodic conditions)
      // We also define a length scale r2 which gives an order of magnitude approximation to the cell size squared
      Scalar xc0[3] = {0}, xc1[3] = {0}, r2 = 0;
      for(int face_point=0; face_point<num_points_per_face; ++face_point){
        Scalar dx2 = 0.;
        for(int dim=0; dim<cell_dim; ++dim){
          xc0[dim] += surface_points(cell_0,lidx_0*num_points_per_face + face_point,dim);
          xc1[dim] += surface_points(cell_1,lidx_1*num_points_per_face + face_point,dim);
          const Scalar dx = surface_points(cell_0,lidx_0*num_points_per_face + face_point,dim) - surface_points(cell_0,lidx_0*num_points_per_face,dim);
          dx2 += dx*dx;
        }

        // Check if the distance squared between these two points is larger than the others (doesn't need to be perfect)
        r2 = (r2 < dx2) ? dx2 : r2;
      }
      for(int dim=0; dim<cell_dim; ++dim){
        xc0[dim] /= double(num_points_per_face);
        xc1[dim] /= double(num_points_per_face);
      }

      // TODO: This needs to be adaptable to having curved faces - for now this will work

      // We need to define a plane with two vectors (transverse and binormal)
      // These will be used with the face centers to construct a local reference frame for aligning points

      // We use the first point on the face to define the normal (may break for curved faces)
      const int example_point_0 = lidx_0*num_points_per_face;
      const int example_point_1 = lidx_1*num_points_per_face;

      // Load the transverse and binormal for the 0 cell (default)
      Scalar t[3] = {surface_rotation_matrices(cell_0,example_point_0,1,0), surface_rotation_matrices(cell_0,example_point_0,1,1), surface_rotation_matrices(cell_0,example_point_0,1,2)};
      Scalar b[3] = {surface_rotation_matrices(cell_0,example_point_0,2,0), surface_rotation_matrices(cell_0,example_point_0,2,1), surface_rotation_matrices(cell_0,example_point_0,2,2)};

      // In case the faces are not antiparallel (e.g. periodic wedge), we may need to change the transverse and binormal
      {

        // Get the normals for each face for constructing one of the plane vectors (transverse)
        const Scalar n0[3] = {surface_rotation_matrices(cell_0,example_point_0,0,0), surface_rotation_matrices(cell_0,example_point_0,0,1), surface_rotation_matrices(cell_0,example_point_0,0,2)};
        const Scalar n1[3] = {surface_rotation_matrices(cell_1,example_point_1,0,0), surface_rotation_matrices(cell_1,example_point_1,0,1), surface_rotation_matrices(cell_1,example_point_1,0,2)};

        // n_0*n_1 == -1 (antiparallel), n_0*n_1 == 1 (parallel - bad), |n_0*n_1| < 1 (other)
        const Scalar n0_dot_n1 = PANZER_DOT(n0,n1);

        // FIXME: Currently virtual cells will set their surface normal along the same direction as the cell they "reflect"
        // This causes a host of issues (e.g. identifying 180 degree periodic wedges), but we have to support virtual cells as a priority
        // Therefore, we will just assume that the ordering is fine (not valid for 180 degree periodic wedges)
        if(Kokkos::fabs(n0_dot_n1 - 1.) < 1.e-8)
          return;

        // Rotated faces case (e.g. periodic wedge) we need to check for non-antiparallel face normals
        if(Kokkos::fabs(n0_dot_n1 + 1.) > 1.e-2){

          // Now we need to define an arbitrary transverse and binormal in the plane across which the faces are anti-symmetric
          // We can do this by setting t = n_0 \times n_1
          PANZER_CROSS(n0,n1,t);

          // Normalize the transverse vector
          const Scalar mag_t = Kokkos::sqrt(PANZER_DOT(t,t));
          t[0] /= mag_t;
          t[1] /= mag_t;
          t[2] /= mag_t;

          // The binormal will be the sum of the normals (does not need to be a right handed system)
          b[0] = n0[0] + n1[0];
          b[1] = n0[1] + n1[1];
          b[2] = n0[2] + n1[2];

          // Normalize the binormal vector
          const Scalar mag_b = Kokkos::sqrt(PANZER_DOT(b,b));
          b[0] /= mag_b;
          b[1] /= mag_b;
          b[2] /= mag_b;

        }
      }

      // Now that we have a reference coordinate plane in which to align our points
      // Point on the transverse/binormal plane
      Scalar p0[2] = {0};
      Scalar p1[2] = {0};

      // Differential position in mesh
      Scalar x0[3] = {0};
      Scalar x1[3] = {0};

      // Iterate through points in cell 1 and find which point they align with in cell 0
      for(int face_point_1=0; face_point_1<num_points_per_face; ++face_point_1){

        // Get the point index in the 1 cell
        const int point_1 = lidx_1*num_points_per_face + face_point_1;

        // Load position shifted by face center
        for(int dim=0; dim<cell_dim; ++dim)
          x1[dim] = surface_points(cell_1,point_1,dim) - xc1[dim];

        // Convert position to transverse/binormal plane
        p1[0] = PANZER_DOT(x1,t);
        p1[1] = PANZER_DOT(x1,b);

        // Compare to points on the other surface
        for(int face_point_0=0; face_point_0<num_points_per_face; ++face_point_0){

          // Get the point index in the 0 cell
          const int point_0 = lidx_0*num_points_per_face + face_point_0;

          // Load position shifted by face center
          for(int dim=0; dim<cell_dim; ++dim)
            x0[dim] = surface_points(cell_0,point_0,dim) - xc0[dim];

          // Convert position to transverse/binormal plane
          p0[0] = PANZER_DOT(x0,t);
          p0[1] = PANZER_DOT(x0,b);

          // Find the distance squared between p0 and p1
          const Scalar p012 = (p0[0]-p1[0])*(p0[0]-p1[0]) + (p0[1]-p1[1])*(p0[1]-p1[1]);

          // TODO: Should this be a constant value, or should we just find the minimum point?
          // If the distance, compared to the size of the cell, is small, we assume these are the same points
          if(p012 / r2 < 1.e-12){
            permutation(cell_1,lidx_1*num_points_per_face+face_point_0) = point_1;
            break;
          }

        }
      }
    });
    PHX::Device::execution_space().fence();
  }

#undef PANZER_DOT
#undef PANZER_CROSS

  return permutation;
}

} // end anonymous namespace

//template<typename DataType>
//using UnmanagedDynRankView = Kokkos::DynRankView<DataType,typename PHX::DevLayout<DataType>::type,PHX::Device,Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

template <typename Scalar>
IntegrationValues2<Scalar>::
IntegrationValues2(const std::string & pre,
                   const bool allocArrays):
  num_cells_(0),
  num_evaluate_cells_(0),
  num_virtual_cells_(-1),
  requires_permutation_(false),
  alloc_arrays_(allocArrays),
  prefix_(pre),
  ddims_(1,0)
{
  resetArrays();
}

template <typename Scalar>
void
IntegrationValues2<Scalar>::
resetArrays()
{
  cub_points_evaluated_ = false;
  side_cub_points_evaluated_ = false;
  cub_weights_evaluated_ = false;
  node_coordinates_evaluated_ = false;
  jac_evaluated_ = false;
  jac_inv_evaluated_ = false;
  jac_det_evaluated_ = false;
  weighted_measure_evaluated_ = false;
  weighted_normals_evaluated_ = false;
  surface_normals_evaluated_ = false;
  surface_rotation_matrices_evaluated_ = false;
  covarient_evaluated_ = false;
  contravarient_evaluated_ = false;
  norm_contravarient_evaluated_ = false;
  ip_coordinates_evaluated_ = false;
  ref_ip_coordinates_evaluated_ = false;

  // TODO: We need to clear the views
}

template <typename Scalar>
void IntegrationValues2<Scalar>::
setupArrays(const Teuchos::RCP<const panzer::IntegrationRule>& ir)
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::setupArrays()",setup_arrays);

  MDFieldArrayFactory af(prefix_,alloc_arrays_);

  typedef panzer::IntegrationDescriptor ID;

  int_rule = ir;

  const int num_nodes = ir->topology->getNodeCount();
  const int num_cells = ir->workset_size;
  const int num_space_dim = ir->topology->getDimension();

  // Specialize content if this is quadrature at a node
  const bool is_node_rule = (num_space_dim==1) and ir->isSide();
  if(not is_node_rule) {
    TEUCHOS_ASSERT(ir->getType() != ID::NONE);
    intrepid_cubature = getIntrepidCubature(*ir);
  }

  const int num_ip = is_node_rule ? 1 : ir->num_points;

  cub_points = af.template buildStaticArray<Scalar,IP,Dim>("cub_points",num_ip, num_space_dim);

  if (ir->isSide() && ir->cv_type == "none")
    side_cub_points = af.template buildStaticArray<Scalar,IP,Dim>("side_cub_points",num_ip,ir->side_topology->getDimension());

  cub_weights = af.template buildStaticArray<Scalar,IP>("cub_weights",num_ip);

  node_coordinates = af.template buildStaticArray<Scalar,Cell,BASIS,Dim>("node_coordinates",num_cells, num_nodes, num_space_dim);

  jac = af.template buildStaticArray<Scalar,Cell,IP,Dim,Dim>("jac",num_cells, num_ip, num_space_dim,num_space_dim);

  jac_inv = af.template buildStaticArray<Scalar,Cell,IP,Dim,Dim>("jac_inv",num_cells, num_ip, num_space_dim,num_space_dim);

  jac_det = af.template buildStaticArray<Scalar,Cell,IP>("jac_det",num_cells, num_ip);

  weighted_measure =  af.template buildStaticArray<Scalar,Cell,IP>("weighted_measure",num_cells, num_ip);

  covarient = af.template buildStaticArray<Scalar,Cell,IP,Dim,Dim>("covarient",num_cells, num_ip, num_space_dim,num_space_dim);

  contravarient = af.template buildStaticArray<Scalar,Cell,IP,Dim,Dim>("contravarient",num_cells, num_ip, num_space_dim,num_space_dim);

  norm_contravarient = af.template buildStaticArray<Scalar,Cell,IP>("norm_contravarient",num_cells, num_ip);

  ip_coordinates = af.template buildStaticArray<Scalar,Cell,IP,Dim>("ip_coordiantes",num_cells, num_ip,num_space_dim);

  ref_ip_coordinates = af.template buildStaticArray<Scalar,Cell,IP,Dim>("ref_ip_coordinates",num_cells, num_ip,num_space_dim);

  weighted_normals = af.template buildStaticArray<Scalar,Cell,IP,Dim>("weighted_normal",num_cells,num_ip,num_space_dim);

  surface_normals = af.template buildStaticArray<Scalar,Cell,IP,Dim>("surface_normals",num_cells, num_ip,num_space_dim);

  surface_rotation_matrices = af.template buildStaticArray<Scalar,Cell,IP,Dim,Dim>("surface_rotation_matrices",num_cells, num_ip,3,3);
}


// ***********************************************************
// * Evaluation of values - NOT specialized
// ***********************************************************
template <typename Scalar>
void IntegrationValues2<Scalar>::
evaluateValues(const PHX::MDField<Scalar,Cell,NODE,Dim> & in_node_coordinates,
               const int in_num_cells,
               const Teuchos::RCP<const SubcellConnectivity> & face_connectivity,
               const int num_virtual_cells)
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::evaluateValues(with virtual cells)",eval_vals_with_virts);

  setup(int_rule, in_node_coordinates, in_num_cells);

  // Setup permutations (only required for surface integrators)
  using ID=panzer::IntegrationDescriptor;
  if((int_rule->getType() == ID::SURFACE) and (face_connectivity != Teuchos::null))
    setupPermutations(face_connectivity, num_virtual_cells);

  // Evaluate everything once permutations are generated
  evaluateEverything();
}

template <typename Scalar>
void
IntegrationValues2<Scalar>::
evaluateValues(const PHX::MDField<Scalar,Cell,NODE,Dim>& in_node_coordinates,
               const PHX::MDField<Scalar,Cell,IP,Dim>& other_ip_coordinates,
               const int in_num_cells)
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::evaluateValues(no virtual cells)",eval_vals_no_virts);

  setup(int_rule, in_node_coordinates, in_num_cells);

  // Setup permutations
  setupPermutations(other_ip_coordinates);

  // Evaluate everything once permutations are generated
  evaluateEverything();
}

template <typename Scalar>
void
IntegrationValues2<Scalar>::
setupPermutations(const Teuchos::RCP<const SubcellConnectivity> & face_connectivity,
                  const int num_virtual_cells)
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::setupPermutations(connectivity)",setup_perms_conn);

  TEUCHOS_ASSERT(not int_rule->isSide());
  TEUCHOS_ASSERT(face_connectivity != Teuchos::null);
  TEUCHOS_TEST_FOR_EXCEPT_MSG(int_rule->getType() != panzer::IntegrationDescriptor::SURFACE,
                              "IntegrationValues2::setupPermutations : Face connectivity based permutations are only required for surface integration schemes");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(num_virtual_cells_ >= 0,
                              "IntegrationValues2::setupPermutations : Number of virtual cells for surface permuations must be >= 0");
  resetArrays();
  side_connectivity_ = face_connectivity;
  num_virtual_cells_ = num_virtual_cells;
  requires_permutation_ = false;
  permutations_ = generateSurfacePermutations<Scalar>(num_evaluate_cells_,*face_connectivity, getCubaturePoints(false,true), getSurfaceRotationMatrices(false,true));
  requires_permutation_ = true;
}

template <typename Scalar>
void
IntegrationValues2<Scalar>::
setupPermutations(const PHX::MDField<Scalar,Cell,IP,Dim> & other_ip_coordinates)
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::setupPermutations(other_coords)",setup_perms_coords);
  resetArrays();
  requires_permutation_ = false;
  permutations_ = generatePermutations<Scalar>(num_evaluate_cells_, getCubaturePoints(false,true), other_ip_coordinates);
  requires_permutation_ = true;
}


template <typename Scalar>
void
IntegrationValues2<Scalar>::
setup(const Teuchos::RCP<const panzer::IntegrationRule>& ir,
      const PHX::MDField<Scalar,Cell,NODE,Dim> & cell_node_coordinates,
      const int num_cells)
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::setup()",setup);

  // Clear arrays just in case we are rebuilding this object
  resetArrays();

  num_cells_ = cell_node_coordinates.extent(0);
  num_evaluate_cells_ = num_cells < 0 ? cell_node_coordinates.extent(0) : num_cells;
  int_rule = ir;

  TEUCHOS_ASSERT(ir->getType() != IntegrationDescriptor::NONE);
  intrepid_cubature = getIntrepidCubature(*ir);

  // Copy node coordinates into owned allocation
  {
    MDFieldArrayFactory af(prefix_,true);

    const int num_space_dim = int_rule->topology->getDimension();
    const int num_nodes = int_rule->topology->getNodeCount();
    TEUCHOS_ASSERT(static_cast<int>(cell_node_coordinates.extent(1)) == num_nodes);

    auto aux = af.template buildStaticArray<Scalar,Cell,NODE,Dim>("node_coordinates",num_cells_,num_nodes, num_space_dim);
    Kokkos::MDRangePolicy<PHX::Device,Kokkos::Rank<3>> policy({0,0,0},{num_evaluate_cells_,num_nodes,num_space_dim});
    Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int & cell, const int & point, const int & dim){
      aux(cell,point,dim) = cell_node_coordinates(cell,point,dim);
    });
    PHX::Device::execution_space().fence();
    node_coordinates = aux;
  }

}

template <typename Scalar>
typename IntegrationValues2<Scalar>::ConstArray_IPDim
IntegrationValues2<Scalar>::
getUniformCubaturePointsRef(const bool cache,
                            const bool force,
                            const bool apply_permutation) const
{
  if(cub_points_evaluated_ and (apply_permutation == requires_permutation_) and not force)
    return cub_points;

  // Only log time if values computed (i.e. don't log if values are already cached)
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::getUniformCubaturePointsRef()",get_uniform_cub_pts_ref);

  Intrepid2::CellTools<PHX::Device::execution_space> cell_tools;
  MDFieldArrayFactory af(prefix_,true);

  int num_space_dim = int_rule->topology->getDimension();
  int num_ip = int_rule->num_points;

  auto aux = af.template buildStaticArray<Scalar,IP,Dim>("cub_points",num_ip, num_space_dim);

  if (int_rule->isSide() && num_space_dim==1) {
    std::cout << "WARNING: 0-D quadrature rule infrastructure does not exist!!! Will not be able to do "
        << "non-natural integration rules.";
    return aux;
  }

  TEUCHOS_TEST_FOR_EXCEPT_MSG(int_rule->cv_type != "none",
                              "IntegrationValues2::getUniformCubaturePointsRef : Cannot build reference cubature points for control volume integration scheme.");

  TEUCHOS_TEST_FOR_EXCEPT_MSG(int_rule->getType() == IntegrationDescriptor::SURFACE,
                              "IntegrationValues2::getUniformCubaturePointsRef : Cannot build reference cubature points for surface integration scheme.");

  auto weights = af.template buildStaticArray<Scalar,IP>("cub_weights",num_ip);

  if (!int_rule->isSide())
    intrepid_cubature->getCubature(aux.get_view(), weights.get_view());
  else {
    auto s_cub_points = af.template buildStaticArray<Scalar,IP,Dim>("side_cub_points",num_ip, num_space_dim-1);

    intrepid_cubature->getCubature(s_cub_points.get_view(), weights.get_view());
    cell_tools.mapToReferenceSubcell(aux.get_view(), s_cub_points.get_view(), num_space_dim-1, int_rule->getSide(), *(int_rule->topology));
  }

  PHX::Device::execution_space().fence();

  if(apply_permutation and requires_permutation_)
    applyBasePermutation(aux, permutations_);

  if(cache and (apply_permutation == requires_permutation_)){
    cub_points = aux;
    cub_points_evaluated_ = true;
  }

  return aux;

}

template <typename Scalar>
typename IntegrationValues2<Scalar>::ConstArray_IPDim
IntegrationValues2<Scalar>::
getUniformSideCubaturePointsRef(const bool cache,
                                const bool force,
                                const bool apply_permutation) const
{
  if(side_cub_points_evaluated_ and (apply_permutation == requires_permutation_) and not force)
    return side_cub_points;

  // Only log time if values computed (i.e. don't log if values are already cached)
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::getUniformSideCubaturePointsRef()",get_uniform_side_cub_pts_ref);

  MDFieldArrayFactory af(prefix_,true);

  int num_space_dim = int_rule->topology->getDimension();
  int num_ip = int_rule->num_points;

  auto aux = af.template buildStaticArray<Scalar,IP,Dim>("side_cub_points",num_ip, num_space_dim-1);

  if (int_rule->isSide() && num_space_dim==1) {
    std::cout << "WARNING: 0-D quadrature rule infrastructure does not exist!!! Will not be able to do "
        << "non-natural integration rules.";
    return aux;
  }

  TEUCHOS_TEST_FOR_EXCEPT_MSG(int_rule->cv_type != "none",
                              "IntegrationValues2::getUniformSideCubaturePointsRef : Cannot build reference cubature points for control volume integration scheme.");

  TEUCHOS_TEST_FOR_EXCEPT_MSG(int_rule->getType() == IntegrationDescriptor::SURFACE,
                              "IntegrationValues2::getUniformSideCubaturePointsRef : Cannot build reference cubature points for surface integration scheme.");

  TEUCHOS_TEST_FOR_EXCEPT_MSG(!int_rule->isSide(),
                              "IntegrationValues2::getUniformSideCubaturePointsRef : Requested side points, which is not supported by integration rule.");

  auto weights = af.template buildStaticArray<Scalar,IP>("cub_weights",num_ip);

  intrepid_cubature->getCubature(aux.get_view(), weights.get_view());

  PHX::Device::execution_space().fence();

  if(apply_permutation and requires_permutation_)
    applyBasePermutation(aux, permutations_);

  if(cache and (apply_permutation == requires_permutation_)){
    side_cub_points = aux;
    side_cub_points_evaluated_ = true;
  }

  return aux;
}

template <typename Scalar>
typename IntegrationValues2<Scalar>::ConstArray_IP
IntegrationValues2<Scalar>::
getUniformCubatureWeightsRef(const bool cache,
                             const bool force,
                             const bool apply_permutation) const
{
  if(cub_weights_evaluated_ and (apply_permutation == requires_permutation_) and not force)
    return cub_weights;

  // Only log time if values computed (i.e. don't log if values are already cached)
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::getUniformCubatureWeightRef()",get_uniform_cub_weights_ref);

  MDFieldArrayFactory af(prefix_,true);

  int num_space_dim = int_rule->topology->getDimension();
  int num_ip = int_rule->num_points;

  auto aux = af.template buildStaticArray<Scalar,IP>("cub_weights",num_ip);

  if (int_rule->isSide() && num_space_dim==1) {
    std::cout << "WARNING: 0-D quadrature rule infrastructure does not exist!!! Will not be able to do "
        << "non-natural integration rules.";
    return aux;
  }

  TEUCHOS_TEST_FOR_EXCEPT_MSG(int_rule->cv_type != "none",
                              "IntegrationValues2::getUniformCubatureWeightsRef : Cannot build reference cubature weights for control volume integration scheme.");

  TEUCHOS_TEST_FOR_EXCEPT_MSG(int_rule->getType() == IntegrationDescriptor::SURFACE,
                              "IntegrationValues2::getUniformCubatureWeightsRef : Cannot build reference cubature weights for surface integration scheme.");

  auto points = af.template buildStaticArray<Scalar,IP,Dim>("cub_points",num_ip,num_space_dim);

  intrepid_cubature->getCubature(points.get_view(), aux.get_view());

  PHX::Device::execution_space().fence();

  if(apply_permutation and requires_permutation_)
    applyBasePermutation(aux, permutations_);

  if(cache and (apply_permutation == requires_permutation_)){
    cub_weights = aux;
    cub_weights_evaluated_ = true;
  }

  return aux;

}

template <typename Scalar>
typename IntegrationValues2<Scalar>::ConstArray_CellBASISDim
IntegrationValues2<Scalar>::
getNodeCoordinates() const
{
  return node_coordinates;
}

template <typename Scalar>
typename IntegrationValues2<Scalar>::ConstArray_CellIPDimDim
IntegrationValues2<Scalar>::
getJacobian(const bool cache,
            const bool force) const
{
  if(jac_evaluated_ and not force)
    return jac;

  // Only log time if values computed (i.e. don't log if values are already cached)
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::getJacobian()",get_jacobian);

  Intrepid2::CellTools<PHX::Device::execution_space> cell_tools;
  MDFieldArrayFactory af(prefix_,true);

  int num_space_dim = int_rule->topology->getDimension();
  int num_ip = int_rule->num_points;

  using ID=panzer::IntegrationDescriptor;
  const bool is_cv = (int_rule->getType() == ID::CV_VOLUME) or (int_rule->getType() == ID::CV_SIDE) or (int_rule->getType() == ID::CV_BOUNDARY);
  const bool is_surface = int_rule->getType() == ID::SURFACE;

  auto aux = af.template buildStaticArray<Scalar,Cell,IP,Dim,Dim>("jac",num_cells_, num_ip, num_space_dim,num_space_dim);

  if(is_cv or is_surface){

    // Don't forget that since we are not caching this, we have to make sure the managed view remains alive while we use the non-const wrapper
    auto const_ref_coord = getCubaturePointsRef(false,force);
    auto ref_coord = PHX::getNonConstDynRankViewFromConstMDField(const_ref_coord);
    auto node_coord = PHX::getNonConstDynRankViewFromConstMDField(getNodeCoordinates());
    
    const auto cell_range = std::make_pair(0,num_evaluate_cells_);
    auto s_ref_coord  = Kokkos::subview(ref_coord,     cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_node_coord = Kokkos::subview(node_coord,    cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_jac        = Kokkos::subview(aux.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());

    cell_tools.setJacobian(s_jac, s_ref_coord, s_node_coord,*(int_rule->topology));

  } else {

    // Don't forget that since we are not caching this, we have to make sure the managed view remains alive while we use the non-const wrapper
    auto const_ref_coord = getUniformCubaturePointsRef(false,force,false);
    auto ref_coord = PHX::getNonConstDynRankViewFromConstMDField(const_ref_coord);
    auto node_coord = PHX::getNonConstDynRankViewFromConstMDField(getNodeCoordinates());
    
    const auto cell_range = std::make_pair(0,num_evaluate_cells_);
    auto s_node_coord = Kokkos::subview(node_coord,    cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_jac        = Kokkos::subview(aux.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());

    cell_tools.setJacobian(s_jac, ref_coord, s_node_coord,*(int_rule->topology));

    if(requires_permutation_)
      applyPermutation(aux, permutations_);

  }

  PHX::Device::execution_space().fence();

  if(cache){
    jac = aux;
    jac_evaluated_ = true;
  }

  return aux;
}

template <typename Scalar>
typename IntegrationValues2<Scalar>::ConstArray_CellIPDimDim
IntegrationValues2<Scalar>::
getJacobianInverse(const bool cache,
                   const bool force) const
{
  if(jac_inv_evaluated_ and not force)
    return jac_inv;

  // Only log time if values computed (i.e. don't log if values are already cached)
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::getJacobianInverse()",get_jacobian_inv);

  Intrepid2::CellTools<PHX::Device::execution_space> cell_tools;
  MDFieldArrayFactory af(prefix_,true);

  const int num_space_dim = int_rule->topology->getDimension();
  const int num_ip = int_rule->num_points;

  auto jacobian = getJacobian(false,force);
  auto aux = af.template buildStaticArray<Scalar,Cell,IP,Dim,Dim>("jac_inv",num_cells_, num_ip, num_space_dim,num_space_dim);

  const auto cell_range = std::make_pair(0,num_evaluate_cells_);
  auto s_jac      = Kokkos::subview(jacobian.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
  auto s_jac_inv  = Kokkos::subview(aux.get_view(),     cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());

  cell_tools.setJacobianInv(s_jac_inv, s_jac);

  PHX::Device::execution_space().fence();

  if(cache){
    jac_inv = aux;
    jac_inv_evaluated_ = true;
  }

  return aux;
}

template <typename Scalar>
typename IntegrationValues2<Scalar>::ConstArray_CellIP
IntegrationValues2<Scalar>::
getJacobianDeterminant(const bool cache,
                       const bool force) const
{
  if(jac_det_evaluated_ and not force)
    return jac_det;

  // Only log time if values computed (i.e. don't log if values are already cached)
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::getJacobianDeterminant()",get_jacobian_det);

  Intrepid2::CellTools<PHX::Device::execution_space> cell_tools;
  MDFieldArrayFactory af(prefix_,true);

  const int num_ip = int_rule->num_points;

  auto jacobian = getJacobian(false,force);
  auto aux = af.template buildStaticArray<Scalar,Cell,IP>("jac_det",num_cells_, num_ip);

  const auto cell_range = std::make_pair(0,num_evaluate_cells_);
  auto s_jac     = Kokkos::subview(jacobian.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
  auto s_jac_det = Kokkos::subview(aux.get_view(),     cell_range,Kokkos::ALL());

  cell_tools.setJacobianDet(s_jac_det, s_jac);

  PHX::Device::execution_space().fence();

  if(cache){
    jac_det = aux;
    jac_det_evaluated_ = true;
  }

  return aux;
}

template <typename Scalar>
typename IntegrationValues2<Scalar>::ConstArray_CellIP
IntegrationValues2<Scalar>::
getWeightedMeasure(const bool cache,
                   const bool force) const
{
  if(weighted_measure_evaluated_ and not force)
    return weighted_measure;

  // Only log time if values computed (i.e. don't log if values are already cached)
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::getWeightedMeasure()",get_wt_meas);

  MDFieldArrayFactory af(prefix_,true);

  const int num_space_dim = int_rule->topology->getDimension();
  const int num_ip = int_rule->num_points;

  auto aux = af.template buildStaticArray<Scalar,Cell,IP>("weighted_measure",num_cells_, num_ip);

  if(int_rule->cv_type != "none"){

    TEUCHOS_TEST_FOR_EXCEPT_MSG(int_rule->cv_type == "side",
                                "IntegrationValues2::getWeightedMeasure : Weighted measure not available for side control volume integrators. Use getWeightedNormals instead.");

    // CV integration uses a single call to map from physical space to the weighted measure - I assume this is slower than what we do with non-cv integration methods

    auto s_cub_points = af.template buildStaticArray<Scalar, Cell, IP, Dim>("cub_points",num_evaluate_cells_,num_ip,num_space_dim);

    auto node_coord = PHX::getNonConstDynRankViewFromConstMDField(getNodeCoordinates());

    const auto cell_range = std::make_pair(0,num_evaluate_cells_);
    auto s_node_coord =       Kokkos::subview(node_coord,    cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_weighted_measure = Kokkos::subview(aux.get_view(),cell_range,Kokkos::ALL());

    intrepid_cubature->getCubature(s_cub_points.get_view(),s_weighted_measure,s_node_coord);

  } else if(int_rule->getType() == IntegrationDescriptor::SURFACE){

    const auto & cell_topology = *int_rule->topology;
    const int cell_dim = cell_topology.getDimension();
    const int num_sides = (cell_dim==1) ? 2 : cell_topology.getSideCount();
    const int cubature_order = int_rule->order();
    const int num_points_on_side = num_ip / num_sides;

    Intrepid2::DefaultCubatureFactory cubature_factory;
    auto jacobian = getJacobian(false,force);

    for(int side=0; side<num_sides; ++side) {

      const int point_offset=side*num_points_on_side;

      // Get the cubature for the side
      Kokkos::DynRankView<double,PHX::Device> side_cub_weights;
      if(cell_dim==1){
        side_cub_weights = Kokkos::DynRankView<double,PHX::Device>("side_cub_weights",num_points_on_side);
        auto tmp_side_cub_weights_host = Kokkos::create_mirror_view(side_cub_weights);
        tmp_side_cub_weights_host(0)=1.;
        Kokkos::deep_copy(side_cub_weights,tmp_side_cub_weights_host);
      } else {

        // Get the face topology from the cell topology
        const shards::CellTopology face_topology(cell_topology.getCellTopologyData(cell_dim-1,side));

        auto ic = cubature_factory.create<PHX::Device::execution_space,double,double>(face_topology,cubature_order);

        side_cub_weights = Kokkos::DynRankView<double,PHX::Device>("side_cub_weights",num_points_on_side);
        auto subcell_cub_points = Kokkos::DynRankView<double,PHX::Device>("subcell_cub_points",num_points_on_side,cell_dim-1);

        // Get the reference face points
        ic->getCubature(subcell_cub_points, side_cub_weights);
      }

      PHX::Device::execution_space().fence();

      // Iterating over face points
      Kokkos::MDRangePolicy<PHX::Device::execution_space,Kokkos::Rank<2>> policy({0,0},{num_evaluate_cells_,num_points_on_side});

      // Calculate measures (quadrature weights in physical space) for this side
      auto side_weighted_measure = Kokkos::DynRankView<Scalar,PHX::Device>("side_weighted_measure",num_evaluate_cells_,num_points_on_side);
      if(cell_dim == 1){
        auto side_cub_weights_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),side_cub_weights);
        Kokkos::deep_copy(side_weighted_measure, side_cub_weights_host(0));
      } else {

        // Copy from complete jacobian to side jacobian
        auto side_jacobian = Kokkos::DynRankView<Scalar,PHX::Device>("side_jac",num_evaluate_cells_,num_points_on_side,cell_dim,cell_dim);

        Kokkos::parallel_for("Copy jacobian to side jacobian",policy,KOKKOS_LAMBDA (const int cell,const int point) {
          for(int dim=0;dim<cell_dim;++dim)
            for(int dim1=0;dim1<cell_dim;++dim1)
              side_jacobian(cell,point,dim,dim1) = jacobian(cell,point_offset+point,dim,dim1);
        });

        auto scratch = af.template buildStaticArray<Scalar,Point>("scratch_for_compute_measure", num_evaluate_cells_*num_points_on_side*num_space_dim*num_space_dim);

        if(cell_dim == 2){
          Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
            computeEdgeMeasure(side_weighted_measure, side_jacobian, side_cub_weights,
                               side,cell_topology,
                               scratch.get_view());
          PHX::Device::execution_space().fence();
        } else if(cell_dim == 3){
          Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
            computeFaceMeasure(side_weighted_measure, side_jacobian, side_cub_weights,
                               side,cell_topology,
                               scratch.get_view());
          PHX::Device::execution_space().fence();
        }
      }


      // Copy to the main array
      Kokkos::parallel_for("copy surface weighted measure values",policy,KOKKOS_LAMBDA (const int cell,const int point) {
        aux(cell,point_offset + point) = side_weighted_measure(cell,point);
      });
      PHX::Device::execution_space().fence();
    }

  } else {

    auto cell_range = std::make_pair(0,num_evaluate_cells_);
    auto s_weighted_measure = Kokkos::subview(aux.get_view(),cell_range,Kokkos::ALL());
    auto cubature_weights = getUniformCubatureWeightsRef(false,force,false);

    if (!int_rule->isSide()) {

      auto s_jac_det = Kokkos::subview(getJacobianDeterminant(false,force).get_view(),cell_range,Kokkos::ALL());
      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
      computeCellMeasure(s_weighted_measure, s_jac_det, cubature_weights.get_view());

    } else if(int_rule->spatial_dimension==3) {

      auto s_jac = Kokkos::subview(getJacobian(false,force).get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
      auto scratch = af.template buildStaticArray<Scalar,Point>("scratch_for_compute_measure", num_evaluate_cells_*num_ip*num_space_dim*num_space_dim);
      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
      computeFaceMeasure(s_weighted_measure, s_jac, cubature_weights.get_view(),
                         int_rule->side, *int_rule->topology,
                         scratch.get_view());

    } else if(int_rule->spatial_dimension==2) {

      auto s_jac = Kokkos::subview(getJacobian(false,force).get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
      auto scratch = af.template buildStaticArray<Scalar,Point>("scratch_for_compute_measure", num_evaluate_cells_*num_ip*num_space_dim*num_space_dim);
      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
      computeEdgeMeasure(s_weighted_measure, s_jac, cubature_weights.get_view(),
                         int_rule->side,*int_rule->topology,
                         scratch.get_view());

    } else {
      TEUCHOS_ASSERT(false);
    }

  }

  PHX::Device::execution_space().fence();

  // Apply permutations if necessary
  if(requires_permutation_)
    applyPermutation(aux, permutations_);

  if(cache){
    weighted_measure = aux;
    weighted_measure_evaluated_ = true;
  }

  return aux;
}

template <typename Scalar>
typename IntegrationValues2<Scalar>::ConstArray_CellIPDim
IntegrationValues2<Scalar>::
getWeightedNormals(const bool cache,
                   const bool force) const
{
  if(weighted_normals_evaluated_ and not force)
    return weighted_normals;

  // Only log time if values computed (i.e. don't log if values are already cached)
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::getWeightedNormals()",get_wt_normals);

  MDFieldArrayFactory af(prefix_,true);

  const int num_space_dim = int_rule->topology->getDimension();
  const int num_ip = int_rule->num_points;

  auto aux = af.template buildStaticArray<Scalar,Cell,IP,Dim>("weighted_normals",num_cells_,num_ip,num_space_dim);

  TEUCHOS_TEST_FOR_EXCEPT_MSG(int_rule->getType() == IntegrationDescriptor::SURFACE,
                              "IntegrationValues2::getWeightedNormals : Cannot build reference weighted normals for surface integration scheme.");

  auto points = af.template buildStaticArray<Scalar,Cell,IP,Dim>("cub_points",num_cells_,num_ip,num_space_dim);

  auto node_coord = PHX::getNonConstDynRankViewFromConstMDField(getNodeCoordinates());

  const auto cell_range = std::make_pair(0,num_evaluate_cells_);
  auto s_cub_points       = Kokkos::subview(points.get_view(),cell_range, Kokkos::ALL(), Kokkos::ALL());
  auto s_weighted_normals = Kokkos::subview(aux.get_view(),   cell_range, Kokkos::ALL(), Kokkos::ALL());
  auto s_node_coord       = Kokkos::subview(node_coord,       cell_range, Kokkos::ALL(), Kokkos::ALL());

  intrepid_cubature->getCubature(s_cub_points,s_weighted_normals,s_node_coord);

  PHX::Device::execution_space().fence();

  // Apply permutations if necessary
  if(requires_permutation_)
    applyPermutation(aux, permutations_);

  if(cache){
    weighted_normals = aux;
    weighted_normals_evaluated_ = true;
  }

  return aux;
}

template <typename Scalar>
typename IntegrationValues2<Scalar>::ConstArray_CellIPDim
IntegrationValues2<Scalar>::
getSurfaceNormals(const bool cache,
                  const bool force) const
{
  if(surface_normals_evaluated_ and not force)
    return surface_normals;

  // Only log time if values computed (i.e. don't log if values are already cached)
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::getSurfaceNormals()",get_surf_normals);

  TEUCHOS_TEST_FOR_EXCEPT_MSG(int_rule->isSide(),
                              "IntegrationValues2::getSurfaceNormals : This call doesn't work with sides (only surfaces).");

  TEUCHOS_TEST_FOR_EXCEPT_MSG(int_rule->cv_type != "none",
                              "IntegrationValues2::getSurfaceNormals : This call does not support control volume integration schemes.");

  TEUCHOS_TEST_FOR_EXCEPT_MSG(int_rule->getType() != IntegrationDescriptor::SURFACE,
                              "IntegrationValues2::getSurfaceNormals : Can only build for surface integrators.");

  Intrepid2::CellTools<PHX::Device::execution_space> cell_tools;
  MDFieldArrayFactory af(prefix_,true);

  const shards::CellTopology & cell_topology = *(int_rule->topology);
  const int cell_dim = cell_topology.getDimension();
  const int subcell_dim = cell_topology.getDimension()-1;
  const int num_subcells = cell_topology.getSubcellCount(subcell_dim);
  const int num_space_dim = int_rule->topology->getDimension();
  const int num_ip = int_rule->num_points;
  const int num_points_on_face = num_ip / num_subcells;

  auto aux = af.template buildStaticArray<Scalar,Cell,IP,Dim>("surface_normals",num_cells_,num_ip,num_space_dim);

  // We only need the jacobian if we're not in 1D
  ConstArray_CellIPDimDim jacobian;
  if(cell_dim != 1)
    jacobian = getJacobian(false,force);

  // We get to build up our surface normals one 'side' at a time
  for(int subcell_index=0; subcell_index<num_subcells; ++subcell_index) {

    const int point_offset = subcell_index * num_points_on_face;;

    // Calculate normals
    auto side_normals = Kokkos::DynRankView<Scalar,PHX::Device>("side_normals",num_evaluate_cells_,num_points_on_face,cell_dim);
    if(cell_dim == 1){

      const int other_subcell_index = (subcell_index==0) ? 1 : 0;
      auto in_node_coordinates_k = getNodeCoordinates().get_view();
      const auto min = std::numeric_limits<typename Sacado::ScalarType<Scalar>::type>::min();

      Kokkos::parallel_for("compute normals 1D",num_evaluate_cells_,KOKKOS_LAMBDA (const int cell) {
        Scalar norm = (in_node_coordinates_k(cell,subcell_index,0) - in_node_coordinates_k(cell,other_subcell_index,0));
        side_normals(cell,0,0) = norm / fabs(norm + min);
      });

    } else {

      // Iterating over side points
      Kokkos::MDRangePolicy<PHX::Device::execution_space,Kokkos::Rank<2>> policy({0,0},{num_evaluate_cells_,num_points_on_face});

      auto side_jacobian = Kokkos::DynRankView<Scalar,PHX::Device>("side_jac",num_evaluate_cells_,num_points_on_face,cell_dim,cell_dim);
      Kokkos::parallel_for("Copy jacobian to side jacobian",policy,KOKKOS_LAMBDA (const int cell,const int point) {
        for(int dim=0;dim<cell_dim;++dim)
          for(int dim1=0;dim1<cell_dim;++dim1)
            side_jacobian(cell,point,dim,dim1) = jacobian(cell,point_offset+point,dim,dim1);
      });

      // Get the "physical side normals"
      cell_tools.getPhysicalSideNormals(side_normals,side_jacobian,subcell_index,cell_topology);

      // Normalize each normal
      Kokkos::parallel_for("Normalize the normals",policy,KOKKOS_LAMBDA (const int cell,const int point) {
        Scalar n = 0.;
        for(int dim=0;dim<cell_dim;++dim){
          n += side_normals(cell,point,dim)*side_normals(cell,point,dim);
        }
        // If n is zero then this is - hopefully - a virtual cell
        if(n > 0.){
          n = Kokkos::sqrt(n);
          for(int dim=0;dim<cell_dim;++dim)
            side_normals(cell,point,dim) /= n;
        }
      });
    }

    PHX::Device::execution_space().fence();

    Kokkos::MDRangePolicy<PHX::Device::execution_space,Kokkos::Rank<2>> policy({0,0},{num_evaluate_cells_,num_points_on_face});
    Kokkos::parallel_for("copy surface normals", policy,KOKKOS_LAMBDA (const int cell,const int point) {
      for(int dim=0;dim<cell_dim;++dim)
        aux(cell,point_offset + point,dim) = side_normals(cell,point,dim);
    });
    PHX::Device::execution_space().fence();
  }

  // Need to correct the virtual cells
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(side_connectivity_ == Teuchos::null,
                                "IntegrationValues2::getSurfaceNormals : Surface normals require a SubcellConnectivity object pass in through the setupPermutations call");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(num_virtual_cells_ < 0,
                                  "IntegrationValues2::getSurfaceNormals : Number of virtual cells (see setupPermutations or evaluateValues) must be greater or equal to zero");

    // Virtual cell normals need to be reversed
    if(num_virtual_cells_ > 0)
      correctVirtualNormals(aux, num_evaluate_cells_ - num_virtual_cells_, num_virtual_cells_, *int_rule->topology, *side_connectivity_);
  }

  if(cache){
    surface_normals = aux;
    surface_normals_evaluated_ = true;
  }

  return aux;

}

template <typename Scalar>
typename IntegrationValues2<Scalar>::ConstArray_CellIPDimDim
IntegrationValues2<Scalar>::
getSurfaceRotationMatrices(const bool cache,
                           const bool force) const
{
  if(surface_rotation_matrices_evaluated_ and not force)
    return surface_rotation_matrices;

  // Only log time if values computed (i.e. don't log if values are already cached)
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::getSurfaceRotationMatrices()",get_surf_rot_mat);

  MDFieldArrayFactory af(prefix_,true);

  const int num_ip = int_rule->num_points;
  const int cell_dim = int_rule->topology->getDimension();

  // This call will handle all the error checking
  auto normals = getSurfaceNormals(false,force).get_static_view();
  auto aux = af.template buildStaticArray<Scalar,Cell,IP,Dim,Dim>("surface_rotation_matrices",num_cells_, num_ip, 3, 3);

  Kokkos::MDRangePolicy<PHX::Device::execution_space,Kokkos::Rank<2>> policy({0,0},{num_evaluate_cells_,num_ip});
  Kokkos::parallel_for("create surface rotation matrices",policy,KOKKOS_LAMBDA (const int cell,const int point) {
    Scalar normal[3];
    for(int i=0;i<3;i++)
      normal[i]=0.;
    for(int dim=0; dim<cell_dim; ++dim)
      normal[dim] = normals(cell,point,dim);

    Scalar transverse[3];
    Scalar binormal[3];
    panzer::convertNormalToRotationMatrix(normal,transverse,binormal);

    for(int dim=0; dim<3; ++dim){
      aux(cell,point,0,dim) = normal[dim];
      aux(cell,point,1,dim) = transverse[dim];
      aux(cell,point,2,dim) = binormal[dim];
    }
  });
  PHX::Device::execution_space().fence();

  // Need to correct the virtual cells
  {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(side_connectivity_ == Teuchos::null,
                                "IntegrationValues2::getSurfaceRotationMatrices : Surface normals require a SubcellConnectivity object pass in through the setupPermutations call");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(num_virtual_cells_ < 0,
                                  "IntegrationValues2::getSurfaceRotationMatrices : Number of virtual cells (see setupPermutations or evaluateValues) must be greater or equal to zero");

    // Virtual cell normals need to be reversed
    if(num_virtual_cells_ > 0)
      correctVirtualRotationMatrices(aux, num_evaluate_cells_ - num_virtual_cells_, num_virtual_cells_, *int_rule->topology, *side_connectivity_);
  }

  if(cache){
    surface_rotation_matrices = aux;
    surface_rotation_matrices_evaluated_ = true;
  }

  return aux;
}

template <typename Scalar>
typename IntegrationValues2<Scalar>::ConstArray_CellIPDimDim
IntegrationValues2<Scalar>::
getCovarientMatrix(const bool cache,
                   const bool force) const
{
  if(covarient_evaluated_ and not force)
    return covarient;

  // Only log time if values computed (i.e. don't log if values are already cached)
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::getCovariantMatrix()",get_cov_mat);

  MDFieldArrayFactory af(prefix_,true);

  const int num_space_dim = int_rule->topology->getDimension();
  const int num_ip = int_rule->num_points;

  auto jacobian = getJacobian(false,force).get_static_view();
  auto aux = af.template buildStaticArray<Scalar,Cell,IP,Dim,Dim>("covarient",num_cells_, num_ip, num_space_dim,num_space_dim);

  Kokkos::MDRangePolicy<PHX::Device,Kokkos::Rank<2>> policy({0,0},{num_evaluate_cells_,num_ip});
  Kokkos::parallel_for("evalaute covarient metric tensor",policy,KOKKOS_LAMBDA (const int cell,const int ip) {
    // g^{ij} = \frac{\parital x_i}{\partial \chi_\alpha}\frac{\parital x_j}{\partial \chi_\alpha}
    for (int i = 0; i < num_space_dim; ++i) {
      for (int j = 0; j < num_space_dim; ++j) {
        Scalar value(0.0);
        for (int alpha = 0; alpha < num_space_dim; ++alpha)
          value += jacobian(cell,ip,i,alpha) * jacobian(cell,ip,j,alpha);
        aux(cell,ip,i,j) = value;
      }
    }
  });
  PHX::Device::execution_space().fence();

  if(cache){
    covarient = aux;
    covarient_evaluated_ = true;
  }

  return aux;
}

template <typename Scalar>
typename IntegrationValues2<Scalar>::ConstArray_CellIPDimDim
IntegrationValues2<Scalar>::
getContravarientMatrix(const bool cache,
                       const bool force) const
{
  if(contravarient_evaluated_ and not force)
    return contravarient;

  // Only log time if values computed (i.e. don't log if values are already cached)
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::getContravarientMatrix()",get_contra_mat);

  MDFieldArrayFactory af(prefix_,true);

  const int num_space_dim = int_rule->topology->getDimension();
  const int num_ip = int_rule->num_points;

  auto cov = getCovarientMatrix(false,force);
  auto aux = af.template buildStaticArray<Scalar,Cell,IP,Dim,Dim>("contravarient",num_cells_, num_ip, num_space_dim,num_space_dim);

  const auto cell_range = std::make_pair(0,num_evaluate_cells_);
  auto s_contravarient = Kokkos::subview(aux.get_view(), cell_range,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
  auto s_covarient = Kokkos::subview(cov.get_view(), cell_range,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);

  Intrepid2::RealSpaceTools<PHX::Device::execution_space>::inverse(s_contravarient, s_covarient);
  PHX::Device::execution_space().fence();

  if(cache){
    contravarient = aux;
    contravarient_evaluated_ = true;
  }

  return aux;
}

template <typename Scalar>
typename IntegrationValues2<Scalar>::ConstArray_CellIP
IntegrationValues2<Scalar>::
getNormContravarientMatrix(const bool cache,
                           const bool force) const
{
  if(norm_contravarient_evaluated_ and not force)
    return norm_contravarient;

  // Only log time if values computed (i.e. don't log if values are already cached)
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::getNormContravarientMatrix()",get_norm_contr_mat);

  MDFieldArrayFactory af(prefix_,true);

  const int num_space_dim = int_rule->topology->getDimension();
  const int num_ip = int_rule->num_points;

  auto con = getContravarientMatrix(false,force).get_static_view();
  auto aux = af.template buildStaticArray<Scalar,Cell,IP>("norm_contravarient",num_cells_, num_ip);

  // norm of g_ij
  Kokkos::MDRangePolicy<PHX::Device,Kokkos::Rank<2>> policy({0,0},{num_evaluate_cells_,num_ip});
  Kokkos::parallel_for("evaluate norm_contravarient",policy,KOKKOS_LAMBDA (const int cell,const int ip) {
    aux(cell,ip) = 0.0;
    for (int i = 0; i < num_space_dim; ++i) {
      for (int j = 0; j < num_space_dim; ++j) {
        aux(cell,ip) += con(cell,ip,i,j) * con(cell,ip,i,j);
      }
    }
    aux(cell,ip) = Kokkos::sqrt(aux(cell,ip));
  });
  PHX::Device::execution_space().fence();

  if(cache){
    norm_contravarient = aux;
    norm_contravarient_evaluated_ = true;
  }

  return aux;
}

template <typename Scalar>
typename IntegrationValues2<Scalar>::ConstArray_CellIPDim
IntegrationValues2<Scalar>::
getCubaturePoints(const bool cache,
                  const bool force) const
{
  if(ip_coordinates_evaluated_ and not force)
    return ip_coordinates;

  // Only log time if values computed (i.e. don't log if values are already cached)
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::getCubaturePoints()",get_cub_pts);

  MDFieldArrayFactory af(prefix_,true);

  const int num_space_dim = int_rule->topology->getDimension();
  const int num_ip = int_rule->num_points;

  auto aux = af.template buildStaticArray<Scalar,Cell,IP,Dim>("ip_coordinates",num_cells_, num_ip, num_space_dim);

  using ID=panzer::IntegrationDescriptor;
  const bool is_cv = (int_rule->getType() == ID::CV_VOLUME) or (int_rule->getType() == ID::CV_SIDE) or (int_rule->getType() == ID::CV_BOUNDARY);
  const bool is_surface = int_rule->getType() == ID::SURFACE;

  auto node_coord = PHX::getNonConstDynRankViewFromConstMDField(getNodeCoordinates());

  if(is_cv){

    // CV integration uses a single call to map from physical space to the weighted measure - I assume this is slower than what we do with non-cv integration methods
    const auto cell_range = std::make_pair(0,num_evaluate_cells_);
    auto s_node_coord = Kokkos::subview(node_coord,    cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_cub_points = Kokkos::subview(aux.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());

    // TODO: We need to pull this apart for control volumes. Right now we calculate both weighted measures/norms and cubature points at the same time
    if(int_rule->cv_type == "side"){
      auto scratch = af.template buildStaticArray<Scalar,Cell,IP,Dim>("scratch",num_evaluate_cells_,num_ip,num_space_dim);
      intrepid_cubature->getCubature(s_cub_points,scratch.get_view(),s_node_coord);
    } else {
      // I think boundary is included as a weighted measure because it has a side embedded in intrepid_cubature
      TEUCHOS_ASSERT((int_rule->getType() == ID::CV_VOLUME) or (int_rule->getType() == ID::CV_BOUNDARY));
      auto scratch = af.template buildStaticArray<Scalar,Cell,IP>("scratch",num_evaluate_cells_,num_ip);
      intrepid_cubature->getCubature(s_cub_points,scratch.get_view(),s_node_coord);
    }

  } else if(is_surface){

    // Don't forget that since we are not caching this, we have to make sure the managed view remains alive while we use the non-const wrapper
    auto const_ref_coord = getCubaturePointsRef(false,force);
    auto ref_coord = PHX::getNonConstDynRankViewFromConstMDField(const_ref_coord);

    const auto cell_range = std::make_pair(0,num_evaluate_cells_);
    auto s_ref_coord  = Kokkos::subview(ref_coord,     cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_coord      = Kokkos::subview(aux.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_node_coord = Kokkos::subview(node_coord,    cell_range,Kokkos::ALL(),Kokkos::ALL());

    Intrepid2::CellTools<PHX::Device::execution_space> cell_tools;
    cell_tools.mapToPhysicalFrame(s_coord, s_ref_coord, s_node_coord, *(int_rule->topology));
  
  } else {

    // Don't forget that since we are not caching this, we have to make sure the managed view remains alive while we use the non-const wrapper
    auto const_ref_coord = getUniformCubaturePointsRef(false,force,false);
    auto ref_coord = PHX::getNonConstDynRankViewFromConstMDField(const_ref_coord);

    const auto cell_range = std::make_pair(0,num_evaluate_cells_);
    auto s_coord      = Kokkos::subview(aux.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_node_coord = Kokkos::subview(node_coord,    cell_range,Kokkos::ALL(),Kokkos::ALL());

    Intrepid2::CellTools<PHX::Device::execution_space> cell_tools;
    cell_tools.mapToPhysicalFrame(s_coord, ref_coord, s_node_coord, *(int_rule->topology));

    if(requires_permutation_)
      applyPermutation(aux, permutations_);

  }

  PHX::Device::execution_space().fence();

  if(cache){
    ip_coordinates = aux;
    ip_coordinates_evaluated_ = true;
  }

  return aux;
}


template <typename Scalar>
typename IntegrationValues2<Scalar>::ConstArray_CellIPDim
IntegrationValues2<Scalar>::
getCubaturePointsRef(const bool cache,
                     const bool force) const
{
  if(ref_ip_coordinates_evaluated_ and not force)
    return ref_ip_coordinates;

  // Only log time if values computed (i.e. don't log if values are already cached)
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::getCubaturePointsRef()",get_cub_pts_ref);

  using ID=panzer::IntegrationDescriptor;
  const bool is_surface = int_rule->getType() == ID::SURFACE;
  const bool is_cv = (int_rule->getType() == ID::CV_VOLUME) or (int_rule->getType() == ID::CV_SIDE) or (int_rule->getType() == ID::CV_BOUNDARY);

  const int num_space_dim = int_rule->topology->getDimension();
  const int num_ip = int_rule->num_points;

  MDFieldArrayFactory af(prefix_,true);

  Intrepid2::CellTools<PHX::Device::execution_space> cell_tools;

  auto aux = af.template buildStaticArray<Scalar,Cell,IP,Dim>("ref_ip_coordinates",num_cells_, num_ip, num_space_dim);

  if(is_cv){

    // Control volume reference points are actually generated from the physical points (i.e. reverse to everything else)

    auto node_coord = PHX::getNonConstDynRankViewFromConstMDField(getNodeCoordinates());

    // Don't forget that since we are not caching this, we have to make sure the managed view remains alive while we use the non-const wrapper
    auto const_coord = getCubaturePoints(false,force);
    auto coord = PHX::getNonConstDynRankViewFromConstMDField(const_coord);

    const auto cell_range = std::make_pair(0,num_evaluate_cells_);
    auto s_ref_coord  = Kokkos::subview(aux.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_coord      = Kokkos::subview(coord,         cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_node_coord = Kokkos::subview(node_coord,    cell_range,Kokkos::ALL(),Kokkos::ALL());

    cell_tools.mapToReferenceFrame(s_ref_coord, s_coord, s_node_coord, *(int_rule->topology));

  } else if(is_surface){

    const auto & cell_topology = *int_rule->topology;
    const int cell_dim = cell_topology.getDimension();
    const int num_sides = (cell_dim==1) ? 2 : cell_topology.getSideCount();
    const int subcell_dim = cell_dim-1;
    const int num_points_on_face = num_ip / num_sides;
    const int order = int_rule->getOrder();

    // Scratch space for storing the points for each side of the cell
    auto side_cub_points2 = af.template buildStaticArray<Scalar,IP,Dim>("side_cub_points",num_points_on_face,cell_dim);

    Intrepid2::DefaultCubatureFactory cubature_factory;

    // We get to build up our cubature one side at a time
    for(int side=0; side<num_sides; ++side) {

      const int point_offset = side*num_points_on_face;

      // Get the cubature for the side
      if(cell_dim==1){
        // In 1D the surface point is either on the left side of the cell, or the right side
        auto side_cub_points_host = Kokkos::create_mirror_view(side_cub_points2.get_view());
        side_cub_points_host(0,0) = (side==0)? -1. : 1.;
        Kokkos::deep_copy(side_cub_points2.get_view(),side_cub_points_host);
      } else {

        // Get the face topology from the cell topology
        const shards::CellTopology face_topology(cell_topology.getCellTopologyData(subcell_dim,side));

        // Create a cubature for the face of the cell
        auto ic = cubature_factory.create<PHX::Device::execution_space,double,double>(face_topology,order);
        auto tmp_side_cub_weights = Kokkos::DynRankView<double,PHX::Device>("tmp_side_cub_weights",num_points_on_face);
        auto tmp_side_cub_points = Kokkos::DynRankView<double,PHX::Device>("tmp_side_cub_points",num_points_on_face,subcell_dim);

        // Get the reference face points
        ic->getCubature(tmp_side_cub_points, tmp_side_cub_weights);

        // Convert from reference face points to reference cell points
        cell_tools.mapToReferenceSubcell(side_cub_points2.get_view(), tmp_side_cub_points, subcell_dim, side, cell_topology);
      }

      PHX::Device::execution_space().fence();

      // Copy from the side allocation to the surface allocation
      Kokkos::MDRangePolicy<PHX::Device::execution_space,Kokkos::Rank<3>> policy({0,0,0},{num_evaluate_cells_,num_points_on_face, num_space_dim});
      Kokkos::parallel_for("copy values",policy,KOKKOS_LAMBDA (const int cell,const int point, const int dim) {
        aux(cell,point_offset + point,dim) = side_cub_points2(point,dim);
      });
      PHX::Device::execution_space().fence();
    }

  } else {

    // Haven't set this up yet
    TEUCHOS_TEST_FOR_EXCEPT_MSG(int_rule->isSide() && num_space_dim==1,
                                "ERROR: 0-D quadrature rule infrastructure does not exist!!! Will not be able to do "
                                 << "non-natural integration rules.");

    auto cub_points2 = getUniformCubaturePointsRef(false,force,false);

    Kokkos::MDRangePolicy<PHX::Device,Kokkos::Rank<3>> policy({0,0,0},{num_evaluate_cells_,num_ip,num_space_dim});
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int & cell, const int & ip, const int & dim){
      aux(cell,ip,dim) = cub_points2(ip,dim);
    });
  }

  PHX::Device::execution_space().fence();

  if(requires_permutation_)
    applyPermutation(aux, permutations_);

  if(cache){
    ref_ip_coordinates = aux;
    ref_ip_coordinates_evaluated_ = true;
  }

  return aux;
}

template <typename Scalar>
void
IntegrationValues2<Scalar>::
evaluateEverything()
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::integrationValues2::evaluateEverything()",eval_everything);

  using ID=panzer::IntegrationDescriptor;
  const bool is_surface = int_rule->getType() == ID::SURFACE;
  const bool is_cv = (int_rule->getType() == ID::CV_VOLUME) or (int_rule->getType() == ID::CV_SIDE) or (int_rule->getType() == ID::CV_BOUNDARY);
  const bool is_side = int_rule->isSide();

  // This will force all values to be re-evaluated
  resetArrays();

  // Base cubature stuff
  if(is_cv){
    getCubaturePoints(true);
    getCubaturePointsRef(true);
  } else {
    if(not is_surface){
      getUniformCubaturePointsRef(true,true,requires_permutation_);
      getUniformCubatureWeightsRef(true,true,requires_permutation_);
      if(is_side)
        getUniformSideCubaturePointsRef(true,true,requires_permutation_);
    }
    getCubaturePointsRef(true);
    getCubaturePoints(true);
  }

  // Measure stuff
  getJacobian(true);
  getJacobianDeterminant(true);
  getJacobianInverse(true);
  if(int_rule->cv_type == "side")
    getWeightedNormals(true);
  else
    getWeightedMeasure(true);

  // Surface stuff
  if(is_surface){
    getSurfaceNormals(true);
    getSurfaceRotationMatrices(true);
  }

  // Stabilization stuff
  if(not (is_surface or is_cv)){
    getContravarientMatrix(true);
    getCovarientMatrix(true);
    getNormContravarientMatrix(true);
  }
}

#define INTEGRATION_VALUES2_INSTANTIATION(SCALAR) \
    template class IntegrationValues2<SCALAR>;

INTEGRATION_VALUES2_INSTANTIATION(panzer::Traits::RealType)

// Disabled FAD support due to long build times on cuda (in debug mode
// it takes multiple hours on some platforms). If we need
// sensitivities wrt coordinates, we can reenable.

// INTEGRATION_VALUES2_INSTANTIATION(panzer::Traits::FadType)

}

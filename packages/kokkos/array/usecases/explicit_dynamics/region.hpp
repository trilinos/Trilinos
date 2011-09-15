#include "BoxMeshFixture.hpp"


template<typename Scalar , class DeviceType>
struct Region;

template<typename Scalar>
struct Region<Scalar ,KOKKOS_MACRO_DEVICE>{

  typedef KOKKOS_MACRO_DEVICE       device_type;
  typedef device_type::size_type    size_type;

  typedef Kokkos::MDArrayView<Scalar,device_type>   scalar_array;
  typedef Kokkos::MDArrayView<int,device_type>      int_array;

  template <class Mesh>
  Region(
      int num_states,
      const Mesh & mesh,
      Scalar arg_lin_bulk_visc,
      Scalar arg_quad_bulk_visc,
      Scalar arg_hg_stiffness,
      Scalar arg_hg_viscosity,
      Scalar youngs_modulus,
      Scalar poissons_ratio,
      Scalar arg_density
      )
    : num_nodes(mesh.nnodes)
    , num_elements(mesh.nelems)
    , spatial_dim(mesh.SpatialDim)
    , lin_bulk_visc(arg_lin_bulk_visc)
    , quad_bulk_visc(arg_quad_bulk_visc)
    , hg_stiffness(arg_hg_stiffness)
    , hg_viscosity(arg_hg_viscosity)
    , two_mu(youngs_modulus/(1.0+poissons_ratio))
    , bulk_modulus(youngs_modulus/(3*(1.0-2.0*poissons_ratio)))
    , density(arg_density)
    // mesh
    , elem_node_connectivity(Kokkos::create_mdarray<int_array>(mesh.elem_node_ids.dimension(0),mesh.elem_node_ids.dimension(1)))
    , node_elem_ids(Kokkos::create_mdarray<int_array>(mesh.node_elem_ids.dimension(0),mesh.node_elem_ids.dimension(1)))
    , node_elem_offset(Kokkos::create_mdarray<int_array>(mesh.node_elem_offset.dimension(0),mesh.node_elem_offset.dimension(1)))
    , model_coords(Kokkos::create_mdarray<scalar_array>(num_nodes,spatial_dim))
    // input/output
    , displacement(Kokkos::create_mdarray<scalar_array>(num_nodes,spatial_dim,num_states))
    , velocity(Kokkos::create_mdarray<scalar_array>(num_nodes,spatial_dim,num_states))
    , acceleration(Kokkos::create_mdarray<scalar_array>(num_nodes,spatial_dim))
    , nodal_mass(Kokkos::create_mdarray<scalar_array>(num_nodes))
    , elem_mass(Kokkos::create_mdarray<scalar_array>(num_elements))
    , stress_new(Kokkos::create_mdarray<scalar_array>(num_elements,6))
    , internal_force(Kokkos::create_mdarray<scalar_array>(num_nodes,spatial_dim))
    , internal_energy(Kokkos::create_mdarray<scalar_array>(num_elements))
    // temporary arrays
    , hg_resist(Kokkos::create_mdarray<scalar_array>(num_elements,12,num_states))
    , rotation(Kokkos::create_mdarray<scalar_array>(num_elements,9,num_states))
    , gradop12(Kokkos::create_mdarray<scalar_array>(num_elements,3,8))
    , element_force(Kokkos::create_mdarray<scalar_array>(num_elements,3,8))
    , hgop(Kokkos::create_mdarray<scalar_array>(num_elements,32,num_states))
    , vel_grad(Kokkos::create_mdarray<scalar_array>(num_elements,9))
    , stretch(Kokkos::create_mdarray<scalar_array>(num_elements,6))
    , vorticity(Kokkos::create_mdarray<scalar_array>(num_elements,3))
    , rot_stretch(Kokkos::create_mdarray<scalar_array>(num_elements,6))
    , rot_stress(Kokkos::create_mdarray<scalar_array>(num_elements,6))
    , mid_vol(Kokkos::create_mdarray<scalar_array>(num_elements))
    , shrmod(Kokkos::create_mdarray<scalar_array>(num_elements))
    , dilmod(Kokkos::create_mdarray<scalar_array>(num_elements))
    , elem_t_step(Kokkos::create_mdarray<scalar_array>(num_elements))
    , hg_energy(Kokkos::create_mdarray<scalar_array>(num_elements))
  {
    Kokkos::deep_copy(elem_node_connectivity, mesh.elem_node_ids);
    Kokkos::deep_copy(node_elem_ids, mesh.node_elem_ids);
    Kokkos::deep_copy(node_elem_offset, mesh.node_elem_offset);
    Kokkos::deep_copy(model_coords, mesh.node_coords);
  }

  const int num_nodes;
  const int num_elements;
  const int spatial_dim;

  const Scalar  lin_bulk_visc;
  const Scalar  quad_bulk_visc;
  const Scalar  hg_stiffness;
  const Scalar  hg_viscosity;
  const Scalar  two_mu;
  const Scalar  bulk_modulus;
  const Scalar  density;

  // mesh connectivity
  int_array elem_node_connectivity;
  int_array node_elem_ids;
  int_array node_elem_offset;

  // input / output
  scalar_array  model_coords;
  scalar_array  displacement;
  scalar_array  velocity;
  scalar_array  acceleration;
  scalar_array  nodal_mass;
  scalar_array  elem_mass;
  scalar_array  stress_new;
  scalar_array  internal_force;
  scalar_array  internal_energy;

  //tempory arrays
  scalar_array  hg_resist;
  scalar_array  rotation;
  scalar_array  gradop12;
  scalar_array  element_force;
  scalar_array  hgop;
  scalar_array  vel_grad;
  scalar_array  stretch;
  scalar_array  vorticity;
  scalar_array  rot_stretch;
  scalar_array  rot_stress;
  scalar_array  mid_vol;
  scalar_array  shrmod;
  scalar_array  dilmod;
  scalar_array  elem_t_step;
  scalar_array  hg_energy;


};

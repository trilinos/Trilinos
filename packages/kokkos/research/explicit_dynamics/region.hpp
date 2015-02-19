/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

template<typename Scalar , class DeviceType>
struct Region;

template<typename Scalar>
struct Region<Scalar ,KOKKOS_MACRO_DEVICE>{

  typedef KOKKOS_MACRO_DEVICE       execution_space;
  typedef execution_space::size_type    size_type;

  typedef Kokkos::MDArray<Scalar,execution_space>   scalar_array;
  typedef Kokkos::MDArray<int,execution_space>      int_array;

  typedef Kokkos::Value<Scalar,execution_space>     scalar;

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
    , dt(Kokkos::create_value<scalar>())
    , prev_dt(Kokkos::create_value<scalar>())
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
    //, elem_t_step(Kokkos::create_mdarray<scalar_array>(num_elements))
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
  scalar_array  model_coords;

  // input / output
  scalar        dt;
  scalar        prev_dt;
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
  //scalar_array  elem_t_step;
  scalar_array  hg_energy;


};

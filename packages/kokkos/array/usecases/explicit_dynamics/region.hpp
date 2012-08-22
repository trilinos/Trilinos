/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

template<typename Scalar , class DeviceType>
struct Region;

template<typename Scalar>
struct Region<Scalar ,KOKKOSARRAY_MACRO_DEVICE>{

  typedef KOKKOSARRAY_MACRO_DEVICE       device_type;
  typedef device_type::size_type    size_type;

  typedef KokkosArray::MDArray<Scalar,device_type>   scalar_array;
  typedef KokkosArray::MDArray<int,device_type>      int_array;

  typedef KokkosArray::Value<Scalar,device_type>     scalar;

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
    , elem_node_connectivity(KokkosArray::create_mdarray<int_array>(mesh.elem_node_ids.dimension(0),mesh.elem_node_ids.dimension(1)))
    , node_elem_ids(KokkosArray::create_mdarray<int_array>(mesh.node_elem_ids.dimension(0),mesh.node_elem_ids.dimension(1)))
    , node_elem_offset(KokkosArray::create_mdarray<int_array>(mesh.node_elem_offset.dimension(0),mesh.node_elem_offset.dimension(1)))
    , model_coords(KokkosArray::create_mdarray<scalar_array>(num_nodes,spatial_dim))
    // input/output
    , dt(KokkosArray::create_value<scalar>())
    , prev_dt(KokkosArray::create_value<scalar>())
    , displacement(KokkosArray::create_mdarray<scalar_array>(num_nodes,spatial_dim,num_states))
    , velocity(KokkosArray::create_mdarray<scalar_array>(num_nodes,spatial_dim,num_states))
    , acceleration(KokkosArray::create_mdarray<scalar_array>(num_nodes,spatial_dim))
    , nodal_mass(KokkosArray::create_mdarray<scalar_array>(num_nodes))
    , elem_mass(KokkosArray::create_mdarray<scalar_array>(num_elements))
    , stress_new(KokkosArray::create_mdarray<scalar_array>(num_elements,6))
    , internal_force(KokkosArray::create_mdarray<scalar_array>(num_nodes,spatial_dim))
    , internal_energy(KokkosArray::create_mdarray<scalar_array>(num_elements))
    // temporary arrays
    , hg_resist(KokkosArray::create_mdarray<scalar_array>(num_elements,12,num_states))
    , rotation(KokkosArray::create_mdarray<scalar_array>(num_elements,9,num_states))
    , gradop12(KokkosArray::create_mdarray<scalar_array>(num_elements,3,8))
    , element_force(KokkosArray::create_mdarray<scalar_array>(num_elements,3,8))
    , hgop(KokkosArray::create_mdarray<scalar_array>(num_elements,32,num_states))
    , vel_grad(KokkosArray::create_mdarray<scalar_array>(num_elements,9))
    , stretch(KokkosArray::create_mdarray<scalar_array>(num_elements,6))
    , vorticity(KokkosArray::create_mdarray<scalar_array>(num_elements,3))
    , rot_stretch(KokkosArray::create_mdarray<scalar_array>(num_elements,6))
    , rot_stress(KokkosArray::create_mdarray<scalar_array>(num_elements,6))
    , mid_vol(KokkosArray::create_mdarray<scalar_array>(num_elements))
    , shrmod(KokkosArray::create_mdarray<scalar_array>(num_elements))
    , dilmod(KokkosArray::create_mdarray<scalar_array>(num_elements))
    //, elem_t_step(KokkosArray::create_mdarray<scalar_array>(num_elements))
    , hg_energy(KokkosArray::create_mdarray<scalar_array>(num_elements))
  {
    KokkosArray::deep_copy(elem_node_connectivity, mesh.elem_node_ids);
    KokkosArray::deep_copy(node_elem_ids, mesh.node_elem_ids);
    KokkosArray::deep_copy(node_elem_offset, mesh.node_elem_offset);
    KokkosArray::deep_copy(model_coords, mesh.node_coords);
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

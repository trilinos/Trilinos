// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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

// ***********************************************************
// * Specializations of setupArrays() for different array types
// ***********************************************************

namespace panzer {

// * Specialization for Kokkos::DynRankView<double,PHX::Device>
template <typename Scalar>
void IntegrationValues2<Scalar>::
setupArraysForNodeRule(const Teuchos::RCP<const panzer::IntegrationRule>& ir)
{
  MDFieldArrayFactory af(prefix,alloc_arrays);

  int num_nodes = ir->topology->getNodeCount();
  int num_cells = ir->workset_size;
  int num_space_dim = ir->topology->getDimension();

  int num_ip = 1;

  dyn_cub_points = af.template buildArray<double,IP,Dim>("cub_points",num_ip, num_space_dim);
  dyn_cub_weights = af.template buildArray<double,IP>("cub_weights",num_ip);

  cub_points = af.template buildStaticArray<Scalar,IP,Dim>("cub_points",num_ip, num_space_dim);

  if (ir->cv_type == "none" && ir->isSide()) {
    dyn_side_cub_points = af.template buildArray<double,IP,Dim>("side_cub_points",num_ip, ir->side_topology->getDimension());
    side_cub_points = af.template buildStaticArray<Scalar,IP,Dim>("side_cub_points",num_ip,ir->side_topology->getDimension());
  }

  if (ir->cv_type != "none") {
    dyn_phys_cub_points = af.template buildArray<double,Cell,IP,Dim>("phys_cub_points",num_cells, num_ip, num_space_dim);
    dyn_phys_cub_weights = af.template buildArray<double,Cell,IP>("phys_cub_weights",num_cells, num_ip);
    if (ir->cv_type == "side") {
      dyn_phys_cub_norms = af.template buildArray<double,Cell,IP,Dim>("phys_cub_norms",num_cells, num_ip, num_space_dim);
    }
  }

  dyn_node_coordinates = af.template buildArray<double,Cell,BASIS,Dim>("node_coordinates",num_cells, num_nodes, num_space_dim);

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

  weighted_normals = af.template buildStaticArray<Scalar,Cell,IP,Dim>("weighted normal",num_cells, num_ip,num_space_dim);

  surface_normals = af.template buildStaticArray<Scalar,Cell,IP,Dim>("surface_normals",num_cells, num_ip,num_space_dim);

  surface_rotation_matrices = af.template buildStaticArray<Scalar,Cell,IP,Dim,Dim>("surface_rotation_matrices",num_cells, num_ip,3,3);

}

template <typename Scalar>
void IntegrationValues2<Scalar>::
setupArrays(const Teuchos::RCP<const panzer::IntegrationRule>& ir)
{
  MDFieldArrayFactory af(prefix,alloc_arrays);

  typedef panzer::IntegrationDescriptor ID;

  int_rule = ir;

  int num_nodes = ir->topology->getNodeCount();
  int num_cells = ir->workset_size;
  int num_space_dim = ir->topology->getDimension();

  // specialize content if this is quadrature at anode
  if(num_space_dim==1 && ir->isSide()) {
    setupArraysForNodeRule(ir);
    return;
  }

  TEUCHOS_ASSERT(ir->getType() != ID::NONE);
  intrepid_cubature = getIntrepidCubature(*ir);

  int num_ip = ir->num_points;

//  Intrepid2::DefaultCubatureFactory cubature_factory;
//
//  if (ir->cv_type == "side")
//    intrepid_cubature = Teuchos::rcp(new Intrepid2::CubatureControlVolumeSide<PHX::Device::execution_space,double,double>(*ir->topology));
//
//  else if (ir->cv_type == "volume")
//    intrepid_cubature = Teuchos::rcp(new Intrepid2::CubatureControlVolume<PHX::Device::execution_space,double,double>(*ir->topology));
//
//  else if (ir->cv_type == "boundary" && ir->isSide())
//    intrepid_cubature = Teuchos::rcp(new Intrepid2::CubatureControlVolumeBoundary<PHX::Device::execution_space,double,double>(*ir->topology,ir->side));
//
//  else if (ir->cv_type == "none" && ir->isSide())
//    intrepid_cubature = cubature_factory.create<PHX::Device::execution_space,double,double>(*(ir->side_topology),
//                                                                                            ir->cubature_degree);
//  else
//    intrepid_cubature = cubature_factory.create<PHX::Device::execution_space,double,double>(*(ir->topology),
//                                                                                            ir->cubature_degree);
//  int num_ip = intrepid_cubature->getNumPoints();

  dyn_cub_points = af.template buildArray<double,IP,Dim>("cub_points",num_ip, num_space_dim);
  dyn_cub_weights = af.template buildArray<double,IP>("cub_weights",num_ip);

  cub_points = af.template buildStaticArray<Scalar,IP,Dim>("cub_points",num_ip, num_space_dim);

  if (ir->isSide() && ir->cv_type == "none") {
    dyn_side_cub_points = af.template buildArray<double,IP,Dim>("side_cub_points",num_ip, ir->side_topology->getDimension());
    side_cub_points = af.template buildStaticArray<Scalar,IP,Dim>("side_cub_points",num_ip,ir->side_topology->getDimension());
  }

  if (ir->cv_type != "none") {
    dyn_phys_cub_points = af.template buildArray<double,Cell,IP,Dim>("phys_cub_points",num_cells, num_ip, num_space_dim);
    dyn_phys_cub_weights = af.template buildArray<double,Cell,IP>("phys_cub_weights",num_cells, num_ip);
    if (ir->cv_type == "side") {
      dyn_phys_cub_norms = af.template buildArray<double,Cell,IP,Dim>("phys_cub_norms",num_cells, num_ip, num_space_dim);
    }
  }

  dyn_node_coordinates = af.template buildArray<double,Cell,BASIS,Dim>("node_coordinates",num_cells, num_nodes, num_space_dim);

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

  scratch_for_compute_side_measure =
      af.template buildStaticArray<Scalar,Point>("scratch_for_compute_side_measure", jac.get_view().span());

}

template <typename Scalar>
Teuchos::RCP<Intrepid2::Cubature<PHX::Device::execution_space,double,double>>
IntegrationValues2<Scalar>::
getIntrepidCubature(const panzer::IntegrationRule & ir) const
{
  typedef panzer::IntegrationDescriptor ID;
  Teuchos::RCP<Intrepid2::Cubature<PHX::Device::execution_space,double,double> > ic;

  Intrepid2::DefaultCubatureFactory cubature_factory;

  if(ir.getType() == ID::CV_SIDE){
    ic = Teuchos::rcp(new Intrepid2::CubatureControlVolumeSide<PHX::Device::execution_space,double,double>(*ir.topology));
  } else if(ir.getType() == ID::CV_VOLUME){
    ic = Teuchos::rcp(new Intrepid2::CubatureControlVolume<PHX::Device::execution_space,double,double>(*ir.topology));
  } else if(ir.getType() == ID::CV_BOUNDARY){
    ic = Teuchos::rcp(new Intrepid2::CubatureControlVolumeBoundary<PHX::Device::execution_space,double,double>(*ir.topology,ir.getSide()));
  } else {
    if(ir.getType() == ID::VOLUME){
      ic = cubature_factory.create<PHX::Device::execution_space,double,double>(*(ir.topology),ir.getOrder());
    } else if(ir.getType() == ID::SIDE){
      ic = cubature_factory.create<PHX::Device::execution_space,double,double>(*(ir.side_topology),ir.getOrder());
    } else if(ir.getType() == ID::SURFACE){
      // closed surface integrals don't exist in intrepid.
    } else {
      TEUCHOS_ASSERT(false);
    }
  }

  return ic;
}


// ***********************************************************
// * Evaluation of values - NOT specialized
// ***********************************************************
template <typename Scalar>
void IntegrationValues2<Scalar>::
evaluateValues(const PHX::MDField<Scalar,Cell,NODE,Dim> & in_node_coordinates,
               const int in_num_cells,
               const Teuchos::RCP<const SubcellConnectivity> & face_connectivity)
{
  typedef panzer::IntegrationDescriptor ID;
  const bool is_surface = int_rule->getType() == ID::SURFACE;
  const bool is_cv = (int_rule->getType() == ID::CV_VOLUME) or (int_rule->getType() == ID::CV_SIDE) or (int_rule->getType() == ID::CV_BOUNDARY);

  TEUCHOS_ASSERT(not (is_surface and is_cv));

  if(is_surface){
    TEUCHOS_TEST_FOR_EXCEPT_MSG(face_connectivity == Teuchos::null,
                                "IntegrationValues2::evaluateValues : Surface integration requires the face connectivity");
    generateSurfaceCubatureValues(in_node_coordinates,in_num_cells,*face_connectivity);
  } else if (is_cv) {
    getCubatureCV(in_node_coordinates, in_num_cells);
    evaluateValuesCV(in_node_coordinates, in_num_cells);
  } else {
    getCubature(in_node_coordinates, in_num_cells);
    evaluateRemainingValues(in_node_coordinates, in_num_cells);
  }
}

template <typename Scalar>
void IntegrationValues2<Scalar>::
getCubature(const PHX::MDField<Scalar,Cell,NODE,Dim>& in_node_coordinates,
            const int in_num_cells)
{

  int num_space_dim = int_rule->topology->getDimension();
  if (int_rule->isSide() && num_space_dim==1) {
    std::cout << "WARNING: 0-D quadrature rule ifrastructure does not exist!!! Will not be able to do "
        << "non-natural integration rules.";
    return;
  }

  Intrepid2::CellTools<PHX::Device::execution_space> cell_tools;

  if (!int_rule->isSide())
    intrepid_cubature->getCubature(dyn_cub_points.get_view(), dyn_cub_weights.get_view());
  else {
    intrepid_cubature->getCubature(dyn_side_cub_points.get_view(), dyn_cub_weights.get_view());

    cell_tools.mapToReferenceSubcell(dyn_cub_points.get_view(),
                                     dyn_side_cub_points.get_view(),
                                     int_rule->spatial_dimension-1,
                                     int_rule->side,
                                     *(int_rule->topology));
  }

  // IP coordinates
  const int num_cells = in_num_cells < 0 ? in_node_coordinates.extent(0) : in_num_cells;
  auto s_ip_coordinates = Kokkos::subview(ip_coordinates.get_view(),std::make_pair(0,num_cells),Kokkos::ALL(),Kokkos::ALL());
  auto s_in_node_coordinates = Kokkos::subview(in_node_coordinates.get_view(),std::make_pair(0,num_cells),Kokkos::ALL(),Kokkos::ALL());
  cell_tools.mapToPhysicalFrame(s_ip_coordinates,
                                dyn_cub_points.get_view(),
                                s_in_node_coordinates,
                                *(int_rule->topology));
}

template <typename Scalar>
void IntegrationValues2<Scalar>::
convertNormalToRotationMatrix(const Scalar normal[3], Scalar transverse[3], Scalar binormal[3])
{
  using T = Scalar;
  
  const T n  = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);

  // If this fails then the geometry for this cell is probably undefined
  if(n > 0.){
    // Make sure transverse is not parallel to normal within some margin of error
    transverse[0]=0.;transverse[1]=1.;transverse[2]=0.;
    if(std::fabs(normal[0]*transverse[0]+normal[1]*transverse[1])>0.9){
      transverse[0]=1.;transverse[1]=0.;
    }

    const T nt = normal[0]*transverse[0]+normal[1]*transverse[1]+normal[2]*transverse[2];

    // Note normal has unit length
    const T mult = nt/(n*n); // = nt

    // Remove normal projection from transverse
    for(int dim=0;dim<3;++dim){
      transverse[dim] = transverse[dim] - mult * normal[dim];
    }

    const T t = sqrt(transverse[0]*transverse[0]+transverse[1]*transverse[1]+transverse[2]*transverse[2]);
    TEUCHOS_ASSERT(t != 0.);
    for(int dim=0;dim<3;++dim){
      transverse[dim] /= t;
    }

    // We assume a right handed system such that b = n \times t
    binormal[0] = (normal[1] * transverse[2] - normal[2] * transverse[1]);
    binormal[1] = (normal[2] * transverse[0] - normal[0] * transverse[2]);
    binormal[2] = (normal[0] * transverse[1] - normal[1] * transverse[0]);

    // Normalize binormal
    const T b = sqrt(binormal[0]*binormal[0]+binormal[1]*binormal[1]+binormal[2]*binormal[2]);
    for(int dim=0;dim<3;++dim){
      binormal[dim] /= b;
    }
  } else {
    transverse[0] = 0.;
    transverse[1] = 0.;
    transverse[2] = 0.;
    binormal[0] = 0.;
    binormal[1] = 0.;
    binormal[2] = 0.;
  }
}
  
template <typename Scalar>
void IntegrationValues2<Scalar>::
swapQuadraturePoints(int cell,
                     int a,
                     int b)
{
  const int new_cell_point = a;
  const int old_cell_point = b;

  const int cell_dim = ref_ip_coordinates.extent(2);

#ifdef PANZER_DEBUG
  TEUCHOS_ASSERT(cell < ip_coordinates.extent_int(0));
  TEUCHOS_ASSERT(a < ip_coordinates.extent_int(1));
  TEUCHOS_ASSERT(b < ip_coordinates.extent_int(1));
  TEUCHOS_ASSERT(cell >= 0);
  TEUCHOS_ASSERT(a >= 0);
  TEUCHOS_ASSERT(b >= 0);
#endif
  
  Scalar hold;

  hold = weighted_measure(cell,new_cell_point);
  weighted_measure(cell,new_cell_point) = weighted_measure(cell,old_cell_point);
  weighted_measure(cell,old_cell_point) = hold;

  hold = jac_det(cell,new_cell_point);
  jac_det(cell,new_cell_point) = jac_det(cell,old_cell_point);
  jac_det(cell,old_cell_point) = hold;

  for(int dim=0;dim<cell_dim;++dim){

    hold = ref_ip_coordinates(cell,new_cell_point,dim);
    ref_ip_coordinates(cell,new_cell_point,dim) = ref_ip_coordinates(cell,old_cell_point,dim);
    ref_ip_coordinates(cell,old_cell_point,dim) = hold;

    hold = ip_coordinates(cell,new_cell_point,dim);
    ip_coordinates(cell,new_cell_point,dim) = ip_coordinates(cell,old_cell_point,dim);
    ip_coordinates(cell,old_cell_point,dim) = hold;

    hold = surface_normals(cell,new_cell_point,dim);
    surface_normals(cell,new_cell_point,dim) = surface_normals(cell,old_cell_point,dim);
    surface_normals(cell,old_cell_point,dim) = hold;

    for(int dim2=0;dim2<cell_dim;++dim2){

      hold = jac(cell,new_cell_point,dim,dim2);
      jac(cell,new_cell_point,dim,dim2) = jac(cell,old_cell_point,dim,dim2);
      jac(cell,old_cell_point,dim,dim2) = hold;

      hold = jac_inv(cell,new_cell_point,dim,dim2);
      jac_inv(cell,new_cell_point,dim,dim2) = jac_inv(cell,old_cell_point,dim,dim2);
      jac_inv(cell,old_cell_point,dim,dim2) = hold;
    }
  }

  // Rotation matrices are always in 3D
  for(int dim=0; dim<3; ++dim){
    for(int dim2=0; dim2<3; ++dim2){
      hold = surface_rotation_matrices(cell,new_cell_point,dim,dim2);
      surface_rotation_matrices(cell,new_cell_point,dim,dim2) = surface_rotation_matrices(cell,old_cell_point,dim,dim2);
      surface_rotation_matrices(cell,old_cell_point,dim,dim2) = hold;
    }
  }
}

template <typename Scalar>
void IntegrationValues2<Scalar>::
generateSurfaceCubatureValues(const PHX::MDField<Scalar,Cell,NODE,Dim>& in_node_coordinates,
                              const int in_num_cells,
                              const SubcellConnectivity & face_connectivity)
{

  TEUCHOS_ASSERT(int_rule->getType() == IntegrationDescriptor::SURFACE);

  Intrepid2::CellTools<PHX::Device::execution_space> cell_tools;

  const shards::CellTopology & cell_topology = *(int_rule->topology);
  const panzer::IntegrationRule & ir = *int_rule;

  const int num_cells = in_num_cells < 0 ? in_node_coordinates.extent(0) : in_num_cells;

  // Copy over coordinates
  {
    const int num_nodes = in_node_coordinates.extent(1);
    const int num_dims = in_node_coordinates.extent(2);

    for(int cell=0; cell<num_cells; ++cell){
      for(int node=0; node<num_nodes; ++node){
        for(int dim=0; dim<num_dims; ++dim){
          node_coordinates(cell,node,dim) = in_node_coordinates(cell,node,dim);
        }
      }
    }
  }

  // NOTE: We are assuming that each face can have a different number of points.
  // Not sure if this is necessary, but it requires a lot of additional allocations

  const int cell_dim = cell_topology.getDimension();
  const int subcell_dim = cell_topology.getDimension()-1;
  const int num_subcells = cell_topology.getSubcellCount(subcell_dim);

  Intrepid2::DefaultCubatureFactory cubature_factory;

  // We get to build up our cubature one face at a time
  {
    int point_offset=0;
    for(int subcell_index=0; subcell_index<num_subcells; ++subcell_index) {

      // Default for 1D
      int num_points_on_face = 1;

      // Get the cubature for the side
      Kokkos::DynRankView<double,PHX::Device> tmp_side_cub_weights;
      Kokkos::DynRankView<double,PHX::Device> tmp_side_cub_points;
      if(cell_dim==1){
        tmp_side_cub_weights = Kokkos::DynRankView<double,PHX::Device>("tmp_side_cub_weights",num_points_on_face);
        tmp_side_cub_points = Kokkos::DynRankView<double,PHX::Device>("cell_tmp_side_cub_points",num_points_on_face,cell_dim);
        tmp_side_cub_weights(0)=1.;
        tmp_side_cub_points(0,0) = (subcell_index==0)? -1. : 1.;
      } else {

        // Get the face topology from the cell topology
        const shards::CellTopology face_topology(cell_topology.getCellTopologyData(subcell_dim,subcell_index));

        auto ic = cubature_factory.create<PHX::Device::execution_space,double,double>(face_topology,ir.getOrder());
        num_points_on_face = ic->getNumPoints();

        tmp_side_cub_weights = Kokkos::DynRankView<double,PHX::Device>("tmp_side_cub_weights",num_points_on_face);
        tmp_side_cub_points = Kokkos::DynRankView<double,PHX::Device>("cell_tmp_side_cub_points",num_points_on_face,cell_dim);

        auto subcell_cub_points = Kokkos::DynRankView<double,PHX::Device>("subcell_cub_points",num_points_on_face,subcell_dim);

        // Get the reference face points
        ic->getCubature(subcell_cub_points, tmp_side_cub_weights);

      // Convert from reference face points to reference cell points
        cell_tools.mapToReferenceSubcell(tmp_side_cub_points,
                                         subcell_cub_points,
                                         subcell_dim,
                                         subcell_index,
                                       cell_topology);
      }


      for(int local_point=0;local_point<num_points_on_face;++local_point){
        const int point = point_offset + local_point;
        for(int dim=0;dim<cell_dim;++dim){
          cub_points(point,dim) = tmp_side_cub_points(local_point,dim);
        }
      }


      // Map from side points to physical points
      auto side_ip_coordinates = Kokkos::DynRankView<Scalar,PHX::Device>("side_ip_coordinates",num_cells,num_points_on_face,cell_dim);
      auto s_node_coordinates = Kokkos::subview(node_coordinates.get_view(),std::make_pair(0,num_cells),Kokkos::ALL,Kokkos::ALL);
      cell_tools.mapToPhysicalFrame(side_ip_coordinates,
                                    tmp_side_cub_points,
                                    s_node_coordinates,
                                    cell_topology);

      // Create a jacobian and his friends for this side
      auto side_jacobian = Kokkos::DynRankView<Scalar,PHX::Device>("side_jac",num_cells,num_points_on_face,cell_dim,cell_dim);
      cell_tools.setJacobian(side_jacobian,
                             tmp_side_cub_points,
                             s_node_coordinates,
                             cell_topology);

      auto side_inverse_jacobian = Kokkos::DynRankView<Scalar,PHX::Device>("side_inv_jac",num_cells,num_points_on_face,cell_dim,cell_dim);
      cell_tools.setJacobianInv(side_inverse_jacobian, side_jacobian);

      auto side_det_jacobian = Kokkos::DynRankView<Scalar,PHX::Device>("side_det_jac",num_cells,num_points_on_face);
      cell_tools.setJacobianDet(side_det_jacobian, side_jacobian);

      // Calculate measures (quadrature weights in physical space) for this side
      auto side_weighted_measure = Kokkos::DynRankView<Scalar,PHX::Device>("side_weighted_measure",num_cells,num_points_on_face);
      if(cell_dim == 1){
        Kokkos::deep_copy(side_weighted_measure, tmp_side_cub_weights(0));
      } else if(cell_dim == 2){
        Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
          computeEdgeMeasure(side_weighted_measure, side_jacobian, tmp_side_cub_weights,
                             subcell_index,cell_topology,
                             scratch_for_compute_side_measure.get_view());
      } else if(cell_dim == 3){
        Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
          computeFaceMeasure(side_weighted_measure, side_jacobian, tmp_side_cub_weights,
                             subcell_index,cell_topology,
                             scratch_for_compute_side_measure.get_view());
      }

      // Calculate normals
      auto side_normals = Kokkos::DynRankView<Scalar,PHX::Device>("side_normals",num_cells,num_points_on_face,cell_dim);
      if(cell_dim == 1){

        int other_subcell_index = (subcell_index==0) ? 1 : 0;

        for(int cell=0;cell<num_cells;++cell){
          Scalar norm = (in_node_coordinates(cell,subcell_index,0) - in_node_coordinates(cell,other_subcell_index,0));
          side_normals(cell,0,0) = norm / fabs(norm+std::numeric_limits<typename Sacado::ScalarType<Scalar>::type>::min());
        }

      } else {

        cell_tools.getPhysicalSideNormals(side_normals,side_jacobian,subcell_index,cell_topology);

        // Normalize each normal
        for(int cell=0;cell<num_cells;++cell){
          for(int point=0;point<num_points_on_face;++point){
            Scalar n = 0.;
            for(int dim=0;dim<cell_dim;++dim){
              n += side_normals(cell,point,dim)*side_normals(cell,point,dim);
            }
            // If n is zero then this is - hopefully - a virtual cell
            if(n > 0.){
              n = std::sqrt(n);
              for(int dim=0;dim<cell_dim;++dim){
                side_normals(cell,point,dim) /= n;
              }
            }
          }
        }

      }

      // Now that we have all these wonderful values, lets copy them to the actual arrays
      for(int cell=0;cell<num_cells;++cell){
        for(int side_point=0; side_point<num_points_on_face;++side_point){
          const int cell_point = point_offset + side_point;

          weighted_measure(cell,cell_point) = side_weighted_measure(cell,side_point);
          jac_det(cell,cell_point) = side_det_jacobian(cell,side_point);
          for(int dim=0;dim<cell_dim;++dim){
            ref_ip_coordinates(cell,cell_point,dim) = cub_points(cell_point,dim);
            ip_coordinates(cell,cell_point,dim)     = side_ip_coordinates(cell,side_point,dim);
            surface_normals(cell,cell_point,dim)    = side_normals(cell,side_point,dim);

            for(int dim2=0;dim2<cell_dim;++dim2){
              jac(cell,cell_point,dim,dim2) = side_jacobian(cell,side_point,dim,dim2);
              jac_inv(cell,cell_point,dim,dim2) = side_inverse_jacobian(cell,side_point,dim,dim2);
            }
          }
        }
      }
      point_offset += num_points_on_face;
    }
  }

  // We also need surface rotation matrices in order to enforce alignment
  {
    const int num_points = ir.getPointOffset(num_subcells);
    Scalar normal[3];
    Scalar transverse[3];
    Scalar binormal[3];
    for(int i=0;i<3;i++){normal[i]=0.;}
    for(int i=0;i<3;i++){transverse[i]=0.;}
    for(int i=0;i<3;i++){binormal[i]=0.;}
    for(int cell=0; cell<num_cells; ++cell){
      for(int subcell_index=0; subcell_index<num_subcells; ++subcell_index){
        for(int point=0; point<num_points; ++point){
          
          for(int dim=0; dim<3; ++dim)
            normal[dim] = 0.0;
          
          for(int dim=0; dim<cell_dim; ++dim){
            normal[dim] = surface_normals(cell,point,dim);
          }
          
          convertNormalToRotationMatrix(normal,transverse,binormal);
          
          for(int dim=0; dim<3; ++dim){
            surface_rotation_matrices(cell,point,0,dim) = normal[dim];
            surface_rotation_matrices(cell,point,1,dim) = transverse[dim];
            surface_rotation_matrices(cell,point,2,dim) = binormal[dim];
          }
        }
      }
    }
  }

  // =========================================================
  // Enforce alignment across surface quadrature points
  
  const int num_points = ip_coordinates.extent_int(1);
  const int num_faces_per_cell = face_connectivity.numSubcellsOnCell(0);
  const int num_points_per_face = num_points / num_faces_per_cell;

  // Now we need to align the cubature points for each face
  // If there is only one point there is no need to align things
  if(num_points_per_face != 1){
    // To enforce that quadrature points on faces are aligned properly we will iterate through faces,
    // map points to a plane shared by the faces, then re-order quadrature points on the "1" face to
    // line up with the "0" face

    // Utility calls
#define PANZER_DOT(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
#define PANZER_CROSS(a,b,c) {c[0] = a[1]*b[2] - a[2]*b[1]; c[1] = a[2]*b[0] - a[0]*b[2]; c[2] = a[0]*b[1] - a[1]*b[0];}

    // Reorder index vector
    std::vector<int> point_order(num_points_per_face);

    // Iterate through faces
    for(int face=0; face<face_connectivity.numSubcells(); ++face){
      // Cells for sides 0 and 1
      const int cell_0 = face_connectivity.cellForSubcell(face,0);
      const int cell_1 = face_connectivity.cellForSubcell(face,1);

      // If this face doesn't connect to anything we don't need to worry about alignment
      if(cell_1 < 0)
        continue;

      // Local face index for sides 0 and 1
      const int lidx_0 = face_connectivity.localSubcellForSubcell(face,0);
      const int lidx_1 = face_connectivity.localSubcellForSubcell(face,1);

      // If the cell exists, then the local face index needs to exist
      TEUCHOS_ASSERT(lidx_1 >= 0);
      
      // To compare points on planes, it is good to center the planes around the origin
      // Find the face center for the left and right cell (may not be the same point - e.g. periodic conditions)
      // We also define a length scale r2 which gives an order of magnitude approximation to the cell size squared
      Scalar xc0[3] = {0}, xc1[3] = {0}, r2 = 0;
      for(int face_point=0; face_point<num_points_per_face; ++face_point){
        Scalar dx2 = 0.;
        for(int dim=0; dim<cell_dim; ++dim){
          xc0[dim] += ip_coordinates(cell_0,lidx_0*num_points_per_face + face_point,dim);
          xc1[dim] += ip_coordinates(cell_1,lidx_1*num_points_per_face + face_point,dim);
          const Scalar dx = ip_coordinates(cell_0,lidx_0*num_points_per_face + face_point,dim) - ip_coordinates(cell_0,lidx_0*num_points_per_face,dim); 
          dx2 += dx*dx;
        }

        // Check if the distance squared between these two points is larger than the others (doesn't need to be perfect)
        r2 = std::max(r2, dx2);
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
        if(std::fabs(n0_dot_n1 - 1.) < 1.e-8)
          continue;

        // Rotated faces case (e.g. periodic wedge) we need to check for non-antiparallel face normals
        if(std::fabs(n0_dot_n1 + 1.) > 1.e-2){

          // Now we need to define an arbitrary transverse and binormal in the plane across which the faces are anti-symmetric
          // We can do this by setting t = n_0 \times n_1
          PANZER_CROSS(n0,n1,t);

          // Normalize the transverse vector
          const Scalar mag_t = std::sqrt(PANZER_DOT(t,t));
          t[0] /= mag_t;
          t[1] /= mag_t;
          t[2] /= mag_t;

          //  The binormal will be the sum of the normals (does not need to be a right handed system)
          b[0] = n0[0] + n1[0];
          b[1] = n0[1] + n1[1];
          b[2] = n0[2] + n1[2];

          // Normalize the binormal vector
          const Scalar mag_b = std::sqrt(PANZER_DOT(b,b));
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
          x1[dim] = ip_coordinates(cell_1,point_1,dim) - xc1[dim];

        // Convert position to transverse/binormal plane
        p1[0] = PANZER_DOT(x1,t);
        p1[1] = PANZER_DOT(x1,b);

        // Set the default for the point order
        point_order[face_point_1] = face_point_1;
        
        // Compare to points on the other surface
        for(int face_point_0=0; face_point_0<num_points_per_face; ++face_point_0){

          // Get the point index in the 0 cell
          const int point_0 = lidx_0*num_points_per_face + face_point_0;

          // Load position shifted by face center
          for(int dim=0; dim<cell_dim; ++dim)
            x0[dim] = ip_coordinates(cell_0,point_0,dim) - xc0[dim];

          // Convert position to transverse/binormal plane
          p0[0] = PANZER_DOT(x0,t);
          p0[1] = PANZER_DOT(x0,b);

          // Find the distance squared between p0 and p1
          const Scalar p012 = (p0[0]-p1[0])*(p0[0]-p1[0]) + (p0[1]-p1[1])*(p0[1]-p1[1]);

          // If the distance, compared to the size of the cell, is small, we assume these are the same points
          if(p012 / r2 < 1.e-12){
            point_order[face_point_1] = face_point_0;
            break;
          }

          // Big problem - there wan't a point that aligned properly
          TEUCHOS_ASSERT(face_point_0 != num_points_per_face-1);

        }
      }
      
      // Now re-order the points on face 1 to correct the alignment
      const int point_offset = lidx_1*num_points_per_face;
      for( int face_point_1 = 0; face_point_1 < num_points_per_face - 1; ++face_point_1 ){
        // While the point is not yet in place - keep swapping until it is in place (N^2 operations - not great)
        while( face_point_1 != point_order[face_point_1] ){
          // We need to swap with the component in this position
          const int face_point_0 = point_order[face_point_1];
          swapQuadraturePoints(cell_1,point_offset+face_point_1,point_offset+face_point_0);
          std::swap( point_order[face_point_1], point_order[face_point_0] );
        }
      }
    }
  }
    
#undef PANZER_DOT
#undef PANZER_CROSS

  // =========================================================
    
  // I'm not sure if these should exist for surface integrals, but here we go!

  // Shakib contravarient metric tensor
  for (int cell = 0; cell < num_cells; ++cell) {
    for (size_type ip = 0; ip < contravarient.extent(1); ++ip) {

      // zero out matrix
      for (size_type i = 0; i < contravarient.extent(2); ++i)
        for (size_type j = 0; j < contravarient.extent(3); ++j)
          covarient(cell,ip,i,j) = 0.0;

      // g^{ij} = \frac{\parital x_i}{\partial \chi_\alpha}\frac{\parital x_j}{\partial \chi_\alpha}
      for (size_type i = 0; i < contravarient.extent(2); ++i) {
        for (size_type j = 0; j < contravarient.extent(3); ++j) {
          for (size_type alpha = 0; alpha < contravarient.extent(2); ++alpha) {
            covarient(cell,ip,i,j) += jac(cell,ip,i,alpha) * jac(cell,ip,j,alpha);
          }
        }
      }

    }
  }

  auto s_contravarient = Kokkos::subview(contravarient.get_view(), std::make_pair(0,num_cells),Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
  auto s_covarient = Kokkos::subview(covarient.get_view(), std::make_pair(0,num_cells),Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
  Intrepid2::RealSpaceTools<PHX::Device::execution_space>::inverse(s_contravarient, s_covarient);

  // norm of g_ij
  for (int cell = 0; cell < num_cells; ++cell) {
    for (size_type ip = 0; ip < contravarient.extent(1); ++ip) {
      norm_contravarient(cell,ip) = 0.0;
      for (size_type i = 0; i < contravarient.extent(2); ++i) {
        for (size_type j = 0; j < contravarient.extent(3); ++j) {
          norm_contravarient(cell,ip) += contravarient(cell,ip,i,j) * contravarient(cell,ip,i,j);
        }
      }
      norm_contravarient(cell,ip) = std::sqrt(norm_contravarient(cell,ip));
    }
  }

}


template <typename Scalar>
void IntegrationValues2<Scalar>::
evaluateRemainingValues(const PHX::MDField<Scalar,Cell,NODE,Dim>& in_node_coordinates,
                        const int in_num_cells)
{
  Intrepid2::CellTools<PHX::Device::execution_space> cell_tools;

  // copy the dynamic data structures into the static data structures
  {
    size_type num_ip = dyn_cub_points.extent(0);
    size_type num_dims = dyn_cub_points.extent(1);

    for (size_type ip = 0; ip < num_ip;  ++ip) {
      cub_weights(ip) = dyn_cub_weights(ip);
      for (size_type dim = 0; dim < num_dims; ++dim)
        cub_points(ip,dim) = dyn_cub_points(ip,dim);
    }
  }

  if (int_rule->isSide()) {
    const size_type num_ip = dyn_cub_points.extent(0), num_side_dims = dyn_side_cub_points.extent(1);
    for (size_type ip = 0; ip < num_ip; ++ip)
      for (size_type dim = 0; dim < num_side_dims; ++dim)
        side_cub_points(ip,dim) = dyn_side_cub_points(ip,dim);
  }

  const int num_cells = in_num_cells < 0 ? in_node_coordinates.extent(0) : in_num_cells;

  {
    size_type num_nodes = in_node_coordinates.extent(1);
    size_type num_dims = in_node_coordinates.extent(2);

    for (int cell = 0; cell < num_cells;  ++cell) {
      for (size_type node = 0; node < num_nodes; ++node) {
        for (size_type dim = 0; dim < num_dims; ++dim) {
          node_coordinates(cell,node,dim) =
              in_node_coordinates(cell,node,dim);
        }
      }
    }
  }

  auto s_in_node_coordinates = Kokkos::subview(in_node_coordinates.get_view(),std::make_pair(0,num_cells),Kokkos::ALL(),Kokkos::ALL());
  auto s_jac = Kokkos::subview(jac.get_view(),std::make_pair(0,num_cells),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
  cell_tools.setJacobian(jac.get_view(),
                         cub_points.get_view(),
                         node_coordinates.get_view(),
                         *(int_rule->topology));

  auto s_jac_inv = Kokkos::subview(jac_inv.get_view(),std::make_pair(0,num_cells),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
  cell_tools.setJacobianInv(s_jac_inv, s_jac);

  auto s_jac_det = Kokkos::subview(jac_det.get_view(),std::make_pair(0,num_cells),Kokkos::ALL());
  cell_tools.setJacobianDet(s_jac_det, s_jac);

  auto s_weighted_measure = Kokkos::subview(weighted_measure.get_view(),std::make_pair(0,num_cells),Kokkos::ALL());
  if (!int_rule->isSide()) {
    Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
    computeCellMeasure(s_weighted_measure, s_jac_det, cub_weights.get_view());
  }
  else if(int_rule->spatial_dimension==3) {
    Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
    computeFaceMeasure(s_weighted_measure, s_jac, cub_weights.get_view(),
                       int_rule->side, *int_rule->topology,
                       scratch_for_compute_side_measure.get_view());
  }
  else if(int_rule->spatial_dimension==2) {
    Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
    computeEdgeMeasure(s_weighted_measure, s_jac, cub_weights.get_view(),
                       int_rule->side,*int_rule->topology,
                       scratch_for_compute_side_measure.get_view());
  }
  else TEUCHOS_ASSERT(false);

  // Shakib contravarient metric tensor
  for (int cell = 0; cell < num_cells; ++cell) {
    for (size_type ip = 0; ip < contravarient.extent(1); ++ip) {

      // zero out matrix
      for (size_type i = 0; i < contravarient.extent(2); ++i)
        for (size_type j = 0; j < contravarient.extent(3); ++j)
          covarient(cell,ip,i,j) = 0.0;

      // g^{ij} = \frac{\parital x_i}{\partial \chi_\alpha}\frac{\parital x_j}{\partial \chi_\alpha}
      for (size_type i = 0; i < contravarient.extent(2); ++i) {
        for (size_type j = 0; j < contravarient.extent(3); ++j) {
          for (size_type alpha = 0; alpha < contravarient.extent(2); ++alpha) {
            covarient(cell,ip,i,j) += jac(cell,ip,i,alpha) * jac(cell,ip,j,alpha);
          }
        }
      }

    }
  }

  auto s_covarient = Kokkos::subview(covarient.get_view(),std::make_pair(0,num_cells),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
  auto s_contravarient = Kokkos::subview(contravarient.get_view(),std::make_pair(0,num_cells),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
  Intrepid2::RealSpaceTools<PHX::Device::execution_space>::inverse(s_contravarient, s_covarient);

  // norm of g_ij
  for (int cell = 0; cell < num_cells; ++cell) {
    for (size_type ip = 0; ip < contravarient.extent(1); ++ip) {
      norm_contravarient(cell,ip) = 0.0;
      for (size_type i = 0; i < contravarient.extent(2); ++i) {
        for (size_type j = 0; j < contravarient.extent(3); ++j) {
          norm_contravarient(cell,ip) += contravarient(cell,ip,i,j) * contravarient(cell,ip,i,j);
        }
      }
      norm_contravarient(cell,ip) = std::sqrt(norm_contravarient(cell,ip));
    }
  }
}

// Find the permutation that maps the set of points coords to other_coords. To
// avoid possible finite precision issues, == is not used, but rather
// min(norm(.)).
template <typename Scalar>
static void
permuteToOther(const PHX::MDField<Scalar,Cell,IP,Dim>& coords,
               const PHX::MDField<Scalar,Cell,IP,Dim>& other_coords,
               std::vector<typename ArrayTraits<Scalar,PHX::MDField<Scalar> >::size_type>& permutation)
{
  typedef typename ArrayTraits<Scalar,PHX::MDField<Scalar> >::size_type size_type;
  // We can safely assume: (1) The permutation is the same for every cell in
  // the workset. (2) The first workset has valid data. Hence we operate only
  // on cell 0.
  const size_type cell = 0;
  const size_type num_ip = coords.extent(1), num_dim = coords.extent(2);
  permutation.resize(num_ip);
  std::vector<char> taken(num_ip, 0);
  for (size_type ip = 0; ip < num_ip; ++ip) {
    // Find an other point to associate with ip.
    size_type i_min = 0;
    Scalar d_min = -1;
    for (size_type other_ip = 0; other_ip < num_ip; ++other_ip) {
      // For speed, skip other points that are already associated.
      if (taken[other_ip]) continue;
      // Compute the distance between the two points.
      Scalar d(0);
      for (size_type dim = 0; dim < num_dim; ++dim) {
        const Scalar diff = coords(cell, ip, dim) - other_coords(cell, other_ip, dim);
        d += diff*diff;
      }
      if (d_min < 0 || d < d_min) {
        d_min = d;
        i_min = other_ip;
      }
    }
    // Record the match.
    permutation[ip] = i_min;
    // This point is now taken.
    taken[i_min] = 1;
  }
}

template <typename Scalar>
void IntegrationValues2<Scalar>::
evaluateValues(const PHX::MDField<Scalar,Cell,NODE,Dim>& in_node_coordinates,
               const PHX::MDField<Scalar,Cell,IP,Dim>& other_ip_coordinates,
               const int in_num_cells)
{
  const int num_cells = in_num_cells < 0 ? in_node_coordinates.extent(0) : in_num_cells;

  if (int_rule->cv_type == "none") {

    getCubature(in_node_coordinates, in_num_cells);

    {
      // Determine the permutation.
      std::vector<size_type> permutation(other_ip_coordinates.extent(1));
      permuteToOther(ip_coordinates, other_ip_coordinates, permutation);
      // Apply the permutation to the cubature arrays.
      MDFieldArrayFactory af(prefix, alloc_arrays);
      {
        const size_type num_ip = dyn_cub_points.extent(0);
        {
          const size_type num_dim = dyn_side_cub_points.extent(1);
          DblArrayDynamic old_dyn_side_cub_points = af.template buildArray<double,IP,Dim>(
            "old_dyn_side_cub_points", num_ip, num_dim);
          old_dyn_side_cub_points.deep_copy(dyn_side_cub_points);
          for (size_type ip = 0; ip < num_ip; ++ip)
            if (ip != permutation[ip])
              for (size_type dim = 0; dim < num_dim; ++dim)
                dyn_side_cub_points(ip, dim) = old_dyn_side_cub_points(permutation[ip], dim);
        }
        {
          const size_type num_dim = dyn_cub_points.extent(1);
          DblArrayDynamic old_dyn_cub_points = af.template buildArray<double,IP,Dim>(
            "old_dyn_cub_points", num_ip, num_dim);
          old_dyn_cub_points.deep_copy(dyn_cub_points);
          for (size_type ip = 0; ip < num_ip; ++ip)
            if (ip != permutation[ip])
              for (size_type dim = 0; dim < num_dim; ++dim)
                dyn_cub_points(ip, dim) = old_dyn_cub_points(permutation[ip], dim);
        }
        {
          DblArrayDynamic old_dyn_cub_weights = af.template buildArray<double,IP>(
            "old_dyn_cub_weights", num_ip);
          old_dyn_cub_weights.deep_copy(dyn_cub_weights);
          for (size_type ip = 0; ip < dyn_cub_weights.extent(0); ++ip)
            if (ip != permutation[ip])
              dyn_cub_weights(ip) = old_dyn_cub_weights(permutation[ip]);
        }
      }
      {
        const size_type num_ip = ip_coordinates.extent(1);
        const size_type num_dim = ip_coordinates.extent(2);
        Array_CellIPDim old_ip_coordinates = af.template buildStaticArray<Scalar,Cell,IP,Dim>(
            "old_ip_coordinates", ip_coordinates.extent(0), num_ip, num_dim);
        Kokkos::deep_copy(old_ip_coordinates.get_static_view(), ip_coordinates.get_static_view());
        for (int cell = 0; cell < num_cells; ++cell)
          for (size_type ip = 0; ip < num_ip; ++ip)
            if (ip != permutation[ip])
              for (size_type dim = 0; dim < num_dim; ++dim)
                ip_coordinates(cell, ip, dim) = old_ip_coordinates(cell, permutation[ip], dim);
      }
      // All subsequent calculations inherit the permutation.
    }

    evaluateRemainingValues(in_node_coordinates, in_num_cells);
  }

  else {

    getCubatureCV(in_node_coordinates, in_num_cells);

    // Determine the permutation.
    std::vector<size_type> permutation(other_ip_coordinates.extent(1));
    permuteToOther(ip_coordinates, other_ip_coordinates, permutation);

    // Apply the permutation to the cubature arrays.
    MDFieldArrayFactory af(prefix, alloc_arrays);
    {
      const size_type workset_size = ip_coordinates.extent(0), num_ip = ip_coordinates.extent(1),
          num_dim = ip_coordinates.extent(2);
      Array_CellIPDim old_ip_coordinates = af.template buildStaticArray<Scalar,Cell,IP,Dim>(
          "old_ip_coordinates", workset_size, num_ip, num_dim);
      Kokkos::deep_copy(old_ip_coordinates.get_static_view(), ip_coordinates.get_static_view());
      Array_CellIPDim old_weighted_normals = af.template buildStaticArray<Scalar,Cell,IP,Dim>(
          "old_weighted_normals", workset_size, num_ip, num_dim);
      Array_CellIP old_weighted_measure = af.template buildStaticArray<Scalar,Cell,IP>(
          "old_weighted_measure", workset_size, num_ip);
      if (int_rule->cv_type == "side")
        Kokkos::deep_copy(old_weighted_normals.get_static_view(), weighted_normals.get_static_view());
      else
        Kokkos::deep_copy(old_weighted_measure.get_static_view(), weighted_measure.get_static_view());
      for (int cell = 0; cell < num_cells; ++cell)
      {
        for (size_type ip = 0; ip < num_ip; ++ip)
        {
          if (ip != permutation[ip]) {
            if (int_rule->cv_type == "boundary" || int_rule->cv_type == "volume")
              weighted_measure(cell, ip) = old_weighted_measure(cell, permutation[ip]);
            for (size_type dim = 0; dim < num_dim; ++dim)
            {
              ip_coordinates(cell, ip, dim) = old_ip_coordinates(cell, permutation[ip], dim);
              if (int_rule->cv_type == "side")
                weighted_normals(cell, ip, dim) = old_weighted_normals(cell, permutation[ip], dim);

            }
          }
        }
      }
    }

    evaluateValuesCV(in_node_coordinates, in_num_cells);
  }
}

template <typename Scalar>
void IntegrationValues2<Scalar>::
getCubatureCV(const PHX::MDField<Scalar,Cell,NODE,Dim>& in_node_coordinates,
              const int in_num_cells)
{
  int num_space_dim = int_rule->topology->getDimension();
  if (int_rule->isSide() && num_space_dim==1) {
    std::cout << "WARNING: 0-D quadrature rule infrastructure does not exist!!! Will not be able to do "
        << "non-natural integration rules.";
    return;
  }

  size_type num_cells = in_num_cells < 0 ? in_node_coordinates.extent(0) : (size_type) in_num_cells;
  std::pair<int,int> cell_range(0,num_cells);
  {
    size_type num_nodes = in_node_coordinates.extent(1);
    size_type num_dims = in_node_coordinates.extent(2);

    for (size_type cell = 0; cell < num_cells;  ++cell) {
      for (size_type node = 0; node < num_nodes; ++node) {
        for (size_type dim = 0; dim < num_dims; ++dim) {
          node_coordinates(cell,node,dim) =
              in_node_coordinates(cell,node,dim);
          dyn_node_coordinates(cell,node,dim) =
              Sacado::ScalarValue<Scalar>::eval(in_node_coordinates(cell,node,dim));
        }
      }
    }
  }

  auto s_dyn_phys_cub_points = Kokkos::subdynrankview(dyn_phys_cub_points.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
  auto s_dyn_node_coordinates = Kokkos::subdynrankview(dyn_node_coordinates.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
  if (int_rule->cv_type == "side") {
    auto s_dyn_phys_cub_norms = Kokkos::subdynrankview(dyn_phys_cub_norms.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    intrepid_cubature->getCubature(s_dyn_phys_cub_points,s_dyn_phys_cub_norms,s_dyn_node_coordinates);
  }
  else {
    auto s_dyn_phys_cub_weights = Kokkos::subdynrankview(dyn_phys_cub_weights.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    intrepid_cubature->getCubature(s_dyn_phys_cub_points,s_dyn_phys_cub_weights,s_dyn_node_coordinates);
  }

  // size_type num_cells = dyn_phys_cub_points.extent(0);
  size_type num_ip =dyn_phys_cub_points.extent(1);
  size_type num_dims = dyn_phys_cub_points.extent(2);

  for (size_type cell = 0; cell < num_cells;  ++cell) {
    for (size_type ip = 0; ip < num_ip;  ++ip) {
      if (int_rule->cv_type != "side")
        weighted_measure(cell,ip) = dyn_phys_cub_weights(cell,ip);
      for (size_type dim = 0; dim < num_dims; ++dim) {
        ip_coordinates(cell,ip,dim) = dyn_phys_cub_points(cell,ip,dim);
        if (int_rule->cv_type == "side")
          weighted_normals(cell,ip,dim) = dyn_phys_cub_norms(cell,ip,dim);
      }
    }
  }

}

template <typename Scalar>
void IntegrationValues2<Scalar>::
evaluateValuesCV(const PHX::MDField<Scalar, Cell, NODE, Dim>& in_node_coordinates,
                 const int in_num_cells)
{

  Intrepid2::CellTools<PHX::Device::execution_space> cell_tools;

  size_type num_cells = in_num_cells < 0 ? in_node_coordinates.extent(0) : (size_type) in_num_cells;

  auto s_ref_ip_coordinates = Kokkos::subview(ref_ip_coordinates.get_view(),std::make_pair(0,(int)num_cells),Kokkos::ALL(),Kokkos::ALL());
  auto s_ip_coordinates = Kokkos::subview(ip_coordinates.get_view(),std::make_pair<int,int>(0,num_cells),Kokkos::ALL(),Kokkos::ALL());
  auto s_node_coordinates = Kokkos::subview(node_coordinates.get_view(),std::make_pair<int,int>(0,num_cells),Kokkos::ALL(),Kokkos::ALL());

  cell_tools.mapToReferenceFrame(s_ref_ip_coordinates,
                                 s_ip_coordinates,
                                 s_node_coordinates,
                                 *(int_rule->topology));

  auto s_jac = Kokkos::subview(jac.get_view(),std::make_pair<int,int>(0,num_cells),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());

  cell_tools.setJacobian(s_jac,
                         s_ref_ip_coordinates,
                         s_node_coordinates,
                         *(int_rule->topology));

  auto s_jac_inv = Kokkos::subview(jac_inv.get_view(),std::make_pair<int,int>(0,num_cells),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());

  cell_tools.setJacobianInv(s_jac_inv, s_jac);

  auto s_jac_det = Kokkos::subview(jac_det.get_view(),std::make_pair<int,int>(0,num_cells),Kokkos::ALL());

  cell_tools.setJacobianDet(s_jac_det, s_jac);
}

#define INTEGRATION_VALUES2_INSTANTIATION(SCALAR) \
    template class IntegrationValues2<SCALAR>;

INTEGRATION_VALUES2_INSTANTIATION(panzer::Traits::RealType)
INTEGRATION_VALUES2_INSTANTIATION(panzer::Traits::FadType)

}

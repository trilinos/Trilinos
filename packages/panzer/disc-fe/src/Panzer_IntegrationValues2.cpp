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
evaluateValues(const PHX::MDField<Scalar,Cell,NODE,Dim> & in_node_coordinates)
{
  typedef panzer::IntegrationDescriptor ID;
  const bool is_surface = int_rule->getType() == ID::SURFACE;
  const bool is_cv = (int_rule->getType() == ID::CV_VOLUME) or (int_rule->getType() == ID::CV_SIDE) or (int_rule->getType() == ID::CV_BOUNDARY);

  TEUCHOS_ASSERT(not (is_surface and is_cv));

  if(is_surface){
    generateSurfaceCubatureValues(in_node_coordinates);
  } else if (is_cv) {
    getCubatureCV(in_node_coordinates);
    evaluateValuesCV(in_node_coordinates);
  } else {
    getCubature(in_node_coordinates);
    evaluateRemainingValues(in_node_coordinates);
  }
}

template <typename Scalar>
void IntegrationValues2<Scalar>::
getCubature(const PHX::MDField<Scalar,Cell,NODE,Dim>& in_node_coordinates)
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
  cell_tools.mapToPhysicalFrame(ip_coordinates.get_view(),
                                dyn_cub_points.get_view(),
                                in_node_coordinates.get_view(),
                                *(int_rule->topology));
}





namespace
{

template <typename array_t, typename scalar_t>
class point_sorter_t
{
public:

  point_sorter_t() = delete;
  point_sorter_t(const array_t & array, const int cell, const int offset):
    _array(array),
    _cell(cell),
    _offset(offset),
    _rel_tol(1.e-12)
  {
    _num_dims=_array.extent(2);
  }


  // This needs to be optimized
  bool operator()(const int & point_a, const int & point_b) const
  {

    if(_num_dims==1){

      const scalar_t & x_a = _array(_cell,_offset+point_a,0);
      const scalar_t & x_b = _array(_cell,_offset+point_b,0);

      const scalar_t rel = std::max(std::fabs(x_a),std::fabs(x_b));

      return test_less(x_a,x_b,rel);

    } else if(_num_dims==2){

      const scalar_t & x_a = _array(_cell,_offset+point_a,0);
      const scalar_t & x_b = _array(_cell,_offset+point_b,0);

      const scalar_t & y_a = _array(_cell,_offset+point_a,1);
      const scalar_t & y_b = _array(_cell,_offset+point_b,1);

      const scalar_t rel_x = std::max(std::fabs(x_a),std::fabs(x_b));
      const scalar_t rel_y = std::max(std::fabs(y_a),std::fabs(y_b));
      const scalar_t rel = std::max(rel_x,rel_y);

      if(test_eq(y_a,y_b,rel)){
        if(test_less(x_a,x_b,rel)){
          // Sort by x
          return true;
        }
      } else if(test_less(y_a,y_b,rel)){
        // Sort by y
        return true;
      }

      // Otherwise b < a
      return false;

    } else if(_num_dims==3){

      const scalar_t & x_a = _array(_cell,_offset+point_a,0);
      const scalar_t & x_b = _array(_cell,_offset+point_b,0);

      const scalar_t & y_a = _array(_cell,_offset+point_a,1);
      const scalar_t & y_b = _array(_cell,_offset+point_b,1);

      const scalar_t & z_a = _array(_cell,_offset+point_a,2);
      const scalar_t & z_b = _array(_cell,_offset+point_b,2);

      const scalar_t rel_x = std::max(std::fabs(x_a),std::fabs(x_b));
      const scalar_t rel_y = std::max(std::fabs(y_a),std::fabs(y_b));
      const scalar_t rel_z = std::max(std::fabs(z_a),std::fabs(z_b));
      const scalar_t rel = std::max(rel_x,std::max(rel_y,rel_z));

      if(test_less(z_a,z_b,rel)){
        // Sort by z
        return true;
      } else if(test_eq(z_a,z_b,rel)){
        if(test_eq(y_a,y_b,rel)){
          if(test_less(x_a,x_b,rel)){
            // Sort by x
            return true;
          }
        } else if(test_less(y_a,y_b,rel)){
          // Sort by y
          return true;
        }
      }
      // Otherwise b < a
      return false;

    } else {
      TEUCHOS_ASSERT(false);
    }
  }

protected:

  bool
  test_eq(const scalar_t & a, const scalar_t & b, const scalar_t & rel) const
  {
    if(rel==0){
      return true;
    }
    return std::fabs(a-b) < _rel_tol * rel;
  }

  bool
  test_less(const scalar_t & a, const scalar_t & b, const scalar_t & rel) const
  {
    if(rel==0){
      return false;
    }
    return (a-b) < -_rel_tol * rel;
  }

  const array_t & _array;
  int _cell;
  int _offset;
  int _num_dims;
  scalar_t _rel_tol;

};

template<typename T>
void
convertNormalToRotationMatrix(const T normal[3], T transverse[3], T binormal[3])
{

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
}

template <typename Scalar>
void IntegrationValues2<Scalar>::
uniqueCoordOrdering(Array_CellIPDim & coords,
                    int cell,
                    int offset,
                    std::vector<int> & order)
{
  for(size_t point_index=0;point_index<order.size();++point_index){
    order[point_index] = point_index;
  }

  // We need to sort the indexes in point_indexes by their ip_coordinate's position in space.
  // We will then use that to sort all of our arrays.

  point_sorter_t<Array_CellIPDim,Scalar> sorter(coords,cell,offset);
  std::sort(order.begin(),order.end(),sorter);
}

template <typename Scalar>
void IntegrationValues2<Scalar>::
generateSurfaceCubatureValues(const PHX::MDField<Scalar,Cell,NODE,Dim>& in_node_coordinates)
{

  TEUCHOS_ASSERT(int_rule->getType() == IntegrationDescriptor::SURFACE);

  Intrepid2::CellTools<PHX::Device::execution_space> cell_tools;

  const shards::CellTopology & cell_topology = *(int_rule->topology);
  const panzer::IntegrationRule & ir = *int_rule;

  // Copy over coordinates
  {
    const int num_cells = in_node_coordinates.extent(0);
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

  const int num_cells = in_node_coordinates.extent(0);
  const int cell_dim = cell_topology.getDimension();
  const int subcell_dim = cell_topology.getDimension()-1;
  const int num_subcells = cell_topology.getSubcellCount(subcell_dim);

  Intrepid2::DefaultCubatureFactory cubature_factory;

  // We get to build up our cubature one face at a time
  int point_offset=0;
  for(int subcell_index=0; subcell_index<num_subcells; ++subcell_index){

    // Default for 1D
    int num_points_on_face = 1;

    // Get the cubature for the side
    Kokkos::DynRankView<double,PHX::Device> side_cub_weights;
    Kokkos::DynRankView<double,PHX::Device> side_cub_points;
    if(cell_dim==1){
      side_cub_weights = Kokkos::DynRankView<double,PHX::Device>("side_cub_weights",num_points_on_face);
      side_cub_points = Kokkos::DynRankView<double,PHX::Device>("cell_side_cub_points",num_points_on_face,cell_dim);
      side_cub_weights(0)=1.;
      side_cub_points(0,0) = (subcell_index==0)? -1. : 1.;
    } else {

      // Get the face topology from the cell topology
      const shards::CellTopology face_topology(cell_topology.getCellTopologyData(subcell_dim,subcell_index));

      auto ic = cubature_factory.create<PHX::Device::execution_space,double,double>(face_topology,ir.getOrder());
      num_points_on_face = ic->getNumPoints();

      side_cub_weights = Kokkos::DynRankView<double,PHX::Device>("side_cub_weights",num_points_on_face);
      side_cub_points = Kokkos::DynRankView<double,PHX::Device>("cell_side_cub_points",num_points_on_face,cell_dim);

      auto subcell_cub_points = Kokkos::DynRankView<double,PHX::Device>("side_cub_points",num_points_on_face,subcell_dim);

      // Get the reference face points
      ic->getCubature(subcell_cub_points, side_cub_weights);

      // Convert from reference face points to reference cell points
      cell_tools.mapToReferenceSubcell(side_cub_points,
                                       subcell_cub_points,
                                       subcell_dim,
                                       subcell_index,
                                       cell_topology);
    }


    for(int local_point=0;local_point<num_points_on_face;++local_point){
      const int point = point_offset + local_point;
      for(int dim=0;dim<cell_dim;++dim){
        cub_points(point,dim) = side_cub_points(local_point,dim);
      }
    }


    // Map from side points to physical points
    auto side_ip_coordinates = Kokkos::DynRankView<Scalar,PHX::Device>("side_ip_coordinates",num_cells,num_points_on_face,cell_dim);
    cell_tools.mapToPhysicalFrame(side_ip_coordinates,
                                  side_cub_points,
                                  node_coordinates.get_view(),
                                  cell_topology);

    // Create a jacobian and his friends for this side
    auto side_jacobian = Kokkos::DynRankView<Scalar,PHX::Device>("side_jac",num_cells,num_points_on_face,cell_dim,cell_dim);
    cell_tools.setJacobian(side_jacobian,
                           side_cub_points,
                           node_coordinates.get_view(),
                           cell_topology);

    auto side_inverse_jacobian = Kokkos::DynRankView<Scalar,PHX::Device>("side_inv_jac",num_cells,num_points_on_face,cell_dim,cell_dim);
    cell_tools.setJacobianInv(side_inverse_jacobian, side_jacobian);

    auto side_det_jacobian = Kokkos::DynRankView<Scalar,PHX::Device>("side_det_jac",num_cells,num_points_on_face);
    cell_tools.setJacobianDet(side_det_jacobian, side_jacobian);

    // Calculate measures (quadrature weights in physical space) for this side
    auto side_weighted_measure = Kokkos::DynRankView<Scalar,PHX::Device>("side_weighted_measure",num_cells,num_points_on_face);
    if(cell_dim == 1){
      Kokkos::deep_copy(side_weighted_measure, side_cub_weights(0));
    } else if(cell_dim == 2){
      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
              computeEdgeMeasure(side_weighted_measure, side_jacobian, side_cub_weights,
                                 subcell_index,cell_topology,
                                 scratch_for_compute_side_measure.get_view());
    } else if(cell_dim == 3){
      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
              computeFaceMeasure(side_weighted_measure, side_jacobian, side_cub_weights,
                                 subcell_index,cell_topology,
                                 scratch_for_compute_side_measure.get_view());
    }

    // Calculate normals
    auto side_normals = Kokkos::DynRankView<Scalar,PHX::Device>("side_normals",num_cells,num_points_on_face,cell_dim);
    if(cell_dim == 1){

      int other_subcell_index = (subcell_index==0) ? 1 : 0;

      for(int cell=0;cell<num_cells;++cell){
        Scalar norm = (in_node_coordinates(cell,subcell_index,0) - in_node_coordinates(cell,other_subcell_index,0));
        side_normals(cell,0,0) = norm / fabs(norm);
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

  // Now we need to sort the cubature points for each face so that they will line up between cells
  {
    for(int subcell_index=0; subcell_index<num_subcells;++subcell_index){

      const int point_offset = ir.getPointOffset(subcell_index);
      const int num_points_on_face = ir.getPointOffset(subcell_index+1) - point_offset;
      std::vector<int> point_indexes(num_points_on_face,-1);

      for(int cell=0; cell<num_cells; ++cell){

        // build a  point index array based on point coordinates
        uniqueCoordOrdering(ip_coordinates,cell,point_offset,point_indexes);

        // Indexes are now sorted, now we swap everything around
        reorder(point_indexes,[=](int a,int b) { swapQuadraturePoints(cell,point_offset+a,point_offset+b); });
      }
    }
  }

  // We also need surface rotation matrices
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
        for(int dim=0; dim<cell_dim; ++dim){
          normal[dim] = surface_normals(cell,point,dim);
        }

        convertNormalToRotationMatrix<Scalar>(normal,transverse,binormal);

        for(int dim=0; dim<3; ++dim){
          surface_rotation_matrices(cell,point,0,dim) = normal[dim];
          surface_rotation_matrices(cell,point,1,dim) = transverse[dim];
          surface_rotation_matrices(cell,point,2,dim) = binormal[dim];
        }
      }
    }
  }


  // I'm not sure if these should exist for surface integrals, but here we go!

  // Shakib contravarient metric tensor
  for (size_type cell = 0; cell < contravarient.extent(0); ++cell) {
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

  Intrepid2::RealSpaceTools<PHX::Device::execution_space>::inverse(contravarient.get_view(), covarient.get_view());

  // norm of g_ij
  for (size_type cell = 0; cell < contravarient.extent(0); ++cell) {
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
evaluateRemainingValues(const PHX::MDField<Scalar,Cell,NODE,Dim>& in_node_coordinates)
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

  {
    size_type num_cells = in_node_coordinates.extent(0);
    size_type num_nodes = in_node_coordinates.extent(1);
    size_type num_dims = in_node_coordinates.extent(2);

    for (size_type cell = 0; cell < num_cells;  ++cell) {
      for (size_type node = 0; node < num_nodes; ++node) {
        for (size_type dim = 0; dim < num_dims; ++dim) {
          node_coordinates(cell,node,dim) =
              in_node_coordinates(cell,node,dim);
        }
      }
    }
  }

  cell_tools.setJacobian(jac.get_view(),
                         cub_points.get_view(),
                         node_coordinates.get_view(),
                         *(int_rule->topology));

  cell_tools.setJacobianInv(jac_inv.get_view(), jac.get_view());

  cell_tools.setJacobianDet(jac_det.get_view(), jac.get_view());

  if (!int_rule->isSide()) {
    Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
    computeCellMeasure(weighted_measure.get_view(), jac_det.get_view(), cub_weights.get_view());
  }
  else if(int_rule->spatial_dimension==3) {
    Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
    computeFaceMeasure(weighted_measure.get_view(), jac.get_view(), cub_weights.get_view(),
                       int_rule->side, *int_rule->topology,
                       scratch_for_compute_side_measure.get_view());
  }
  else if(int_rule->spatial_dimension==2) {
    Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
    computeEdgeMeasure(weighted_measure.get_view(), jac.get_view(), cub_weights.get_view(),
                       int_rule->side,*int_rule->topology,
                       scratch_for_compute_side_measure.get_view());
  }
  else TEUCHOS_ASSERT(false);

  // Shakib contravarient metric tensor
  for (size_type cell = 0; cell < contravarient.extent(0); ++cell) {
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

  Intrepid2::RealSpaceTools<PHX::Device::execution_space>::inverse(contravarient.get_view(), covarient.get_view());

  // norm of g_ij
  for (size_type cell = 0; cell < contravarient.extent(0); ++cell) {
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
               const PHX::MDField<Scalar,Cell,IP,Dim>& other_ip_coordinates)
               {
  if (int_rule->cv_type == "none") {

    getCubature(in_node_coordinates);

    {
      // Determine the permutation.
      std::vector<size_type> permutation(other_ip_coordinates.extent(1));
      permuteToOther(ip_coordinates, other_ip_coordinates, permutation);
      // Apply the permutation to the cubature arrays.
      MDFieldArrayFactory af(prefix, alloc_arrays);
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
      {
        const size_type num_cells = ip_coordinates.extent(0), num_ip = ip_coordinates.extent(1),
            num_dim = ip_coordinates.extent(2);
        Array_CellIPDim old_ip_coordinates = af.template buildStaticArray<Scalar,Cell,IP,Dim>(
            "old_ip_coordinates", num_cells, num_ip, num_dim);
        Kokkos::deep_copy(old_ip_coordinates.get_static_view(), ip_coordinates.get_static_view());
        for (size_type cell = 0; cell < num_cells; ++cell)
          for (size_type ip = 0; ip < num_ip; ++ip)
            if (ip != permutation[ip])
              for (size_type dim = 0; dim < num_dim; ++dim)
                ip_coordinates(cell, ip, dim) = old_ip_coordinates(cell, permutation[ip], dim);
      }
      // All subsequent calculations inherit the permutation.
    }

    evaluateRemainingValues(in_node_coordinates);
  }

  else {

    getCubatureCV(in_node_coordinates);

    // Determine the permutation.
    std::vector<size_type> permutation(other_ip_coordinates.extent(1));
    permuteToOther(ip_coordinates, other_ip_coordinates, permutation);

    // Apply the permutation to the cubature arrays.
    MDFieldArrayFactory af(prefix, alloc_arrays);
    {
      const size_type num_cells = ip_coordinates.extent(0), num_ip = ip_coordinates.extent(1),
          num_dim = ip_coordinates.extent(2);
      Array_CellIPDim old_ip_coordinates = af.template buildStaticArray<Scalar,Cell,IP,Dim>(
          "old_ip_coordinates", num_cells, num_ip, num_dim);
      Kokkos::deep_copy(old_ip_coordinates.get_static_view(), ip_coordinates.get_static_view());
      Array_CellIPDim old_weighted_normals = af.template buildStaticArray<Scalar,Cell,IP,Dim>(
          "old_weighted_normals", num_cells, num_ip, num_dim);
      Array_CellIP old_weighted_measure = af.template buildStaticArray<Scalar,Cell,IP>(
          "old_weighted_measure", num_cells, num_ip);
      if (int_rule->cv_type == "side")
        Kokkos::deep_copy(old_weighted_normals.get_static_view(), weighted_normals.get_static_view());
      else
        Kokkos::deep_copy(old_weighted_measure.get_static_view(), weighted_measure.get_static_view());
      for (size_type cell = 0; cell < num_cells; ++cell)
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

    evaluateValuesCV(in_node_coordinates);
  }
               }

template <typename Scalar>
void IntegrationValues2<Scalar>::
getCubatureCV(const PHX::MDField<Scalar,Cell,NODE,Dim>& in_node_coordinates)
{
  int num_space_dim = int_rule->topology->getDimension();
  if (int_rule->isSide() && num_space_dim==1) {
    std::cout << "WARNING: 0-D quadrature rule infrastructure does not exist!!! Will not be able to do "
        << "non-natural integration rules.";
    return;
  }
  {
    size_type num_cells = in_node_coordinates.extent(0);
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

  if (int_rule->cv_type == "side")
    intrepid_cubature->getCubature(dyn_phys_cub_points.get_view(),dyn_phys_cub_norms.get_view(),dyn_node_coordinates.get_view());
  else
    intrepid_cubature->getCubature(dyn_phys_cub_points.get_view(),dyn_phys_cub_weights.get_view(),dyn_node_coordinates.get_view());

  size_type num_cells = dyn_phys_cub_points.extent(0);
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
evaluateValuesCV(const PHX::MDField<Scalar, Cell, NODE, Dim>& /* in_node_coordinates */)
{

  Intrepid2::CellTools<PHX::Device::execution_space> cell_tools;

  cell_tools.mapToReferenceFrame(ref_ip_coordinates.get_view(),
                                 ip_coordinates.get_view(),
                                 node_coordinates.get_view(),
                                 *(int_rule->topology));

  cell_tools.setJacobian(jac.get_view(),
                         ref_ip_coordinates.get_view(),
                         node_coordinates.get_view(),
                         *(int_rule->topology));

  cell_tools.setJacobianInv(jac_inv.get_view(), jac.get_view());

  cell_tools.setJacobianDet(jac_det.get_view(), jac.get_view());

}

#define INTEGRATION_VALUES2_INSTANTIATION(SCALAR) \
    template class IntegrationValues2<SCALAR>;

INTEGRATION_VALUES2_INSTANTIATION(panzer::Traits::RealType)
INTEGRATION_VALUES2_INSTANTIATION(panzer::Traits::FadType)

}

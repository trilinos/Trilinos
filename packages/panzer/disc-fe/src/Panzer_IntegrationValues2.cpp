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

#include "Shards_CellTopology.hpp"

#include "Kokkos_DynRankView.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_RealSpaceTools.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_ArrayTools.hpp"
#include "Intrepid2_CubatureControlVolume.hpp"
#include "Intrepid2_CubatureControlVolumeSide.hpp"

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

    if (ir->isSide()) {
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
    
    dyn_node_coordinates = af.template buildArray<double,Cell,IP,Dim>("node_coordinates",num_cells, num_ip, num_space_dim);

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

  }
  
  template <typename Scalar>
  void IntegrationValues2<Scalar>::
  setupArrays(const Teuchos::RCP<const panzer::IntegrationRule>& ir)
  {
    MDFieldArrayFactory af(prefix,alloc_arrays);

    int_rule = ir;
    
    int num_nodes = ir->topology->getNodeCount();
    int num_cells = ir->workset_size;
    int num_space_dim = ir->topology->getDimension();

    // specialize content if this is quadrature at anode
    if(num_space_dim==1 && ir->isSide()) {
       setupArraysForNodeRule(ir); 
       return;
    }

    Intrepid2::DefaultCubatureFactory<double,DblArrayDynamic>
      cubature_factory;
    
    if (ir->cv_type == "side")
       intrepid_cubature = Teuchos::rcp(new Intrepid2::CubatureControlVolumeSide<double,DblArrayDynamic,DblArrayDynamic>(ir->topology));

    else if (ir->cv_type == "volume")
       intrepid_cubature = Teuchos::rcp(new Intrepid2::CubatureControlVolume<double,DblArrayDynamic,DblArrayDynamic>(ir->topology));

    else if (ir->cv_type == "none" && ir->isSide())
       intrepid_cubature = cubature_factory.create(*(ir->side_topology),
                                                  ir->cubature_degree);
    else
       intrepid_cubature = cubature_factory.create(*(ir->topology),
                                                   ir->cubature_degree);
    

    int num_ip = intrepid_cubature->getNumPoints();

    dyn_cub_points = af.template buildArray<double,IP,Dim>("cub_points",num_ip, num_space_dim);
    dyn_cub_weights = af.template buildArray<double,IP>("cub_weights",num_ip);

    cub_points = af.template buildStaticArray<Scalar,IP,Dim>("cub_points",num_ip, num_space_dim);

    if (ir->isSide()) {
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

    dyn_node_coordinates = af.template buildArray<double,Cell,IP,Dim>("node_coordinates",num_cells, num_ip, num_space_dim);
    
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

  }

// ***********************************************************
// * Evaluation of values - NOT specialized
// ***********************************************************
  template <typename Scalar>
  void IntegrationValues2<Scalar>::
  evaluateValues(const PHX::MDField<Scalar,Cell,NODE,Dim> & in_node_coordinates)
  {
    if (int_rule->cv_type != "none") {
       evaluateValuesCV(in_node_coordinates);
    }
    else {
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
    
    Intrepid2::CellTools<Scalar> cell_tools;
    
    if (!int_rule->isSide())
      intrepid_cubature->getCubature(dyn_cub_points, dyn_cub_weights);
    else {
      intrepid_cubature->getCubature(dyn_side_cub_points, dyn_cub_weights);

      cell_tools.mapToReferenceSubcell(dyn_cub_points, 
				       dyn_side_cub_points,
				       int_rule->spatial_dimension-1,
				       int_rule->side, 
				       *(int_rule->topology));
    }

    // IP coordinates
    cell_tools.mapToPhysicalFrame(ip_coordinates, dyn_cub_points, in_node_coordinates, *(int_rule->topology));
  }


  template <typename Scalar>
  void IntegrationValues2<Scalar>::
  evaluateRemainingValues(const PHX::MDField<Scalar,Cell,NODE,Dim>& in_node_coordinates)
  {
    Intrepid2::CellTools<Scalar> cell_tools;

    // copy the dynamic data structures into the static data structures
    {
      size_type num_ip = dyn_cub_points.dimension(0);
      size_type num_dims = dyn_cub_points.dimension(1);
     
      for (size_type ip = 0; ip < num_ip;  ++ip) {
        cub_weights(ip) = dyn_cub_weights(ip);
        for (size_type dim = 0; dim < num_dims; ++dim)
          cub_points(ip,dim) = dyn_cub_points(ip,dim);
      }
    }

    if (int_rule->isSide()) {
      const size_type num_ip = dyn_cub_points.dimension(0), num_side_dims = dyn_side_cub_points.dimension(1);
      for (size_type ip = 0; ip < num_ip; ++ip)
        for (size_type dim = 0; dim < num_side_dims; ++dim)
          side_cub_points(ip,dim) = dyn_side_cub_points(ip,dim);
    }
    
    {
      size_type num_cells = in_node_coordinates.dimension(0);
      size_type num_nodes = in_node_coordinates.dimension(1);
      size_type num_dims = in_node_coordinates.dimension(2);
     
      for (size_type cell = 0; cell < num_cells;  ++cell) {
	for (size_type node = 0; node < num_nodes; ++node) {
	  for (size_type dim = 0; dim < num_dims; ++dim) {
	    node_coordinates(cell,node,dim) = 
	      in_node_coordinates(cell,node,dim);
	  }
	}
      }
    }

    cell_tools.setJacobian(jac, cub_points, node_coordinates, 
			   *(int_rule->topology));
    
    cell_tools.setJacobianInv(jac_inv, jac);
    
    cell_tools.setJacobianDet(jac_det, jac);
    
    if (!int_rule->isSide()) {
       Intrepid2::FunctionSpaceTools::
         computeCellMeasure<Scalar>(weighted_measure, jac_det, cub_weights);
    }
    else if(int_rule->spatial_dimension==3) {
       Intrepid2::FunctionSpaceTools::
         computeFaceMeasure<Scalar>(weighted_measure, jac, cub_weights,int_rule->side,*int_rule->topology);
    }
    else if(int_rule->spatial_dimension==2) {
       Intrepid2::FunctionSpaceTools::
         computeEdgeMeasure<Scalar>(weighted_measure, jac, cub_weights,int_rule->side,*int_rule->topology);
    }
    else TEUCHOS_ASSERT(false);
    
    // Shakib contravarient metric tensor
    for (size_type cell = 0; cell < contravarient.dimension(0); ++cell) {
      for (size_type ip = 0; ip < contravarient.dimension(1); ++ip) {

	// zero out matrix
	for (size_type i = 0; i < contravarient.dimension(2); ++i)
	  for (size_type j = 0; j < contravarient.dimension(3); ++j)
	    covarient(cell,ip,i,j) = 0.0;
	   
	// g^{ij} = \frac{\parital x_i}{\partial \chi_\alpha}\frac{\parital x_j}{\partial \chi_\alpha}
	for (size_type i = 0; i < contravarient.dimension(2); ++i) {
	  for (size_type j = 0; j < contravarient.dimension(3); ++j) {
	    for (size_type alpha = 0; alpha < contravarient.dimension(2); ++alpha) {
	      covarient(cell,ip,i,j) += jac(cell,ip,i,alpha) * jac(cell,ip,j,alpha);
	    }
	  }
	}	

      }
    }

    Intrepid2::RealSpaceTools<Scalar>::inverse(contravarient, covarient);

    // norm of g_ij
    for (size_type cell = 0; cell < contravarient.dimension(0); ++cell) {
      for (size_type ip = 0; ip < contravarient.dimension(1); ++ip) {
	norm_contravarient(cell,ip) = 0.0;
	for (size_type i = 0; i < contravarient.dimension(2); ++i) {
	  for (size_type j = 0; j < contravarient.dimension(3); ++j) {
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
    const size_type num_ip = coords.dimension(1), num_dim = coords.dimension(2);
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
    getCubature(in_node_coordinates);

    {
      // Determine the permutation.
      std::vector<size_type> permutation(other_ip_coordinates.dimension(1));
      permuteToOther(ip_coordinates, other_ip_coordinates, permutation);
      // Apply the permutation to the cubature arrays.
      MDFieldArrayFactory af(prefix, alloc_arrays);
      const size_type num_ip = dyn_cub_points.dimension(0);
      {
        const size_type num_dim = dyn_side_cub_points.dimension(1);
        DblArrayDynamic old_dyn_side_cub_points = af.template buildArray<double,IP,Dim>(
          "old_dyn_side_cub_points", num_ip, num_dim);
        old_dyn_side_cub_points.deep_copy(dyn_side_cub_points);
        for (size_type ip = 0; ip < num_ip; ++ip)
          if (ip != permutation[ip])
            for (size_type dim = 0; dim < num_dim; ++dim)
              dyn_side_cub_points(ip, dim) = old_dyn_side_cub_points(permutation[ip], dim);
      }
      {
        const size_type num_dim = dyn_cub_points.dimension(1);
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
        for (size_type ip = 0; ip < dyn_cub_weights.dimension(0); ++ip)
          if (ip != permutation[ip])
            dyn_cub_weights(ip) = old_dyn_cub_weights(permutation[ip]);
      }
      {
        const size_type num_cells = ip_coordinates.dimension(0), num_ip = ip_coordinates.dimension(1),
          num_dim = ip_coordinates.dimension(2);
        Array_CellIPDim old_ip_coordinates = af.template buildStaticArray<Scalar,Cell,IP,Dim>(
          "old_ip_coordinates", num_cells, num_ip, num_dim);
        Kokkos::deep_copy(old_ip_coordinates.get_kokkos_view(), ip_coordinates.get_kokkos_view());
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

  template <typename Scalar>
  void IntegrationValues2<Scalar>::
  evaluateValuesCV(const PHX::MDField<Scalar,Cell,NODE,Dim> & in_node_coordinates)
  {
  
      Intrepid2::CellTools<Scalar> cell_tools;

     {
      size_type num_cells = in_node_coordinates.dimension(0);
      size_type num_nodes = in_node_coordinates.dimension(1);
      size_type num_dims = in_node_coordinates.dimension(2);

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

    if (int_rule->cv_type == "volume")
      intrepid_cubature->getCubature(dyn_phys_cub_points,dyn_phys_cub_weights,dyn_node_coordinates);

    else if (int_rule->cv_type == "side")
      intrepid_cubature->getCubature(dyn_phys_cub_points,dyn_phys_cub_norms,dyn_node_coordinates);

    if (int_rule->cv_type == "volume")
    {
        size_type num_cells = dyn_phys_cub_points.dimension(0);
        size_type num_ip = dyn_phys_cub_points.dimension(1);
        size_type num_dims = dyn_phys_cub_points.dimension(2);

        for (size_type cell = 0; cell < num_cells;  ++cell) {
          for (size_type ip = 0; ip < num_ip;  ++ip) {
             weighted_measure(cell,ip) = dyn_phys_cub_weights(cell,ip);
             for (size_type dim = 0; dim < num_dims; ++dim)
               ip_coordinates(cell,ip,dim) = dyn_phys_cub_points(cell,ip,dim);
           }
        }
        cell_tools.mapToReferenceFrame(ref_ip_coordinates, ip_coordinates, node_coordinates,
                                    *(int_rule->topology),-1);

        cell_tools.setJacobian(jac, ref_ip_coordinates, node_coordinates,
                           *(int_rule->topology));

    }
    else if (int_rule->cv_type == "side")
    {
        size_type num_cells = dyn_phys_cub_points.dimension(0);
        size_type num_ip = dyn_phys_cub_points.dimension(1);
        size_type num_dims = dyn_phys_cub_points.dimension(2);

        for (size_type cell = 0; cell < num_cells;  ++cell) {
          for (size_type ip = 0; ip < num_ip;  ++ip) {
             for (size_type dim = 0; dim < num_dims; ++dim) {
               ip_coordinates(cell,ip,dim) = dyn_phys_cub_points(cell,ip,dim);
               weighted_normals(cell,ip,dim) = dyn_phys_cub_norms(cell,ip,dim);
             }
           }
        }

        cell_tools.mapToReferenceFrame(ref_ip_coordinates, ip_coordinates, node_coordinates,
                                       *(int_rule->topology),-1);
        cell_tools.setJacobian(jac, ref_ip_coordinates, node_coordinates,
                               *(int_rule->topology));
     }

     cell_tools.setJacobianInv(jac_inv, jac);
 
     cell_tools.setJacobianDet(jac_det, jac);

}

#define INTEGRATION_VALUES2_INSTANTIATION(SCALAR) \
template class IntegrationValues2<SCALAR>;

INTEGRATION_VALUES2_INSTANTIATION(panzer::Traits::RealType)
INTEGRATION_VALUES2_INSTANTIATION(panzer::Traits::FadType)

}

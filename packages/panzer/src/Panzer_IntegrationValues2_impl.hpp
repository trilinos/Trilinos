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

#ifndef PANZER_INTEGRATION_VALUES2_IMPL_HPP
#define PANZER_INTEGRATION_VALUES2_IMPL_HPP

#include "Shards_CellTopology.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"

// ***********************************************************
// * Specializations of setupArrays() for different array types
// ***********************************************************

namespace panzer {
  
  // * Specialization for Intrepid::FieldContainer<double>
  template <typename Scalar,
            template <typename DataT,
               typename Tag0, typename Tag1, typename Tag2,
               typename Tag3, typename Tag4, typename Tag5,
               typename Tag6, typename Tag7> class Array >
  template <typename ArrayFactory>
  inline
  void IntegrationValues2<Scalar,Array>::
  setupArraysForNodeRule(const Teuchos::RCP<const panzer::IntegrationRule>& ir,const ArrayFactory & af)
  {
    int num_nodes = ir->topology->getNodeCount();
    int num_cells = ir->workset_size;
    int num_space_dim = ir->topology->getDimension();

    int num_ip = 1;

    dyn_cub_points = af.template buildArray<Scalar,IP,Dim>("cub_points",num_ip, num_space_dim);
    dyn_cub_weights = af.template buildArray<Scalar,IP>("cub_weights",num_ip);

    cub_points = af.template buildStaticArray<Scalar,IP,Dim>("cub_points",num_ip, num_space_dim);

    if (ir->isSide())
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
  }
  
  template <typename Scalar,
            template <typename DataT,
               typename Tag0, typename Tag1, typename Tag2,
               typename Tag3, typename Tag4, typename Tag5,
               typename Tag6, typename Tag7> class Array >
  template <typename ArrayFactory>
  inline
  void IntegrationValues2<Scalar,Array>::
  setupArrays(const Teuchos::RCP<const panzer::IntegrationRule>& ir,const ArrayFactory & af)
  {
    int_rule = ir;
    
    int num_nodes = ir->topology->getNodeCount();
    int num_cells = ir->workset_size;
    int num_space_dim = ir->topology->getDimension();

    // specialize content if this is quadrature at anode
    if(num_space_dim==1 && ir->isSide()) {
       setupArraysForNodeRule(ir,af); 
       return;
    }

    Intrepid::DefaultCubatureFactory<Scalar,Array<Scalar,void,void,void,void,void,void,void,void> >
      cubature_factory;
    
    if (ir->isSide())
      intrepid_cubature = cubature_factory.create(*(ir->side_topology), 
						  ir->cubature_degree);
    else
      intrepid_cubature = cubature_factory.create(*(ir->topology), 
						  ir->cubature_degree);

    int num_ip = intrepid_cubature->getNumPoints();

    dyn_cub_points = af.template buildArray<Scalar,IP,Dim>("cub_points",num_ip, num_space_dim);
    dyn_cub_weights = af.template buildArray<Scalar,IP>("cub_weights",num_ip);

    cub_points = af.template buildStaticArray<Scalar,IP,Dim>("cub_points",num_ip, num_space_dim);

    if (ir->isSide())
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
  }

// ***********************************************************
// * Evaluation of values - NOT specialized
// ***********************************************************

  template <typename Scalar,
            template <typename DataT,
               typename Tag0, typename Tag1, typename Tag2,
               typename Tag3, typename Tag4, typename Tag5,
               typename Tag6, typename Tag7> class Array >
  template<typename NodeCoordinateArray>
  inline
  void IntegrationValues2<Scalar,Array>::
    evaluateValues(const NodeCoordinateArray& in_node_coordinates)
  {
    int num_space_dim = int_rule->topology->getDimension();
    if (int_rule->isSide() && num_space_dim==1) {
       std::cout << "WARNING: 0-D quadrature rule ifrastructure does not exist!!! Will not be able to do "
                 << "non-natural integration rules.";
       return; 
    }
    
    Intrepid::CellTools<Scalar> cell_tools;
    
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

    cell_tools.setJacobianTemp(jac, cub_points, node_coordinates, 
			   *(int_rule->topology));
    
    cell_tools.setJacobianInvTemp(jac_inv, jac);
    
    cell_tools.setJacobianDetTemp(jac_det, jac);
    
#ifdef INTREPID_KOKKOS_PROBLEMS
    if (!int_rule->isSide()) {
       Intrepid::FunctionSpaceTools::
         computeCellMeasure<Scalar>(weighted_measure, jac_det, cub_weights);
    }
    else if(int_rule->spatial_dimension==3) {
       Intrepid::FunctionSpaceTools::
         computeFaceMeasure<Scalar>(weighted_measure, jac, cub_weights,int_rule->side,*int_rule->topology);
    }
    else if(int_rule->spatial_dimension==2) {
       Intrepid::FunctionSpaceTools::
         computeEdgeMeasure<Scalar>(weighted_measure, jac, cub_weights,int_rule->side,*int_rule->topology);
    }
    else TEUCHOS_ASSERT(false);
#endif
    
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

    Intrepid::RealSpaceTools<Scalar>::inverseTemp(contravarient, covarient);

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

    // IP coordinates
    {
      cell_tools.mapToPhysicalFrameTemp(ip_coordinates, cub_points, node_coordinates, *(int_rule->topology));
    }

  }
}

#endif

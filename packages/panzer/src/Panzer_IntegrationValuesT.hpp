
#include "Shards_CellTopology.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"
#include "Panzer_ArrayTraits.hpp"

// ***********************************************************
// * Specializations of setupArrays() for different array types
// ***********************************************************

namespace panzer {
  
  // * Specialization for Intrepid::FieldContainer<double>
  template<>
  inline
  void panzer::IntegrationValues<double,Intrepid::FieldContainer<double> >::
  setupArrays(const Teuchos::RCP<panzer::IntegrationRule>& ir)
  {
    int_rule = ir;
    
    int num_nodes = ir->topology->getNodeCount();
    int num_cells = ir->workset_size;
    int num_space_dim = ir->topology->getDimension();

    Intrepid::DefaultCubatureFactory<double,Intrepid::FieldContainer<double> > 
      cubature_factory;
    
    if (ir->isSide())
      intrepid_cubature = cubature_factory.create(*(ir->side_topology), 
						  ir->cubature_degree);
    else
      intrepid_cubature = cubature_factory.create(*(ir->topology), 
						  ir->cubature_degree);

    int num_ip = intrepid_cubature->getNumPoints();

    cub_points = Intrepid::FieldContainer<double>(num_ip, num_space_dim);

    if (ir->isSide())
      side_cub_points = 
	Intrepid::FieldContainer<double>(num_ip, 
					 ir->side_topology->getDimension());
    
    cub_weights = Intrepid::FieldContainer<double>(num_ip);
    
    node_coordinates = 
      Intrepid::FieldContainer<double>(num_cells, num_nodes, num_space_dim);
    
    jac = Intrepid::FieldContainer<double>(num_cells, num_ip, num_space_dim,
					   num_space_dim);
    
    jac_inv = Intrepid::FieldContainer<double>(num_cells, num_ip, 
					       num_space_dim,
					       num_space_dim);
    
    jac_det = Intrepid::FieldContainer<double>(num_cells, num_ip);
    
    weighted_measure = 
      Intrepid::FieldContainer<double>(num_cells, num_ip);
    
    covarient = 
      Intrepid::FieldContainer<double>(num_cells, num_ip, 
				       num_space_dim,
				       num_space_dim);

    contravarient = 
      Intrepid::FieldContainer<double>(num_cells, num_ip, 
				       num_space_dim,
				       num_space_dim);

    norm_contravarient = 
      Intrepid::FieldContainer<double>(num_cells, num_ip);

  }
  
// ***********************************************************
// * Evaluation of values - NOT specialized
// ***********************************************************

  template<typename Scalar, typename Array>
  template<typename NodeCoordinateArray>
  inline
  void panzer::IntegrationValues<Scalar,Array>::
    evaluateValues(const NodeCoordinateArray& in_node_coordinates)
  {
    
    Intrepid::CellTools<Scalar> cell_tools;
    
    if (!int_rule->isSide())
      intrepid_cubature->getCubature(cub_points, cub_weights);
    else {
      intrepid_cubature->getCubature(side_cub_points, cub_weights);
      
      cell_tools.mapToReferenceSubcell(cub_points, 
				       side_cub_points,
				       int_rule->spatial_dimension-1,
				       int_rule->side, 
				       *(int_rule->topology));
    }


    {
      typedef typename 
	panzer::ArrayTraits<Scalar,NodeCoordinateArray>::size_type size_type;

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
    
    Intrepid::FunctionSpaceTools::
      computeCellMeasure<Scalar>(weighted_measure, jac_det, cub_weights);
    
    // Shakib contravarient metric tensor
    typedef typename 
      panzer::ArrayTraits<Scalar,Array>::size_type size_type;

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

    Intrepid::RealSpaceTools<Scalar>::inverse(contravarient, covarient);

    /*
    for (std::size_t cell = 0; cell < contravarient.dimension(0); ++cell) {
      std::cout << "cell = " << cell << std::endl;
      for (std::size_t ip = 0; ip < contravarient.dimension(1); ++ip) {
      std::cout << "  ip = " << ip << std::endl;
	for (std::size_t i = 0; i < contravarient.dimension(2); ++i) {
	  for (std::size_t j = 0; j < contravarient.dimension(2); ++j) {
	    std::cout << "contravarient(" << i << "," << j << ") = " << contravarient(cell,ip,i,j) << std::endl;
	  }
	}
      }
    }
    */

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
}

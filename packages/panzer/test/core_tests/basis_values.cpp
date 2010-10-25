#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Panzer_CellData.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_IntegrationValues.hpp"
#include "Panzer_ArrayTraits.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Panzer_BasisValues.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using panzer::IntegrationRule;
using Intrepid::FieldContainer;

namespace panzer {

  TEUCHOS_UNIT_TEST(integration_values, volume)
  {
    const int num_cells = 20;
    const int base_cell_dimension = 2;
    const panzer::CellData cell_data(num_cells, base_cell_dimension);

    const int cubature_degree = 2;    
    RCP<IntegrationRule> int_rule = 
      rcp(new IntegrationRule(cubature_degree, cell_data));
    
    panzer::IntegrationValues<double,Intrepid::FieldContainer<double> > 
      int_values;

    int_values.setupArrays(int_rule);

    const int num_vertices = int_rule->topology->getNodeCount();
    FieldContainer<double> node_coordinates(num_cells, num_vertices,
					    base_cell_dimension);



    // Set up node coordinates.  Here we assume the following
    // ordering.  This needs to be consistent with shards topology,
    // otherwise we will get negative determinates

    // 3(0,1)---2(1,1)
    //   |    0  |
    //   |       |
    // 0(0,0)---1(1,0)

    typedef panzer::ArrayTraits<double,FieldContainer<double> >::size_type size_type;
    const size_type x = 0;
    const size_type y = 1;
    for (size_type cell = 0; cell < node_coordinates.dimension(0); ++cell) {
      node_coordinates(cell,0,x) = 0.0;
      node_coordinates(cell,0,y) = 0.0;
      node_coordinates(cell,1,x) = 1.0;
      node_coordinates(cell,1,y) = 0.0;
      node_coordinates(cell,2,x) = 1.0;
      node_coordinates(cell,2,y) = 1.0;
      node_coordinates(cell,3,x) = 0.0;
      node_coordinates(cell,3,y) = 1.0;
    }

    int_values.evaluateValues(node_coordinates);
    
    const std::string basis_type = "Q2";
  
    RCP<panzer::Basis> basis = rcp(new panzer::Basis(basis_type, *int_rule));

    panzer::BasisValues<double,Intrepid::FieldContainer<double> > basis_values;

    basis_values.setupArrays(basis);
    
    basis_values.evaluateValues(int_values.cub_points,
				int_values.jac_inv,
				int_values.weighted_measure);

  }

}

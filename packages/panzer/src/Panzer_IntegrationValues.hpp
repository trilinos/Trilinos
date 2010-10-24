
#ifndef PANZER_INTEGRATION_VALUES_HPP
#define PANZER_INTEGRATION_VALUES_HPP

#include "Teuchos_RCP.hpp"
#include "Panzer_IntegrationRule.hpp"

namespace panzer {

  template<typename Scalar,typename Array>
  struct IntegrationValues {
    
    //! Sizes/allocates memory for arrays
    void setupArrays(const Teuchos::RCP<panzer::IntegrationRule>& ir);

    template<typename NodeCoordinateArray>
    void evaluateValues(const NodeCoordinateArray& node_coordinates);

    Array cub_points;          // <IP,Dim>
    Array side_cub_points;     // <IP,SideDim> points on face topology (dim-1)
    Array cub_weights;         // <IP>
    Array node_coordinates;    // <Cell,NODE,Dim>
    Array jac;                 // <Cell,IP,Dim,Dim>
    Array jac_inv;             // <Cell,IP,Dim,Dim>
    Array jac_det;             // <Cell,IP>
    Array weighted_measure;    // <Cell,IP>

    Teuchos::RCP<panzer::IntegrationRule> int_rule;

    Teuchos::RCP< Intrepid::Cubature<Scalar,Array> > intrepid_cubature;

    // for Shakib stabilization <Cell,IP,Dim,Dim>
    Array covarient; 
    Array contravarient; 
    Array norm_contravarient; 

  };

} // namespace panzer

#include "Panzer_IntegrationValuesT.hpp"

#endif

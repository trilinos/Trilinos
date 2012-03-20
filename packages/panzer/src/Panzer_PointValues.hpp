#ifndef PANZER_POINT_VALUES_HPP
#define PANZER_POINT_VALUES_HPP

#include "Teuchos_RCP.hpp"
#include "Panzer_PointRule.hpp"

namespace panzer {

  template <typename Scalar,typename Array>
  struct PointValues {
    
    //! Sizes/allocates memory for arrays
    template <typename ArrayFactory>
    void setupArrays(const Teuchos::RCP<const panzer::PointRule>& pr,const ArrayFactory & af);

    template <typename NodeCoordinateArray,typename PointCoordinateArray>
    void evaluateValues(const NodeCoordinateArray & node_coordinates,const PointCoordinateArray & point_coordinates);

    Array coords_ref;          // <Point,Dim>
    Array node_coordinates;    // <Cell,NODE,Dim>
    Array jac;                 // <Cell,Point,Dim,Dim>
    Array jac_inv;             // <Cell,Point,Dim,Dim>
    Array jac_det;             // <Cell,Point>

    // cell points
    Array point_coords;      // <Cell,Point,Dim>

    Teuchos::RCP<const panzer::PointRule> point_rule;
  };

} // namespace panzer

#include "Panzer_PointValues_impl.hpp"

#endif

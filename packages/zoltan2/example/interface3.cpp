// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// Questions? Contact Lee Ann Riesen (lriesen@sandia.gov)
//
// ***********************************************************************
// @HEADER


// ***********************************************************************
// This is pseudocode showing the a Zoltan2 interface where the central
// Zoltan2 object is a Mapping.  It has user input and objectives and 
// it can compute a mapping from input objects to parts, colors or an
// ordering.
//
// Use case:
//
// Coordinates are stored in an Epetra_MultiVector.
// Partition the coordinates using RCB.
// Evaluate the imbalance.
// Partition again using RIB.
// Evaluate the imbalance again and keep the best result.
//
// ***********************************************************************

#include <Zoltan2.hpp>
#include <Epetra_MultiVector.h>
#include <Teuchos_ParameterList.hpp>

int main()
{

  Epetra::MultiVector &mv = get_my_coordinates();
  
  // Epetra_MultiVector uses
  //    float    for scalar values
  //    int      for local ids
  //    int      for global ids
  //
  //    4th template parameter is application GID type and we're
  //    just using the Epetra global id.
  //
  // The adapter provide a uniform interface to the input.

  Z2::Epetra_MultiVectorAdapter<float, int, int, int> input(mv);

  Z2_Interface3::PartMapping<float, int, int, int> mapping(input);

  // A Mapping has an Objective and we set the objective here.

  mapping.partition_vertices();
  mapping.unit_weights();
  mapping.rectilinear_regions();        // will do RCB
  mapping.uniform_part_sizes();
  mapping.set_number_of_parts(num_procs);

  // Parameters - In interface 3 it's not as important to separate objectives
  //   from parameters.  But it may still be a good idea.

  Teuchos::ParameterList params;
  params.set("Debug Level", "3");

  mapping.set_parameters(params);

  float imbalance_before = mapping.imbalance();

  mapping.solve();

  std::vector<int> &gids = my_gids(mv);
  std::vector<int> rcb_parts;

  mapping.get_mapping(gids, rcb_parts);

  float rcb_imbalance= mapping.imbalance();

  // Try RIB

  mapping.unset_rectilinear_regions();
  params.set("METHOD", "RIB");

  mapping.solve();

  std::vector<int> rib_parts;

  mapping.get_mapping(gids, rib_parts);

  float rib_imbalance= mapping.imbalance();

  // Pick the best

  if ((rib_imbalance < rcb_imbalance) &&
      (rib_imbalance < imbalance_before) ){

    redistribute_mv(mv, rib_parts);
  }
  else if ((rcb_imbalance < rib_imbalance) &&
           (rcb_imbalance < imbalance_before) ){

    redistribute_mv(mv, rcb_parts);
  }
}

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
// This is pseudocode showing the a Zoltan2 interface where Zoltan methods
// are functions in the Zoltan2 namespace.
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

#include <Zoltan2_Functions.hpp>
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

  Z2::PartitioningObjective<float, int, int, int> objective(mv);

  // The objective defines the problem to be solved without
  // telling Zoltan how to solve it.

  objective.partition_vertices();
  objective.unit_weights();
  objective.rectilinear_regions();
  objective.uniform_part_sizes();
  objective.set_number_of_parts(num_procs);

  Z2::PartitioningResult<float, int> result();

  // A parameter list can be used to set non-problem-related
  // parameters and also to give specific directions to the
  // function on how to solve the problem.

  Teuchos::ParameterList params;
  params.set("Debug Level", "3");

  Z2_Interface2::partition(objective, params, result);

  // Try another approach

  objective.unset_rectilinear_regions();

  params.set("Method", "RIB");

  Z2::PartitioningResult<float, int> second_result();
  
  Z2_Interface2::partition(objective, params, second_result);

  // Compare the balance achieved in the two results

  Z2::evaluate(objective, result);
  Z2::evaluate(objective, second_result);

  if (second_result.imbalance < result.imbalance){
    result = second_result;
  }

  // Application uses result queries to build a new Epetra_Map and
  //  redistribute the vector.

  redistribute_my_vector(mv, result);
}

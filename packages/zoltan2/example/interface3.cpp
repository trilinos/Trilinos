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
// Zoltan2 object is a Problem.  It has user input and objectives and 
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

  // Objective describes the problem to be solved

  Z2::Objective<float, int, int, int> objective(input)

  objective.partitionVertices();
  objective.unitWeights();
  objective.rectilinearRegions();        // will do RCB
  objective.uniformPartSizes();
  objective.setNumberOfParts(numprocs);

  // Parameters - govern the behavior of the Problem, but not
  //    really related to computing the solution.

  Teuchos::ParameterList params;
  params.set("Debug Level", "3");

  Z2::Problem<float, int, int,int>  problem(input, params, objective);

  float imbalance_before = problem.imbalance();

  problem.solve();

  std::map<GNO, int> rcb_mapping;

  problem.getMapping(rcb_mapping);

  float rcb_imbalance= problem.imbalance();

  // Try RIB

  objective.unsetRectilinearRegions();

  problem.setObjective(objective)

  params.set("METHOD", "RIB");

  problem.solve();

  std::map<GNO, int> rib_mapping;

  problem.getMapping(rib_mapping);

  float rib_imbalance= problem.imbalance();

  // Pick the best

  if ((rib_imbalance < rcb_imbalance) &&
      (rib_imbalance < imbalance_before) ){

    Z2:Migrate(input, rib_mapping);
  }
  else if ((rcb_imbalance < rib_imbalance) &&
           (rcb_imbalance < imbalance_before) ){

    Z2:Migrate(input, rcb_mapping);
  }
}

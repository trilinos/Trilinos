#ifndef MLAPI_DEFAULTS_H
#define MLAPI_DEFAULTS_H

#include "Teuchos_ParameterList.hpp"

namespace MLAPI 
{

// ====================================================================== 
//! Sets default values in input \c List.
// ====================================================================== 

void SetDefaults(Teuchos::ParameterList& List)
{
  // general defaults for multilevel
  List.set("max levels", 10);
  List.set("PDE equations", 1);

  // defaults for aggregation
  List.set("aggregation: type", "Uncoupled");
  List.set("aggregation: threshold", 0.0);
  List.set("aggregation: damping factor", 1.333);
  List.set("coarse: max size", 32);
  List.set("eigen-analysis: type", "Anorm");

  // defaults for smoother
  List.set("smoother: type", "symmetric Gauss-Seidel");
  List.set("smoother: sweeps", 1);
  List.set("smoother: damping factor", 0.67);

  // defaults for coarse solver
  List.set("coarse: type", "Amesos-KLU");

}

} // namespace MLAPI

#endif // MLAPI_DEFAULTS_H

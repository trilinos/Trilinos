// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_AdaptivePreconditionerFactory_hpp__
#define __Teko_AdaptivePreconditionerFactory_hpp__

#include "Teuchos_RCP.hpp"
#include "Teko_PreconditionerFactory.hpp"
#include "Teko_Utilities.hpp"

namespace Teko {

struct AdaptiveSolverSettings {
  int numAppliesCycle            = 100;
  double targetResidualReduction = 0.1;
};

/** Adaptive sub-block solver.
  * Given a user-provided schedule of sub-block solvers, progressively try more sub-block solvers
  * until sufficient residual reduction has been reached.
  *
  \code
       <Parameter name="Type" type="string" value="Adaptive"/>

       <!-- Target relative residual reduction. Default: 1e-1. -->
       <Parameter name="Target Residual Reduction" type="double" value="1e-2"/>

       <!-- Specify schedule of inverses:
         Robustness/cost should increase, with, e.g., 3 more expensive and robust than 2, 2 more
  expensive and robust than 1, etc.
       -->
       <Parameter name="Inverse Type 1" type="string" value="GMRES"/>
       <Parameter name="Preconditioner Type 1" type="string" value="Jacobi"/>
       <Parameter name="Inverse Type 2" type="string" value="GMRES"/>
       <Parameter name="Preconditioner Type 2" type="string" value="MueLu"/>
       <Parameter name="Inverse Type 3" type="string" value="Amesos2"/>

       <!--
         If a system has more rows than the given maximum size, avoid constructing the
  preconditioner/solver at this position. For example, this can be used to prevent a direct solver
  from being used on a sufficiently large sub-block. Default:
  std::numeric_limits<SizeOrdinalType>::max().
       -->
       <Parameter name="Maximum Size 3" type="long long int" value="10000"/>

       <!-- Number of times to try a solver before resetting to the first solver. Default: 100.-->
       <Parameter name="Number of Successful Applications Before Resetting" type="int" value="100"/>
  \endcode
  */
class AdaptivePreconditionerFactory : public PreconditionerFactory {
 public:
  using SizeOrdinalType = Thyra::Ordinal;

  AdaptivePreconditionerFactory() = default;

  ~AdaptivePreconditionerFactory() override = default;

  LinearOp buildPreconditionerOperator(LinearOp& lo, PreconditionerState& state) const override;

  void initializeFromParameterList(const Teuchos::ParameterList& pl) override;

  const std::vector<Teuchos::RCP<InverseFactory>>& get_inverses() const { return inverses; }

  const std::vector<Teuchos::RCP<InverseFactory>>& get_preconditioners() const {
    return preconditioners;
  }

 private:
  std::vector<Teuchos::RCP<InverseFactory>> inverses;
  std::vector<Teuchos::RCP<InverseFactory>> preconditioners;
  std::vector<SizeOrdinalType> maximumSizes;
  AdaptiveSolverSettings settings{};
};

}  // end namespace Teko

#endif
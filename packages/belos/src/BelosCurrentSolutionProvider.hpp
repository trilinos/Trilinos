// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_CURRENT_SOLUTION_PROVIDER_HPP
#define BELOS_CURRENT_SOLUTION_PROVIDER_HPP

/*! \file BelosCurrentSolutionProvider.hpp
    \brief Optional interface for iterations that can provide a current solution estimate.
*/

#include "BelosConfigDefs.hpp"
#include "BelosIteration.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosTypes.hpp"

#include "Teuchos_RCP.hpp"

namespace Belos {

template<class ScalarType, class MV, class OP, class DM = DefaultDenseMatrix<int,ScalarType>>
class CurrentSolutionProvider : virtual public Iteration<ScalarType,MV,OP,DM> {
private:
  typedef MultiVecTraits<ScalarType,MV,DM> MVT;

public:
  virtual ~CurrentSolutionProvider() {}

  //! Whether this iteration can currently provide a solution estimate.
  virtual bool hasCurrentSolution() const {
    return !this->getProblem().updateSolution(Teuchos::null).is_null();
  }

  //! Get the current update in solution space, if available separately.
  virtual Teuchos::RCP<const MV> getCurrentSolutionUpdate() const {
    return Teuchos::null;
  }

  //! Get the current solution estimate in solution space.
  virtual Teuchos::RCP<MV> getCurrentSolution() const {
    Teuchos::RCP<MV> curSoln = this->getProblem().updateSolution(Teuchos::null);
    if (curSoln.is_null()) {
      return Teuchos::null;
    }
    return MVT::CloneCopy(*curSoln);
  }
};

} // end Belos namespace

#endif /* BELOS_CURRENT_SOLUTION_PROVIDER_HPP */

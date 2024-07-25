// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_LINESEARCH_CURVILINEAR_H
#define NOX_LINESEARCH_CURVILINEAR_H

#include "NOX_Common.H"

#ifdef WITH_PRERELEASE

#include "NOX_LineSearch_Generic.H" // base class
#include "NOX_LineSearch_Utils_Printing.H"  // class data member
#include "NOX_LineSearch_Utils_Slope.H"     // class data member

// Forward declarations
namespace NOX {
  namespace Direction {
    class Tensor;
  }
  class LineSearchCounters;
}

namespace NOX {
namespace LineSearch {

  /*! \brief A curvilinear line search algorithm developed for tensor methods.

    This object implements the curvilinear line search developed for
    tensor methods.  Additionally it implements a number of other
    lines searches for direct comparison.
  */
class Tensor : public Generic {

public:

  //! Constructor
  Tensor(const Teuchos::RCP<NOX::GlobalData>& gd,
     Teuchos::ParameterList& params);

  //! Destructor
  ~Tensor();

  // derived
  bool reset(const Teuchos::RCP<NOX::GlobalData>& gd,
         Teuchos::ParameterList& lsParams);

  // derived
  bool compute(NOX::Abstract::Group& newgrp, double& step,
           const NOX::Abstract::Vector& dir,
           const NOX::Solver::Generic& s);

protected:

  //! Performs a standard tensor linesearch (tensor or Newton direction)
  bool performLinesearch(NOX::Abstract::Group& newsoln,
             double& step,
             const NOX::Abstract::Vector& dir,
             const NOX::Solver::Generic& s,
             const NOX::Direction::Tensor& direction);

  //! Returns true if converged.
  bool checkConvergence(double newValue, double oldValue, double oldSlope, double step,
            double eta, int nIters, int nNonlinearIters) const;

  //! Select lambda for linesearch (quadratic or halving)
  double selectLambda(double newf, double oldf,
              double oldfprime, double lambda);

  //! Prints a warning message saying that the slope is negative
  void printBadSlopeWarning(double slope) const;

protected:

  //! Type of tensor linesearch
  enum LineSearchType {Curvilinear, Standard, Dual, Newton};

  //! Algorithms used to determine convergence of the line search
  enum ConvergenceCriteriaType {
    //! Sufficient decrease condition
    ArmijoGoldstein,
    //! Ared/Pred condition
    AredPred,
    //! Just accept the first step
    None
  };

  //! Types of lambda selection
  enum LambdaSelectionType {Halving, Quadratic};

  //! Flag for lambda selection (Halving/Quadratic)
  LambdaSelectionType lambdaSelection;

  //! Choice of convergence criteria
  ConvergenceCriteriaType suffDecrCond;

  //! Choice of line search
  LineSearchType lsType;

  //! Flag for Newton direction
  bool isNewtonDirection;

  //! Minimum step length (i.e., when we give up)
  double minStep;

  //! Default step
  double defaultStep;

  //! Default step for linesearch failure
  double recoveryStep;

  //! Maximum iterations
  int maxIters;

  //! Scaling factor for the Armijo-Goldstein condition
  double alpha;

  //! Global data pointer.  Keep this so the parameter list remains valid.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  //! Pointer to the input parameter list.  We need this to create an "Output" sublist to store output parameters from the counter object.
  Teuchos::ParameterList* paramsPtr;

  //! Printing Utils
  NOX::Utils utils;

  //! Common line search printing utilities.
  NOX::LineSearch::Utils::Printing print;

  //! Common common counters for line searches.
  NOX::LineSearchCounters* counter;

  //! Common slope calculations for line searches.
  NOX::LineSearch::Utils::Slope slopeObj;

  //! Counter for number of Jacobian-vector products
  mutable int multsJv;

};
} // namespace LineSearch
} // namespace NOX

#endif

#endif

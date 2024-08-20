// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_LINESEARCH_UTILS_PRINTING_H
#define NOX_LINESEARCH_UTILS_PRINTING_H

#include "NOX_LineSearch_Generic.H"
#include "NOX_Utils.H"
#include "NOX_Common.H"

// Forward declaration
namespace NOX {
  class Utils;
}

namespace NOX {

namespace LineSearch {

namespace Utils {

//! %Common line search utilites for printing line search information to the screen.
/*!
  All line searches should print output results in a similar format.
  This utility provides common output routines.
*/

class Printing : public NOX::Utils {

public:

  //! Default constructor
  Printing(const Teuchos::RCP<NOX::Utils>& u);

  //! Destructor
  virtual ~Printing();

  void reset(const Teuchos::RCP<NOX::Utils>& u);

  //! Prints the opening information.
  void printOpeningRemarks(const std::string& lineSearchName) const;

  //! Print out step information for the inner iterations of a line search algorithm.
  /*!
    Example of output from the inner iterations of a Polynomial line search:

\verbatim
************************************************************************
-- Polynomial Line Search --
  1: step = 1.000e+00 oldf = 2.403e+00 newf = 1.076e+03
  2: step = 1.000e-01 oldf = 2.403e+00 newf = 4.440e+00
  3: step = 1.000e-02 oldf = 2.403e+00 newf = 2.394e+00 (STEP ACCEPTED!)
************************************************************************
\endverbatim

    \param unscaleF - If this is true (the default), than the values
    printed are \f$ \sqrt{2 * {\rm oldf}} \f$ and \f$ \sqrt{2 * {\rm
    newf}} \f$. This is to accomodate the standard merit function, \f$
    \phi(x) = \frac{1}{2} \|F(x)\|^2 \f$.

  */
  void printStep(int n, double step, double oldf, double newf,
         const std::string s = "", bool unscaleF = true) const;

};
} // namespace Utils
} // namespace LineSearch
} // namespace NOX

#endif

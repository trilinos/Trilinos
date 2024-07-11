// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   NOX_EXCEPTIONS_H
#define   NOX_EXCEPTIONS_H

namespace NOX {
namespace Exceptions {

  /*!
    \brief A generic exception class denoting a solver failure.

    This class is a member of namespace NOX::Exceptions.
  */
  class SolverFailure : public std::logic_error {

  public:

    //! Constructor.
    SolverFailure(const std::string& what) : std::logic_error(what) {};

  }; // end of class SolverFailure

} // end of namespace Exceptions
} // end of namespace NOX

#endif // NOX_EXCEPTIONS_H

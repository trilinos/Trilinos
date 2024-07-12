// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_MULTIPHYSICS_DATAEXCHANGE_INTERFACE_H
#define NOX_MULTIPHYSICS_DATAEXCHANGE_INTERFACE_H

#include "NOX_Common.H"

namespace NOX {
namespace Multiphysics {
namespace DataExchange {

  /*!
    \brief Provides a set of interfaces for users to provide information about exchanging data between registered NOX solvers
  */

class Interface {

public:

  //! Constructor
  Interface() {};

  //! Destructor
  virtual ~Interface() {};

  //! Exchange data for all registered problems
  virtual bool exchangeAllData() = 0;

  //! Exchange data for a specified problem - brings needed data from others to this problem
  virtual bool exchangeDataTo(int solverId) = 0;

};
} // namespace DataExchange
} // namespace Multiphysics
} // namespace NOX

#endif

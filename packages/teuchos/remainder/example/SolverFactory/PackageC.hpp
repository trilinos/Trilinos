// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGE_C_HPP
#define PACKAGE_C_HPP

//
// Header file for Package C.
//

#include "Common.hpp"

namespace C {

  //
  // This solver is independent of other solvers.
  //
  template<class MV, class OP, class NormType>
  class Solver5 : public Common::LinearSolverTestBase<MV, OP, NormType> {
  protected:
    std::string name () const {
      return "Solver5";
    }

  public:
    virtual ~Solver5 () {}

    void solve (MV& /* X */, const MV& /* B */ ) {
      std::cout << this->name () << "::solve START" << std::endl;
      std::cout << this->name () << "::solve END" << std::endl;
    }
  };

  //
  // This solver is independent of other solvers.
  //
  template<class MV, class OP, class NormType>
  class Solver6 : public Common::LinearSolverTestBase<MV, OP, NormType> {
  protected:
    std::string name () const {
      return "Solver6";
    }

  public:
    virtual ~Solver6 () {}

    void solve (MV& /* X */, const MV& /* B */ ) {
      std::cout << this->name () << "::solve START" << std::endl;
      std::cout << this->name () << "::solve END" << std::endl;
    }
  };

  //
  // Package C's solver factory.
  //
  template<class MV, class OP, class NormType>
  class FactoryC : public Trilinos::Details::LinearSolverFactory<MV, OP, NormType> {
  public:
    Teuchos::RCP<Trilinos::Details::LinearSolver<MV, OP, NormType> >
    getLinearSolver (const std::string& solverName)
    {
      typedef Trilinos::Details::LinearSolver<MV, OP, NormType> solver_type;

      if (solverName == "5") {
        return Teuchos::RCP<solver_type> (new Solver5<MV, OP, NormType> ());
      }
      else if (solverName == "6") {
        return Teuchos::RCP<solver_type> (new Solver6<MV, OP, NormType> ());
      }
      else {
        std::ostringstream err;
        err << "C::FactoryC::getLinearSolver: Invalid solver name \""
            << solverName << "\"";
        throw std::invalid_argument (err.str ());
      }
    }
  };

} // namespace C

#endif // PACKAGE_C_HPP

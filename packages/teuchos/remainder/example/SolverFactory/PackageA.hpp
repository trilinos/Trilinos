// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGE_A_HPP
#define PACKAGE_A_HPP

//
// Header file for Package A.
//

#include "Common.hpp"

namespace A {

  //
  // This solver is independent of other solvers.
  //
  template<class MV, class OP, class NormType>
  class Solver1 : public Common::LinearSolverTestBase<MV, OP, NormType> {
  protected:
    std::string name () const {
      return "Solver1";
    }

  public:
    virtual ~Solver1 () {}

    void solve (MV& /* X */, const MV& /* Y */) {
      std::cout << this->name () << "::solve START" << std::endl;
      std::cout << this->name () << "::solve END" << std::endl;
    }
  };

  //
  // This solver uses Solver4 from Package B.
  //
  template<class MV, class OP, class NormType>
  class Solver2 : public Common::LinearSolverTestBase<MV, OP, NormType> {
  protected:
    std::string name () const {
      return "Solver2";
    }

  public:
    virtual ~Solver2 () {}

    void solve (MV& X, const MV& Y) {
      std::cout << this->name () << "::solve START" << std::endl;

      Teuchos::RCP<Trilinos::Details::LinearSolver<MV, OP, NormType> > solverB4 =
        Trilinos::Details::getLinearSolver<MV, OP, NormType> ("B", "4");
      if (solverB4.get () == NULL) {
        throw std::runtime_error ("Solver4 from package B has not been registered!");
      }
      // A real implementation would probably do something to X and Y
      // before or after calling the "inner" solver.
      solverB4->solve (X, Y);
      std::cout << this->name () << "::solve END" << std::endl;
    }
  };

  //
  // Package A's solver factory.
  //
  template<class MV, class OP, class NormType>
  class FactoryA : public Trilinos::Details::LinearSolverFactory<MV, OP, NormType> {
  public:
    // Get an instance of a solver from a particular package
    Teuchos::RCP<Trilinos::Details::LinearSolver<MV, OP, NormType> >
    getLinearSolver (const std::string& solverName)
    {
      typedef Trilinos::Details::LinearSolver<MV, OP, NormType> solver_type;

      if (solverName == "1") {
        return Teuchos::RCP<solver_type> (new Solver1<MV, OP, NormType> ());
      }
      else if (solverName == "2") {
        return Teuchos::RCP<solver_type> (new Solver2<MV, OP, NormType> ());
      }
      else {
        std::ostringstream err;
        err << "A::FactoryA::getLinearSolver: Invalid solver name \"" << solverName << "\"";
        throw std::invalid_argument (err.str ());
      }
    }
  };

} // namespace A

#endif // PACKAGE_A_HPP

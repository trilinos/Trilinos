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
  template<class MV, class OP>
  class Solver1 : public Common::SolverTestBase<MV, OP> {
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
  template<class MV, class OP>
  class Solver2 : public Common::SolverTestBase<MV, OP> {
  protected:
    std::string name () const {
      return "Solver2";
    }

  public:
    virtual ~Solver2 () {}

    void solve (MV& X, const MV& Y) {
      std::cout << this->name () << "::solve START" << std::endl;

      Teuchos::RCP<Trilinos::Details::Solver<MV, OP> > solverB4 = Trilinos::Details::getSolver<MV, OP> ("B", "4");
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
  template<class MV, class OP>
  class FactoryA : public Trilinos::Details::SolverFactory<MV, OP> {
  public:
    // Get an instance of a solver from a particular package
    Teuchos::RCP<Trilinos::Details::Solver<MV, OP> > getSolver (const std::string& solverName) {
      if (solverName == "1") {
        return Teuchos::RCP<Trilinos::Details::Solver<MV, OP> > (new Solver1<MV, OP> ());
      }
      else if (solverName == "2") {
        return Teuchos::RCP<Trilinos::Details::Solver<MV, OP> > (new Solver2<MV, OP> ());
      }
      else {
        std::ostringstream err;
        err << "A::FactoryA::getSolver: Invalid solver name \"" << solverName << "\"";
        throw std::invalid_argument (err.str ());
      }
    }
  };

} // namespace A

#endif // PACKAGE_A_HPP

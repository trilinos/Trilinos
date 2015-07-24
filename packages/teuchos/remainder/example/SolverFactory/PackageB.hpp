#ifndef PACKAGE_B_HPP
#define PACKAGE_B_HPP

//
// Header file for Package B.
//

#include "Common.hpp"

namespace B {

  //
  // This solver is independent of other solvers.
  //
  template<class MV, class OP>
  class Solver3 : public Common::SolverTestBase<MV, OP> {
  protected:
    std::string name () const {
      return "Solver3";
    }

  public:
    virtual ~Solver3 () {}

    void solve (MV& /* X */, const MV& /* Y */ ) {
      std::cout << this->name () << "::solve START" << std::endl;
      std::cout << this->name () << "::solve END" << std::endl;
    }
  };

  //
  // This solver uses Solver1 from package A.
  //
  template<class MV, class OP>
  class Solver4 : public Common::SolverTestBase<MV, OP> {
  protected:
    std::string name () const {
      return "Solver3";
    }

  public:
    virtual ~Solver4 () {}

    void solve (MV& X, const MV& B) {
      std::cout << this->name () << "::solve START" << std::endl;

      Teuchos::RCP<Trilinos::Details::Solver<MV, OP> > solverA1 =
        Trilinos::Details::getSolver<MV, OP> ("A", "1");
      if (solverA1.get () == NULL) {
        std::runtime_error ("Solver1 from package A has not been registered!");
      }
      solverA1->solve (X, B);

      std::cout << this->name () << "::solve END" << std::endl;
    }
  };

  //
  // Package B's solver factory.
  //
  template<class MV, class OP>
  class FactoryB : public Trilinos::Details::SolverFactory<MV, OP> {
  public:
    Teuchos::RCP<Trilinos::Details::Solver<MV, OP> > getSolver (const std::string& solverName) {
      if (solverName == "3") {
        return Teuchos::RCP<Trilinos::Details::Solver<MV, OP> > (new Solver3<MV, OP> ());
      }
      else if (solverName == "4") {
        return Teuchos::RCP<Trilinos::Details::Solver<MV, OP> > (new Solver4<MV, OP> ());
      }
      else {
        std::ostringstream err;
        err << "B::FactoryB::getSolver: Invalid solver name \"" << solverName << "\"";
        throw std::invalid_argument (err.str ());
      }
    }
  };

} // namespace B

#endif // PACKAGE_B_HPP


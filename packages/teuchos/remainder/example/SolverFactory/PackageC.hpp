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
  template<class MV, class OP>
  class Solver5 : public Common::SolverTestBase<MV, OP> {
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
  template<class MV, class OP>
  class Solver6 : public Common::SolverTestBase<MV, OP> {
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
  template<class MV, class OP>
  class FactoryC : public Trilinos::Details::SolverFactory<MV, OP> {
  public:
    Teuchos::RCP<Trilinos::Details::Solver<MV, OP> >
    getSolver (const std::string& solverName)
    {
      if (solverName == "5") {
        return Teuchos::RCP<Trilinos::Details::Solver<MV, OP> > (new Solver5<MV, OP> ());
      }
      else if (solverName == "6") {
        return Teuchos::RCP<Trilinos::Details::Solver<MV, OP> > (new Solver6<MV, OP> ());
      }
      else {
        std::ostringstream err;
        err << "C::FactoryC::getSolver: Invalid solver name \"" << solverName << "\"";
        throw std::invalid_argument (err.str ());
      }
    }
  };

} // namespace C

#endif // PACKAGE_C_HPP

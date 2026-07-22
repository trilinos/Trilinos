// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef COMMON_HPP
#define COMMON_HPP

//
// Header file for classes common to all "packages" in this example.
//

#include "Trilinos_Details_LinearSolver.hpp"
#include "Trilinos_Details_LinearSolverFactory.hpp"
#include <iostream>
#include <sstream>

// Namespace for classes common to all "packages" in this example.
namespace Common {

  // Stub of a MultiVector (MV) class, templated on Scalar type (the
  // type of its entries).
  template<class Scalar>
  class MultiVector {};

  // Stub of an Operator (OP) class, templated on Scalar type (the
  // template parameter of the MultiVector specialization that it
  // uses).
  template<class Scalar>
  class Operator {
  public:
    typedef MultiVector<Scalar> MV;

    void apply (MV& /* Y */, const MV& /* X */) {
      std::cout << "Operator<" << typeid (Scalar).name () << ">::apply" << std::endl;
    }
  };

  // Base classes of Trilinos::Details::LinearSolver must implement
  // all the pure virtual methods of that interface.  This base class
  // only exists to make the example more concise.  Its subclasses
  // must implement solve(), name(), and the virtual destructor.
  template<class MV, class OP, class NormType>
  class LinearSolverTestBase : public Trilinos::Details::LinearSolver<MV, OP, NormType> {
  protected:
    virtual std::string name () const = 0;

  public:
    virtual ~LinearSolverTestBase () {}

    void setMatrix (const Teuchos::RCP<const OP>& A) {
      std::cout << this->name () << "::setMatrix" << std::endl;
      A_ = A;
    }

    Teuchos::RCP<const OP> getMatrix () const {
      std::cout << this->name () << "::getMatrix" << std::endl;
      return A_; // this could be null if setMatrix wasn't called
    }

    void setParameters (const Teuchos::RCP<Teuchos::ParameterList>& /* params */ ) {
      std::cout << this->name () << "::setParameters" << std::endl;
    }

    void symbolic () {
      std::cout << this->name () << "::symbolic" << std::endl;
    }

    void numeric () {
      std::cout << this->name () << "::numeric" << std::endl;
    }

  private:
    Teuchos::RCP<const OP> A_; // the matrix given to setMatrix
  };

} // namespace Common

#endif // COMMON_HPP

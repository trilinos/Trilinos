// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "PackageB.hpp"

namespace B {

  // Creating an instance of this "object" registers B::FactoryB<MV,
  // OP> with the central registry of packages' factories.  That lets
  // B::getLinearSolver create solvers from package B.
  template<class MV, class OP, class NormType>
  class RegisterFactoryB {
  public:
    RegisterFactoryB () {
#ifdef HAVE_TEUCHOSCORE_CXX11
      typedef std::shared_ptr<Trilinos::Details::LinearSolverFactory<MV, OP, NormType> > ptr_type;
#else
      typedef Teuchos::RCP<Trilinos::Details::LinearSolverFactory<MV, OP, NormType> > ptr_type;
#endif // HAVE_TEUCHOSCORE_CXX11

      ptr_type factory (new FactoryB<MV, OP, NormType> ());
      Trilinos::Details::registerLinearSolverFactory<MV, OP, NormType> ("B", factory);
    }
  };

} // namespace B

namespace { // (anonymous)
  //
  // See PackageA.cpp for an explanation of the macro and its use.
  //
#define INSTMACRO( SCALAR ) \
  B::RegisterFactoryB< Common::MultiVector< SCALAR >, Common::Operator< SCALAR >, SCALAR > registerer_##SCALAR;

  //B::RegisterFactoryB< Common::MultiVector<double>, Common::Operator<double>, double > registerer_double;
  INSTMACRO( double )

  //B::RegisterFactoryB< Common::MultiVector<float>, Common::Operator<float>, float > registerer_float;
  INSTMACRO( float )

} // namespace (anonymous)


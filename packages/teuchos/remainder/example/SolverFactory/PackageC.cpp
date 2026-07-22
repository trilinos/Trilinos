// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "PackageC.hpp"

namespace C {

  // Creating an instance of this "object" registers C::FactoryC<MV,
  // OP> with the central registry of packages' factories.  That lets
  // C::getLinearSolver create solvers from package C.
  template<class MV, class OP, class NormType>
  class RegisterFactoryC {
  public:
    RegisterFactoryC () {
#ifdef HAVE_TEUCHOSCORE_CXX11
      typedef std::shared_ptr<Trilinos::Details::LinearSolverFactory<MV, OP, NormType> > ptr_type;
#else
      typedef Teuchos::RCP<Trilinos::Details::LinearSolverFactory<MV, OP, NormType> > ptr_type;
#endif // HAVE_TEUCHOSCORE_CXX11

      ptr_type factory (new FactoryC<MV, OP, NormType> ());
      Trilinos::Details::registerLinearSolverFactory<MV, OP, NormType> ("C", factory);
    }
  };

} // namespace C

namespace { // (anonymous)
  //
  // See PackageA.cpp for an explanation of the macro and its use.
  //
#define INSTMACRO( SCALAR ) \
  C::RegisterFactoryC< Common::MultiVector< SCALAR >, Common::Operator< SCALAR >, SCALAR > registerer_##SCALAR;

  //C::RegisterFactoryC< Common::MultiVector<double>, Common::Operator<double>, double > registerer_double;
  INSTMACRO( double )

  //C::RegisterFactoryC< Common::MultiVector<float>, Common::Operator<float>, float > registerer_float;
  INSTMACRO( float )

} // namespace (anonymous)


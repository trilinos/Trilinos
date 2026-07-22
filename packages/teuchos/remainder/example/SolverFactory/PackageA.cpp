// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "PackageA.hpp"

namespace A {

  // Creating an instance of this "object" registers A::FactoryA<MV,
  // OP> with the central registry of packages' factories.  That lets
  // getLinearSolver create solvers from package A.
  template<class MV, class OP, class NormType>
  class RegisterFactoryA {
  public:
    RegisterFactoryA () {
#ifdef HAVE_TEUCHOSCORE_CXX11
      typedef std::shared_ptr<Trilinos::Details::LinearSolverFactory<MV, OP, NormType> > ptr_type;
#else
      typedef Teuchos::RCP<Trilinos::Details::LinearSolverFactory<MV, OP, NormType> > ptr_type;
#endif // HAVE_TEUCHOSCORE_CXX11

      ptr_type factory (new FactoryA<MV, OP, NormType> ());
      Trilinos::Details::registerLinearSolverFactory<MV, OP, NormType> ("A", factory);
    }
  };

} // namespace A

namespace { // (anonymous)

  // For each triple of types (MV, OP, NormType) of interest, register
  // A::FactoryA<MV, OP, NormType>.  We use MV =
  // Common::MultiVector<Scalar> and OP = Common::Operator<Scalar>
  // here, for various Scalar types.  NormType = Scalar simulates the
  // case where Scalar is real (not complex).
  //
  // This is a stub of what you likely will want to do with Tpetra and
  // its downstream solver packages.  See the public documentation of
  // Trilinos::Details::LinearSolverFactory for details.
  //
  // The ## operator in a macro appends two things.  For example, with
  // INSTMACRO( float ), registerer_##SCALAR becomes registerer_float.
  // This ensures that the different instances of RegisterFactoryA
  // have different names.

#define INSTMACRO( SCALAR ) \
  A::RegisterFactoryA< Common::MultiVector< SCALAR >, Common::Operator< SCALAR >, SCALAR > registerer_##SCALAR;

  //A::RegisterFactoryA< Common::MultiVector<double>, Common::Operator<double>, double > registerer_double;
  INSTMACRO( double )

  //A::RegisterFactoryA< Common::MultiVector<float>, Common::Operator<float>, float > registerer_float;
  INSTMACRO( float )

} // namespace (anonymous)

#include "PackageC.hpp"

namespace C {

  // Creating an instance of this "object" registers C::FactoryC<MV,
  // OP> with the central registry of packages' factories.  That lets
  // C::getSolver create solvers from package C.
  template<class MV, class OP>
  class RegisterFactoryC {
  public:
    RegisterFactoryC () {
#ifdef HAVE_TEUCHOSCORE_CXX11
      typedef std::shared_ptr<Trilinos::Details::SolverFactory<MV, OP> > ptr_type;
#else
      typedef Teuchos::RCP<Trilinos::Details::SolverFactory<MV, OP> > ptr_type;
#endif // HAVE_TEUCHOSCORE_CXX11

      ptr_type factory (new FactoryC<MV, OP> ());
      Trilinos::Details::registerFactory<MV, OP> ("C", factory);
    }
  };

} // namespace C

namespace { // (anonymous)
  //
  // See PackageA.cpp for an explanation of the macro and its use.
  //
#define INSTMACRO( SCALAR ) \
  C::RegisterFactoryC< Common::MultiVector< SCALAR >, Common::Operator< SCALAR > > registerer_##SCALAR;

  //C::RegisterFactoryC< Common::MultiVector<double>, Common::Operator<double> > registerer_double;
  INSTMACRO( double )

  //C::RegisterFactoryC< Common::MultiVector<float>, Common::Operator<float> > registerer_float;
  INSTMACRO( float )

} // namespace (anonymous)


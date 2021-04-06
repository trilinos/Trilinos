#include "ROL_KokkosVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

int main( int argc, char* argv[] ) {

  using RealT  = float;
  using Device = Kokkos::Cuda;
  using Vector = ROL::KokkosVector<RealT,Device>;

  Kokkos::initialize(argc,argv); {

    Teuchos::GlobalMPISession mpiSession(&argc,&argv);

    int iprint = argc - 1;
    ROL::Ptr<std::ostream> os;
    ROL::nullstream bhs;
    if( iprint > 0 ) os = ROL::makePtrFromRef(std::cout);
    else             os = ROL::makePtrFromRef(bhs);

    int errorFlag = 0;
    auto errtol = ROL::ROL_THRESHOLD<RealT>();

    int n = 1000;

    try {

      Vector x(n), y(n), z(n);
      x.randomize();
      y.randomize();
      z.randomize();

      // Standard Tests
      auto consistency = x.checkVector(y, z, true, *os );
      ROL::StdVector<RealT> checkvec( ROL::makePtrFromRef(consistency) );

    
    } catch( std::exception& err ) {
      *os << err.what() << std::endl;
      errorFlag = -1000;
    }; // end try

    if( errorFlag != 0 ) std::cout << "End Result: TEST FAILED\n";
    else                 std::cout << "End Result: TEST PASSED\n";
  } Kokkos::finalize();

  return 0;
}

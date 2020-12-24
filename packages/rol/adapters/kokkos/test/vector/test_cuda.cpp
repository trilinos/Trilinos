#include "ROL_KokkosVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_DualView.hpp"

int main( int argc, char* argv[] ) {

  using RealT   = float;
  using VectorT = ROL::KokkosVector<RealT,Kokkos::Cuda>;

  Kokkos::initialize(argc,argv); {

    Teuchos::GlobalMPISession mpiSession(&argc,&argv);

    int iprint = argc - 1;
    ROL::Ptr<std::ostream> os;
    ROL::nullstream bhs;
    if( iprint > 0 ) os = ROL::makePtrFromRef(std::cout);
    else             os = ROL::makePtrFromRef(bhs);

    int errorFlag = 0;
    auto errtol = ROL::ROL_THRESHOLD<RealT>();

    int n = 100;

    try {

      VectorT x(n), y(n), z(n);
      x.randomize();
      y.randomize();
      z.randomize();

      *os << "\n\nCreated three random vectors of length " << n << std::endl;
      *os << "||x|| = " << x.norm() << std::endl;
      *os << "||y|| = " << y.norm() << std::endl;
      *os << "||z|| = " << z.norm() << std::endl;

      // Standard Tests
      auto consistency = x.checkVector(y, z, true, *os );
      ROL::StdVector<RealT> checkvec( ROL::makePtrFromRef(consistency) );
  
      if( checkvec.norm() > std::sqrt(ROL::ROL_EPSILON<RealT>()) ) ++errorFlag;    
      
      // Basis tests
      // Set to first basis vector
      auto zp = x.basis(0);
      auto znorm = zp->norm();
      *os << "Norm or ROL::Vector z (first basis vector): " << znorm << "\n";
      if( std::abs(znorm-1.0) > errtol ) {
        *os << "---> POSSIBLE ERROR ABOVE!\n";
        ++errorFlag;
      }
      
      // set to middle basis vector
      zp = x.basis(n/2);
      znorm = zp->norm();
      *os << "\nNorm of ROL::Vector z ('middle' basis vector): " << znorm << "\n";
      if( std::abs(znorm-1.0) > errtol ) {
        *os << "---> POSSIBLE ERROR ABOVE!\n";
        ++errorFlag;
      }
 
      // set to last basis vector
      zp = x.basis(n-1);
      znorm = zp->norm();
      *os << "\nNorm of ROL::Vector z ('middle' basis vector): " << znorm << "\n";
      if( std::abs(znorm-1.0) > errtol ) {
        *os << "---> POSSIBLE ERROR ABOVE!\n";
        ++errorFlag;
      }
       
      // Repeat the checkVector tests with a zero vector
      x.scale(0.0);
      consistency = x.checkVector(x,x,true,*os);
      if( checkvec.norm() > 0.0 ) ++errorFlag;

    } catch( std::exception& err ) {
      *os << err.what() << std::endl;
      errorFlag = -1000;
    }; // end try

    if( errorFlag != 0 ) std::cout << "End Result: TEST FAILED\n";
    else                 std::cout << "End Result: TEST PASSED\n";
  } Kokkos::finalize();

  return 0;
}

#include "ROL_Types.hpp"
#include "ROL_KokkosVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using RealT   = float;
using VectorT = ROL::KokkosVector<RealT,Kokkos::Cuda>;

// 1) Create two vectors with the same number of elements:
// 
//    - A ROL::KokkosVector on a device (GPU)
//    - A ROL::StdVector on the host (CPU)
//
// 2) Fill both vectors with the same randomly-generated element values
//
// 3) Apply a ROL::Elementwise::UnaryFunction to the StdVector and the corresponding 
//    Visitor-Factory-generated ROL::Elementwise::KokkosUnaryFunction to the KokkosVector
//    
// 4) Compare the elemental sums and norms of the two vectors

template<typename Generator, typename Distribution> 
int test_unary_function( const ROL::Elementwise::UnaryFunction<RealT>& uf, int n, 
                         Generator& gen, Distribution& dist ) {

  auto errtol = std::sqrt(ROL::ROL_EPSILON<RealT>()*n);

  ROL::KokkosVector<RealT,Kokkos::Cuda> xkv(n), ekv(n);
  ROL::StdVector<RealT> xsv(n), esv(n);

  ekv.setScalar(1);
  esv.setScalar(1);
  
  auto xkv_host = xkv.view_host();
  auto xsv_ptr  = xsv.getVector();

  for( int i=0; i<n; ++i ) {
    auto rval = dist(gen);  
    xkv_host(i) = rval;
    (*xsv_ptr)[i] = rval;
  }
  
  xkv.modify_host();
  xkv.sync_device();

  xkv.applyUnary(uf);
  xsv.applyUnary(uf);

  auto sum_error  = std::abs(ekv.dot(xkv) - esv.dot(xsv));
  auto norm_error = std::abs(xkv.norm() - xsv.norm());

  int errorFlag = (sum_error > errtol) || (norm_error > errtol);

  return errorFlag;
}



int main( int argc, char* argv[] ) {


  Kokkos::initialize(argc,argv); {

    Teuchos::GlobalMPISession mpiSession(&argc,&argv);

    int iprint = argc - 1;
    ROL::Ptr<std::ostream> os;
    ROL::nullstream bhs;
    if( iprint > 0 ) os = ROL::makePtrFromRef(std::cout);
    else             os = ROL::makePtrFromRef(bhs);

    int errorFlag = 0;
    int n = 1000;

    std::mt19937_64 gen;
    std::uniform_real_distribution<RealT> udist(1,2);
    std::normal_distribution<RealT> ndist(0,1);

    try {

      ROL::Elementwise::AbsoluteValue<RealT> absval;
      errorFlag += test_unary_function(absval,n,gen,ndist);

      ROL::Elementwise::Heaviside<RealT> theta;
      errorFlag += 2 * test_unary_function(theta,n,gen,ndist);

      ROL::Elementwise::Logarithm<RealT> log_e;
      errorFlag += 4 * test_unary_function(log_e,n,gen,udist);

      ROL::Elementwise::Power<RealT> cube(3);
      errorFlag += 8 * test_unary_function(cube,n,gen,ndist);

      ROL::Elementwise::Reciprocal<RealT> recip;
      errorFlag += 16 * test_unary_function(recip,n,gen,udist);

      ROL::Elementwise::Round<RealT> round;
      errorFlag += 32 * test_unary_function(round,n,gen,ndist);

      ROL::Elementwise::Shift<RealT> shift(3);
      errorFlag += 64 * test_unary_function(shift,n,gen,ndist);

      ROL::Elementwise::Sign<RealT> sign;
      errorFlag += 128 * test_unary_function(sign,n,gen,ndist);

      ROL::Elementwise::SquareRoot<RealT> sqrt;
      errorFlag += 256 * test_unary_function(sqrt,n,gen,udist);

      ROL::Elementwise::ThresholdLower<RealT> threshl(0);
      errorFlag += 512 * test_unary_function(threshl,n,gen,ndist);

      ROL::Elementwise::ThresholdLower<RealT> threshu(0);
      errorFlag += 1024 * test_unary_function(threshu,n,gen,ndist);
      
    } catch( std::exception& err ) {
      *os << err.what() << std::endl;
      errorFlag = -2048;
    }; // end try

    if( errorFlag != 0 ) {
      std::cout << "End Result: TEST FAILED\n";
      *os << "Error Code: " << errorFlag << std::endl; 
    }
    else                 std::cout << "End Result: TEST PASSED\n";
  } Kokkos::finalize();

  return 0;
}

#include "ROL_Types.hpp"
#include "ROL_KokkosVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using RealT   = float;
using VectorT = ROL::KokkosVector<RealT,Kokkos::Cuda>;

template<typename Generator, typename Distribution> 
int test_binary_function( const ROL::Elementwise::BinaryFunction<RealT>& bf, int n, 
                          Generator& gen, Distribution& dist ) {

  auto errtol = std::sqrt(ROL::ROL_EPSILON<RealT>()*n);

  ROL::KokkosVector<RealT,Kokkos::Cuda> xkv(n), ykv(n), ekv(n);
  ROL::StdVector<RealT> xsv(n), ysv(n), esv(n);

  ekv.setScalar(1);
  esv.setScalar(1);
  
  auto xkv_host = xkv.view_host();
  auto xsv_ptr  = xsv.getVector();

  auto ykv_host = ykv.view_host();
  auto ysv_ptr  = ysv.getVector();

  for( int i=0; i<n; ++i ) {
    auto xval = dist(gen);  
    auto yval = dist(gen);  
    xkv_host(i) = xval;
    ykv_host(i) = yval;
    (*xsv_ptr)[i] = xval;
    (*ysv_ptr)[i] = yval;
  }
  
  xkv.modify_host();
  xkv.sync_device();
  ykv.modify_host();
  ykv.sync_device();

  ykv.applyBinary(bf,xkv);
  ysv.applyBinary(bf,xsv);

  auto sum_error = std::abs(ekv.dot(ykv) - esv.dot(ysv));
  auto norm_error = std::abs(ykv.norm() - ysv.norm());
  
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

      ROL::Elementwise::Axpy<RealT> axpy(std::sqrt(2));
      errorFlag += test_binary_function(axpy,n,gen,ndist);

      ROL::Elementwise::Aypx<RealT> aypx(std::sqrt(3));
      errorFlag += 2 * test_binary_function(aypx,n,gen,ndist);
 
      ROL::Elementwise::Divide<RealT> div;
      errorFlag += 4 * test_binary_function(div,n,gen,udist);

      ROL::Elementwise::DivideAndInvert<RealT> idiv;
      errorFlag += 8 * test_binary_function(idiv,n,gen,udist);

      ROL::Elementwise::Greater<RealT> grtr;
      errorFlag += 16 * test_binary_function(grtr,n,gen,ndist);

      ROL::Elementwise::Lesser<RealT> lsr;
      errorFlag += 32 * test_binary_function(lsr,n,gen,ndist);

      ROL::Elementwise::Max<RealT> max;
      errorFlag += 64 * test_binary_function(max,n,gen,ndist);

      ROL::Elementwise::Min<RealT> min;
      errorFlag += 128 * test_binary_function(min,n,gen,ndist);

      ROL::Elementwise::Plus<RealT> pls;
      errorFlag += 256 * test_binary_function(pls,n,gen,ndist);

      ROL::Elementwise::Set<RealT> st;
      errorFlag += 512 * test_binary_function(st,n,gen,ndist);

    } catch( std::exception& err ) {
      *os << err.what() << std::endl;
      errorFlag = -1024;
    }; // end try

    if( errorFlag != 0 ) {
      std::cout << "End Result: TEST FAILED\n";
      *os << "Error Code: " << errorFlag << std::endl; 
    }
    else                 std::cout << "End Result: TEST PASSED\n";
  } Kokkos::finalize();

  return 0;
}

#include "ROL_Types.hpp"
#include "ROL_KokkosVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using IntT    = int;
using RealT   = float;
using VectorT = ROL::KokkosVector<RealT,Kokkos::Cuda>;


template<typename ScalarT, typename Generator, typename Distribution> 
int test_reduction_op( ScalarT, const ROL::Elementwise::ReductionOp<ScalarT>& r, int n, 
                       Generator& gen, Distribution& dist ) {

  auto errtol = std::sqrt(ROL::ROL_EPSILON<ScalarT>()*n);

  ROL::KokkosVector<ScalarT,Kokkos::Cuda> xkv(n);
  ROL::StdVector<ScalarT> xsv(n);

  auto xkv_host = xkv.view_host();
  auto xsv_ptr  = xsv.getVector();

  for( int i=0; i<n; ++i ) {
    auto rval = dist(gen);  
    xkv_host(i) = rval;
    (*xsv_ptr)[i] = rval;
  }
  
  xkv.modify_host();
  xkv.sync_device();

  ScalarT result_kv = xkv.reduce(r);
  ScalarT result_sv = xsv.reduce(r);

  auto error  = std::abs(result_kv - result_sv);

  int errorFlag = error > errtol;

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
    std::normal_distribution<RealT> ndist(0,1);
    std::uniform_int_distribution<IntT> idist(-10,10);

    try {

      ROL::Elementwise::EuclideanNormSquared<RealT> enorm2;
      errorFlag += test_reduction_op(RealT{}, enorm2,n,gen,ndist);

//      ROL::Elementwise::ReductionAnd<IntT> rand;
//      errorFlag += 2 * test_reduction_op(IntT{}, rand,n,gen,idist);

      ROL::Elementwise::ReductionMax<RealT> rmax;
      errorFlag += 4 * test_reduction_op(RealT{}, rmax,n,gen,ndist);

      ROL::Elementwise::ReductionMin<RealT> rmin;
      errorFlag += 8 * test_reduction_op(RealT{}, rmin,n,gen,ndist);

      ROL::Elementwise::ReductionSum<RealT> rsum;
      errorFlag += 16 * test_reduction_op(RealT{}, rsum,n,gen,ndist);
     
    } catch( std::exception& err ) {
      *os << err.what() << std::endl;
      errorFlag = -32;
    }; // end try

    if( errorFlag != 0 ) {
      std::cout << "End Result: TEST FAILED\n";
      *os << "Error Code: " << errorFlag << std::endl; 
    }
    else                 std::cout << "End Result: TEST PASSED\n";
  } Kokkos::finalize();

  return 0;
}


#include <sstream>
#include <stdexcept>

template< class Device >
struct UnitTestReduceSum {

  typedef Device device_type ;
  typedef size_t value_type ;

  KOKKOS_DEVICE_FUNCTION
  static void init( value_type & result )
  { result = 0 ; }

  KOKKOS_DEVICE_FUNCTION
  static void join( value_type & result , const value_type & source )
  { result += source ; }

  KOKKOS_DEVICE_FUNCTION
  void operator()( const size_t iwork , value_type & result ) const
  { result += 1 ; }


  UnitTestReduceSum() {}

  static void test( const size_t work_count )
  {
    const size_t result =
      Kokkos::parallel_reduce( work_count , UnitTestReduceSum() );

    if ( result != work_count ) {
      std::ostringstream msg ;
      msg << "UnitTestReduceSum::test( "
          << work_count 
          << " ) FAILED : result = "
          << result ;
      throw std::runtime_error( msg.str() );
    }
  }
};



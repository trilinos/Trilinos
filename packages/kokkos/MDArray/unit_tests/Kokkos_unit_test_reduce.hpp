
#include <sstream>
#include <stdexcept>

template< class Device >
struct UnitTestReduceSum {

  typedef Device device_type ;

  typedef struct ValueType {
    size_t count ;
    size_t weighted ;
  } value_type ;

  KOKKOS_DEVICE_FUNCTION
  static void init( value_type & result )
  { result.count = 0 ; result.weighted = 0 ; }

  KOKKOS_DEVICE_FUNCTION
  static void join( value_type & result , const value_type & source )
  { result.count += source.count ; result.weighted += source.weighted ; }

  KOKKOS_DEVICE_FUNCTION
  void operator()( const size_t iwork , value_type & result ) const
  { result.count += 1 ; result.weighted += iwork + 1 ; }


  UnitTestReduceSum() {}

  static void test( const size_t work_count )
  {
    const value_type result =
      Kokkos::parallel_reduce( work_count , UnitTestReduceSum() );
    const size_t weighted = ( work_count * ( work_count + 1 ) ) / 2 ;

    if ( result.count != work_count || result.weighted != weighted ) {
      std::ostringstream msg ;
      msg << "UnitTestReduceSum::test( "
          << work_count 
          << " weighted[" << weighted
          << "] ) FAILED : result = "
          << result.count << " , " << result.weighted ;
      throw std::runtime_error( msg.str() );
    }
  }
};



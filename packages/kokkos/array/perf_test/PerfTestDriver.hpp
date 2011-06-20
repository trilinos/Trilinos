
#include <impl/Kokkos_Preprocessing_macros.hpp>

//------------------------------------------------------------------------

namespace Test {

enum { NUMBER_OF_TRIALS = 5 };

template< class DeviceType >
void run_test_hexgrad( int exp_beg , int exp_end );

template< class DeviceType >
void run_test_gramschmidt( int exp_beg , int exp_end );

template<>
void run_test_hexgrad< KOKKOS_MACRO_DEVICE >( int exp_beg , int exp_end )
{
  std::string label_hexgrad ;
  label_hexgrad.append( "\"HexGrad< double , " );
  label_hexgrad.append( KOKKOS_MACRO_TO_STRING( KOKKOS_MACRO_DEVICE ) );
  label_hexgrad.append( " >\"" );

  try {
    for (int i = exp_beg ; i < exp_end ; ++i) {
      double min_seconds = 0.0 ;
      double max_seconds = 0.0 ;
      double avg_seconds = 0.0 ;

      const int parallel_work_length = 1<<i;

      for ( int j = 0 ; j < NUMBER_OF_TRIALS ; ++j ) {
        const double seconds = HexGrad< double , KOKKOS_MACRO_DEVICE >::test(parallel_work_length) ;

        if ( 0 == j ) {
          min_seconds = seconds ;
          max_seconds = seconds ;
        }
        else {
          if ( seconds < min_seconds ) min_seconds = seconds ;
          if ( seconds > max_seconds ) max_seconds = seconds ;
        }
        avg_seconds += seconds ;
      }
      avg_seconds /= NUMBER_OF_TRIALS ;

      std::cout << label_hexgrad
                << " , " << parallel_work_length
                << " , " << min_seconds
                << " , " << ( min_seconds / parallel_work_length )
                << std::endl ;
    }
    std::cout << "PASSED : " << label_hexgrad << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "FAILED : " << label_hexgrad << " : " << x.what() << std::endl ;
  }
}

template<>
void run_test_gramschmidt< KOKKOS_MACRO_DEVICE >( int exp_beg , int exp_end )
{
  std::string label_gramschmidt ;
  label_gramschmidt.append( "\"GramSchmidt< double , " );
  label_gramschmidt.append( KOKKOS_MACRO_TO_STRING( KOKKOS_MACRO_DEVICE ) );
  label_gramschmidt.append( " >\"" );

  try {
    for (int i = exp_beg ; i < exp_end ; ++i) {
      double min_seconds = 0.0 ;
      double max_seconds = 0.0 ;
      double avg_seconds = 0.0 ;

      const int parallel_work_length = 1<<i;

      for ( int j = 0 ; j < NUMBER_OF_TRIALS ; ++j ) {
        const double seconds = ModifiedGramSchmidt< double , KOKKOS_MACRO_DEVICE >::test(parallel_work_length, 32 ) ;

        if ( 0 == j ) {
          min_seconds = seconds ;
          max_seconds = seconds ;
        }
        else {
          if ( seconds < min_seconds ) min_seconds = seconds ;
          if ( seconds > max_seconds ) max_seconds = seconds ;
        }
        avg_seconds += seconds ;
      }
      avg_seconds /= NUMBER_OF_TRIALS ;

      std::cout << label_gramschmidt
                << " , " << parallel_work_length
                << " , " << min_seconds
                << " , " << ( min_seconds / parallel_work_length )
                << std::endl ;
    }
    std::cout << "PASSED : " << label_gramschmidt << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "FAILED : " << label_gramschmidt << " : " << x.what() << std::endl ;
  }
}

}


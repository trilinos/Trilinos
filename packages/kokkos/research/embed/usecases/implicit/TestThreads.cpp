#include <Kokkos_Macros.hpp>
#include <ParallelComm.hpp>

#include <sstream>
#include <Kokkos_Core.hpp>

#include <TestImplicit.hpp>

namespace Test {

int test_host( comm::Machine machine , std::istream & input )
{
  const unsigned numa_count       = Kokkos::hwloc::get_available_numa_count();
  const unsigned cores_per_numa   = Kokkos::hwloc::get_available_cores_per_numa();
  const unsigned threads_per_core = Kokkos::hwloc::get_available_threads_per_core();

  unsigned team_count = numa_count ;
  unsigned threads_per_team = cores_per_numa * threads_per_core ;

  unsigned elem_beg = 3 ;
  unsigned elem_end = 4 ;
  unsigned run = 1 ;

  while ( ! input.eof() ) {
    std::string which ;

    input >> which ;

    if ( which == std::string("team") ) {
      input >> team_count ;
      input >> threads_per_team ;
    }
    else if ( which == std::string("implicit") ) {
      input >> elem_beg ;
      input >> elem_end ;
      input >> run ;
    }
    else {
      std::cerr << "Expected \"gang #Gang #Worker\" OR \"implicit #ElemBeg #ElemEnd #Run\""
                << std::endl ;
      return -1 ;
    }
  }

  if ( 0 == comm::rank( machine ) ) {
    std::cout << "\"P" << comm::rank( machine )
              << ": hwloc[ " << numa_count
              << " x " << cores_per_numa
              << " x " << threads_per_core
              << " ] Threads[ " << team_count
              << " x " << threads_per_team
              << " ]\"" << std::endl ;
  }

  Kokkos::Threads::initialize( team_count * threads_per_team , numa_count );

  {
    std::ostringstream label ;

    label << "Scalar, Threads[" << team_count << "x" << threads_per_team << "]" ;

    implicit_driver<double,Kokkos::Threads>(
      label.str().c_str() , machine , team_count , elem_beg , elem_end , run );
  }

  {
    std::ostringstream label ;

    label << "Ensemble[32], Threads[" << team_count << "x" << threads_per_team << "]" ;

    implicit_driver< Kokkos::Array<double,32> , Kokkos::Threads>(
      label.str().c_str() , machine , team_count , elem_beg , elem_end , run );
  }

  Kokkos::Threads::finalize();

  return 0 ;
}

}


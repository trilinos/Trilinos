
#include <sstream>
#include <KokkosArray_Host.hpp>
#include <KokkosArray_hwloc.hpp>

#include <ParallelComm.hpp>
#include <TestImplicit.hpp>

namespace Test {

int test_host( comm::Machine machine , std::istream & input )
{
  const std::pair<unsigned,unsigned> core_top  = KokkosArray::hwloc::get_core_topology();
  const unsigned                     core_size = KokkosArray::hwloc::get_core_capacity();

  std::pair<unsigned,unsigned> gang_top( core_top.first , core_top.second * core_size );

  unsigned elem_beg = 3 ;
  unsigned elem_end = 4 ;
  unsigned run = 1 ;

  while ( ! input.eof() ) {
    std::string which ;

    input >> which ;

    if ( which == std::string("gang") ) {
      input >> gang_top.first ;
      input >> gang_top.second ;
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
              << ": hwloc[ " << core_top.first
              << " x " << core_top.second
              << " x " << core_size
              << " ] Host[ " << gang_top.first
              << " x " << gang_top.second
              << " ]\"" << std::endl ;
  }

  KokkosArray::Host::initialize( gang_top , core_top );

  {
    std::ostringstream label ;

    label << "Scalar, Host[" << gang_top.first << "x" << gang_top.second << "]" ;

    implicit_driver<double,KokkosArray::Host>(
      label.str().c_str() , machine , gang_top.first , elem_beg , elem_end , run );
  }

  {
    std::ostringstream label ;

    label << "Ensemble[32], Host[" << gang_top.first << "x" << gang_top.second << "]" ;

    implicit_driver< KokkosArray::Array<double,32> , KokkosArray::Host>(
      label.str().c_str() , machine , gang_top.first , elem_beg , elem_end , run );
  }

  KokkosArray::Host::finalize();

  return 0 ;
}

}


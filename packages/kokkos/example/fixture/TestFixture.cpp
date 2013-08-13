
#include <utility>
#include <iostream>

#include <Kokkos_Threads.hpp>
#include <Kokkos_Cuda.hpp>
#include <BoxElemPart.hpp>
#include <BoxElemFixture.hpp>

namespace Kokkos {
namespace Example {

template< class Device > void test_fixture();


#if defined( __CUDACC__ )
typedef Kokkos::Cuda UseDevice ;
#else
typedef Kokkos::Threads UseDevice ;
#endif

template< class Device >
struct FixtureVerifyElemNodeCoord
{
  typedef Device device_type ;

  typedef struct { size_t success , error ; } value_type ;

  typedef Kokkos::Example::BoxElemFixture< Device , Kokkos::Example::BoxElemPart::ElemLinear > FixtureType ;

  FixtureType m_fixture ;

  KOKKOS_INLINE_FUNCTION
  void init( value_type & update ) const { update.success = update.error = 0 ; }

  KOKKOS_INLINE_FUNCTION
  void join( volatile       value_type & update ,
             volatile const value_type & input ) const
    {
      update.success += input.success ;
      update.error += input.error ;
    }
  

  KOKKOS_INLINE_FUNCTION
  void operator()( size_t ielem , value_type & update ) const
  {
    unsigned node_coord[ FixtureType::ElemNode ][3] ;

    for ( unsigned i = 0 ; i < FixtureType::ElemNode ; ++i ) {
      const unsigned node_id = m_fixture.elem_node(ielem,i);
      node_coord[i][0] = m_fixture.node_grid(node_id,0);
      node_coord[i][1] = m_fixture.node_grid(node_id,1);
      node_coord[i][2] = m_fixture.node_grid(node_id,2);
    }

    int error = 0 ;
    for ( unsigned i = 1 ; i < FixtureType::ElemNode ; ++i ) {
      if ( node_coord[0][0] + m_fixture.elem_node_local(i,0) != node_coord[i][0] ||
           node_coord[0][1] + m_fixture.elem_node_local(i,1) != node_coord[i][1] ||
           node_coord[0][2] + m_fixture.elem_node_local(i,2) != node_coord[i][2] ) {
        error = 1 ;
      }
    }

    if ( error ) {
      ++update.error ;
    }
    else {
      ++update.success ;
    }
  }

  FixtureVerifyElemNodeCoord( const FixtureType & f ) : m_fixture(f) {}
};


template<>
void test_fixture< UseDevice >()
{
  typedef Kokkos::Example::BoxElemFixture< UseDevice , Kokkos::Example::BoxElemPart::ElemLinear > FixtureType ;

  const Kokkos::Example::BoxElemPart::Decompose
    decompose = Kokkos::Example::BoxElemPart:: DecomposeElem ; // DecomposeElem | DecomposeNode ;

  const unsigned global_size = 256 ;
  const unsigned global_nx = 100 ;
  const unsigned global_ny = 120 ;
  const unsigned global_nz = 140 ;

  for ( unsigned my_rank = 0 ; my_rank < global_size ; ++my_rank ) {

    const FixtureType fixture( decompose , global_size , my_rank , global_nx , global_ny , global_nz );

    // Verify grid coordinates of element's nodes
    
    FixtureVerifyElemNodeCoord<UseDevice>::value_type result = { 0 , 0 };

    Kokkos::parallel_reduce( fixture.elem_node().dimension_0() , FixtureVerifyElemNodeCoord<UseDevice>( fixture ) , result );

    if ( result.error ) {
      std::cout << "P[" << my_rank << ":" << global_size
                << "] Fixture elem_node_coord"
                << " success(" << result.success << ")"
                << " error(" << result.error << ")"
                << std::endl ;
    }

    // Check send/recv alignment


  }
}


} /* namespace Example */
} /* namespace Kokkos */


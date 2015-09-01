#pragma once
#ifndef __TEAM_VIEW_HPP__
#define __TEAM_VIEW_HPP__

/// \file team_view.hpp
/// \brief Team view is inherited from matrix view and typedef of team policy type
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;

  template<typename MatViewType,
           typename TeamFactoryType>
  class TeamView : public MatViewType {
  public:
    typedef typename MatViewType::value_type   value_type;
    typedef typename MatViewType::ordinal_type ordinal_type;

    typedef TeamFactoryType team_factory_type;
    typedef typename team_factory_type::policy_type policy_type;
  };
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
  namespace Impl {

    //  The Kokkos::View allocation will by default assign each allocated datum to zero.
    //  This is not the required initialization behavior when
    //  non-trivial objects are used within a Kokkos::View.
    //  Create a partial specialization of the Kokkos::Impl::AViewDefaultConstruct
    //  to replace the assignment initialization with placement new initialization.
    //
    //  This work-around is necessary until a TBD design refactorization of Kokkos::View.

    template< class ExecSpace , typename T1, typename T2 >
    struct ViewDefaultConstruct< ExecSpace , Example::TeamView<T1,T2> , true >
    {
      typedef Example::TeamView<T1,T2> type ;
      type * const m_ptr ;

      KOKKOS_FORCEINLINE_FUNCTION
      void operator()( const typename ExecSpace::size_type& i ) const
      { new(m_ptr+i) type(); }

      ViewDefaultConstruct( type * pointer , size_t capacity )
        : m_ptr( pointer )
      {
        Kokkos::RangePolicy< ExecSpace > range( 0 , capacity );
        parallel_for( range , *this );
        ExecSpace::fence();
      }
    };

  } // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

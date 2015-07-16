#pragma once
#ifndef __EXAMPLE_KOKKOS_TEAM_BARRIER_HPP__
#define __EXAMPLE_KOKKOS_TEAM_BARRIER_HPP__

#define CACHELINE_SIZE 64
#define BARRIER_MASK 0x01010101

namespace Kokkos { 

  namespace Experimental {
      
    // CoreBarrier : team_base
    // =======================
    class SimpleCoreBarrier {
    private:
      volatile int m_arrive;
      volatile int m_depart;
      int m_padding[CACHELINE_SIZE/sizeof(int)-2];

    public:
      SimpleCoreBarrier() 
        : m_arrive(0),
          m_depart(0) {}

      KOKKOS_INLINE_FUNCTION 
      void 
      reset() { m_arrive = 0; m_depart = 0; }

      KOKKOS_INLINE_FUNCTION 
      void 
      set_arrive(const int team_rank) { ((char*)&m_arrive)[team_rank] = 1; }

      KOKKOS_INLINE_FUNCTION 
      void 
      set_depart(const int team_rank) { ((char*)&m_depart)[team_rank] = 1; }
        
      KOKKOS_INLINE_FUNCTION 
      volatile int& 
      arrive()  { return m_arrive; }

      KOKKOS_INLINE_FUNCTION 
      volatile int& 
      depart()  { return m_depart; }

    };
    typedef class SimpleCoreBarrier SimpleCoreBarrierType;

    // TeamMember : test wrapper
    // =========================
    class SimpleTeamMember {
    private:
      int m_team_size;
      int m_team_rank;
      int m_team_rank_rev;

      SimpleCoreBarrierType *m_core_barrier;

    public:
      SimpleTeamMember(const int arg_team_size,
                       const int arg_team_rank,
                       const int arg_team_rank_rev,
                       SimpleCoreBarrierType *arg_core_barrier) 
        : m_team_size(arg_team_size),
          m_team_rank(arg_team_rank),
          m_team_rank_rev(arg_team_rank_rev),
          m_core_barrier(arg_core_barrier) {}

    public:

      KOKKOS_INLINE_FUNCTION 
      bool 
      team_fan_in() const {
        if (m_team_size != 1) {
          m_core_barrier->set_arrive(m_team_rank);
          Kokkos::Impl::spinwait(m_core_barrier->arrive(), BARRIER_MASK);
            
          m_core_barrier->set_depart(m_team_rank);          
          Kokkos::Impl::spinwait(m_core_barrier->depart(), BARRIER_MASK);
        }
        return !m_team_rank_rev;
      }
        
      KOKKOS_INLINE_FUNCTION void team_fan_out() const {
        if (m_team_size != 1) {          
          m_core_barrier->reset();
        }
      } 
      KOKKOS_INLINE_FUNCTION void team_barrier() const {
        team_fan_in();
        team_fan_out();
      }
        
    };
    typedef class SimpleTeamMember SimpleTeamMemberType;
  }
}

#endif

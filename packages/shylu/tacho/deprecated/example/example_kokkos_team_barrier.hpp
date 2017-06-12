#pragma once
#ifndef __EXAMPLE_KOKKOS_TEAM_BARRIER_HPP__
#define __EXAMPLE_KOKKOS_TEAM_BARRIER_HPP__

#define CACHELINE_SIZE 64
#define BARRIER_MASK 0x0101010101010101

namespace Kokkos { 

  namespace Experimental {
      
    // CoreBarrier : team_base
    // =======================
    class SimpleCoreBarrier {
    private:
      volatile int64_t m_arrive;
      volatile int64_t m_depart;
      int64_t m_padding[CACHELINE_SIZE/sizeof(int64_t)-2];

    public:
      SimpleCoreBarrier() 
        : m_arrive(0),
          m_depart(0) {}

      KOKKOS_INLINE_FUNCTION 
      int
      set_arrive(const int team_rank) { 
        const int flip = !((char*)&m_arrive)[team_rank];
        ((char*)&m_arrive)[team_rank] = flip;
        return flip;
      }

      KOKKOS_INLINE_FUNCTION 
      int
      set_depart(const int team_rank) { 
        const int flip = !((char*)&m_depart)[team_rank];
        ((char*)&m_depart)[team_rank] = flip; 
        return flip;
      }
        
      KOKKOS_INLINE_FUNCTION 
      int64_t
      arrive() const { return m_arrive; }

      KOKKOS_INLINE_FUNCTION 
      int64_t
      depart() const { return m_depart; }

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
          const int flip = m_core_barrier->set_arrive(m_team_rank);
          const int64_t mask = static_cast<int64_t>(flip*BARRIER_MASK) 
            >> CHAR_BIT*(sizeof(int64_t) - m_team_size);
          while (m_core_barrier->arrive() != mask); 
        }
        return !m_team_rank_rev;
      }
        
      KOKKOS_INLINE_FUNCTION void team_fan_out() const {
        if (m_team_size != 1) {  
          const int flip = m_core_barrier->set_depart(m_team_rank);
          const int64_t mask = static_cast<int64_t>(flip*BARRIER_MASK) 
            >> CHAR_BIT*(sizeof(int64_t) - m_team_size);
          while (m_core_barrier->depart() != mask);  
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

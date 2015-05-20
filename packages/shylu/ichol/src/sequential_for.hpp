#pragma once
#ifndef __SEQUENTIAL_FOR_HPP__
#define __SEQUENTIAL_FOR_HPP__

/// \file sequential_for.hpp
/// \brief A wrapper of "for" loop.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;

  // Team thread communicator
  // ========================
  class TeamThreadMember {
  public:
    TeamThreadMember() { }
    int team_rank() const { return 0; } 
    int team_size() const { return 1; } 
    void team_barrier() const { }
  };

  // Team policy
  // ===========
  class TeamPolicy {
  public:
    typedef class TeamThreadMember member_type;
    static member_type member_null() { return member_type(); }
  };

  // Team thread loop region for sequential loop
  // ===========================================
  template<typename OrdinalType, typename MemberType>
  class TeamThreadLoopRegion {
  public:
    typedef OrdinalType ordinal_type;
    typedef MemberType member_type;
    
  private:
    const member_type &_thread;
    const ordinal_type _begin, _end;

  public:
    TeamThreadLoopRegion(const member_type  &thread,
                         const ordinal_type &begin, 
                         const ordinal_type &end) 
      : _thread(thread),
        _begin(begin),
        _end(end) { }

    ordinal_type Begin() const { return _begin; }
    ordinal_type End()   const { return _end;   }
  };

  // Sequential loop region generator
  // ================================
  template<typename OrdinalType, typename MemberType>
  TeamThreadLoopRegion<OrdinalType,MemberType> TeamThreadLoop(const MemberType &thread,
                                                              const OrdinalType &begin,
                                                              const OrdinalType &end) {
    return TeamThreadLoopRegion<OrdinalType,MemberType>(thread, begin, end);
  }

  // Sequential for loop (debugging use only)
  // =======================================
  class SequentialFor {
  public:
    template<typename TeamThreadLoopRegionType, typename LambdaType>
    SequentialFor(const TeamThreadLoopRegionType &loop, 
                  const LambdaType &lambda) {
      for (auto i=loop.Begin();i<loop.End();++i) 
        lambda(i);
    }
  };
  
}

#endif

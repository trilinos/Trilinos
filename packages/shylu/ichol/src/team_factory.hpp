#pragma once
#ifndef __TEAM_FACTORY_HPP__
#define __TEAM_FACTORY_HPP__

/// \file task_factory.hpp
/// \brief A wrapper for task policy and future with a provided space type.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;

  template<typename PolicyType,
           template<typename,typename> class ThreadLoopRegionType>
  class TeamFactory {
  public:
    typedef PolicyType policy_type;
    
    template<typename OrdinalType, typename MemberType>
    static
    ThreadLoopRegionType<OrdinalType,MemberType> 
    createThreadLoopRegion(const MemberType  &thread, 
                           const OrdinalType &begin,
                           const OrdinalType &end) {
      return ThreadLoopRegionType<OrdinalType,MemberType>(thread, begin, end);
    }

  };

}

#endif

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

#endif

#pragma once
#ifndef __CRS_TEAM_VIEW_HPP__
#define __CRS_TEAM_VIEW_HPP__

/// \file crs_team_view.hpp
/// \brief Team view is inherited from matrix view and typedef of team policy type
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;

  template<typename CrsMatBaseType,
           typename TeamFactoryType>
  class CrsTeamView : public CrsMatrixView<CrsMatBaseType> {
  public:
    typedef typename CrsMatBaseType::value_type   value_type;
    typedef typename CrsMatBaseType::ordinal_type ordinal_type;

    typedef TeamFactoryType team_factory_type;
    typedef typename team_factory_type::policy_type policy_type;

    using CrsMatrixView<CrsMatBaseType>::CrsMatrixView;
  };
}

#endif

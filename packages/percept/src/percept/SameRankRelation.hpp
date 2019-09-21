// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_SameRankRelation_hpp
#define percept_SameRankRelation_hpp

//#include <percept/PerceptMesh.hpp>

#include <vector>
//#include <algorithm>
//#include <iostream>

/// define only one of these to be 1
#define SAME_RANK_RELATION_MAP_TYPE_BOOST 1
#define SAME_RANK_RELATION_MAP_TYPE_TR1 0
#define SAME_RANK_RELATION_MAP_TYPE_STD 0

#if SAME_RANK_RELATION_MAP_TYPE_BOOST
#include <boost/unordered_map.hpp>
#endif

#if SAME_RANK_RELATION_MAP_TYPE_STD
#include <map>
#endif

#if SAME_RANK_RELATION_MAP_TYPE_TR1
#include <tr1/unordered_map>
#endif

  namespace percept {


    typedef std::vector<stk::mesh::Entity> PerceptEntityVector;
    typedef stk::mesh::Entity SameRankRelationKey;
    typedef PerceptEntityVector SameRankRelationValue;
    typedef boost::unordered_map<SameRankRelationKey, SameRankRelationValue > SameRankRelation;


  }

#endif

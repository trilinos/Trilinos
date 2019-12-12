// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_UniformRefiner_hpp
#define adapt_UniformRefiner_hpp

#include <adapt/Refiner.hpp>

  namespace percept {

    //========================================================================================================================
    //========================================================================================================================
    //========================================================================================================================
    //template<class UniformRefinerPattern>
    class UniformRefiner : public Refiner
    {
    public:
      UniformRefiner(percept::PerceptMesh& eMesh, UniformRefinerPatternBase & bp, stk::mesh::FieldBase *proc_rank_field=0);

    protected:
    };



  }

#endif

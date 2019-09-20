// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_ElementRefinePredicate_hpp
#define adapt_ElementRefinePredicate_hpp

#include <adapt/IElementBasedAdapterPredicate.hpp>

  namespace percept {

    // Can be instantiated by the user, or used to define your own
    struct ElementRefinePredicate : public IElementBasedAdapterPredicate {

      ElementRefinePredicate(PerceptMesh& eMesh, stk::mesh::Selector* selector=0, stk::mesh::FieldBase *field=0, double tolerance=0.0) :
        IElementBasedAdapterPredicate(eMesh, selector, field, tolerance)
      {
      }

      /// Return DO_REFINE, DO_UNREFINE, DO_NOTHING
      int operator()(const stk::mesh::Entity entity);
    };

  }

#endif

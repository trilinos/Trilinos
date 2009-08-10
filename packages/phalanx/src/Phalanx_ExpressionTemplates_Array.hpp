// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef PHALANX_EXPRESSION_TEMPLATES_ARRAY_HPP
#define PHALANX_EXPRESSION_TEMPLATES_ARRAY_HPP

#include <cassert>
#include <vector>
#include "Phalanx_ConfigDefs.hpp"

#include "Phalanx_ExpressionTemplates_ExprScalar.hpp"

namespace PHX {

  template<typename Ordinal, typename Scalar, 
	   typename Rep = std::vector<Scalar> >
  class ExprArray {

  private:

    Rep expr_rep;

  public:

    explicit ExprArray(Ordinal size) : expr_rep(size) {}

    ExprArray(Rep const& r) : expr_rep(r) {}

    ExprArray& operator=(ExprArray const& b) {
      assert(this->size() == b.size());
      for (Ordinal i=0; i < b.size(); ++i)
	expr_rep[i] = b[i];
      return *this;
    }

    template<typename Ordinal2, typename Scalar2, typename Rep2>
    ExprArray& operator=(ExprArray<Ordinal2,Scalar2,Rep2> const& b) {
      assert(this->size() == b.size());
      for (Ordinal i=0; i < b.size(); ++i)
	expr_rep[i] = b[i];
      return *this;
    }

    Ordinal size() const { return expr_rep.size(); }

    Scalar operator[](Ordinal idx) const {
      assert(idx<this->size());
      return expr_rep[idx];
    }

    Scalar& operator[](Ordinal idx) {
      assert(idx<this->size());
      return expr_rep[idx];
    }

    Rep const& rep() const {
      return expr_rep;
    }

    Rep rep() {
      return expr_rep;
    }

  };

}

#endif

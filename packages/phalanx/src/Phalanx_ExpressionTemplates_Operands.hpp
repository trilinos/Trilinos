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

#ifndef PHALANX_EXPRESSION_TEMPLATES_OPERANDS_HPP
#define PHALANX_EXPRESSION_TEMPLATES_OPERANDS_HPP

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_ExpressionTemplates_Traits.hpp"

namespace PHX {

  //! Addition
  template<typename Ordinal, typename Scalar, typename OP1, typename OP2>
  class ExprAdd {
    
  private:

    typename PHX::ExprTraits<Ordinal,OP1>::ExprRef op1;
    typename PHX::ExprTraits<Ordinal,OP2>::ExprRef op2;
    
  public:
    
    ExprAdd(OP1 const& a, OP2 const& b) : op1(a), op2(b) {};
    
    Scalar operator[] (Ordinal idx) const { return op1[idx] + op2[idx]; }
    
    Ordinal size() const 
    {
      return op1.size() != 0 ? op1.size() : op2.size(); 
    } 
    
  };

  //! Subtraction
  template<typename Ordinal, typename Scalar, typename OP1, typename OP2>
  class ExprSubtr {
    
  private:

    typename PHX::ExprTraits<Ordinal,OP1>::ExprRef op1;
    typename PHX::ExprTraits<Ordinal,OP2>::ExprRef op2;
    
  public:
    
    ExprSubtr(OP1 const& a, OP2 const& b) : op1(a), op2(b) {};
    
    Scalar operator[] (Ordinal idx) const { return op1[idx] - op2[idx]; }
    
    Ordinal size() const 
    {
      return op1.size() != 0 ? op1.size() : op2.size(); 
    } 
    
  };

  //! Multiplication
  template<typename Ordinal, typename Scalar, typename OP1, typename OP2>
  class ExprMult {


  private:
    typename PHX::ExprTraits<Ordinal,OP1>::ExprRef op1;
    typename PHX::ExprTraits<Ordinal,OP2>::ExprRef op2;

  public:
    
    ExprMult(OP1 const& a, OP2 const& b) : op1(a), op2(b) {}
    
    Scalar operator[] (Ordinal idx) const { return op1[idx] * op2[idx]; }

    Ordinal size() const 
    {
      return op1.size() != 0 ? op1.size() : op2.size(); 
    }   
  };
  
  //! Division
  template<typename Ordinal, typename Scalar, typename OP1, typename OP2>
  class ExprDiv {


  private:
    typename PHX::ExprTraits<Ordinal,OP1>::ExprRef op1;
    typename PHX::ExprTraits<Ordinal,OP2>::ExprRef op2;

  public:
    
    ExprDiv(OP1 const& a, OP2 const& b) : op1(a), op2(b) {}
    
    Scalar operator[] (Ordinal idx) const { return op1[idx] / op2[idx]; }

    Ordinal size() const 
    {
      return op1.size() != 0 ? op1.size() : op2.size(); 
    }   
  };
  

}

#endif

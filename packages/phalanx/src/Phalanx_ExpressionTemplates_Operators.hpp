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

#ifndef PHALANX_EXPRESSION_TEMPLATES_OPERATORS_HPP
#define PHALANX_EXPRESSION_TEMPLATES_OPERATORS_HPP

#include <cassert>
#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_ExpressionTemplates_Operands.hpp"
#include "Phalanx_ExpressionTemplates_Array.hpp"

namespace PHX {

  // ***********
  // Addition
  // ***********

  //! ExprArray + ExprArray
  template<typename Ordinal, typename Scalar, typename R1, typename R2>
  PHX::ExprArray<Ordinal, Scalar, PHX::ExprAdd<Ordinal,Scalar,R1,R2> > 
  operator+(PHX::ExprArray<Ordinal,Scalar,R1> const& a, 
	    PHX::ExprArray<Ordinal,Scalar,R2> const& b) {
    return 
      PHX::ExprArray<Ordinal, Scalar, PHX::ExprAdd<Ordinal,Scalar,R1,R2> > 
      (PHX::ExprAdd<Ordinal,Scalar,R1,R2>(a.rep(), b.rep()));
  }    

  //! Scalar + ExprArray
  template<typename Ordinal, typename Scalar, typename R2>
  PHX::ExprArray<Ordinal, Scalar, PHX::ExprAdd<Ordinal,Scalar,PHX::ExprScalar<Ordinal,Scalar>,R2> > 
  operator+(Scalar const& s, PHX::ExprArray<Ordinal,Scalar,R2> const& b) {
    return 
      PHX::ExprArray<Ordinal, Scalar, PHX::ExprAdd<Ordinal,Scalar,PHX::ExprScalar<Ordinal,Scalar>,R2> > 
      (PHX::ExprAdd<Ordinal,Scalar,PHX::ExprScalar<Ordinal,Scalar>,R2>(ExprScalar<Ordinal,Scalar>(s), b.rep()));
  }    

  //! ExprArray + Scalar
  template<typename Ordinal, typename Scalar, typename R1>
  PHX::ExprArray<Ordinal, Scalar, PHX::ExprAdd<Ordinal,Scalar,R1,PHX::ExprScalar<Ordinal,Scalar> > > 
  operator+(PHX::ExprArray<Ordinal,Scalar,R1> const& a, Scalar const& s) {
    return 
      PHX::ExprArray<Ordinal, Scalar, PHX::ExprAdd<Ordinal,Scalar,R1,PHX::ExprScalar<Ordinal,Scalar> > > 
      (PHX::ExprAdd<Ordinal,Scalar,R1,PHX::ExprScalar<Ordinal,Scalar> >(a.rep(), ExprScalar<Ordinal,Scalar>(s)));
  }    

  // ***********
  // Subtration
  // ***********

  //! ExprArray - ExprArray
  template<typename Ordinal, typename Scalar, typename R1, typename R2>
  PHX::ExprArray<Ordinal, Scalar, PHX::ExprSubtr<Ordinal,Scalar,R1,R2> > 
  operator-(PHX::ExprArray<Ordinal,Scalar,R1> const& a, 
	    PHX::ExprArray<Ordinal,Scalar,R2> const& b) {
    return 
      PHX::ExprArray<Ordinal, Scalar, PHX::ExprSubtr<Ordinal,Scalar,R1,R2> > 
      (PHX::ExprSubtr<Ordinal,Scalar,R1,R2>(a.rep(), b.rep()));
  }    

  //! Scalar - ExprArray
  template<typename Ordinal, typename Scalar, typename R2>
  PHX::ExprArray<Ordinal, Scalar, PHX::ExprSubtr<Ordinal,Scalar,PHX::ExprScalar<Ordinal,Scalar>,R2> > 
  operator-(Scalar const& s, PHX::ExprArray<Ordinal,Scalar,R2> const& b) {
    return 
      PHX::ExprArray<Ordinal, Scalar, PHX::ExprSubtr<Ordinal,Scalar,PHX::ExprScalar<Ordinal,Scalar>,R2> > 
      (PHX::ExprSubtr<Ordinal,Scalar,PHX::ExprScalar<Ordinal,Scalar>,R2>(ExprScalar<Ordinal,Scalar>(s), b.rep()));
  }    

  //! ExprArray - Scalar
  template<typename Ordinal, typename Scalar, typename R1>
  PHX::ExprArray<Ordinal, Scalar, PHX::ExprSubtr<Ordinal,Scalar,R1,PHX::ExprScalar<Ordinal,Scalar> > > 
  operator-(PHX::ExprArray<Ordinal,Scalar,R1> const& a, Scalar const& s) {
    return 
      PHX::ExprArray<Ordinal, Scalar, PHX::ExprSubtr<Ordinal,Scalar,R1,PHX::ExprScalar<Ordinal,Scalar> > > 
      (PHX::ExprSubtr<Ordinal,Scalar,R1,PHX::ExprScalar<Ordinal,Scalar> >(a.rep(), ExprScalar<Ordinal,Scalar>(s)));
  }    

  // ***********
  // Multiplication
  // ***********

  //! ExprArray * ExprArray
  template<typename Ordinal, typename Scalar, typename R1, typename R2>
  PHX::ExprArray<Ordinal, Scalar, PHX::ExprMult<Ordinal,Scalar,R1,R2> > 
  operator*(PHX::ExprArray<Ordinal,Scalar,R1> const& a, 
	    PHX::ExprArray<Ordinal,Scalar,R2> const& b) {
    return 
      PHX::ExprArray<Ordinal, Scalar, PHX::ExprMult<Ordinal,Scalar,R1,R2> > 
      (PHX::ExprMult<Ordinal,Scalar,R1,R2>(a.rep(), b.rep()));
  }    

  //! Scalar * ExprArray
  template<typename Ordinal, typename Scalar, typename R2>
  PHX::ExprArray<Ordinal, Scalar, PHX::ExprMult<Ordinal,Scalar,PHX::ExprScalar<Ordinal,Scalar>,R2> > 
  operator*(Scalar const& s, PHX::ExprArray<Ordinal,Scalar,R2> const& b) {
    return 
      PHX::ExprArray<Ordinal, Scalar, PHX::ExprMult<Ordinal,Scalar,PHX::ExprScalar<Ordinal,Scalar>,R2> > 
      (PHX::ExprMult<Ordinal,Scalar,PHX::ExprScalar<Ordinal,Scalar>,R2>(ExprScalar<Ordinal,Scalar>(s), b.rep()));
  }    

  //! ExprArray * Scalar
  template<typename Ordinal, typename Scalar, typename R1>
  PHX::ExprArray<Ordinal, Scalar, PHX::ExprMult<Ordinal,Scalar,R1,PHX::ExprScalar<Ordinal,Scalar> > > 
  operator*(PHX::ExprArray<Ordinal,Scalar,R1> const& a, Scalar const& s) {
    return 
      PHX::ExprArray<Ordinal, Scalar, PHX::ExprMult<Ordinal,Scalar,R1,PHX::ExprScalar<Ordinal,Scalar> > > 
      (PHX::ExprMult<Ordinal,Scalar,R1,PHX::ExprScalar<Ordinal,Scalar> >(a.rep(), ExprScalar<Ordinal,Scalar>(s)));
  }    

  // ***********
  // Division
  // ***********

  //! ExprArray / ExprArray
  template<typename Ordinal, typename Scalar, typename R1, typename R2>
  PHX::ExprArray<Ordinal, Scalar, PHX::ExprDiv<Ordinal,Scalar,R1,R2> > 
  operator/(PHX::ExprArray<Ordinal,Scalar,R1> const& a, 
	    PHX::ExprArray<Ordinal,Scalar,R2> const& b) {
    return 
      PHX::ExprArray<Ordinal, Scalar, PHX::ExprDiv<Ordinal,Scalar,R1,R2> > 
      (PHX::ExprDiv<Ordinal,Scalar,R1,R2>(a.rep(), b.rep()));
  }    

  //! Scalar / ExprArray
  template<typename Ordinal, typename Scalar, typename R2>
  PHX::ExprArray<Ordinal, Scalar, PHX::ExprDiv<Ordinal,Scalar,PHX::ExprScalar<Ordinal,Scalar>,R2> > 
  operator/(Scalar const& s, PHX::ExprArray<Ordinal,Scalar,R2> const& b) {
    return 
      PHX::ExprArray<Ordinal, Scalar, PHX::ExprDiv<Ordinal,Scalar,PHX::ExprScalar<Ordinal,Scalar>,R2> > 
      (PHX::ExprDiv<Ordinal,Scalar,PHX::ExprScalar<Ordinal,Scalar>,R2>(ExprScalar<Ordinal,Scalar>(s), b.rep()));
  }    

  //! ExprArray / Scalar
  template<typename Ordinal, typename Scalar, typename R1>
  PHX::ExprArray<Ordinal, Scalar, PHX::ExprDiv<Ordinal,Scalar,R1,PHX::ExprScalar<Ordinal,Scalar> > > 
  operator/(PHX::ExprArray<Ordinal,Scalar,R1> const& a, Scalar const& s) {
    return 
      PHX::ExprArray<Ordinal, Scalar, PHX::ExprDiv<Ordinal,Scalar,R1,PHX::ExprScalar<Ordinal,Scalar> > > 
      (PHX::ExprDiv<Ordinal,Scalar,R1,PHX::ExprScalar<Ordinal,Scalar> >(a.rep(), ExprScalar<Ordinal,Scalar>(s)));
  }    

}

#endif

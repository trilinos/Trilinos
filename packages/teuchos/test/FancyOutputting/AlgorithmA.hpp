/*
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_VerboseObject.hpp"

#ifndef TEUCHOS_ALGORITHM_A_HPP
#define TEUCHOS_ALGORITHM_A_HPP


void someDumbFunction( std::ostream &out, const std::string &indentSpacer );


void someLessDumbFunction( std::ostream &out_arg );


// This is a typical numerical class that derives from VerboseObject and does
// outputting.  Note that the use of the OSTab class requires initialization
// using VerboseObject::getOSTab(...) which takes care of the hassles and is
// easy to use.
//
// This class also derives from ParameterListAcceptor and uses helper
// functio  ns to read options for VerboseObject from a parameter sublist.
class AlgorithmA
  : public Teuchos::VerboseObject<AlgorithmA>,
    public Teuchos::ParameterListAcceptor
{
public:

  // Constructor(s)

  AlgorithmA();

  // Overridden from ParameterListAccpetor

  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);

  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();

  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();

  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  // Other functions

  void doAlgorithm();

private:

  enum EAlgoType { ALGO_BOB, ALGO_JOHN, ALGO_HARRY };

  static const std::string toString( AlgorithmA::EAlgoType algoType );

  Teuchos::RCP<Teuchos::ParameterList> paramList_;
  EAlgoType algoType_;
  double algoTol_;
  
};


#endif // TEUCHOS_ALGORITHM_A_HPP

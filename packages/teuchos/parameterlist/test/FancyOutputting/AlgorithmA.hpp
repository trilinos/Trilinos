// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

  // Overridden from ParameterListAcceptor

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

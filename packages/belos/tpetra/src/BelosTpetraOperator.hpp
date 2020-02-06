//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

/*! \file BelosTpetraOperator.h
    \brief This file provides a Tpetra::Operator interface so Belos can be integrated into
     other codes as an abstract operator.
*/

#ifndef BELOS_TPETRA_OPERATOR_H
#define BELOS_TPETRA_OPERATOR_H

#include "Tpetra_MultiVector_decl.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_Map.hpp"

#include "BelosLinearProblem.hpp"
#include "BelosStatusTest.hpp"
#include "BelosOutputManager.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosSolverFactory.hpp"

#include "Teuchos_ParameterList.hpp"

/*! \class Belos::TpetraOperator
    \brief This class provides an interface to the Tpetra::Operator class, so Belos can be 
    integrated into other codes as an abstract operator.
*/

namespace Belos {

///////////////////////////////////////////////////////////////
//-------- class BelosTpetraOperator --------------------
//
// This class will allow a Belos solver to be called as an Tpetra_Operator.
// Thus, it can use itself as a preconditioner if need be.  It can
// also be used as the inner iteration of Anasazi :)
//
///////////////////////////////////////////////////////////////
template<class Scalar = ::Tpetra::Details::DefaultTypes::scalar_type,
         class LocalOrdinal = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
         class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
         class Node = ::Tpetra::Details::DefaultTypes::node_type>
class TpetraOperator : public virtual Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
public:
  //! @name Typedefs
  //@{
  //! The Tpetra Multivector 
  typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> Tpetra_MultiVector;
  //! The Tpetra Operator
  typedef Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> Tpetra_Operator;
  //! The Tpetra Map
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Tpetra_Map;

  //! @name Constructor / Destructor
  //@{ 
  
  //! Constructor
  TpetraOperator( const Teuchos::RCP<LinearProblem<Scalar,Tpetra_MultiVector,Tpetra_Operator> >& lp, 
		  const Teuchos::RCP<Teuchos::ParameterList>& plist,
                  bool initSolnVec = false )  
  : lp_(lp), 
    plist_(plist),
    initSolnVec_(initSolnVec)
  {
  // Get the solver's name from the parameter list, use block Gmres by default.
  std::string solver = plist_->get("Solver", "BlockGmres");

  // Create a label for this Tpetra_Operator.
  std::string solver_name = "Belos " + solver + " Solver";
 
  // Copy std::string to character array.  
  // Not using conversion routine copy() because it's not supported by RW on Janus. (HKT 11/13/2003) 
  Solver.resize(solver_name.length()+1);
  for (int i=0; i<(int)solver_name.length(); i++) {
    Solver[i] = solver_name[i];
  }   
  Solver[solver_name.length()] = 0;

  //  
  // Create solver and solve problem.  This is inefficient, an instance of the solver should
  // exist already and just be reset with a new RHS.
  //  
  if (solver == "BlockGmres") {
    solver = "Block Gmres";
  }   
  else if (solver == "PseudoBlockGmres") {
    solver = "Pseudo Block Gmres";
  }
  else if (solver == "BlockCG") {
    solver = "Block CG";
  }
  else if (solver == "PseudoBlockCG") {
    solver = "Pseudo Block CG";
  }

  plist_->remove( "Solver" );
  
  Belos::SolverFactory<Scalar,Tpetra_MultiVector,Tpetra_Operator> factory;
  solver_ = factory.create( solver, plist );
  solver_->setProblem( lp_ );

}

  
  //! Destructor
  virtual ~TpetraOperator() {}
  //@}
  
  //! @name Operator application methods
  //@{ 
  
  //! Apply the operator.
  void apply( const Tpetra_MultiVector &X, Tpetra_MultiVector &Y, 
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
           Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const
  {
    Teuchos::RCP<const Tpetra_MultiVector> vec_X = Teuchos::rcp( &X, false );
    Teuchos::RCP<Tpetra_MultiVector> vec_Y = Teuchos::rcp( &Y, false );
    if (initSolnVec_)
    {
      vec_Y->putScalar( 0.0 );
      lp_->setInitResVec( vec_X );
    }
    lp_->setProblem( vec_Y, vec_X );
    Belos::ReturnType ret = solver_->solve();
  };
  //@}
  
  //! @name Attribute access functions
  //@{ 

  //! Return the label of the operator.
  const char* Label() const { return(&Solver[0]); };
  
  //! Return whether the operator supports applying its transpose or conjugate transpose.
  bool hasTransposeApply() const { return(false); };
  
  //! Return the domain map for this operator.
  Teuchos::RCP<const Tpetra_Map> getDomainMap() const
   { return (lp_->getOperator()->getDomainMap()); };
  
  //! Return the range map for this operator.
  Teuchos::RCP<const Tpetra_Map> getRangeMap() const	
   { return (lp_->getOperator()->getRangeMap()); };
  //@}	   

private:

  Teuchos::RCP<SolverManager<Scalar,Tpetra_MultiVector,Tpetra_Operator> > solver_;
  Teuchos::RCP<LinearProblem<Scalar,Tpetra_MultiVector,Tpetra_Operator> > lp_;
  Teuchos::RCP<Teuchos::ParameterList> plist_;

  std::vector<char> Solver;
  bool initSolnVec_;
};

} //end namespace Belos

// end of file BELOS_TPETRA_OPERATOR_H
#endif 


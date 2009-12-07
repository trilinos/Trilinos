// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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

#ifndef TPETRA_FLIPOP_HPP
#define TPETRA_FLIPOP_HPP

#include "Tpetra_InverseOperator.hpp"

namespace Tpetra {

  //! \brief Flip an InverseOperator to make it behave like an Operator.
  /*!
     FlipOp wraps an InverseOperator in the Operator interface, meaning that
     calling FlipOp::apply internally calls the wrapped object's applyInverse
     method.

     This is intended to be used for wrapping Tifpack::Preconditioner objects
     and passing them to Belos::LinearProblem. Belos views preconditioners as
     Operators, and applies them using an Apply implemented in an OperatorTraits
     template (see BelosTpetraAdaptor.hpp in the belos package).

     Tifpack::Preconditioner objects are InverseOperators, and the action of the
     preconditioner is applied using 'applyInverse'. The mis-match between
     Belos' and Tifpack's view of preconditioners motivates this FlipOp wrapper.
   */
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class FlipOp : public Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
public:
  FlipOp(const Teuchos::RCP<const Tifpack::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& InvOp)
   : m_invop(InvOp) {}
  ~FlipOp(){}

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const
  { return m_invop->getDomainMap(); }

  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const
  { return m_invop->getRangeMap(); }

  void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
               Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
               Teuchos::ETransp mode = Teuchos::NO_TRANS) const
  { m_invop->applyInverse(X, Y, mode); }

  bool hasTransposeApply() const
  { return m_invop->hasTransposeApply(); }

private:
  Teuchos::RCP<const Tifpack::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > m_invop;
};//class FlipOp

} // Tpetra namespace

#endif // TPETRA_FLIPOP_HPP

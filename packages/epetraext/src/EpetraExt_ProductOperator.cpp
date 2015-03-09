//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#include "EpetraExt_ProductOperator.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

namespace EpetraExt {

// Constructors / initializers / accessors

ProductOperator::ProductOperator()
{}

ProductOperator::ProductOperator(
        const int                                            numOp
        ,const Teuchos::RefCountPtr<const Epetra_Operator>   op[]
        ,const Teuchos::ETransp                              opTrans[]
        ,const EApplyMode                                    opInverse[]
        )
{
        initialize(numOp,op,opTrans,opInverse);
}

void ProductOperator::initialize(
        const int                                  numOp
        ,const Teuchos::RCP<const Epetra_Operator> op[]
        ,const Teuchos::ETransp                    opTrans[]
        ,const EApplyMode                          opInverse[]
        )
{
#ifdef _DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION(
                numOp < 1, std::invalid_argument
                ,"ProductOperator::initialize(...): Error!"
                );
        // ToDo: Validate maps for operators!
#endif // _DEBUG
        Op_.resize(numOp);
        Op_trans_.resize(numOp);
        Op_inverse_.resize(numOp);
        std::copy( op, op + numOp, Op_.begin() );
        std::copy( opTrans, opTrans + numOp, Op_trans_.begin() );
        std::copy( opInverse, opInverse + numOp, Op_inverse_.begin() );
        UseTranspose_ = false;
        // Wipe cache vectors so that they will be reset just to be safe
        range_vecs_.resize(0);
        domain_vecs_.resize(0);
}

void ProductOperator::uninitialize(
        int                                  *numOp
        ,Teuchos::RCP<const Epetra_Operator> op[]
        ,Teuchos::ETransp                    opTrans[]
        ,EApplyMode                          opInverse[]
        )
{
#ifdef _DEBUG
        TEUCHOS_TEST_FOR_EXCEPTION(
                (op != NULL || opTrans != NULL || opInverse!=NULL) && numOp==NULL
                ,std::invalid_argument
                ,"ProductOperator::uninitialize(...): Error!"
                );
#endif // _DEBUG
        if(numOp) {
                *numOp = Op_.size();
                if(op) std::copy( Op_.begin(), Op_.end(), op );
                if(opTrans) std::copy( Op_trans_.begin(), Op_trans_.end(), opTrans );
                if(opInverse) std::copy( Op_inverse_.begin(), Op_inverse_.end(), opInverse );
        }
        UseTranspose_ = false;
        Op_.resize(0);
        Op_trans_.resize(0);
        Op_inverse_.resize(0);
        range_vecs_.resize(0);
        domain_vecs_.resize(0);
}

void ProductOperator::applyConstituent(
        const int                  k
        ,Teuchos::ETransp          opTrans
        ,EApplyMode                opInverse
        ,const Epetra_MultiVector  &X_k
        ,Epetra_MultiVector        *Y_k
        ) const
{
  Epetra_Operator &Op_k = const_cast<Epetra_Operator&>(*Op_[k]); // Okay since we put back UseTranspose!
  bool oldUseTranspose = Op_k.UseTranspose();
  Op_k.SetUseTranspose((opTrans==Teuchos::NO_TRANS) != (Op_trans_[k]==Teuchos::NO_TRANS));
  const bool applyInverse_k = (opInverse==APPLY_MODE_APPLY) != (Op_inverse_[k]==APPLY_MODE_APPLY);
  const int err = ! applyInverse_k ? Op_[k]->Apply(X_k,*Y_k) : Op_[k]->ApplyInverse(X_k,*Y_k);
  Op_k.SetUseTranspose(oldUseTranspose);
  TEUCHOS_TEST_FOR_EXCEPTION(
    err!=0, std::runtime_error, "ProductOperator::applyConstituent(...): Error,"
    " Op["<<k<<"]." << (!applyInverse_k?"Apply":"ApplyInverse") << "(...) "
    "returned err = " << err << " with Op["<<k<<"].UseTranspose() = "<<
    Op_[k]->UseTranspose() << "!");
}

// Overridden from Epetra_Operator

int ProductOperator::SetUseTranspose(bool useTranspose)
{
        assertInitialized();
        UseTranspose_ = useTranspose;
        return 0;
}

int ProductOperator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  assertInitialized();
  const int numOp = this->num_Op();
  // Setup the temporary vectors
  initializeTempVecs(false);
  // Apply the constituent operators one at a time!
  if( !UseTranspose_ ) {
    //
    // Forward Mat-vec: Y = M * X (See main documenation)
    //
    for( int k = numOp-1; k >= 0; --k ) {
      const Epetra_MultiVector  &X_k = ( k==numOp-1 ? X : *range_vecs_[k]   );
      Epetra_MultiVector        &Y_k = ( k==0        ? Y : *range_vecs_[k-1] );
      applyConstituent(k,Teuchos::NO_TRANS,APPLY_MODE_APPLY,X_k,&Y_k);
    }
  }
  else if( UseTranspose_ ) {
    //
    // Adjoint Mat-vec: Y = M' * X (See main documentation)
    //
    for( int k = 0; k <= numOp-1; ++k ) {
      const Epetra_MultiVector  &X_k = ( k==0         ? X : *domain_vecs_[k-1] );
      Epetra_MultiVector        &Y_k = ( k==numOp-1  ? Y : *domain_vecs_[k]   );
      applyConstituent(k,Teuchos::TRANS,APPLY_MODE_APPLY,X_k,&Y_k);
    }
  }
  return 0;
}

int ProductOperator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  assertInitialized();
  const int numOp = this->num_Op();
  // Setup the temporary vectors
  initializeTempVecs(true);
  // Apply the constituent operators one at a time!
  if( !UseTranspose_ ) {
    //
    // Forward Inverse Mat-vec: Y = inv(M) * X (See main documenation)
    //
    for( int k = 0; k <= numOp-1; ++k ) {
      const Epetra_MultiVector  &X_k = ( k==0         ? X : *domain_vecs_[k-1] );
      Epetra_MultiVector        &Y_k = ( k==numOp-1  ? Y : *domain_vecs_[k]   );
      applyConstituent(k,Teuchos::NO_TRANS,APPLY_MODE_APPLY_INVERSE,X_k,&Y_k);
    }
  }
  else if( UseTranspose_ ) {
    //
    // Adjoint Invese Mat-vec: Y = inv(M') * X (See main documentation)
    //
    for( int k = numOp-1; k >= 0; --k ) {
      const Epetra_MultiVector  &X_k = ( k==numOp-1 ? X : *range_vecs_[k]   );
      Epetra_MultiVector        &Y_k = ( k==0        ? Y : *range_vecs_[k-1] );
      applyConstituent(k,Teuchos::TRANS,APPLY_MODE_APPLY_INVERSE,X_k,&Y_k);
    }
  }
  return 0;
}

double ProductOperator::NormInf() const
{
        assertInitialized();
        return -1.0;
}

const char* ProductOperator::Label() const
{
        assertInitialized();
        return NULL;
}

bool ProductOperator::UseTranspose() const
{
        assertInitialized();
        return UseTranspose_;
}

bool ProductOperator::HasNormInf() const
{
        assertInitialized();
        return false;
}

const Epetra_Comm&
ProductOperator::Comm() const
{
        assertInitialized();
        return Op_.front()->OperatorRangeMap().Comm();
}

const Epetra_Map&
ProductOperator::OperatorDomainMap() const
{
        assertInitialized();
        return ( Op_trans_.back()==Teuchos::NO_TRANS
                                         ? Op_.back()->OperatorDomainMap()
                                         : Op_.back()->OperatorRangeMap()
                                         );
}

const Epetra_Map&
ProductOperator::OperatorRangeMap() const
{
        assertInitialized();
        return ( Op_trans_.front()==Teuchos::NO_TRANS
                                         ? Op_.front()->OperatorRangeMap()
                                         : Op_.front()->OperatorDomainMap()
                                         );
}

// private

void ProductOperator::initializeTempVecs(bool applyInverse) const
{
  const int numOp = this->num_Op ();
  if (numOp > 0) {
    // FIXME (mfh 24 Mar 2014): I added the parentheses around the ||
    // below to silence a compiler warning.  I'm concerned that the
    // original author of that code didn't understand that && takes
    // precedence over ||, but I didn't want to change the meaning of
    // the original code.
    if (((! UseTranspose_ && ! applyInverse) || (UseTranspose_ && applyInverse))
        && range_vecs_.size () == 0) {
      //
      // Forward Mat-vec
      //
      // We need to create storage to hold:
      //
      //  T[k-1] = M[k]*T[k]
      //
      //    for k = numOp-1...1
      //
      //      where: T[numOp-1] = X (input vector)
      //
      range_vecs_.resize (numOp - 1);
      for (int k = numOp-1; k >= 1; --k) {
        range_vecs_[k-1] = Teuchos::rcp (new Epetra_Vector ((Op_trans_[k]==Teuchos::NO_TRANS) == (Op_inverse_[k]==APPLY_MODE_APPLY)
                                                            ? Op_[k]->OperatorRangeMap ()
                                                            : Op_[k]->OperatorDomainMap ()));
      }
    }
    // FIXME (mfh 24 Mar 2014): I added the parentheses around the ||
    // below to silence a compiler warning.  I'm concerned that the
    // original author of that code didn't understand that && takes
    // precedence over ||, but I didn't want to change the meaning of
    // the original code.
    else if (((UseTranspose_ && ! applyInverse) || (! UseTranspose_ && applyInverse))
             && domain_vecs_.size () == 0) {
      //
      // Adjoint Mat-vec
      //
      // We need to create storage to hold:
      //
      //   T[k] = M[k]'*T[k-1]
      //
      //     for k = 0...numOp-2
      //
      //       where: T[-1]       = X (input vector)
      //
      domain_vecs_.resize (numOp - 1);
      for (int k = 0; k <= numOp - 2; ++k) {
        domain_vecs_[k] = Teuchos::rcp (new Epetra_Vector ((Op_trans_[k]==Teuchos::NO_TRANS) == (Op_inverse_[k]==APPLY_MODE_APPLY)
                                                           ? Op_[k]->OperatorDomainMap ()
                                                           : Op_[k]->OperatorRangeMap ()));
      }
    }
  }
}

} // namespace EpetraExt

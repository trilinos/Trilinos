// @HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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

#include "EpetraExt_ProductOperator.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

namespace EpetraExt {

// Constructors / initializers / accessors

ProductOperator::ProductOperator()
{}

ProductOperator::ProductOperator(
	const int                                            num_Op
	,const Teuchos::RefCountPtr<const Epetra_Operator>   Op[]
	,const Teuchos::ETransp                              Op_trans[]
	,const EApplyMode                                    Op_inverse[]
	)
{
	initialize(num_Op,Op,Op_trans,Op_inverse);
}

void ProductOperator::initialize(
	const int                                      num_Op
	,const Teuchos::RefCountPtr<const Epetra_Operator>   Op[]
	,const Teuchos::ETransp                        Op_trans[]
	,const EApplyMode                              Op_inverse[]
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		num_Op < 1, std::invalid_argument
		,"ProductOperator::initialize(...): Error!"
		);
	// ToDo: Validate maps for operators!
#endif // _DEBUG
	Op_.resize(num_Op);
	Op_trans_.resize(num_Op);
	Op_inverse_.resize(num_Op);
	std::copy( Op, Op + num_Op, Op_.begin() );
	std::copy( Op_trans, Op_trans + num_Op, Op_trans_.begin() );
	std::copy( Op_inverse, Op_inverse + num_Op, Op_inverse_.begin() );
	UseTranspose_ = false;
	// Wipe cache vectors so that they will be reset just to be safe
	range_vecs_.resize(0);
	domain_vecs_.resize(0);
}

void ProductOperator::uninitialize(
	int                                      *num_Op
	,Teuchos::RefCountPtr<const Epetra_Operator>   Op[]
	,Teuchos::ETransp                        Op_trans[]
	,EApplyMode                              Op_inverse[]
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		(Op!=NULL || Op_trans!=NULL || Op_inverse!=NULL) && num_Op==NULL
		,std::invalid_argument
		,"ProductOperator::uninitialize(...): Error!"
		);
#endif // _DEBUG
	if(num_Op) {
		*num_Op = Op_.size();
		if(Op) std::copy( Op_.begin(), Op_.end(), Op );
		if(Op_trans) std::copy( Op_trans_.begin(), Op_trans_.end(), Op_trans );
		if(Op_inverse) std::copy( Op_inverse_.begin(), Op_inverse_.end(), Op_inverse );
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
	,Teuchos::ETransp          Op_trans
	,EApplyMode                Op_inverse
	,const Epetra_MultiVector  &X_k
	,Epetra_MultiVector        *Y_k
	) const
{
  Epetra_Operator &Op_k = const_cast<Epetra_Operator&>(*Op_[k]); // Okay since we put back UseTranspose!
  bool oldUseTranspose = Op_k.UseTranspose();
	Op_k.SetUseTranspose((Op_trans==Teuchos::NO_TRANS)!=(Op_trans_[k]==Teuchos::NO_TRANS));
	const bool applyInverse_k = (Op_inverse==APPLY_MODE_APPLY)!=(Op_inverse_[k]==APPLY_MODE_APPLY);
	const int err = !applyInverse_k ? Op_[k]->Apply(X_k,*Y_k) :  Op_[k]->ApplyInverse(X_k,*Y_k);
  Op_k.SetUseTranspose(oldUseTranspose);
	TEST_FOR_EXCEPTION(
		err!=0, std::runtime_error
		,"ProductOperator::applyConstituent(...): Error, Op["<<k<<"]."
		<<(!applyInverse_k?"Apply":"ApplyInverse")<<"(...) returned "
		"err = "<<err<<" with Op["<<k<<"].UseTranspose() = "<<Op_[k]->UseTranspose()<<"!"
		);
}

// Overridden from Epetra_Operator

int ProductOperator::SetUseTranspose(bool UseTranspose)
{
	assertInitialized();
	UseTranspose_ = UseTranspose;
	return 0;
}

int ProductOperator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
	assertInitialized();
	const int num_Op = this->num_Op();
	// Setup the temporary vectors
	initializeTempVecs(false);
	// Apply the constituent operators one at a time!
	if( !UseTranspose_ ) {
		//
		// Forward Mat-vec: Y = M * X (See main documenation)
		//
		for( int k = num_Op-1; k >= 0; --k ) {
			const Epetra_MultiVector  &X_k = ( k==num_Op-1 ? X : *range_vecs_[k]   );
			Epetra_MultiVector        &Y_k = ( k==0        ? Y : *range_vecs_[k-1] );
			applyConstituent(k,Teuchos::NO_TRANS,APPLY_MODE_APPLY,X_k,&Y_k);
		}
	}
	else if( UseTranspose_ ) {
		//
		// Adjoint Mat-vec: Y = M' * X (See main documentation)
		//
		for( int k = 0; k <= num_Op-1; ++k ) {
			const Epetra_MultiVector  &X_k = ( k==0         ? X : *domain_vecs_[k-1] );
			Epetra_MultiVector        &Y_k = ( k==num_Op-1  ? Y : *domain_vecs_[k]   );
			applyConstituent(k,Teuchos::TRANS,APPLY_MODE_APPLY,X_k,&Y_k);
		}
	}
	return 0;
}

int ProductOperator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
	assertInitialized();
	const int num_Op = this->num_Op();
	// Setup the temporary vectors
	initializeTempVecs(true);
	// Apply the constituent operators one at a time!
	if( !UseTranspose_ ) {
		//
		// Forward Inverse Mat-vec: Y = inv(M) * X (See main documenation)
		//
		for( int k = 0; k <= num_Op-1; ++k ) {
			const Epetra_MultiVector  &X_k = ( k==0         ? X : *domain_vecs_[k-1] );
			Epetra_MultiVector        &Y_k = ( k==num_Op-1  ? Y : *domain_vecs_[k]   );
			applyConstituent(k,Teuchos::NO_TRANS,APPLY_MODE_APPLY_INVERSE,X_k,&Y_k);
		}
	}
	else if( UseTranspose_ ) {
		//
		// Adjoint Invese Mat-vec: Y = inv(M') * X (See main documentation)
		//
		for( int k = num_Op-1; k >= 0; --k ) {
			const Epetra_MultiVector  &X_k = ( k==num_Op-1 ? X : *range_vecs_[k]   );
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
	const int num_Op = this->num_Op();
	if( num_Op > 0 ) {
		if( ( !UseTranspose_ && !applyInverse ) || ( UseTranspose_ && applyInverse )
				&& range_vecs_.size()==0
			)
		{
			//
			// Forward Mat-vec
			//
			// We need to create storage to hold:
			//
			//  T[k-1] = M[k]*T[k]
			//
			//    for k = num_Op-1...1
			//
			//      where: T[num_Op-1] = X (input vector)
			//
			range_vecs_.resize(num_Op-1);
			for( int k = num_Op-1; k >= 1; --k ) {
				range_vecs_[k-1] = Teuchos::rcp(
					new Epetra_Vector(
						(Op_trans_[k]==Teuchos::NO_TRANS) == (Op_inverse_[k]==APPLY_MODE_APPLY)
						? Op_[k]->OperatorRangeMap()
						: Op_[k]->OperatorDomainMap()
						)
					);
			}
		}
		else if( ( UseTranspose_ && !applyInverse ) || ( !UseTranspose_ && applyInverse )
						 && domain_vecs_.size()==0
			)
		{
			//
			// Adjoint Mat-vec
			//
			// We need to create storage to hold:
			//
			//   T[k] = M[k]'*T[k-1]
			//
			//     for k = 0...num_Op-2
			//
			//       where: T[-1]       = X (input vector)
			//
			domain_vecs_.resize(num_Op-1);
			for( int k = 0; k <= num_Op-2; ++k ) {
				domain_vecs_[k] = Teuchos::rcp(
					new Epetra_Vector(
						(Op_trans_[k]==Teuchos::NO_TRANS) == (Op_inverse_[k]==APPLY_MODE_APPLY)
						? Op_[k]->OperatorDomainMap()
						: Op_[k]->OperatorRangeMap()
						)
					);
			}
		}
	}
}

} // namespace EpetraExt

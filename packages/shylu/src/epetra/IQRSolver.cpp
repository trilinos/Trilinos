//@HEADER
// ************************************************************************
//
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
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

/** \file IQRSolver.h

    \brief Encapsulates the IQR inexact solver functionality

    \author Radu Popescu

*/

#include <IQRSolver.h>

#include <cmath>

#include <Ifpack_config.h>

#ifdef HAVE_IFPACK_DYNAMIC_FACTORY
#include "Ifpack_DynamicFactory.h"
#else
#include "Ifpack.h"
#endif

namespace IQR {

IQRSolver::IQRSolver(const Teuchos::ParameterList& pList)
	: gmresStateManager_(),
	  prec_(0),
	  isComputed_(false),
	  krylovDim_(0.0),
	  numIter_(0),
	  doScaling_(false),
	  useFullIQR_(true),
	  precType_("Amesos"),
	  precAmesosType_("Amesos_Klu")
{
	// We shouldn't modify the external parameter list
	Teuchos::ParameterList params = pList;

	// Default mode, G or IQR?
    useFullIQR_ = params.get<bool>("Use full IQR", true);

	// Default parameters for IQR
	krylovDim_ = params.get<double>("IQR Krylov Dim", 0.5);
	numIter_ = params.get<int>("IQR Number Iterations", 0);
    doScaling_ = params.get<bool>("IQR Scaling", false);
    precType_ = params.get<string>("IQR Initial Prec Type", "Amesos");
    precAmesosType_ = params.get<string>("IQR Initial Prec Amesos Type", "Amesos_Klu");
}

IQRSolver::~IQRSolver()
{
}

int IQRSolver::Compute(const ShyLU_Probing_Operator& S,
						const Epetra_MultiVector& B,
						Epetra_MultiVector & X)
{
	// The preconditioner is used regardless of the value ofuseFullIQR_
	Epetra_CrsMatrix* G = S.G_;
#ifdef HAVE_IFPACK_DYNAMIC_FACTORY
    Ifpack_DynamicFactory Factory;
#else
    Ifpack Factory;
#endif
	prec_ = Teuchos::rcp(Factory.Create(precType_, G, 1));
	Teuchos::ParameterList pList;
	pList.set("amesos: solver type", precAmesosType_, "");
	pList.set("Reindex", true);
	prec_->SetParameters(pList);

	IFPACK_CHK_ERR(prec_->Initialize());
	IFPACK_CHK_ERR(prec_->Compute());

	// Only for IQR
	if (useFullIQR_) {
		// Computation phase - only during the first outer GMRES iteration
		int sSize = S.OperatorDomainMap().NumGlobalElements();
		int kSize = std::floor(krylovDim_ * sSize);

		gmresStateManager_ = Teuchos::rcp(
				new GMRESStateManager(S.OperatorDomainMap(),
									  kSize, false, doScaling_));
		gmresStateManager_->isFirst = true;

		double tol = 1e-10;
		IQR::IdPreconditioner L;
		IQR::GMRES<Epetra_Operator,
				   Epetra_MultiVector,
				   IQR::IdPreconditioner,
				   Ifpack_Preconditioner,
				   GMRESStateManager,
				   std::vector<double>, double>
			(S, X, B, &L, &(*prec_), *(gmresStateManager_), kSize, tol);
		gmresStateManager_->P2 = &(*prec_);

//			if (! Xs.Comm().MyPID()) {
//				std::cout << "KSIZE: " << kSize
//						  << ", SSIZE: " << sSize
//						  << ", TOL: " << tol << std::endl;
//			}
	}
	isComputed_ = true;

	return 0;
}

int IQRSolver::Solve(const ShyLU_Probing_Operator& S,
		   	   	   	  const Epetra_MultiVector& B,
		   	   	   	  Epetra_MultiVector & X)
{
	if (! isComputed_) {
		IFPACK_CHK_ERR(Compute(S, B, X));
	}
	// Solve phase
	if (! useFullIQR_) {
		IFPACK_CHK_ERR(prec_->ApplyInverse(B, X));
	} else {
		if (numIter_ > 0) {
			GMRESStateManager newGmresManager(S.OperatorDomainMap(),
											  numIter_+1, false);
		double tol=1e-10;
		IQR::IdPreconditioner L;
		IQR::GMRES<Epetra_Operator, Epetra_MultiVector,
				   IQR::IdPreconditioner, GMRESStateManager, GMRESStateManager,
				   std::vector<double>, double>
			(S, X, B, &L, &*(gmresStateManager_), newGmresManager, numIter_, tol);
		} else {
			IFPACK_CHK_ERR(gmresStateManager_->ApplyInverse(B, X));
		}
	}

	return 0;
}

} /* namespace IQR */

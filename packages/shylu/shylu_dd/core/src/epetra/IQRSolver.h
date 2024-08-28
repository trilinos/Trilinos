// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file IQRSolver.h

    \brief Encapsulates the IQR inexact solver functionality

    \author Radu Popescu

*/

#ifndef IQRSOLVER_H_
#define IQRSOLVER_H_

#if defined(ShyLU_DDCore_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ShyLU_DDCore package is deprecated"
#endif
#endif

#include <vector>

#include <Epetra_BlockMap.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>

#include <shylu_internal_gmres_tools.h>
#include <shylu_internal_gmres.h>

#include <Ifpack_Preconditioner.h>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <shylu_probing_operator.h>

namespace IQR {

class IQRSolver {
public:
	// Public type definitions
	typedef IQR::GMRESManager<Epetra_BlockMap,
							  Epetra_MultiVector,
							  std::vector<std::vector<double> >,
							  std::vector<double> > GMRESStateManager;

	// Public constructor and destructor
	IQRSolver(const Teuchos::ParameterList& pList);
	virtual ~IQRSolver();

	// Public methods
	int Solve(const ShyLU_Probing_Operator& S,
			   const Epetra_MultiVector& B,
			   Epetra_MultiVector & X);

private:
	// Private methods
	int Compute(const ShyLU_Probing_Operator& S,
			     const Epetra_MultiVector& B,
			     Epetra_MultiVector & X);

	// Private data
	Teuchos::RCP<GMRESStateManager> gmresStateManager_;
	Teuchos::RCP<Ifpack_Preconditioner> prec_;
	bool isComputed_;
    double krylovDim_;
    int numIter_;
    bool doScaling_;
    bool useFullIQR_;
    std::string precType_;
    std::string precAmesosType_;
};

} /* namespace IQR */

#endif /* IQRSOLVER_H_ */

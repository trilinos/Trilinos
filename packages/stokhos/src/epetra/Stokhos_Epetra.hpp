// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_EPETRA_HPP
#define STOKHOS_EPETRA_HPP

#include "Stokhos.hpp"

#ifdef HAVE_STOKHOS_EPETRAEXT

#include "Stokhos_ParallelData.hpp"
#include "Stokhos_VectorOrthogPoly.hpp"
#include "Stokhos_VectorOrthogPolyTraitsEpetra.hpp"
#include "Stokhos_Sparse3TensorUtilities.hpp"

// SG Operators
#include "Stokhos_MLPrecOp.hpp"
#include "Stokhos_MatrixFreeOperator.hpp"
#include "Stokhos_KLMatrixFreeOperator.hpp"
#include "Stokhos_KLReducedMatrixFreeOperator.hpp"
#include "Stokhos_FullyAssembledOperator.hpp"
#include "Stokhos_SGOperatorFactory.hpp"

// SG Preconditioners
#include "Stokhos_MeanBasedPreconditioner.hpp"
#include "Stokhos_GaussSeidelPreconditioner.hpp"
#include "Stokhos_ApproxGaussSeidelPreconditioner.hpp"
#include "Stokhos_ApproxJacobiPreconditioner.hpp"
#include "Stokhos_KroneckerProductPreconditioner.hpp"
#include "Stokhos_FullyAssembledPreconditioner.hpp"
#include "Stokhos_SGPreconditionerFactory.hpp"

// SG Model Evaluators
#include "Stokhos_SGModelEvaluator.hpp"
#include "Stokhos_SGModelEvaluator_Interlaced.hpp"
#include "Stokhos_SGModelEvaluator_Adaptive.hpp"
#include "Stokhos_SGQuadModelEvaluator.hpp"
#include "Stokhos_SGInverseModelEvaluator.hpp"
#include "Stokhos_ResponseStatisticModelEvaluator.hpp"

// MP Operators
#include "Stokhos_BlockDiagonalOperator.hpp"

// MP Preconditioners
#include "Stokhos_MPBlockDiagonalPreconditioner.hpp"
#include "Stokhos_MPPreconditionerFactory.hpp"

// MP Model Evaluators
#include "Stokhos_MPModelEvaluator.hpp"
#include "Stokhos_MPModelEvaluatorAdapter.hpp"
#include "Stokhos_SGQuadMPModelEvaluator.hpp"
#include "Stokhos_MPInverseModelEvaluator.hpp"

#endif // HAVE_STOKHOS_EPETRAEXT

#endif // STOKHOS_EPETRA_HPP 

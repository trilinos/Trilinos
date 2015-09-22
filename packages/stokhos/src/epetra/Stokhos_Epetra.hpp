// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
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

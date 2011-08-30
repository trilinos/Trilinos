// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_EPETRA_HPP
#define STOKHOS_EPETRA_HPP

#include "Stokhos.hpp"

#include "Stokhos_ParallelData.hpp"
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

#endif // STOKHOS_EPETRA_HPP 

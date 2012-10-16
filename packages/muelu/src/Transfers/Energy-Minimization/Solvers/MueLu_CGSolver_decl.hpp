// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_CGSOLVER_DECL_HPP
#define MUELU_CGSOLVER_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include <Xpetra_CrsMatrixFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SolverBase.hpp"
#include "MueLu_Constraint_fwd.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

  /*!
    @class CGSolver class.
    @brief Implements conjugate gradient algorithm for energy-minimization

    This is CG applied to the problem
    \f[ \min_P  \frac12 \sum_i (p_i)^T A p_i \f]
    subject to
    \f[ P(i,j) = 0 \quad \mbox{if} \; (i,j) \mbox{ is not in SparsityPattern} \f]
    and
    \f[ P cnull =  fnull. \f]

    \note
        For now we assume, that matrices have real values. The case of complex values is not explored.
    */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class CGSolver : public SolverBase<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> {
#undef MUELU_CGSOLVER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors / destructors
    //@{

    /*! Constructor
      \param Its -- Number of performed iterations
      */
    CGSolver(size_t Its);

    //@}

    //! @name Iterate methods.
    //@{

    //! Iterate
    void Iterate(const Matrix& A, const Constraint& C, const Matrix& P0, const MultiVector& B, RCP<Matrix>& P) const;

    //@}
  private:
    SC Frobenius(const Matrix& A, const Matrix& B) const;

  private:
    size_t nIts_;           //!< Number of performed iterations
  };

} // namespace MueLu

#define MUELU_CGSOLVER_SHORT
#endif // MUELU_CGSOLVER_DECL_HPP

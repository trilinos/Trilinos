// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_LU2x2InverseOp.hpp"

#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

using namespace Thyra;
using Teuchos::ArrayRCP;
using Teuchos::dyn_cast;
using Teuchos::RCP;

namespace Teko {

// Thyra::LinearOpBase requirements
////////////////////////////////////////////////////////////////////////
LU2x2InverseOp::LU2x2InverseOp(const BlockedLinearOp& A, const LinearOp& invA00,
                               const LinearOp& invS)
    : A_(A),
      hatInvA00_(invA00),
      tildeInvA00_(invA00),
      invS_(invS),
      A10_(A->getBlock(1, 0)),
      A01_(A->getBlock(0, 1)) {
  using Teuchos::tuple;
  using Thyra::productVectorSpace;

  // create and store range space
  productRange_ = productVectorSpace<double>(tuple(hatInvA00_->range(), invS_->range())());

  // create and store domain space
  productDomain_ = productVectorSpace<double>(tuple(hatInvA00_->domain(), invS_->domain())());
}

/** \brief This constructor explicitly takes the parts of \f$ A \f$ required to
 *        build the inverse operator.
 *
 * This constructor explicitly takes the parts of \f$ A \f$ required to build
 * the inverse operator.
 *
 * \param[in] A The block \f$ 2 \times 2 \f$ \f$A\f$ operator.
 * \param[in] hatInvA00  An approximate inverse of \f$ \hat{A}_{00} \f$
 * \param[in] tildeInvA00  An approximate inverse of \f$ \tilde{A}_{00} \f$
 * \param[in] invS  An approximate inverse of \f$ S = -A_{11} + A_{10} A_{00}^{-1} A_{01} \f$.
 */
LU2x2InverseOp::LU2x2InverseOp(const BlockedLinearOp& A, const LinearOp& hatInvA00,
                               const LinearOp& tildeInvA00, const LinearOp& invS)
    : A_(A),
      hatInvA00_(hatInvA00),
      tildeInvA00_(tildeInvA00),
      invS_(invS),
      A10_(A->getBlock(1, 0)),
      A01_(A->getBlock(0, 1)) {
  using Teuchos::tuple;
  using Thyra::productVectorSpace;

  // create and store range space
  productRange_ = productVectorSpace<double>(tuple(hatInvA00_->range(), invS_->range())());

  // create and store domain space
  productDomain_ = productVectorSpace<double>(tuple(hatInvA00_->domain(), invS_->domain()));
}

void LU2x2InverseOp::implicitApply(const BlockedMultiVector& x, BlockedMultiVector& y,
                                   const double alpha, const double beta) const {
  // get src blocks
  MultiVector f = getBlock(0, x);  // f
  MultiVector g = getBlock(1, x);  // g

  // get extra storage
  MultiVector ps = deepcopy(g);  // this is need b/c of p = -inv(S)*p

  // get destination blocks
  MultiVector u = getBlock(0, y);  // u (u^)
  MultiVector p = getBlock(1, y);  // p (p^)

  // for efficiency make copies of u and p
  MultiVector uc, pc;  // copies of u and p
  if (beta != 0) {
    // perform a deep copy
    uc = deepcopy(u);
    pc = deepcopy(p);
  } else {
    // perform a shallow copy

    // uc and pc point
    // to the data of u and p
    uc = u;
    pc = p;
  }

  // set temporary operator for performing inv(A_00)*A_01
  LinearOp invA00_A01 = Thyra::multiply(tildeInvA00_, A01_);

  // compute actual product
  applyOp(hatInvA00_, f, uc);              // u   = inv(A_00) * f
  applyOp(A10_, uc, ps, -1.0, 1.0);        // ps += -A_10*u
  applyOp(invS_, ps, pc, -1.0);            // p   = -inv(S)*ps
  applyOp(invA00_A01, pc, uc, -1.0, 1.0);  // u  += -inv(A_00)*A_01*p

  // scale result by alpha
  if (beta != 0) {
    update(alpha, uc, beta, u);  // u = alpha * uc + beta * u
    update(alpha, pc, beta, p);  // p = alpha * pc + beta * p
  } else if (alpha != 1.0) {
    scale(alpha, u);  // u = alpha * u
    scale(alpha, p);  // p = alpha * p
  }
}

void LU2x2InverseOp::describe(Teuchos::FancyOStream& out_arg,
                              const Teuchos::EVerbosityLevel verbLevel) const {
  using Teuchos::OSTab;

  RCP<Teuchos::FancyOStream> out = rcp(&out_arg, false);
  OSTab tab(out);
  switch (verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW: *out << this->description() << std::endl; break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME: {
      *out << Teuchos::Describable::description() << "{"
           << "rangeDim=" << this->range()->dim() << ",domainDim=" << this->domain()->dim()
           << "}\n";
      {
        *out << "[invS]:\n";
        *out << Teuchos::describe(*invS_, verbLevel);

        *out << "[hatInvA00]:\n";
        *out << Teuchos::describe(*hatInvA00_, verbLevel);

        *out << "[tildeInvA00]:\n";
        *out << Teuchos::describe(*tildeInvA00_, verbLevel);

        *out << "[A_10]:\n";
        *out << Teuchos::describe(*A10_, verbLevel);

        *out << "[A_01]:\n";
        *out << Teuchos::describe(*A01_, verbLevel);
      }
      break;
    }
    default: TEUCHOS_TEST_FOR_EXCEPT(true);  // Should never get here!
  }
}

}  // end namespace Teko

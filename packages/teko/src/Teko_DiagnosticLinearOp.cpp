// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_DiagnosticLinearOp.hpp"

#include "Teuchos_TimeMonitor.hpp"

#include "Thyra_MultiVectorStdOps.hpp"

namespace Teko {

/** \brief This constructor explicitly takes the linear operator
 *        that needs to be wrapped and a string for output that describes
 *        the diagnostics.
 */
DiagnosticLinearOp::DiagnosticLinearOp(const Teuchos::RCP<std::ostream>& ostrm,
                                       const ModifiableLinearOp& A,
                                       const std::string& diagnosticString)
    : outputStream_(ostrm),
      wrapOpA_(A),
      wrapOpA_lo_(A),
      diagString_(diagnosticString),
      timer_(diagnosticString) {}

DiagnosticLinearOp::DiagnosticLinearOp(const Teuchos::RCP<std::ostream>& ostrm, const LinearOp& A,
                                       const std::string& diagnosticString)
    : outputStream_(ostrm),
      wrapOpA_(Teuchos::null),
      wrapOpA_lo_(A),
      diagString_(diagnosticString),
      timer_(diagnosticString) {}

/** \brief This constructor explicitly takes the linear operator
 *        that needs to be wrapped and a string for output that describes
 *        the diagnostics.
 */
DiagnosticLinearOp::DiagnosticLinearOp(const Teuchos::RCP<std::ostream>& ostrm,
                                       const LinearOp& fwdOp, const ModifiableLinearOp& A,
                                       const std::string& diagnosticString)
    : outputStream_(ostrm),
      wrapOpA_(A),
      wrapOpA_lo_(A),
      fwdOp_(fwdOp),
      diagString_(diagnosticString),
      timer_(diagnosticString) {}

DiagnosticLinearOp::~DiagnosticLinearOp() {
  double elapsedTime = totalTime();
  int applications   = numApplications();

  (*outputStream_) << "DiagnosticLinearOp \"" << diagString_ << "\": "
                   << "elapsed = " << elapsedTime << ", "
                   << "applications = " << applications << ", ";
  if (applications > 0)
    (*outputStream_) << "timer/app = " << elapsedTime / double(applications) << std::endl;
  else
    (*outputStream_) << "timer/app = "
                     << "none" << std::endl;
}

/** @brief Perform a matrix vector multiply with this operator.
 *
 * The <code>apply</code> function takes one vector as input
 * and applies the inverse \f$ LDU \f$ decomposition. The result
 * is returned in \f$y\f$. If this operator is reprsented as \f$M\f$ then
 * \f$ y = \alpha M x + \beta y \f$ (ignoring conjugation!).
 *
 * @param[in]     x
 * @param[in,out] y
 * @param[in]     alpha (default=1)
 * @param[in]     beta  (default=0)
 */
void DiagnosticLinearOp::implicitApply(const MultiVector& x, MultiVector& y, const double alpha,
                                       const double beta) const {
  Teko_DEBUG_SCOPE("DiagnosticLinearOp::implicityApply", 10);

  // start timer on construction, end on destruction
  Teuchos::TimeMonitor monitor(timer_, false);

  MultiVector z;  // for temporary storage dealing with nozero beta
  if (beta != 0.0) z = deepcopy(y);

  wrapOpA_lo_->apply(Thyra::NOTRANS, *x, y.ptr(), alpha, beta);

  // print residual if there is a fwd Op
  bool printResidual = (fwdOp_ != Teuchos::null);
  if (printResidual) {
    // compute residual
    MultiVector residual = Teko::deepcopy(x);
    // fwdOp_->apply(Thyra::NOTRANS,*y,residual.ptr(),-1.0,1.0);

    fwdOp_->apply(Thyra::NOTRANS, *y, residual.ptr(), -1.0, alpha);
    if (beta != 0.0) fwdOp_->apply(Thyra::NOTRANS, *z, residual.ptr(), beta, 1.0);

    // calculate norms
    std::vector<double> norms(y->domain()->dim());     // size of column count
    std::vector<double> rhsNorms(x->domain()->dim());  // size of column count
    Thyra::norms_2<double>(*residual, Teuchos::arrayViewFromVector(norms));
    Thyra::norms_2<double>(*x, Teuchos::arrayViewFromVector(rhsNorms));

    // print out residual norms
    (*outputStream_) << "DiagnosticLinearOp \"" << diagString_ << "\": residual = [";
    for (std::size_t i = 0; i < norms.size(); ++i)
      (*outputStream_) << " " << std::scientific << std::setprecision(4)
                       << norms[i] / rhsNorms[i];  // << " (" <<rhsNorms[i]<<") ";
    (*outputStream_) << " ]" << std::endl;

    residualNorm_ = norms[0];
  }
}

}  // end namespace Teko

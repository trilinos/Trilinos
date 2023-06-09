/*@HEADER
// ***********************************************************************
//
//         Stratimikos: Thyra-based strategies for linear solvers
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
// Questions? Contact Jennifer A. Loe (jloe@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/
#ifndef THYRA_BELOS_TPETRA_PRECONDITIONERFACTORY_DECL_HPP
#define THYRA_BELOS_TPETRA_PRECONDITIONERFACTORY_DECL_HPP

#include "Thyra_PreconditionerFactoryBase.hpp"

namespace Thyra {

/** \brief Concrete preconditioner factory subclass based on Belos.
 * (Yes, Belos solvers can also be used as preconditioners!)
 */
template <typename MatrixType>
class BelosTpetraPreconditionerFactory :
  public PreconditionerFactoryBase<typename MatrixType::scalar_type> {
public:
  /** \brief . */
  typedef typename MatrixType::scalar_type scalar_type;

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  BelosTpetraPreconditionerFactory();
  //@}

  /** @name Overridden from PreconditionerFactoryBase */
  //@{

  /** \brief . */
  bool isCompatible(const LinearOpSourceBase<scalar_type> &fwdOp) const;

  /** \brief . */
  Teuchos::RCP<PreconditionerBase<scalar_type> > createPrec() const;

  /** \brief . */
  void initializePrec(
    const Teuchos::RCP<const LinearOpSourceBase<scalar_type> > &fwdOp,
    PreconditionerBase<scalar_type> *prec,
    const ESupportSolveUse supportSolveUse
    ) const;

  /** \brief . */
  void uninitializePrec(
    PreconditionerBase<scalar_type> *prec,
    Teuchos::RCP<const LinearOpSourceBase<scalar_type> > *fwdOp,
    ESupportSolveUse *supportSolveUse
    ) const;

  //@}

  /** @name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> &paramList);
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

  /** \name Public functions overridden from Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  Teuchos::RCP<Teuchos::ParameterList> paramList_;

};

} // namespace Thyra

#endif // THYRA_BELOS_TPETRA_PRECONDITIONERFACTORY_DECL_HPP

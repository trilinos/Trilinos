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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef THYRA_MUELU_REFMAXWELL_PRECONDITIONER_FACTORY_DECL_HPP
#define THYRA_MUELU_REFMAXWELL_PRECONDITIONER_FACTORY_DECL_HPP

#include <MueLu_ConfigDefs.hpp>

#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)

// Stratimikos needs Thyra, so we don't need special guards for Thyra here
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_DiagonalLinearOpBase.hpp"
#include "Thyra_XpetraLinearOp.hpp"
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#ifdef HAVE_MUELU_EPETRA
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#endif

#include "Teuchos_Ptr.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Time.hpp"

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_ThyraUtils.hpp>

#include <MueLu_XpetraOperator_decl.hpp>  // todo fix me
#include <MueLu_RefMaxwell.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <Xpetra_TpetraOperator.hpp>
#include <Xpetra_TpetraHalfPrecisionOperator.hpp>
#ifdef HAVE_MUELU_EPETRA
#include <MueLu_EpetraOperator.hpp>
#include <Xpetra_EpetraOperator.hpp>
#endif

#include "Thyra_PreconditionerFactoryBase.hpp"

#include <Tpetra_KokkosCompat_DefaultNode.hpp>

#include <list>

namespace Thyra {

/** @brief Concrete preconditioner factory subclass for Thyra based on MueLu.
    @ingroup MueLuAdapters
    Add support for MueLu preconditioners in Thyra. This class provides an interface both
    for Epetra and Tpetra.

    The general implementation only handles Tpetra. For Epetra there is a specialization
    on SC=double, LO=int, GO=int and NO=EpetraNode.
*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class MueLuRefMaxwellPreconditionerFactory : public PreconditionerFactoryBase<Scalar> {
 public:
  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  MueLuRefMaxwellPreconditionerFactory();
  //@}

  /** @name Overridden from PreconditionerFactoryBase */
  //@{

  /** \brief . */
  bool isCompatible(const LinearOpSourceBase<Scalar>& fwdOp) const;
  /** \brief . */
  Teuchos::RCP<PreconditionerBase<Scalar> > createPrec() const;
  /** \brief . */
  void initializePrec(const Teuchos::RCP<const LinearOpSourceBase<Scalar> >& fwdOp,
                      PreconditionerBase<Scalar>* prec,
                      const ESupportSolveUse supportSolveUse) const;
  /** \brief . */
  void uninitializePrec(PreconditionerBase<Scalar>* prec,
                        Teuchos::RCP<const LinearOpSourceBase<Scalar> >* fwdOp,
                        ESupportSolveUse* supportSolveUse) const;

  //@}

  /** @name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList);
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

  /** \name Public functions overridden from Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  // ToDo: Add an override of describe(...) to give more detail!

  //@}

 private:
  Teuchos::RCP<Teuchos::ParameterList> paramList_;
};

}  // namespace Thyra

#endif  // #ifdef HAVE_MUELU_STRATIMIKOS

#endif  // THYRA_MUELU_REFMAXWELL_PRECONDITIONER_FACTORY_DECL_HPP

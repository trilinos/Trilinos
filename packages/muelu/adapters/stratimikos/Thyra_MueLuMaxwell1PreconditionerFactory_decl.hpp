// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_MUELU_MAXWELL1_PRECONDITIONER_FACTORY_DECL_HPP
#define THYRA_MUELU_MAXWELL1_PRECONDITIONER_FACTORY_DECL_HPP

#include <MueLu_ConfigDefs.hpp>

#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)

#include <Tpetra_KokkosCompat_DefaultNode.hpp>
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_XpetraLinearOp.hpp"
#include <Thyra_MueLuPreconditionerFactory.hpp>

namespace Thyra {

/** @brief Concrete preconditioner factory subclass for Thyra based on MueLu.
    @ingroup MueLuAdapters
    Add support for MueLu's Maxwell1 preconditioner in Thyra.
*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class MueLuMaxwell1PreconditionerFactory : public PreconditionerFactoryBase<Scalar> {
 public:
  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  MueLuMaxwell1PreconditionerFactory();
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

#endif  // THYRA_MUELU_MAXWELL1_PRECONDITIONER_FACTORY_DECL_HPP

// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_MUELU_PRECONDITIONER_FACTORY_DECL_HPP
#define THYRA_MUELU_PRECONDITIONER_FACTORY_DECL_HPP

#include <MueLu_ConfigDefs.hpp>

#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)

// Stratimikos needs Thyra, so we don't need special guards for Thyra here
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_DiagonalLinearOpBase.hpp"
#include "Thyra_XpetraLinearOp.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#ifdef HAVE_MUELU_EPETRA
#include "Thyra_EpetraLinearOp.hpp"
#endif

#include "Teuchos_Ptr.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Time.hpp"

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_ThyraUtils.hpp>

#include <MueLu_Hierarchy.hpp>
#include <MueLu_HierarchyManager.hpp>
#include <MueLu_HierarchyUtils.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_MasterList.hpp>
#include <MueLu_XpetraOperator_decl.hpp>  // todo fix me
#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <Xpetra_TpetraHalfPrecisionOperator.hpp>
#ifdef HAVE_MUELU_EPETRA
#include <MueLu_EpetraOperator.hpp>
#endif

#include "Thyra_PreconditionerFactoryBase.hpp"

#include <Tpetra_KokkosCompat_DefaultNode.hpp>

#include <list>

namespace Thyra {

using Teuchos::RCP;
using Teuchos::rcp;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
struct Converters {
  static bool replaceWithXpetra(ParameterList& paramList, std::string parameterName);
};

#ifdef HAVE_MUELU_EPETRA
template <class GlobalOrdinal>
struct Converters<double, int, GlobalOrdinal, Tpetra::KokkosCompat::KokkosSerialWrapperNode> {
  static bool replaceWithXpetra(ParameterList& paramList, std::string parameterName);
};
#endif

/** @brief Concrete preconditioner factory subclass for Thyra based on MueLu.
    @ingroup MueLuAdapters
    Add support for MueLu preconditioners in Thyra. This class provides an interface both
    for Epetra and Tpetra.
*/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class MueLuPreconditionerFactory : public PreconditionerFactoryBase<Scalar> {
 public:
  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  MueLuPreconditionerFactory();

  ~MueLuPreconditionerFactory();
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
  // Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> > CreateXpetraPreconditioner(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > op, const Teuchos::ParameterList& paramList, Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> > coords, Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > nullspace) const;

  Teuchos::RCP<Teuchos::ParameterList> paramList_;
};

}  // namespace Thyra

#endif  // #ifdef HAVE_MUELU_STRATIMIKOS

#endif  // THYRA_MUELU_PRECONDITIONER_FACTORY_DECL_HPP

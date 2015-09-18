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
#ifndef THYRA_MUELU_PRECONDITIONER_FACTORY_DECL_HPP
#define THYRA_MUELU_PRECONDITIONER_FACTORY_DECL_HPP

#include <MueLu_ConfigDefs.hpp>

#ifdef HAVE_MUELU_STRATIMIKOS

#include <MueLu_Hierarchy_fwd.hpp>
#include <MueLu_Utilities.hpp>

#include "Thyra_PreconditionerFactoryBase.hpp"

#include "Kokkos_DefaultNode.hpp"


namespace Thyra {

  /** \brief Concrete preconditioner factory subclass for Thyra based on MueLu.
   *
   * Add support for MueLu preconditioners in Thyra. This class provides an interface both
   * for Epetra and Tpetra
   */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class MueLuPreconditionerFactory : public PreconditionerFactoryBase<Scalar> {
  public:

    /** @name Constructors/initializers/accessors */
    //@{

    /** \brief . */
    MueLuPreconditionerFactory();
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
                        const ESupportSolveUse supportSolveUse
                       ) const;
    /** \brief . */
    void uninitializePrec(PreconditionerBase<Scalar>* prec,
                          Teuchos::RCP<const LinearOpSourceBase<Scalar> >* fwdOp,
                          ESupportSolveUse* supportSolveUse
                         ) const;

    //@}

    /** @name Overridden from Teuchos::ParameterListAcceptor */
    //@{

    /** \brief . */
    void                                          setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList);
    /** \brief . */
    Teuchos::RCP<Teuchos::ParameterList>          unsetParameterList();
    /** \brief . */
    Teuchos::RCP<Teuchos::ParameterList>          getNonconstParameterList();
    /** \brief . */
    Teuchos::RCP<const Teuchos::ParameterList>    getParameterList() const;
    /** \brief . */
    Teuchos::RCP<const Teuchos::ParameterList>    getValidParameters() const;
    //@}

    /** \name Public functions overridden from Describable. */
    //@{

    /** \brief . */
    std::string description() const;

    // ToDo: Add an override of describe(...) to give more detail!

    //@}

  private:

    Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> > CreateXpetraPreconditioner(Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > op, const Teuchos::ParameterList& paramList, Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > coords, Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > nullspace) const;

    Teuchos::RCP<Teuchos::ParameterList> paramList_;

  };

} // namespace Thyra

#endif // #ifdef HAVE_MUELU_STRATIMIKOS

#endif // THYRA_MUELU_PRECONDITIONER_FACTORY_DECL_HPP

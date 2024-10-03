// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_Interpolation_hpp__
#define __Panzer_Interpolation_hpp__

#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#endif

#ifdef PANZER_HAVE_EPETRA_STACK
#include "Thyra_EpetraThyraWrappers.hpp"
#endif

#include <Tpetra_Operator.hpp>
#include <Tpetra_Import.hpp>


namespace panzer {

  Teuchos::RCP<Thyra::LinearOpBase<double>> buildInterpolation(const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > &linObjFactory,
                                                               const std::string &domain_basis_name, const std::string &range_basis_name,
                                                               Intrepid2::EOperator op = Intrepid2::OPERATOR_VALUE,
                                                               size_t worksetSize = 1000,
                                                               const bool matrixFree = false);


  Teuchos::RCP<Thyra::LinearOpBase<double>> buildInterpolation(const Teuchos::RCP<const panzer::ConnManager> &conn,
                                                               const Teuchos::RCP<panzer::DOFManager> &domain_ugi,
                                                               const Teuchos::RCP<panzer::DOFManager> &range_ugi,
                                                               const std::string &domain_basis_name, const std::string &range_basis_name,
                                                               Intrepid2::EOperator op = Intrepid2::OPERATOR_VALUE,
                                                               size_t worksetSize = 1000, const bool forceVectorial = false,
                                                               const bool useTpetra = true,
                                                               const bool matrixFree = false);


  template <class Scalar = Tpetra::Operator<>::scalar_type,
            class LocalOrdinal = typename Tpetra::Operator<Scalar>::local_ordinal_type,
            class GlobalOrdinal = typename Tpetra::Operator<Scalar, LocalOrdinal>::global_ordinal_type,
            class Node = typename Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
  class MatrixFreeInterpolationOp : public Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

  public:

    typedef PHX::Device DeviceSpace;

    MatrixFreeInterpolationOp(const Teuchos::RCP<const panzer::ConnManager> &conn,
                              const Teuchos::RCP<panzer::DOFManager> &_domain_ugi,
                              const Teuchos::RCP<panzer::DOFManager> &_range_ugi,
                              const std::string& _domain_basis_name,
                              const std::string& _range_basis_name,
                              Intrepid2::EOperator _op=Intrepid2::OPERATOR_VALUE,
                              size_t _worksetSize=1000);

    void allocateColumnMapVector(size_t numVectors);

    // Pre-compute elements that own a DoF
    void precomputeOwnersAndOrientations(const Teuchos::RCP<const panzer::ConnManager> &conn);

    //! Returns the Tpetra::Map object associated with the domain of this operator.
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
    getDomainMap() const {
      return domainMap_;
    }

    //! Returns the Tpetra::Map object associated with the range of this operator.
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
    getRangeMap() const {
      return rangeMap_;
    }

    void
    apply (const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
           Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
           Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

    void
    applyNonTransposed (const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                        Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
                        Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                        Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

    void
    applyTransposed (const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                     Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
                     Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                     Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

    bool
    hasTransposeApply() const {
      return true;
    };

    void setName(std::string &_name) {
      name = _name;
    };

  private:

    std::string name;

    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > domainMap_;
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rangeMap_;
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > columnMap_;
    Teuchos::RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > import_;
    Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > colmapMV_;

    std::string domain_basis_name;
    std::string range_basis_name;
    Intrepid2::EOperator op;
    size_t worksetSize;

    Teuchos::RCP<panzer::DOFManager> domain_ugi;
    Teuchos::RCP<panzer::DOFManager> range_ugi;

    Kokkos::View<LocalOrdinal*,DeviceSpace> owner_d_;
    std::map<std::string, Kokkos::DynRankView<Intrepid2::Orientation,DeviceSpace> > orientations_;
  };
}
#endif

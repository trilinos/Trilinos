#ifndef _MiniEM_MatrixFreeInterpolationOp_hpp_
#define _MiniEM_MatrixFreeInterpolationOp_hpp_

#include <Tpetra_Operator.hpp>
#include <Tpetra_Import.hpp>

#include "Panzer_Traits.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Panzer_STKConnManager.hpp"


namespace mini_em {

  template <class Scalar = Tpetra::Operator<>::scalar_type,
            class LocalOrdinal = typename Tpetra::Operator<Scalar>::local_ordinal_type,
            class GlobalOrdinal = typename Tpetra::Operator<Scalar, LocalOrdinal>::global_ordinal_type,
            class Node = typename Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
  class MatrixFreeInterpolationOp : public Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

  public:

    typedef PHX::Device DeviceSpace;

    MatrixFreeInterpolationOp(const std::string& _name,
                              Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > _linObjFactory,
                              const std::string& _lo_basis_name,
                              const std::string& _ho_basis_name,
                              Intrepid2::EOperator _op=Intrepid2::OPERATOR_VALUE,
                              size_t _worksetSize=1000);

    void allocateColumnMapVector(size_t numVectors);

    // Pre-compute elements that own a DoF
    void precomputeOwners();

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

    void
    setupNodeOnlyConnManager();

  private:

    std::string name;

    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > domainMap_;
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rangeMap_;
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > columnMap_;
    Teuchos::RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > import_;
    Teuchos::RCP<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > colmapMV_;

    Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linObjFactory;
    std::string lo_basis_name;
    std::string ho_basis_name;
    Intrepid2::EOperator op;
    size_t worksetSize;

    Teuchos::RCP<const panzer::BlockedDOFManager> blockedDOFMngr;
    Teuchos::RCP<panzer::DOFManager> lo_ugi;
    Teuchos::RCP<panzer::DOFManager> ho_ugi;

    Kokkos::View<LocalOrdinal*,DeviceSpace> owner_d_;

    Teuchos::RCP<panzer_stk::STKConnManager> node_conn_;
  };

}
#endif

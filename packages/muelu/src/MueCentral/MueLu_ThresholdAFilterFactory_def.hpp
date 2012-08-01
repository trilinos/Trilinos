#ifndef MUELU_THRESHOLDAFILTERFACTORY_DEF_HPP
#define MUELU_THRESHOLDAFILTERFACTORY_DEF_HPP

#include <Xpetra_Operator.hpp>
#include <Xpetra_CrsOperator.hpp>

#include "MueLu_ThresholdAFilterFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ThresholdAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ThresholdAFilterFactory(const std::string& ename, const FactoryBase* fac, const Scalar threshold)
    : varName_(ename), factory_(fac), threshold_(threshold)
  { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  ThresholdAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~ThresholdAFilterFactory() {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ThresholdAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput(varName_,factory_,this);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void ThresholdAFilterFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & currentLevel) const {
    FactoryMonitor m(*this, "A filter (thresholding)", currentLevel);

    typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> OOperator; //TODO
    typedef Xpetra::CrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsOOperator; //TODO

    RCP<OOperator> Ain = currentLevel.Get< RCP<OOperator> >(varName_, factory_);

    // create new empty Operator
    RCP<CrsOOperator> Aout = rcp(new CrsOOperator(Ain->getRowMap(),Ain->getGlobalMaxNumRowEntries(),Xpetra::StaticProfile)); //FIXME

    // loop over local rows
    for(size_t row=0; row<Ain->getNodeNumRows(); row++)
      {
        size_t nnz = Ain->getNumEntriesInLocalRow(row);

        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        Ain->getLocalRowView(row, indices, vals);

        TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::ThresholdAFilterFactory::Build: number of nonzeros not equal to number of indices? Error.");

        Teuchos::ArrayRCP<GlobalOrdinal> indout(indices.size(),Teuchos::ScalarTraits<GlobalOrdinal>::zero());
        Teuchos::ArrayRCP<Scalar> valout(indices.size(),Teuchos::ScalarTraits<Scalar>::zero());
        size_t nNonzeros = 0;
        for(size_t i=0; i<(size_t)indices.size(); i++) {
          if(std::abs(vals[i]) > std::abs(threshold_) || indices[i]==(LocalOrdinal)row) {
            indout[nNonzeros] = Ain->getColMap()->getGlobalElement(indices[i]); // LID -> GID (column)
            valout[nNonzeros] = vals[i];
            nNonzeros++;
          }
        }

        indout.resize(nNonzeros);
        valout.resize(nNonzeros);

        Aout->insertGlobalValues(Ain->getRowMap()->getGlobalElement(row), indout.view(0,indout.size()), valout.view(0,valout.size()));
      }

    Aout->fillComplete(Ain->getDomainMap(), Ain->getRangeMap());

    GetOStream(Statistics0, 0) << "Nonzeros in " << varName_ << "(input): " << Ain->getGlobalNumEntries() << ", Nonzeros after filtering " << varName_ << " (parameter: " << threshold_ << "): " << Aout->getGlobalNumEntries() << std::endl;

    currentLevel.Set(varName_, Teuchos::rcp_dynamic_cast<OOperator>(Aout), this);
  }

} // namespace MueLu

#endif // MUELU_THRESHOLDAFILTERFACTORY_DEF_HPP

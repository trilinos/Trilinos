/*
 * MueLu_ThresholdAFilterFactory.hpp
 *
 *  Created on: 14.10.2011
 *      Author: tobias
 */

#ifndef MUELU_THRESHOLDAFILTERFACTORY_HPP_
#define MUELU_THRESHOLDAFILTERFACTORY_HPP_

#include <Teuchos_TestForException.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"



namespace MueLu {

  /*!
    @class ThresholdAFilterFactory class.
    @brief Factory for building a thresholded operator.

  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class ThresholdAFilterFactory : public SingleLevelFactoryBase {

    #include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    ThresholdAFilterFactory(const std::string& ename, const FactoryBase* fac, const Scalar threshold)
      : varName_(ename), factory_(fac), threshold_(threshold)
    { }

    //! Destructor.
    virtual ~ThresholdAFilterFactory() {}
    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const {
      currentLevel.DeclareInput(varName_,factory_);
    }

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    void Build(Level & currentLevel) const {
      RCP<Operator> Ain = currentLevel.Get< RCP<Operator> >(varName_, factory_);

      Monitor m(*this, "A filter (thresholding)");

      // create new empty Operator
      RCP<CrsOperator> Aout = rcp(new CrsOperator(Ain->getRowMap(),Ain->getGlobalMaxNumRowEntries(),Xpetra::StaticProfile));

      // loop over local rows
      for(size_t row=0; row<Ain->getNodeNumRows(); row++)
      {
        size_t nnz = Ain->getNumEntriesInLocalRow(row);

        Teuchos::ArrayView<const LocalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        Ain->getLocalRowView(row, indices, vals);

        TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::ThresholdAFilterFactory::Build: number of nonzeros not equal to number of indices? Error.");

        Teuchos::ArrayRCP<LocalOrdinal> indout(indices.size(),Teuchos::ScalarTraits<LocalOrdinal>::zero());
        Teuchos::ArrayRCP<Scalar> valout(indices.size(),Teuchos::ScalarTraits<Scalar>::zero());
        size_t nNonzeros = 0;
        for(size_t i=0; i<(size_t)indices.size(); i++) {
          if((Scalar)abs(vals[i]) > threshold_ || indices[i]==(LocalOrdinal)row) {
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

      currentLevel.Set(varName_, Teuchos::rcp_dynamic_cast<Operator>(Aout), this);
    }

    //@}

  private:
    std::string         varName_;   ///< name of input and output variable
    const FactoryBase*  factory_;   ///< generating factory of input variable
    const Scalar        threshold_; ///< threshold parameter


  }; // class ThresholdAFilterFactory

} //#ifndef MUELU_THRESHOLDAFILTERFACTORY_HPP_namespace MueLu

#define MUELU_THRESHOLDAFILTERFACTORY_SHORT

#endif /* MUELU_THRESHOLDAFILTERFACTORY_HPP_ */

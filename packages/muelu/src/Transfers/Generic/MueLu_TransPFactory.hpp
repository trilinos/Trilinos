#ifndef MUELU_TRANSPFACTORY_HPP
#define MUELU_TRANSPFACTORY_HPP

#include <iostream>

#include "Xpetra_CrsOperator.hpp"

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_RFactory.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  /*!
    @class TransPFactory class.
    @brief Factory for building restriction operators.

    This factory currently depends on an underlying matrix-matrix multiply with the identity
    matrix to do the transpose.  This should probably be fixed at some point.
  */

template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class TransPFactory : public RFactory {

#include "MueLu_UseShortNames.hpp"

    template<class AA, class BB, class CC, class DD, class EE>
    inline friend std::ostream& operator<<(std::ostream& os, TransPFactory<AA,BB,CC,DD,EE> &factory);

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    TransPFactory() {
      //Teuchos::OSTab tab(this->out_);
      //MueLu_cout(Teuchos::VERB_HIGH) << "TransPFactory: Instantiating a new factory" << std::endl;
    }

    //! Destructor.
    virtual ~TransPFactory() {}
    //@}

    //! Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const { }

    //@}

    //! @name Build methods.
    //@{
/*
   FIXME this uses the Tpetra RowMatrixTransposer.  This has revealed a bug somewhere.
   FIXME so disabling it right now.
    bool BuildR(Level & fineLevel, Level & coarseLevel) {

      RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("TransPFactory::BuildR"));
      timer->start(true);

      RCP<Operator> P = coarseLevel.Get< RCP<Operator> >("P");
      RCP<Operator> R = Utils::Transpose(P,true); //true indicated optimize storage
      coarseLevel.Set("R", R);

      timer->stop();
      MemUtils::ReportTimeAndMemory(*timer, *(P->getRowMap()->getComm()));

      return true;
    } //BuildR
*/

    bool Build(Level & fineLevel, Level & coarseLevel) const {
      return BuildR(fineLevel,coarseLevel);
    }

    bool BuildR(Level & fineLevel, Level & coarseLevel) const {

      std::ostringstream buf; buf << coarseLevel.GetLevelID();
      RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("TransPFactory::OldBuildR_"+buf.str()));
      timer->start(true);

      Teuchos::OSTab tab(this->out_);
      Teuchos::ParameterList matrixList;
      RCP<Operator> P = coarseLevel.Get< RCP<Operator> >("P");
      //doesn't work -- bug in EpetraExt?
      //RCP<Operator> I = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, CrsOperator>("Identity",P->getRangeMap(),matrixList);
      //      RCP<CrsOperator> I = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, CrsOperator>("Identity",P->getDomainMap(),matrixList);
      //RCP<Operator> R = Utils::TwoMatrixMultiply(P,true,I,false); //doesn't work -- bug in EpetraExt?
      //      RCP<Operator> R = Utils::TwoMatrixMultiply(I,false,P,true);

      RCP<Operator> R= Utils2<SC,LO,GO>::Transpose(P,true);

      coarseLevel.Set("R", R, this);

      timer->stop();
      MemUtils::ReportTimeAndMemory(*timer, *(P->getRowMap()->getComm()));

      return true;
    } //BuildR

    //@}

    //! @name Set methods.
    //@{
    void UsePtent(bool ToF) {
      throw(Exceptions::NotImplemented("TransPFactory.UsePtent()")); //TODO
    }
    //@}



  }; //class TransPFactory

  //! Friend print function.
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::ostream& operator<<(std::ostream& os, TransPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> &factory) {
    os << "Printing a TransPFactory object" << std::endl;
    return os;
  }

} //namespace MueLu

#define MUELU_TRANSPFACTORY_SHORT

#endif //ifndef MUELU_TRANSPFACTORY_HPP

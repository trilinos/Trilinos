#ifndef MUELU_RAPFACTORY_HPP
#define MUELU_RAPFACTORY_HPP

#include <iostream>
#include "MueLu_OperatorFactory.hpp"

namespace MueLu {
/*!
  @class RAPFactory class.
  @brief Factory for building coarse matrices.
*/
  template<class ScalarType, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
class RAPFactory : public OperatorFactory<ScalarType,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> {

#include "MueLu_UseShortNames.hpp"


  //JG to JJH: use Teuchos::Describable instead ?
  template<class AA, class BB, class CC, class DD, class EE>
  inline friend std::ostream& operator<<(std::ostream& os, RAPFactory<AA,BB,CC,DD,EE> &factory);

  public:
    //@{ Constructors/Destructors.
    RAPFactory() {}

    virtual ~RAPFactory() {}
    //@}

    //@{ Build methods.
    bool Build(Level &fineLevel, Level &coarseLevel) {
      Teuchos::OSTab tab(this->getOStream());
      MueLu_cout(Teuchos::VERB_HIGH) << "RAPFactory: Building a coarse operator" << std::endl;
#ifdef CTHULHU_USE_EPETRA
      //CTHULHU_RCP_DYNAMIC_CAST(Epetra_CrsMatrix, A, epA, "Problem casting A to Epetra_CrsMatrix");
      MueLu_cout(Teuchos::VERB_LOW) << "call the Epetra matrix-matrix multiply here" << std::endl;
#else
      MueLu_cout(Teuchos::VERB_LOW) << "Tpetra has no matrix-matrix multiply yet" << std::endl;
#endif

      return true;
    }
    //@}


}; //class RAPFactory

//! Friend print method.
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
std::ostream& operator<<(std::ostream& os, RAPFactory<Scalar,LocalOrdinal,GO,Node,LocalMatOps> &factory) {
  os << "Printing RAPFactory object" << std::endl;
  return os;
}

} //namespace MueLu

#define MUELU_RAPFACTORY_SHORT

#endif //ifndef MUELU_RAPFACTORY_HPP

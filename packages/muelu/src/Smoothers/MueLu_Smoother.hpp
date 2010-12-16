#ifndef MUELU_SMOOTHER_HPP
#define MUELU_SMOOTHER_HPP

#include "Teuchos_VerboseObject.hpp"

namespace MueLu {

/*!
  @class Smoother class.
  @brief Smoother class. Just a stub right now.
*/

  template<class ScalarType,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class Smoother : public Teuchos::VerboseObject<Smoother<ScalarType,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > {

#include "MueLu_UseShortNames.hpp"

  private:

  public:

    //! Constructor
    Smoother() {std::cout << "Instantiating a new smoother" << std::endl;}

    virtual ~Smoother() {}

}; //class Smoother

} //namespace MueLu

#define MUELU_SMOOTHER_SHORT

#endif //ifndef MUELU_SMOOTHER_HPP

#ifndef MUELU_SMOOTHER_HPP
#define MUELU_SMOOTHER_HPP

#include "Teuchos_VerboseObject.hpp"

namespace MueLu {

/*!
  @class Smoother class.
  @brief Smoother class. Just a stub right now.
*/

template<class Scalar,class LO, class GO, class Node>
class Smoother :  public Teuchos::VerboseObject<Smoother<Scalar,LO,GO,Node> > {

  private:

  public:

    //! Constructor
    Smoother() {std::cout << "Instantiating a new smoother" << std::endl;}

    virtual ~Smoother() {}

}; //class Smoother

} //namespace MueLu

#endif //ifndef MUELU_SMOOTHER_HPP

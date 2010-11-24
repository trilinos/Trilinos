#ifndef MUELU_SMOOTHER_HPP
#define MUELU_SMOOTHER_HPP

#include "Teuchos_VerboseObject.hpp"

namespace MueLu {

/*!
  @class Smoother class.
  @brief Smoother class. Just a stub right now.
*/

  template<class Scalar,class LO, class GO, class NO, class LMO>
  class Smoother : public Teuchos::VerboseObject<Smoother<Scalar,LO,GO,NO,LMO> > {

  private:

  public:

    //! Constructor
    Smoother() {std::cout << "Instantiating a new smoother" << std::endl;}

    virtual ~Smoother() {}

}; //class Smoother

} //namespace MueLu

#endif //ifndef MUELU_SMOOTHER_HPP

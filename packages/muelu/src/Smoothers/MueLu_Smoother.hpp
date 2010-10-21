#ifndef MUELU_SMOOTHER_HPP
#define MUELU_SMOOTHER_HPP

namespace MueLu {

template<class Scalar,class LO, class GO, class Node>
class Smoother {

  private:

  public:

    Smoother() {std::cout << "Instantiating a new smoother" << std::endl;}

    virtual ~Smoother() {}

}; //class Smoother

} //namespace MueLu

#endif //ifndef MUELU_SMOOTHER_HPP

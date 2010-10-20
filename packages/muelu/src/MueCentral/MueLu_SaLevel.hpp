#ifndef MUELU_SALEVEL_HPP
#define MUELU_SALEVEL_HPP

#include "Teuchos_RefCountPtr.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

template<class SC,class LO, class GO, class NO>
class SaLevel : public Level<SC,LO,GO,NO> {

  typedef Tpetra::CrsMatrix<SC,LO,GO,NO> Operator;
  typedef Tpetra::Vector<SC,LO,GO,NO> Vector;

  template<class AA, class BB, class CC, class DD>
  friend inline std::ostream& operator<<(std::ostream &os, SaLevel<AA,BB,CC,DD> const & level);

  private:
    Teuchos::RCP<Vector>   nullSpace_;
    //aggInfo_;
    Teuchos::RCP<Operator> Ptent_;
    Teuchos::RCP<Operator> Rtent_;

  public:

  //! Constructor
  SaLevel() : Ptent_(Teuchos::null), Rtent_(Teuchos::null) {}

  //! Destructor
  virtual ~SaLevel() {}

  //! Set nullspace.
  void SetNullSpace(Teuchos::RCP<Vector> const &nullspace) {
    nullSpace_ = nullspace;
  }

}; //class SaLevel

template<class SC,class LO, class GO, class NO>
std::ostream& operator<<(std::ostream &os, SaLevel<SC,LO,GO,NO> const & level)
{
  os << level.Ptent_;
  os << level.Rtent_;
  return os;
}

} //namespace MueLu


#endif //ifndef MUELU_SALEVEL_HPP

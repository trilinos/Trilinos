#ifndef MUELU_FAKESMOOTHERPROTOTYPE_HPP
#define MUELU_FAKESMOOTHERPROTOTYPE_HPP

// This is a minimal implementation of the SmootherPrototype interface. Used by unit tests. Do not use elsewhere.

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherPrototype.hpp"

namespace MueLu {

  class Level;

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class FakeSmootherPrototype : public SmootherPrototype<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> {

#include "MueLu_UseShortNames.hpp"

  public:

    FakeSmootherPrototype(int param=0) : param_(param) {}

    virtual ~FakeSmootherPrototype() {}

    virtual RCP<SmootherPrototype> Copy() const { return rcp(new FakeSmootherPrototype(*this)); }

    void Setup(Level &) { SmootherPrototype::IsSetup(true); }

    void Apply(MultiVector &x, MultiVector const &rhs, bool const &InitialGuessIsZero) { 
      TEST_FOR_EXCEPTION(1, Exceptions::NotImplemented, "MueLu::FakeSmootherPrototype()::Apply(): this class is for test purpose only.")
    }
    
    void SetParam(int param) { param_ = param; }

    int GetParam(int param) const { return param_; }

  private:
    int param_;
  
  };

}

#define MUELU_FAKESMOOTHERPROTOTYPE_SHORT

#endif //ifndef MUELU_FAKESMOOTHERPROTOTYPE_HPP

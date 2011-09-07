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

    FakeSmootherPrototype(int param=0) : param_(param), numOfSetup_(0), numOfSetupCall_(0) {}

    virtual ~FakeSmootherPrototype() {}

    virtual RCP<SmootherPrototype> Copy() const { 
      TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == true, Exceptions::RuntimeError, "Not a prototype. Do not copy"); // test not mandatory, but it is the only use case that we need.
      return rcp(new FakeSmootherPrototype(*this)); 
    }

    void Setup(Level &) { 
      numOfSetupCall_++; 
      if (SmootherPrototype::IsSetup()) return;

      numOfSetup_++;

      SmootherPrototype::IsSetup(true); 
    }

    void Apply(MultiVector &x, MultiVector const &rhs, bool const &InitialGuessIsZero) { 
      TEST_FOR_EXCEPTION(1, Exceptions::NotImplemented, "MueLu::FakeSmootherPrototype()::Apply(): this class is for test purpose only.")
    }
    
    void SetParam(int param) { param_ = param; }

    int GetParam() const { return param_; }

    int GetNumOfSetup() const { return numOfSetup_; }
    int GetNumOfSetupCall() const { return numOfSetupCall_; }

  private:
    int param_;
  
    int numOfSetup_;
    int numOfSetupCall_;
  };

}

#define MUELU_FAKESMOOTHERPROTOTYPE_SHORT

#endif //ifndef MUELU_FAKESMOOTHERPROTOTYPE_HPP

#ifndef MUELU_SMOOTHERFACTORY_HPP
#define MUELU_SMOOTHERFACTORY_HPP

#include <iostream>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherFactoryBase.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

  /*!
    @class SmootherFactory
    @brief Generic Smoother Factory for generating the smoothers of the MG hierarchy

    This factory is generic and can produce any kind of Smoother.
 
    The design of the Smoother factory is based on the prototype design
    pattern.  The smoother factory uses prototypical instances of
    smoothers, and prototypes are cloned to produce the smoothers of the MG
    hierarchy.
 
    See also:
    - http://en.wikipedia.org/wiki/Factory_method_pattern
    - http://en.wikipedia.org/wiki/Prototype_pattern
    - http://www.oodesign.com/prototype-pattern.html
 
    There is one prototype for the pre-smoother and one for the
    post-smoother. Thus, this factory can produce two different smoothers
    for pre and post-smoothing. Prototypes are stored in this factory.
 
    The type of smoother to create is determined by the prototypical
    instances. Prototypes store their own parameters (nits, omega), and the
    factory doesn't have any knowledge about that.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class SmootherFactory : public SmootherFactoryBase<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> {

#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    /*!
      @brief Constructor

      TODO update the documentation from Matlab syntax

      The constructor takes as arguments pre-configured smoother
      object(s). They are used as prototypes by the factory.
      
      SYNTAX

      SmootherFactory(PreSmootherPrototype_, PostSmootherPrototype_);
      
      PreSmootherPrototype_  - prototype for pre-smoothers  (SmootherPrototype)
      PostSmootherPrototype_ - prototype for post-smoothers
      (SmootherPrototype,optional,default=PreSmootherPrototype_)

      EXAMPLES:

      SmooFactory = SmootherFactory(ChebySmoother(5,1/30))

      or 

      nIts        = 5;
      lambdaRatio = 1/30;
      Smoo        = ChebySmoother()
      Smoo        = Smoo.SetIts(nIts);
      Smoo        = Smoo.SetLambdaRatio(1/30); 
      SmooFactory = SmootherFactory(Smoo);

      To use different smoothers for pre and post smoothing, two prototypes can be passed in as argument:
      PreSmoo     = ChebySmoother(2, 1/30)
      PostSmoo    = ChebySmoother(10,1/30)
      SmooFactory = SmootherFactory(PreSmoo, PostSmoo);

      [] is also a valid smoother which do nothing:
      PostSmoo    = ChebySmoother(10,1/30)
      SmooFactory = SmootherFactory([], PostSmoo);
    */
    SmootherFactory(RCP<SmootherPrototype> preProto=Teuchos::null, RCP<SmootherPrototype> postProto = Teuchos::null)
    {
      if (preProto == Teuchos::null)
        throw(Exceptions::RuntimeError("Presmoother prototype cannot be null")); //TODO
      PreSmootherPrototype_ = preProto;
      if (postProto == Teuchos::null)
        PostSmootherPrototype_ = preProto;
      else
        PostSmootherPrototype_ = postProto;
    }

    //! Copy constructor.
    SmootherFactory(SmootherFactory &smooFact)
    {
      throw(Exceptions::NotImplemented("copy constructor"));
    }

    virtual ~SmootherFactory() {}
    //@}

    //! @name Build methods.
    //@{

    /*! @brief Creates pre and post smoothers.

    Factory.Build() clones its prototypes and calls Setup() on the
    new objects to create fully functional smoothers for the current
    hierarchy level. Note that smoothers do their own setup
    (ie: ILUSmoother.Setup() computes the ILU factorization). 
       
    If pre and post smoother are identical, the Setup() phase
    is done only once: to create the post-smoother, the
    pre-smoother is duplicated and its parameters are changed
    according to the parameters of the post-smoother prototype.
       
    If the parameters of pre and post smoothers are not the same,
    the Setup() phase is also done only once when parameters
    don't change the result of the setup computation.
    */
    bool Build(Level& currentLevel) const {
      return BuildSmoother(currentLevel, BOTH);
    }
    
    bool BuildSmoother(Level & currentLevel, PreOrPost const &pop = BOTH) const {
      RCP<SmootherPrototype> preSmoo;
      RCP<SmootherPrototype> postSmoo;
      
      if ((pop == BOTH || pop == PRE) && (PreSmootherPrototype_ != Teuchos::null)) {
        preSmoo = PreSmootherPrototype_->Copy();
        //preSmoo = rcp( new SmootherPrototype(PreSmootherPrototype_) );
        //TODO if outputlevel high enough
        //TODO preSmoo.Print();
        preSmoo->Setup(currentLevel);
        
        // Level Set
        currentLevel.NewSet("PreSmoother", preSmoo, this);
      }
      
      if ((pop == BOTH || pop == POST) && (PostSmootherPrototype_ != Teuchos::null))
          {
            if (pop == BOTH && PreSmootherPrototype_ == PostSmootherPrototype_) {
              
              // Very simple reuse. TODO: should be done in MueMat too
              postSmoo = preSmoo;
              
            } else if (pop == BOTH &&
                       PreSmootherPrototype_ != Teuchos::null &&
                       PreSmootherPrototype_->GetType() == PostSmootherPrototype_->GetType()) {
              
              // More complex reuse case: need implementation of CopyParameters() and a smoothers smart enough to know when parameters affect the setup phase.
              
              // YES: post-smoother == pre-smoother 
              // => copy the pre-smoother to avoid the setup phase of the post-smoother.
              postSmoo = preSmoo->Copy();
              // If the post-smoother parameters are different from
              // pre-smoother, the parameters stored in the post-smoother
              // prototype are copied in the new post-smoother object.
              postSmoo->CopyParameters(PostSmootherPrototype_);
              // If parameters don't influence the Setup phase (it is the case
              // for Jacobi, Chebyshev...), PostSmoo is already setup. Nothing
              // more to do. In the case of ILU, parameters of the smoother
              // are in fact the parameters of the Setup phase. The call to
              // CopyParameters resets the smoother (only if parameters are
              // different) and we must call Setup() again.
              postSmoo->Setup(currentLevel);

              // TODO: if CopyParameters do not exist, do setup twice.

            } else {
              
              // NO reuse: pop==POST or post-smoother != pre-smoother
              // Copy the prototype and run the setup phase.
              postSmoo = PostSmootherPrototype_->Copy();
              postSmoo->Setup(currentLevel);
              
            }
            
            // Level Set
            std::cout << "POST" << std::endl;
            currentLevel.NewSet("PostSmoother", postSmoo, this);
          }
        
        return true;//?
        
    } //Build()
    
    //@}

    //! @name Set/Get methods.
    //@{

    //! Set smoother prototypes.
    void SetSmootherPrototypes(RCP<SmootherPrototype> &preProto, RCP<SmootherPrototype> &postProto)
    {
      PreSmootherPrototype_ = preProto;
      PostSmootherPrototype_ = postProto;
    }

    //! Get smoother prototypes.
    void GetSmootherPrototypes(RCP<SmootherPrototype> &preProto, RCP<SmootherPrototype> &postProto)
    {
      preProto = PreSmootherPrototype_;
      postProto = PostSmootherPrototype_;
    }

  private:
    RCP<SmootherPrototype> PreSmootherPrototype_;
    RCP<SmootherPrototype> PostSmootherPrototype_;

    //@}

  }; //class SmootherFactory

} //namespace MueLu

#define MUELU_SMOOTHERFACTORY_SHORT

#endif //ifndef MUELU_SMOOTHERFACTORY_HPP

// TODO: add a simpler test (PreSmoo==PostSmoo) to try to reuse directly presmoo for postsmoo even if CopyParameter is not implemented

// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SMOOTHERFACTORY_DECL_HPP
#define MUELU_SMOOTHERFACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"

#include "MueLu_SmootherFactoryBase.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_SmootherPrototype_fwd.hpp"
#include "MueLu_SmootherBase_fwd.hpp"

namespace MueLu {

/*!
  @class SmootherFactory
  @ingroup MueLuSmootherClasses
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

template <class Scalar        = double,
          class LocalOrdinal  = int,
          class GlobalOrdinal = LocalOrdinal,
          class Node          = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class SmootherFactory : public SmootherFactoryBase {
 public:
  typedef Scalar scalar_type;
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node node_type;

 private:
#undef MUELU_SMOOTHERFACTORY_SHORT
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

    SmootherFactory(preSmootherPrototype_, postSmootherPrototype_);

    preSmootherPrototype_  - prototype for pre-smoothers  (SmootherPrototype)
    postSmootherPrototype_ - prototype for post-smoothers
    (SmootherPrototype,optional,default=preSmootherPrototype_)

    EXAMPLES:

    SmootherFactory = SmootherFactory(ChebySmoother(5,1/30))

    or

    nIts        = 5;
    lambdaRatio = 1/30;
    Smoother        = ChebySmoother()
    Smoother        = Smoother.SetIts(nIts);
    Smoother        = Smoother.SetLambdaRatio(1/30);
    SmootherFactory = SmootherFactory(Smoother);

    To use different smoothers for pre and post smoothing, two prototypes can be passed in as argument:
    PreSmoother     = ChebySmoother(2, 1/30)
    PostSmoother    = ChebySmoother(10,1/30)
    SmootherFactory = SmootherFactory(PreSmoother, PostSmoother);

    [] is also a valid smoother which do nothing:
    PostSmoother    = ChebySmoother(10,1/30)
    SmootherFactory = SmootherFactory([], PostSmoother);
  */
  // Note: Teuchos::null as parameter allowed (= no smoother)
  // Note: precondition: input smoother must not be Setup()
  SmootherFactory(RCP<SmootherPrototype> preAndPostSmootherPrototype = Teuchos::null);

  SmootherFactory(RCP<SmootherPrototype> preSmootherPrototype, RCP<SmootherPrototype> postSmootherPrototype);

  virtual ~SmootherFactory();
  //@}

  //! @name Set/Get methods.
  //@{

  //! Set smoother prototypes.
  void SetSmootherPrototypes(RCP<SmootherPrototype> preAndPostSmootherPrototype);

  //! Set smoother prototypes.
  void SetSmootherPrototypes(RCP<SmootherPrototype> preSmootherPrototype, RCP<SmootherPrototype> postSmootherPrototype);

  //! Get smoother prototypes.
  void GetSmootherPrototypes(RCP<SmootherPrototype>& preSmootherPrototype, RCP<SmootherPrototype>& postSmootherPrototype) const;

  //! Input
  //@{

  RCP<const ParameterList> GetValidParameterList() const;

  void DeclareInput(Level& currentLevel) const;

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
  void Build(Level& currentLevel) const;

  void BuildSmoother(Level& currentLevel, const PreOrPost preOrPost = BOTH) const;  // Build()

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const;

  //! Print the object with some verbosity level to an FancyOStream object.
  using MueLu::Describable::describe;  // overloading, not hiding
  void describe(Teuchos::FancyOStream& out, const VerbLevel verbLevel = Default) const;

  //@}

 private:
  RCP<SmootherPrototype> preSmootherPrototype_;
  RCP<SmootherPrototype> postSmootherPrototype_;

  void CheckPrototypes() const;

  //@}

};  // class SmootherFactory

}  // namespace MueLu

// TODO: doc: setup done twice if PostSmoother object != PreSmoother object and no adv. reused capability
#define MUELU_SMOOTHERFACTORY_SHORT
#endif  // MUELU_SMOOTHERFACTORY_DECL_HPP

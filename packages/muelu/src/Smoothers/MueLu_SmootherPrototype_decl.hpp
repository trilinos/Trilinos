// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_SMOOTHERPROTOTYPE_DECL_HPP
#define MUELU_SMOOTHERPROTOTYPE_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherPrototype_fwd.hpp"
#include "MueLu_SmootherBase.hpp"
#include "MueLu_Factory.hpp"

namespace MueLu {

class Level;

/*!
  @class SmootherPrototype
  @ingroup MueLuSmootherClasses
  @brief Base class for smoother prototypes

  A smoother prototype is a smoother which can be in two states:
  - ready to be duplicated (parameters defined)
  - ready to be used (setup phase completed)

  'Smoother prototypes' can be fully copied using the Copy() method.
*/

template <class Scalar        = SmootherBase<>::scalar_type,
          class LocalOrdinal  = typename SmootherBase<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename SmootherBase<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class SmootherPrototype : public SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node>,
                          public Factory {
 public:
  typedef Scalar scalar_type;
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node node_type;

 private:
#undef MUELU_SMOOTHERPROTOTYPE_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //!@nameConstructors/Destructors.
  //@{

  SmootherPrototype();

  virtual ~SmootherPrototype();

  //@}

  //! Input
  //@{

  virtual void DeclareInput(Level &currentLevel) const = 0;

  //@}

  //! @name Build methods.
  //@{

  virtual void Setup(Level &) = 0;

  virtual RCP<SmootherPrototype> Copy() const = 0;

  //@}

  //! @name Get/Set methods.
  //@{

  //! Get the state of a smoother prototype.
  bool IsSetup() const;

  //! Set the state of a smoother prototype.
  // Developpers: this method must be called by your Setup() method.
  void IsSetup(bool const &ToF);

  //@}

  //! @name Implements FactoryBase interface
  //@{
  virtual void CallBuild(Level & /* requestedLevel */) const {
    TEUCHOS_TEST_FOR_EXCEPT(true);
  }

  //!
  virtual void CallDeclareInput(Level &requestedLevel) const {
    DeclareInput(requestedLevel);
  }

  //@}

 private:
  bool isSetup_;

};  // class SmootherPrototype

}  // namespace MueLu

// TODO: private copy constructor
// TODO: update comments

#define MUELU_SMOOTHERPROTOTYPE_SHORT
#endif  // MUELU_SMOOTHERPROTOTYPE_DECL_HPP

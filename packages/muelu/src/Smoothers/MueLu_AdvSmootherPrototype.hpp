// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_ADVSMOOTHERPROTOTYPE_HPP
#define MUELU_ADVSMOOTHERPROTOTYPE_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SmootherPrototype.hpp"

namespace MueLu {

class Level;

/*!
  @class AdvSmootherPrototype

  'Advanced Smoother prototypes' can be fully copied using the Copy() method.
  They can also copy the parameters of
  another smoother object of the same type (CopyParameters()). Both
  capabilities are used by the SmootherFactory.
*/

template <class Scalar        = SmootherPrototype<>::scalar_type,
          class LocalOrdinal  = typename SmootherPrototype<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename SmootherPrototype<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class AdvSmootherPrototype : public SmootherPrototypex<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_ADVSMOOTHERPROTOTYPE_HPP
#include "MueLu_UseShortNames.hpp"

 public:
  //!@nameConstructors/Destructors.
  //@{
  AdvSmootherPrototype()
    : type_("undefined") {}

  virtual ~AdvSmootherPrototype() {}
  //@}

  //! @name Build methods.
  //@{

  virtual void CopyParameters(const AdvSmootherPrototype& smootherPrototype) = 0;

  //@}

  //! @name Get/Set methods.
  //@{

  //! Get the smoother type.
  std::string GetType() const { return type_; }

  /*! @brief Set the smoother type.
    This method must be called by constructors of derived classes.
  */
  // TODO: remove, type_ should be const
  void SetType(std::string& type) { type_ = type; }

  //@}

 private:
  std::string type_;

};  // class AdvSmootherPrototype

}  // namespace MueLu

#define MUELU_ADVSMOOTHERPROTOTYPE_SHORT

#endif  // ifndef MUELU_ADVSMOOTHERPROTOTYPE_HPP

// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_DESCRIBABLE_DECL_HPP
#define MUELU_DESCRIBABLE_DECL_HPP

#include <string>                      // for string
#include "Teuchos_FancyOStream.hpp"    // for FancyOStream
#include "Teuchos_VerbosityLevel.hpp"  // for EVerbosityLevel
#include "Teuchos_Describable.hpp"

#include "MueLu_VerbosityLevel.hpp"

namespace MueLu {

/*!
   @class Describable
   @brief Base class for MueLu classes

   @ingroup MueLuBaseClasses
*/
class Describable
  : public Teuchos::Describable {
  mutable std::string shortClassName_ = "";  // cached so that we don't have to call demangleName() every time; mutable so that ShortClassName() can initialize lazily while remaining const

 public:
  //! Destructor.
  virtual ~Describable();

  //! @name MueLu Describe
  //@{

  virtual void describe(Teuchos::FancyOStream &out_arg, const VerbLevel verbLevel = Default) const;

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  virtual std::string description() const;

  //! Print the object with some verbosity level to an FancyOStream object.
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const;

  //@}

  //! Return the class name of the object, without template parameters and without namespace
  virtual std::string ShortClassName() const;

};  // class Describable

}  // namespace MueLu

#define MUELU_DESCRIBABLE_SHORT
#endif  // MUELU_DESCRIBABLE_DECL_HPP

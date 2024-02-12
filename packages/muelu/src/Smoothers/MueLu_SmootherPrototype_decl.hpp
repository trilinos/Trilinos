// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
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

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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
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

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class AdvSmootherPrototype : public SmootherPrototypex<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> {
#undef MUELU_ADVSMOOTHERPROTOTYPE_HPP
#include "MueLu_UseShortNames.hpp"

  public:
    //!@nameConstructors/Destructors.
    //@{
    AdvSmootherPrototype()
      : type_("undefined")
    {}

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
    //TODO: remove, type_ should be const
    void SetType(std::string & type) { type_ = type; }

    //@}

  private:
    std::string type_;

  }; //class AdvSmootherPrototype

} //namespace MueLu

#define MUELU_ADVSMOOTHERPROTOTYPE_SHORT

#endif //ifndef MUELU_ADVSMOOTHERPROTOTYPE_HPP

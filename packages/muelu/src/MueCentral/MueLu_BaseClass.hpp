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
#ifndef MUELU_BASECLASS_HPP
#define MUELU_BASECLASS_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_VerboseObject.hpp"
#include "MueLu_Describable.hpp"

namespace MueLu {

  /*!
     @class BaseClass class.
     @brief Base class for MueLu classes

     @ingroup MueLuBaseClasses
  */
  class BaseClass
    : public VerboseObject, public Describable
  {

  public:

    //! @name Constructors/Destructors
    //@{

    //! Destructor.
    virtual ~BaseClass() {}

    //@}

  }; // class BaseClass

} // namespace MueLu

//! Helper macro for implementing Describable::describe() for BaseClass objects.
//  This macro defines ostream out0 that print only on root node. It print description() and indent the ostream.
//  Note: Runtime1 displays basic parameter information when Parameters0 is not enabled.
#define MUELU_DESCRIBE                                                  \
  using std::endl;                                                      \
  Teuchos::FancyOStream& out0 = (VerboseObject::GetProcRankVerbose() == 0) ? out : VerboseObject::GetBlackHole(); \
                                                                        \
  if ((verbLevel & Runtime1) && (!(verbLevel & Parameters0)))           \
    out << description() << std::endl;                                  \
  else if (verbLevel & Runtime0)                                        \
    out << BaseClass::description() << std::endl;                       \
                                                                        \
  Teuchos::OSTab tab1(out);                                             \
  //

#define MUELU_BASECLASS_SHORT
#endif // ifndef MUELU_BASECLASS_HPP

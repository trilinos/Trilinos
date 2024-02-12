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
#ifndef MUELU_DESCRIBABLE_HPP
#define MUELU_DESCRIBABLE_HPP

#include "MueLu_Describable.hpp"

namespace MueLu {

Describable::~Describable() {}

void Describable::describe(Teuchos::FancyOStream &out_arg, const VerbLevel /* verbLevel */) const {
  Teuchos::RCP<Teuchos::FancyOStream> out = rcp(&out_arg, false);  // JG: no idea why we have to do that, but it's how Teuchos::Describable::describe() is implemented
  Teuchos::OSTab tab(out);
  *out << this->description() << std::endl;
}

std::string Describable::description() const {
  std::string str = Teuchos::Describable::description();

  // remove template parameters
  size_t found = str.find_first_of("<");
  if (found != std::string::npos)
    return str.substr(0, found);

  return str;
}

void Describable::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const { describe(out, toMueLuVerbLevel(verbLevel)); }

std::string Describable::ShortClassName() const {
  if (shortClassName_.empty()) {
    std::string str = Teuchos::Describable::description();

    // remove template parameters
    {
      size_t found = str.find_first_of("<");
      if (found != std::string::npos)
        str = str.substr(0, found);
    }

    // remove namespace
    {
      size_t found = str.find_last_of(":");
      if (found != std::string::npos)
        str = str.substr(found + 1);
    }
    shortClassName_ = str;
  }
  return shortClassName_;
}

}  // namespace MueLu

#define MUELU_DESCRIBABLE_SHORT
#endif  // MUELU_DESCRIBABLE_HPP

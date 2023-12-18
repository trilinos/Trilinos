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
// Helper to get ride of template parameters

// This file can be use for two purpose:
// 1) As an header of a user program.
//    In this case, this file must be included *after* other headers
//    and the types Scalar, LocalOrdinal, GlobalOrdinal, Node must be defined.
//    Note also that there is no #ifndef/#endif to protect again the multiple inclusion of this file.
//    User should create is own header file including this one:
//
//    Example:
//     #ifndef MY_HEADER
//     #define MY_HEADER
//     #include <MueLu_UseDefaultTypes.hpp>
//     #include <MueLu_UseShortNames.hpp>
//     #endif
//
// 2) Inside of MueLu to enhance the readability.
//
// template <class Scalar = Xpetra::MultiVector<>::scalar_type,
//           class LocalOrdinal = typename Xpetra::MultiVector<Scalar>::local_ordinal_type,
//           class GlobalOrdinal = typename Xpetra::MultiVector<Scalar, LocalOrdinal>::global_ordinal_type,
//           class Node = typename Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
//  class TpetraMultiVector : public virtual Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
//
//  #include <MueLu_UseShortNames.hpp>
//
//  myMethod(RCP<const Map> & map) { [...] } // instead of myMethod(RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map)
//
//  [...]
//
// }
//

#include "MueLu_UseShortNamesOrdinal.hpp"
#include "MueLu_UseShortNamesScalar.hpp"

//! @file MueLu_UseShortNamesOrdinal.hpp

// TODO / NOTE: This file should not be included at the global scope (to avoid name collision)

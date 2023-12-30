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
#ifndef MUELU_KEEPTYPE_HPP
#define MUELU_KEEPTYPE_HPP

namespace MueLu {

//! Keep status of a variable of Level.
// Several keep status can be set at the same time
enum KeepEnum {
  UserData = 0x1,  //!< User data are always kept. This flag is set automatically when Level::Set("data", data) is used. The keep status of the variable is not propagated to coarser level (if you use Level::Build()).
  Keep     = 0x2,  //!< Always keep data, even accross run. This flag is set by Level::Keep(). This flag is propagated to coarser level by Level::Build().
  Final    = 0x4,  //!< Keep data only for this run. Used to keep data useful for Hierarchy::Iterate(). Data will be deleted if setup phase is re-run. This flag is set by default for A, P, R, PreSmoother and PostSmoother of NoFactory by Hierarchy::Setup(). Not propagated by Level::Build().

  NextRun = UserData | Keep,  //!< Both UserData and Keep flags force data to be kept and reused for the next run. Do not use MueLu::NextRun in AddKeepFlag. Use it only for testing keep == UserData || keep == Keep.
  All     = UserData | Keep | Final
};

//!
typedef short KeepType;  // TODO: name it KeepFlag?

}  // namespace MueLu

#endif  // MUELU_KEEPTYPE_HPP

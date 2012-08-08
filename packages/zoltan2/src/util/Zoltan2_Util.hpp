// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_Util.hpp
 *  \brief A gathering of useful namespace methods.
 *  \todo Should each class of utility functions be in a separate source file
 *         instead of having a source file with the unhelpful name of Util?
 */

#ifndef ZOLTAN2_UTIL_HPP
#define ZOLTAN2_UTIL_HPP

#include <Zoltan2_Standards.hpp>
#include <Teuchos_DefaultComm.hpp>

namespace Zoltan2{

long getProcessKilobytes();

template <typename scalar_t>
  inline bool outsideRegion(scalar_t val, scalar_t mark, double epsilon){
    return ((val < mark-epsilon) || (val > mark+epsilon));
}

#ifdef HAVE_ZOLTAN2_MPI

RCP<Teuchos::MpiComm<int> > MPI2Teuchos(const MPI_Comm &comm);

RCP<const Teuchos::MpiComm<int> > MPI2TeuchosConst(const MPI_Comm &comm);

MPI_Comm Teuchos2MPI(const RCP<Comm<int> > &comm);

MPI_Comm TeuchosConst2MPI(const RCP<const Comm<int> > &comm);

#endif

} // namespace Zoltan2

#endif

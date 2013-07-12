// @HEADER
// ***********************************************************************
//
//            Domi: Multidimensional Datastructures Package
//                 Copyright (2013) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef DOMI_MDCOMM_HPP
#define DOMI_MDCOMM_HPP

// Teuchos includes
#include "Teuchos_Comm.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"

// Domi includes
#include "Domi_ConfigDefs.hpp"
#include "Domi_Slice.hpp"

namespace Domi
{

typedef Teuchos::RCP< const Teuchos::Comm<int> > TeuchosComm;

class MDComm
{
public:

  typedef Teuchos::Array< int >::size_type size_type;

  MDComm(const TeuchosComm teuchosComm,
         const Teuchos::ArrayView< int > & axisSizes);
  MDComm(const TeuchosComm teuchosComm,
         int numDims);
  MDComm(const TeuchosComm teuchosComm,
         int numDims,
         const Teuchos::ArrayView< int > & axisSizes);
  MDComm(const Teuchos::RCP< MDComm > parent,
         const Teuchos::ArrayView< Slice > & slices);
  ~MDComm();
  TeuchosComm getTeuchosComm() const;
  int getNumDims() const;
  int getAxisSize(int axis) const;
  int getAxisRank(int axis) const;
  int getLowerNeighbor(int axis) const;
  int getUpperNeighbor(int axis) const;

protected:

  MDComm();

private:

  TeuchosComm           _teuchosComm;
  Teuchos::Array< int > _axisSizes;
  Teuchos::Array< int > _axisRanks;
  Teuchos::Array< int > _axisStrides;

};

}  // namespace Domi

#endif

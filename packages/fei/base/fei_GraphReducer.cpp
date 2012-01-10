/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#include <fei_macros.hpp>

#include <fei_GraphReducer.hpp>
#include <fei_TemplateUtils.hpp>
#include <fei_VectorSpace.hpp>

#undef fei_file
#define fei_file "fei_GraphReducer.cpp"
#include <fei_ErrMacros.hpp>

//----------------------------------------------------------------------------
fei::GraphReducer::GraphReducer(fei::SharedPtr<fei::Reducer> reducer,
                              fei::SharedPtr<fei::Graph> target)
  : reducer_(reducer),
    target_(target)
{
}

//----------------------------------------------------------------------------
fei::GraphReducer::~GraphReducer()
{
}

//----------------------------------------------------------------------------
int fei::GraphReducer::addIndices(int row, int len, const int* indices)
{
  reducer_->addGraphIndices(1, &row, len, indices, *target_);
  return(0);
}

//----------------------------------------------------------------------------
int fei::GraphReducer::addSymmetricIndices(int numIndices, int* indices,
                                         bool diagonal)
{
  reducer_->addSymmetricGraphIndices(numIndices, indices, diagonal, *target_);
  return(0);
}

//----------------------------------------------------------------------------
int fei::GraphReducer::writeLocalGraph(FEI_OSTREAM& os, bool debug,
				    bool prefixLinesWithPoundSign)
{
  return(target_->writeLocalGraph(os, debug, prefixLinesWithPoundSign));
}

//----------------------------------------------------------------------------
int fei::GraphReducer::writeRemoteGraph(FEI_OSTREAM& os)
{
  return(target_->writeRemoteGraph(os));
}

//----------------------------------------------------------------------------
int fei::GraphReducer::gatherFromOverlap()
{
  reducer_->assembleReducedGraph(target_.get(), false);
  return(target_->gatherFromOverlap());
}


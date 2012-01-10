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

#include <snl_fei_PointBlockMap.hpp>
#include <fei_ctg_set.hpp>
#undef fei_file
#define fei_file "snl_fei_PointBlockMap.cpp"
#include <fei_ErrMacros.hpp>

//----------------------------------------------------------------------------
snl_fei::PointBlockMap::PointBlockMap()
  : ptEqns_(NULL),
    blkEqns_(NULL),
    maxSize_(0),
    ptEqualBlk_(false)
{
  ptEqns_ = new std::map<int,int>;
  blkEqns_ = new std::map<int,std::pair<int,int> >;
}

//----------------------------------------------------------------------------
snl_fei::PointBlockMap::~PointBlockMap()
{
  delete ptEqns_;
  delete blkEqns_;
}

//----------------------------------------------------------------------------
void snl_fei::PointBlockMap::setPtEqualBlk()
{
  ptEqualBlk_ = true;
}

//----------------------------------------------------------------------------
int snl_fei::PointBlockMap::setEqn(int ptEqn, int blkEqn)
{
  return( setEqn(ptEqn, blkEqn, 1) );
}

//----------------------------------------------------------------------------
int snl_fei::PointBlockMap::setEqn(int ptEqn, int blkEqn, int blkSize)
{
  if (ptEqualBlk_ == true) {
    if (ptEqn != blkEqn) return(-1);
    else return(0);
  }

  ptEqns_->insert(std::pair<int,int>(ptEqn, blkEqn));

  //check whether blkEqn is already stored in blkEqns_.
  //if it is not, then insert it along with the pair ptEqn,blkSize.
  //if it is, check whether the already-associated ptEqn is greater than
  //the incoming ptEqn, and replace if so.
  //We want to have blkEqn mapped to the lower ptEqn.

  std::pair<int, int> newpair;
  std::map<int,std::pair<int,int> >::iterator
    b_iter = blkEqns_->find(blkEqn);

  if (b_iter == blkEqns_->end()) {
    newpair.first = ptEqn;
    newpair.second = blkSize;
    blkEqns_->insert(std::pair<int,std::pair<int,int> >(blkEqn, newpair));
  }
  else {
    newpair = (*b_iter).second;
    if (newpair.first > ptEqn) {
      newpair.first = ptEqn;
      newpair.second = blkSize;
      (*b_iter).second = newpair;
    }
  }

  return(0);
}

//----------------------------------------------------------------------------
int snl_fei::PointBlockMap::setBlkEqnSize(int blkEqn, int size)
{
  if (ptEqualBlk_ == true) return(0);

  std::pair<int,int> newpair;
  std::map<int,std::pair<int,int> >::iterator
    b_iter = blkEqns_->find(blkEqn);
  if (b_iter == blkEqns_->end()) {
    return(-1);
  }

  newpair = (*b_iter).second;
  newpair.second = size;
  (*b_iter).second = newpair;

  if (maxSize_ < size) maxSize_ = size;

  return(0);
}

//----------------------------------------------------------------------------
int snl_fei::PointBlockMap::getBlkEqnSize(int blkEqn)
{
  if (ptEqualBlk_ == true) return(1);

  std::map<int,std::pair<int,int> >::iterator
    b_iter = blkEqns_->find(blkEqn);

  if (b_iter != blkEqns_->end()) {
    return((*b_iter).second.second);
  }

  return(-1);
}

//----------------------------------------------------------------------------
int snl_fei::PointBlockMap::eqnToBlkEqn(int eqn) const
{
  if (ptEqualBlk_ == true) return(eqn);

  int blkEqn = -1;
  std::map<int,int>::iterator p_iter = ptEqns_->find(eqn);
  if (p_iter != ptEqns_->end()) blkEqn = (*p_iter).second;

  return(blkEqn);
}

//----------------------------------------------------------------------------
int snl_fei::PointBlockMap::blkEqnToPtEqn(int blkEqn) const
{
  if (ptEqualBlk_ == true) return(blkEqn);

  
  std::map<int,std::pair<int,int> >::iterator
    b_iter = blkEqns_->find(blkEqn);
  if (b_iter == blkEqns_->end()) {
    return(-1);
  }

  return((*b_iter).second.first);
}

//----------------------------------------------------------------------------
int snl_fei::PointBlockMap::getBlkEqnInfo(int blkEqn, int& ptEqn, int& blkSize)
{
  if (ptEqualBlk_ == true) {
    ptEqn = blkEqn;
    blkSize = 1;
    return(0);
  }

  std::map<int,std::pair<int,int> >::iterator
    b_iter = blkEqns_->find(blkEqn);
  if (b_iter == blkEqns_->end()) {
    return(-1);
  }

  ptEqn = (*b_iter).second.first;
  blkSize = (*b_iter).second.second;

  return(0);
}

//----------------------------------------------------------------------------
int snl_fei::PointBlockMap::getPtEqnInfo(int ptEqn,
					 int& blkEqn,
					 int& blkOffset)
{
  if (ptEqualBlk_ == true) {
    blkEqn = ptEqn;
    blkOffset = 0;
    return(0);
  }

  std::map<int,int>::iterator
    p_iter = ptEqns_->find(ptEqn);
  if (p_iter == ptEqns_->end()) {
    return(-1);
  }

  blkEqn = (*p_iter).second;

  std::map<int,std::pair<int,int> >::iterator
    b_iter = blkEqns_->find(blkEqn);

  std::pair<int,int> bpair = (*b_iter).second;

  blkOffset = ptEqn - bpair.first;

  return(0);
}

//----------------------------------------------------------------------------
int snl_fei::PointBlockMap::getBlkEqnOffset(int blkEqn, int eqn)
{
  if (ptEqualBlk_ == true) return(0);

  int blkOffset = 0;
  int err = getPtEqnInfo(eqn, blkEqn, blkOffset);
  if (err != 0) return(err);

  return(blkOffset);
}

//----------------------------------------------------------------------------
bool snl_fei::PointBlockMap::isExactlyBlkEqn(int ptEqn)
{
  if (ptEqualBlk_==true) return(true);

  std::map<int,int>::iterator
    p_iter = ptEqns_->find(ptEqn);
  if (p_iter == ptEqns_->end()) {
    return(false);
  }

  return( getBlkEqnOffset((*p_iter).first, ptEqn) == 0 );
}

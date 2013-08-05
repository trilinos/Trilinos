/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#include "mrtr_pnode.H"
#include "mrtr_interface.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 07/05|
 *----------------------------------------------------------------------*/
MOERTEL::ProjectedNode::ProjectedNode(const MOERTEL::Node& basenode, 
                                   const double* xi, 
                                   MOERTEL::Segment* pseg) :
MOERTEL::Node(basenode),
orthseg_(-1)
{
  pseg_ = pseg;
  if (xi)
  {
    xi_[0] = xi[0];
    xi_[1] = xi[1];
  }
  else
  {
    xi_[0] = 999.0;
    xi_[1] = 999.0;
  }
}

/*----------------------------------------------------------------------*
 |  ctor for orthogonal projection (public)                  mwgee 08/05|
 *----------------------------------------------------------------------*/
MOERTEL::ProjectedNode::ProjectedNode(const MOERTEL::Node& basenode, 
                                   const double* xi, 
                                   MOERTEL::Segment* pseg,
                                   int orthseg) :
MOERTEL::Node(basenode),
orthseg_(orthseg)
{
  pseg_ = pseg;
  if (xi)
  {
    xi_[0] = xi[0];
    xi_[1] = xi[1];
  }
  else
  {
    xi_[0] = 999.0;
    xi_[1] = 999.0;
  }
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 07/05|
 *----------------------------------------------------------------------*/
MOERTEL::ProjectedNode::ProjectedNode(MOERTEL::ProjectedNode& old) :
MOERTEL::Node(old)
{
  pseg_    = old.pseg_;
  xi_[0]   = old.xi_[0];
  xi_[1]   = old.xi_[1];
  orthseg_ = old.orthseg_;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 06/05|
 *----------------------------------------------------------------------*/
MOERTEL::ProjectedNode::~ProjectedNode()
{
  pseg_ = NULL; // this is just a 'referencing' ptr, not in charge of destroying
}

/*----------------------------------------------------------------------*
 |  print node                                               mwgee 06/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::ProjectedNode::Print() const
{ 
  std::cout << "Projected ";
  const MOERTEL::Node& basenode = dynamic_cast<const MOERTEL::Node&>(*this);
  std::cout << basenode;
  if (pseg_)
  {
    std::cout << "is on ";
    std::cout << *pseg_;
    std::cout << "at xi[0]/[1] = " << xi_[0] << "/" << xi_[1];
  }
  else
  {
    std::cout << "on Segment !!!!!NULL!!!!! at xi[0]/[1] = " << xi_[0] << "/" << xi_[1];
  }
  std::cout << "orth to seg " << orthseg_ << std::endl;
  return true;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 06/05|
 *----------------------------------------------------------------------*/
std::ostream& operator << (std::ostream& os, const MOERTEL::ProjectedNode& pnode)
{ 
  pnode.Print();
  return (os);
}

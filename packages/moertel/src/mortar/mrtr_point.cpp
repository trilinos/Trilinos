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
#include "mrtr_point.H"
#include "mrtr_overlap.hpp"
#include "mrtr_utils.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/05|
 *----------------------------------------------------------------------*/
MOERTEL::Point::Point(const int id, const double* xi, int out) :
id_(id),
outputlevel_(out)
{
  xi_[0] = xi[0];
  xi_[1] = xi[1];
  node_ = Teuchos::null;
  vals_[0].clear();
  vals_[1].clear();
  vals_[2].clear();
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 10/05|
 *----------------------------------------------------------------------*/
MOERTEL::Point::~Point()
{
  vals_[0].clear();
  vals_[1].clear();
  vals_[2].clear();
}
/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/05|
 *----------------------------------------------------------------------*/
std::ostream& operator << (std::ostream& os, const MOERTEL::Point& point)
{ 
  point.Print();
  return (os);
}
/*----------------------------------------------------------------------*
 |  print (public)                                           mwgee 10/05|
 *----------------------------------------------------------------------*/
void MOERTEL::Point::Print() const
{
  std::cout << "Point " << id_ << " xi[0]/[1] = " << xi_[0] << " / " << xi_[1] << std::endl;
  if (node_ != Teuchos::null)
    std::cout << *node_;
  return;
}

/*----------------------------------------------------------------------*
 |  store shape function values (public)                     mwgee 10/05|
 *----------------------------------------------------------------------*/
void MOERTEL::Point::StoreFunctionValues(int place, double* val, int valdim)
{
  vals_[place].resize(valdim);
  for (int i=0; i<valdim; ++i)
    vals_[place][i] = val[i];
  return;
}

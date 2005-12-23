/*
#@HEADER
# ************************************************************************
#
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#include "mrtr_point.H"
#include "mrtr_overlap.H"
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
  node_ = null;
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
ostream& operator << (ostream& os, const MOERTEL::Point& point)
{ 
  point.Print();
  return (os);
}
/*----------------------------------------------------------------------*
 |  print (public)                                           mwgee 10/05|
 *----------------------------------------------------------------------*/
void MOERTEL::Point::Print() const
{
  cout << "Point " << id_ << " xi[0]/[1] = " << xi_[0] << " / " << xi_[1] << endl;
  if (node_ != null)
    cout << *node_;
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

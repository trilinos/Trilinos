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
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Glen Hansen (Glen.Hansen@inl.gov)
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
  if (node_ != Teuchos::null)
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

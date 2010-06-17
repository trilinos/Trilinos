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
  cout << "Projected ";
  const MOERTEL::Node& basenode = dynamic_cast<const MOERTEL::Node&>(*this);
  cout << basenode;
  if (pseg_)
  {
    cout << "is on ";
    cout << *pseg_;
    cout << "at xi[0]/[1] = " << xi_[0] << "/" << xi_[1];
  }
  else
  {
    cout << "on Segment !!!!!NULL!!!!! at xi[0]/[1] = " << xi_[0] << "/" << xi_[1];
  }
  cout << "orth to seg " << orthseg_ << endl;
  return true;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 06/05|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const MOERTEL::ProjectedNode& pnode)
{ 
  pnode.Print();
  return (os);
}

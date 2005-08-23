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
#ifdef TRILINOS_PACKAGE

#include "mrtr_utils.H"
#include "mrtr_segment.H"
#include "mrtr_segment_linear1D.H"
#include "mrtr_node.H"

/*----------------------------------------------------------------------*
 | allocate a segment depending on the type                 mwgee 07/05|
 *----------------------------------------------------------------------*/
MRTR::Segment* MRTR::AllocateSegment(int type)
{
  switch (type)
  {
    case MRTR::Segment::seg_Linear1D:
      {
        MRTR::Segment_Linear1D* tmp = new MRTR::Segment_Linear1D();
        return tmp;
      }
    break;
    case MRTR::Segment::seg_none:
      cout << "***ERR*** MRTR::AllocateSegment:\n"
           << "***ERR*** type is func_none, cannot allocate.\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    break;
    default:
      cout << "***ERR*** MRTR::AllocateSegment:\n"
           << "***ERR*** type is unknown, cannot allocate new Segment\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    break;
  }

  return NULL;
}


/*----------------------------------------------------------------------*
 | allocate a function depending on the type                 mwgee 07/05|
 *----------------------------------------------------------------------*/
MRTR::Function* MRTR::AllocateFunction(MRTR::Function::FunctionType type)
{
  switch (type)
  {
    case MRTR::Function::func_Constant1D:
      {
        MRTR::Function_Constant1D* tmp = new MRTR::Function_Constant1D();
        return tmp;
      }
    break;
    case MRTR::Function::func_Linear1D:
      {
        MRTR::Function_Linear1D* tmp = new MRTR::Function_Linear1D();
        return tmp;
      }
    break;
    case MRTR::Function::func_DualLinear1D:
      {
        MRTR::Function_DualLinear1D* tmp = new MRTR::Function_DualLinear1D();
        return tmp;
      }
    break;
    case MRTR::Function::func_none:
      cout << "***ERR*** MRTR::AllocateFunction:\n"
           << "***ERR*** type is func_none, cannot allocate.\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    break;
    default:
      cout << "***ERR*** MRTR::AllocateFunction:\n"
           << "***ERR*** type is unknown, cannot allocate new Function\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    break;
  }

  return NULL;
}


/*----------------------------------------------------------------------*
 | destroy a segment map                                     mwgee 07/05|
 *----------------------------------------------------------------------*/
bool MRTR::DestroyMap(map<int,MRTR::Segment*>& m)
{
  map<int,MRTR::Segment*>::iterator curr;
  for (curr=m.begin(); curr != m.end(); ++curr)
    if (curr->second) delete curr->second;
  m.clear();  
  return true;
}

/*----------------------------------------------------------------------*
 | destroy a node map                                        mwgee 07/05|
 *----------------------------------------------------------------------*/
bool MRTR::DestroyMap(map<int,MRTR::Node*>& m)
{
  map<int,MRTR::Node*>::iterator curr;
  for (curr=m.begin(); curr != m.end(); ++curr)
    if (curr->second) delete curr->second;
  m.clear();  
  return true;
}

#endif // TRILINOS_PACKAGE

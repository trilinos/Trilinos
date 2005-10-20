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
#include "mrtr_segment_bilineartri.H"
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
    case MRTR::Segment::seg_BiLinearTri:
      {
        MRTR::Segment_BiLinearTri* tmp = new MRTR::Segment_BiLinearTri();
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
    case MRTR::Function::func_LinearTri:
      {
        MRTR::Function_LinearTri* tmp = new MRTR::Function_LinearTri();
        return tmp;
      }
    break;
    case MRTR::Function::func_DualLinearTri:
      {
        MRTR::Function_DualLinearTri* tmp = new MRTR::Function_DualLinearTri();
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

/*----------------------------------------------------------------------*
 | do cross product                                          mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::cross(double* out, const double* g1, const double* g2)
{
  out[0] = g1[1]*g2[2] - g1[2]*g2[1];
  out[1] = g1[2]*g2[0] - g1[0]*g2[2];
  out[2] = g1[0]*g2[1] - g1[1]*g2[0];
  return true;
}

/*----------------------------------------------------------------------*
 | do dot product                                            mwgee 10/05|
 *----------------------------------------------------------------------*/
double MRTR::dot(const double* g1, const double* g2, const int dim)
{
  double result=0.0;
  for (int i=0; i<dim; ++i) result+=g1[i]*g2[i];
  return result;
}

/*----------------------------------------------------------------------*
 | do 2x2 solve                                              mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::solve22(const double A[][2], double* x, const double* b)
{
  double det = A[0][0]*A[1][1]-A[0][1]*A[1][0];
  if (abs(det)<1.0e-10)
  {
    cout << "***ERR*** MRTR::solve22:\n"
         << "***ERR*** Determinant is zero\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  det = 1/det;
  x[0] = det*A[1][1]*b[0]-det*A[0][1]*b[1];
  x[1] = det*A[0][0]*b[1]-det*A[1][0]*b[0];
  return true;
}

/*----------------------------------------------------------------------*
 | do 3x3 solve                                              mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::solve33(const double A[][3], double* x, const double* b)
{
  double B[3][3];
  
  B[0][0] =   A[0][0]*A[2][2] - A[2][1]*A[1][2];
  B[1][0] = - A[1][0]*A[2][2] + A[2][0]*A[1][2];
  B[2][0] =   A[1][0]*A[2][1] - A[2][0]*A[1][1];
  B[0][1] = - A[0][1]*A[2][2] + A[2][1]*A[0][2];
  B[1][1] =   A[0][0]*A[2][2] - A[2][0]*A[0][2];
  B[2][1] = - A[0][0]*A[2][1] + A[2][0]*A[0][1];
  B[0][2] =   A[0][1]*A[1][2] - A[1][1]*A[0][2];
  B[1][2] = - A[0][0]*A[1][2] + A[1][0]*A[0][2];
  B[2][2] =   A[0][0]*A[1][1] - A[1][0]*A[0][1];
  
  const double detinv = 1./(A[0][0]*B[0][0]+A[0][1]*B[1][0]+A[0][2]*B[2][0]);

  for (int i=0; i<3; ++i)
  {
    x[i] = 0.;
    for (int j=0; j<3; ++j)
      x[i] += B[i][j]*b[j]*detinv;
  }

  return true;
}

/*----------------------------------------------------------------------*
 | get the '10' digit from a pos. int                        mwgee 10/05|
 *----------------------------------------------------------------------*/
int MRTR::digit_ten(int number)
{
  number = abs(number);
  if (number<10) return 0;
  number /= 10; 
  number = number%10;
  return number;
}

/*----------------------------------------------------------------------*
 | swap 2 kinds                                              mwgee 10/05|
 | this template is given in mrtr_utils.H                               |
 *----------------------------------------------------------------------*/
// template<typename kind> void swap(kind& a, kind& b);


#endif // TRILINOS_PACKAGE

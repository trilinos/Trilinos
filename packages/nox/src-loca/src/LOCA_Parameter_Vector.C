// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "LOCA_Parameter_Vector.H" // class definition

#include "Teuchos_TestForException.hpp" // for errors

LOCA::ParameterVector::ParameterVector() : 
  x(),
  l()
{
}

LOCA::ParameterVector::ParameterVector(const LOCA::ParameterVector& source) :
  x(source.x),
  l(source.l)
{
}

LOCA::ParameterVector* 
LOCA::ParameterVector::clone() const
{
  LOCA::ParameterVector* y = new LOCA::ParameterVector(*this);
  return y;
}

LOCA::ParameterVector::~ParameterVector()
{
}

int 
LOCA::ParameterVector::addParameter(string label, double value)
{
  unsigned int size = x.size();
  x.push_back(value);
  l.push_back(label);
  return (size);
}

bool 
LOCA::ParameterVector::init(double value)
{
  for (unsigned int i = 0; i < x.size(); i ++)
    x[i] = value;
  return true;
}

bool 
LOCA::ParameterVector::scale(double value)
{
  for (unsigned int i = 0; i < x.size(); i ++)
    x[i] *= value;
  return true;
}

bool 
LOCA::ParameterVector::scale(const LOCA::ParameterVector& p)
{
  //make sure vectors are of compatible size
  if (this->x.size() != p.x.size())
    return false;

  for (unsigned int i = 0; i < x.size(); i ++)
    x[i] *= p[i];
  return true;
}

bool 
LOCA::ParameterVector::update(double alpha, 
			      const LOCA::ParameterVector& alphaVector, 
			      double b)
{  //make sure vectors are of compatible size
  if (x.size() != alphaVector.x.size())
    return false;

  for (unsigned int i = 0; i < x.size(); i ++) {
    x[i] *= b;
    x[i] += alpha*alphaVector[i];
  }
  return true;
}

LOCA::ParameterVector& 
LOCA::ParameterVector::operator=(const LOCA::ParameterVector& source)
{
  x = source.x;
  l = source.l;
  return *this;
}

double& 
LOCA::ParameterVector::operator[] (unsigned int i)
{
  TEST_FOR_EXCEPTION(i >= x.size(), 
		     std::out_of_range,
		     "Error:  LOCA::ParameterVector::operator[]:  " << 
		     " Index " << i << " is out of range!");
  return x[i];
}

const double& 
LOCA::ParameterVector::operator[] (unsigned int i) const
{
  TEST_FOR_EXCEPTION(i >= x.size(), 
		     std::out_of_range,
		     "Error:  LOCA::ParameterVector::operator[]:  " << 
		     " Index " << i << " is out of range!");
  return x[i];
}

void 
LOCA::ParameterVector::setValue(unsigned int i, double value)
{
  TEST_FOR_EXCEPTION(i >= x.size(), 
		     std::out_of_range,
		     "Error:  LOCA::ParameterVector::setValue():  " << 
		     " Index " << i << " is out of range!");

  x[i] = value;
  return;
}

void 
LOCA::ParameterVector::setValue(string label, double value)
{
  for (unsigned int i = 0; i < x.size(); i++) {
    if (l[i] == label) {
      x[i] = value;
      return;
    }
  }

  TEST_FOR_EXCEPTION(true, 
		     std::invalid_argument,
		     "Error:  LOCA::ParameterVector::setValue():  " << 
		     " Label " << label << " is not valid!");
}

double 
LOCA::ParameterVector::getValue(unsigned int i) const
{
  TEST_FOR_EXCEPTION(i >= x.size(), 
		     std::out_of_range,
		     "Error:  LOCA::ParameterVector::getValue():  " << 
		     " Index " << i << " is out of range!");
  return x[i];
}

double 
LOCA::ParameterVector::getValue(string label) const
{
  for (unsigned int i = 0; i < x.size(); i++) {
    if (l[i] == label)
      return x[i];
  }

  TEST_FOR_EXCEPTION(true, 
		     std::invalid_argument,
		     "Error:  LOCA::ParameterVector::getValue():  " << 
		     " Label " << label << " is not valid!");
  return 0.0;
}

int 
LOCA::ParameterVector::getIndex(string label) const
{
  for (unsigned int i = 0; i < x.size(); i++) {
    if (l[i] == label)
      return i;
  }

  TEST_FOR_EXCEPTION(true, 
		     std::invalid_argument,
		     "Error:  LOCA::ParameterVector::getIndex():  " << 
		     " Label " << label << " is not valid!");
  return -1;
}

double* 
LOCA::ParameterVector::getDoubleArrayPointer()
{
  if (x.size() == 0)
    return NULL;
  return &x[0];
}

bool 
LOCA::ParameterVector::isParameter(string label) const
{
  for (unsigned int i = 0; i < x.size(); i++) {
    if (l[i] == label)
      return true;
  }
  return false;
}

string 
LOCA::ParameterVector::getLabel(unsigned int i) const
{
  return l[i];
}

int 
LOCA::ParameterVector::length() const
{
  return x.size();
}

void 
LOCA::ParameterVector::print(ostream& stream) const
{
  stream << "LOCA::ParameterVector \n(size = " << x.size() << ")";
  for (unsigned int i = 0; i < x.size(); i++) {
    stream << "\n    " << i << "    " << l[i] << " = " << x[i];
  }
  stream << std::endl;
}

ostream& 
operator<<(ostream& stream, const LOCA::ParameterVector& p)
{
  p.print(stream);
  return stream;
}

const vector<double>& 
LOCA::ParameterVector::getValuesVector() const
{
  return x;
}

const vector<string>& 
LOCA::ParameterVector::getNamesVector() const
{
  return l;
}

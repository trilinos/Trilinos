// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "LOCA_Parameter_Vector.H" // class definition

#include "LOCA_Utils.H"            // print utilities

using namespace LOCA;

ParameterVector::ParameterVector()
{
}

ParameterVector::ParameterVector(const ParameterVector& source) :
  x(source.x),
  l(source.l)
{
}

ParameterVector* ParameterVector::clone() const
{
  ParameterVector* y = new ParameterVector(*this);
  return y;
}

ParameterVector::~ParameterVector()
{
}

int ParameterVector::addParameter(string label, double value)
{
  unsigned int size = x.size();
  x.push_back(value);
  l.push_back(label);
  return (size);
}

bool ParameterVector::init(double value)
{
  for (unsigned int i = 0; i < x.size(); i ++)
    x[i] = value;
  return true;
}

bool ParameterVector::scale(double value)
{
  for (unsigned int i = 0; i < x.size(); i ++)
    x[i] *= value;
  return true;
}

bool ParameterVector::scale(const ParameterVector& p)
{
  //make sure vectors are of compatible size
  if (this->x.size() != p.x.size())
    return false;

  for (unsigned int i = 0; i < x.size(); i ++)
    x[i] *= p[i];
  return true;
}

bool ParameterVector::update(double alpha, 
			     const ParameterVector& alphaVector, 
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

ParameterVector& ParameterVector::operator=(const ParameterVector& source)
{
  x = source.x;
  l = source.l;
  return *this;
}

double& ParameterVector::operator[] (int i)
{
  if ((i < 0) || (i >= x.size())) {
    if (Utils::doPrint(Utils::Error)) {
      cout << "ERROR: LOCA::Parameter::Vector::operator[] - index is out "
	   << "of range!" << endl;
    }
    throw "NOX Error";
  }
  return x[i];
}

const double& ParameterVector::operator[] (int i) const
{
  if ((i < 0) || (i >= x.size())) {
    if (Utils::doPrint(Utils::Error)) {
      cout << "ERROR: LOCA::Parameter::Vector::operator[] const - index is "
	   << "out of range!" << endl;
    }
    throw "NOX Error";
  }
  return x[i];
}

void ParameterVector::setValue(int i, double value)
{
  if ((i < 0) || (i >= x.size())) {
    if (Utils::doPrint(Utils::Error)) {
      cout << "ERROR: LOCA::Parameter::Vector::setValue() - index is "
	   << "out of range!" << endl;
    }
    throw "NOX Error";
  }

  x[i] = value;
  return;
}

void ParameterVector::setValue(string label, double value)
{
  if (!isParameter(label)) {
    if (Utils::doPrint(Utils::Error)) {
      cout << "ERROR: LOCA::Parameter::Vector::setValue() - label is "
	   << "not valid!" << endl;
    }
    throw "NOX Error";
  }
  
  for (unsigned int i = 0; i < x.size(); i++) {
    if (l[i] == label)
      x[i] = value;
  }
  
  return;
}

double ParameterVector::getValue(int i) const
{
  if ((i < 0) || (i >= x.size())) {
    if (Utils::doPrint(Utils::Error)) {
      cout << "ERROR: LOCA::Parameter::Vector::getValue(int) - index is "
	   << "out of range!" << endl;
    }
    throw "NOX Error";
  }
  return x[i];
}

double ParameterVector::getValue(string label) const
{
  for (unsigned int i = 0; i < x.size(); i++) {
    if (l[i] == label)
      return x[i];
  }

  if (Utils::doPrint(Utils::Error)) {
    cout << "ERROR: LOCA::Parameter::Vector::getValue(string) - label is "
	 << "not valid!" << endl;
  }
  throw "NOX Error";

  return 0.0;
}

int ParameterVector::getIndex(string label) const
{
  for (unsigned int i = 0; i < x.size(); i++) {
    if (l[i] == label)
      return i;
  }

  if (Utils::doPrint(Utils::Warning)) {
    cout << "Warning: LOCA::ParameterVector::getIndex() - The string \"" 
	 << label << "\" was not found!" << endl;
  }

  return -1;
}

double* ParameterVector::getDoubleArrayPointer()
{
  return &x[0];
}

bool ParameterVector::isParameter(string label) const
{
  for (unsigned int i = 0; i < x.size(); i++) {
    if (l[i] == label)
      return true;
  }
  return false;
}

string ParameterVector::getLabel(int i) const
{
  return l[i];
}

int ParameterVector::length() const
{
  return x.size();
}

void ParameterVector::print(ostream& stream) const
{
  stream << "LOCA::ParameterVector \n(size = " << x.size() << ")";
  for (unsigned int i = 0; i < x.size(); i++) {
    stream << "\n    " << i << "    " << l[i] << " = " << x[i];
  }
  cout << endl;
}

ostream& operator<<(ostream& stream, const ParameterVector& p)
{
  p.print(stream);
  return stream;
}

// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#include "Teuchos_Assert.hpp" // for errors

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
LOCA::ParameterVector::addParameter(std::string label, double value)
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
  TEUCHOS_TEST_FOR_EXCEPTION(i >= x.size(), 
		     std::out_of_range,
		     "Error:  LOCA::ParameterVector::operator[]:  " << 
		     " Index " << i << " is out of range!");
  return x[i];
}

const double& 
LOCA::ParameterVector::operator[] (unsigned int i) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(i >= x.size(), 
		     std::out_of_range,
		     "Error:  LOCA::ParameterVector::operator[]:  " << 
		     " Index " << i << " is out of range!");
  return x[i];
}

void 
LOCA::ParameterVector::setValue(unsigned int i, double value)
{
  TEUCHOS_TEST_FOR_EXCEPTION(i >= x.size(), 
		     std::out_of_range,
		     "Error:  LOCA::ParameterVector::setValue():  " << 
		     " Index " << i << " is out of range!");

  x[i] = value;
  return;
}

void 
LOCA::ParameterVector::setValue(std::string label, double value)
{
  for (unsigned int i = 0; i < x.size(); i++) {
    if (l[i] == label) {
      x[i] = value;
      return;
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION(true, 
		     std::invalid_argument,
		     "Error:  LOCA::ParameterVector::setValue():  " << 
		     " Label " << label << " is not valid!");
}

double 
LOCA::ParameterVector::getValue(unsigned int i) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(i >= x.size(), 
		     std::out_of_range,
		     "Error:  LOCA::ParameterVector::getValue():  " << 
		     " Index " << i << " is out of range!");
  return x[i];
}

double 
LOCA::ParameterVector::getValue(std::string label) const
{
  for (unsigned int i = 0; i < x.size(); i++) {
    if (l[i] == label)
      return x[i];
  }

  TEUCHOS_TEST_FOR_EXCEPTION(true, 
		     std::invalid_argument,
		     "Error:  LOCA::ParameterVector::getValue():  " << 
		     " Label " << label << " is not valid!");
  return 0.0;
}

int 
LOCA::ParameterVector::getIndex(std::string label) const
{
  for (unsigned int i = 0; i < x.size(); i++) {
    if (l[i] == label)
      return i;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(true, 
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
LOCA::ParameterVector::isParameter(std::string label) const
{
  for (unsigned int i = 0; i < x.size(); i++) {
    if (l[i] == label)
      return true;
  }
  return false;
}

std::string 
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
LOCA::ParameterVector::print(std::ostream& stream) const
{
  stream << "LOCA::ParameterVector \n(size = " << x.size() << ")";
  for (unsigned int i = 0; i < x.size(); i++) {
    stream << "\n    " << i << "    " << l[i] << " = " << x[i];
  }
  stream << std::endl;
}

std::ostream& 
LOCA::operator<<(std::ostream& stream, const LOCA::ParameterVector& p)
{
  p.print(stream);
  return stream;
}

const std::vector<double>& 
LOCA::ParameterVector::getValuesVector() const
{
  return x;
}

const std::vector<std::string>& 
LOCA::ParameterVector::getNamesVector() const
{
  return l;
}

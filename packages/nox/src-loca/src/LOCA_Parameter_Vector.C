// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
#ifdef __NVCC__
#pragma diag_suppress code_is_unreachable
#endif
  return 0.0;
#ifdef __NVCC__
#pragma diag_warning code_is_unreachable
#endif
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
#ifdef __NVCC__
#pragma diag_suppress code_is_unreachable
#endif
  return -1;
#ifdef __NVCC__
#pragma diag_warning code_is_unreachable
#endif
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

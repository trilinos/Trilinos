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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "NOX_Common.H"
#include "NOX_TSF_Vector.H"
#include "NOX_Utils.H"
#include "TSFCoreVectorStdOps.hpp"
#include "NOX_Random.H" // for Random class


NOX::TSF::Vector::Vector(const NOX::TSF::Vector& source, 
			    NOX::CopyType type)
{
 switch (type) 
 {
    
  case NOX::DeepCopy:
    x = source.x.copy();
    break;

  case NOX::ShapeCopy:
    x = ((source.x).space()).createMember();
    break;

  default:
    cerr << "NOX:TSF::Vector - invalid CopyType for copy constructor." << endl;
    throw "NOX TSF Error";
  }
}

NOX::TSF::Vector::Vector(const TSFExtended::Vector<double>& source, 
			    NOX::CopyType type)
{
  switch (type) 
 {
    
  case NOX::DeepCopy:
    x = source.copy();
    break;

  case NOX::ShapeCopy:
    x = ((source).space()).createMember();
    break;

  default:
    cerr << "NOX:TSF::Vector - invalid CopyType for copy constructor." << endl;
    throw "NOX TSF Error";
  }
}


TSFExtended::Vector<double>& NOX::TSF::Vector::getTSFVector()
{
  return x;
}
 
const TSFExtended::Vector<double>& NOX::TSF::Vector::getTSFVector() const
{
  return x;
}


NOX::Abstract::Vector& NOX::TSF::Vector::operator=(
					   const NOX::Abstract::Vector& source)
{
  return operator=(dynamic_cast<const NOX::TSF::Vector&>(source));
}

NOX::Abstract::Vector& NOX::TSF::Vector::operator=(
					   const NOX::TSF::Vector& source)
{
  // in TSFExtended operator= results in a shallow copy while 
  // acceptCopyOf(source.x) provides the deep copy we want
  x = source.getTSFVector().copy();
  return *this;
}

// void NOX::TSF::Vector::setElement(int i, const double& value)
// {
//   return x.setElement(i,value);
// }

// const double& NOX::TSF::Vector::getElement(int i) const
// {
//   return x.getElement(i);
// }   

//const double& NOX::TSF::Vector::operator() (int i) const
//{
//  return x.getElement(i);
//}  

NOX::Abstract::Vector& NOX::TSF::Vector::init(double value)
{
  x.setToConstant(value);
  return *this;
}


NOX::Abstract::Vector& NOX::TSF::Vector::abs(
					     const NOX::Abstract::Vector& base)
{
  return abs(dynamic_cast<const NOX::TSF::Vector&>(base));
}

NOX::Abstract::Vector& NOX::TSF::Vector::abs(
					     const NOX::TSF::Vector& base)
{
  x.acceptCopyOf(base.x);
  x.abs();
  return *this;
}

NOX::Abstract::Vector& NOX::TSF::Vector::reciprocal(
					    const NOX::Abstract::Vector& base)
{
  return reciprocal(dynamic_cast<const NOX::TSF::Vector&>(base));
}

NOX::Abstract::Vector& NOX::TSF::Vector::reciprocal(
					    const NOX::TSF::Vector& base)
{
  x.acceptCopyOf(base.x);
  x.reciprocal();
  return *this;
}

NOX::Abstract::Vector& NOX::TSF::Vector::scale(double alpha)
{
  x.scale(alpha);
  return *this;
}

NOX::Abstract::Vector& NOX::TSF::Vector::update(
					       double alpha, 
					       const NOX::Abstract::Vector& a, 
					       double gamma)
{
  return update( alpha, dynamic_cast<const NOX::TSF::Vector&>(a), gamma);
}

NOX::Abstract::Vector& NOX::TSF::Vector::update(
						 double alpha, 
						 const NOX::TSF::Vector& a, 
						 double gamma)
{
  x.update(alpha,a.x,gamma);
  return *this;
}

NOX::Abstract::Vector& NOX::TSF::Vector::update(
					      double alpha, 
					      const NOX::Abstract::Vector& a, 
					      double beta, 
					      const NOX::Abstract::Vector& b,
					      double gamma)
{
  return update(alpha, dynamic_cast<const NOX::TSF::Vector&>(a), 
		beta, dynamic_cast<const NOX::TSF::Vector&>(b), gamma);
}

NOX::Abstract::Vector& NOX::TSF::Vector::update(
					       double alpha, 
					       const NOX::TSF::Vector& a, 
					       double beta, 
					       const NOX::TSF::Vector& b,
					       double gamma)
{
  x.update(alpha,a.x,beta,b.x,gamma);
  return *this;
}

NOX::Abstract::Vector& NOX::TSF::Vector::scale(
					      const NOX::Abstract::Vector& a)
{  
  return scale(dynamic_cast<const NOX::TSF::Vector&>(a));
}

NOX::Abstract::Vector& NOX::TSF::Vector::scale(const NOX::TSF::Vector& a)
{  
  x.dotStar(a.x);
  return *this;
}

NOX::Abstract::Vector* NOX::TSF::Vector::clone(NOX::CopyType type) const
{
  return new NOX::TSF::Vector(*this, type);
}

double NOX::TSF::Vector::norm(NOX::Abstract::Vector::NormType type) const
{
  
  if (this->length() == 0)
    return 0.0;

  int i;			// counter
  double value;			// final answer

  switch (type) 
  {
  case MaxNorm:
    value = x.normInf();
    break;
  case OneNorm:
    value = x.norm1();
    break;
  case TwoNorm:
  default:
    value = x.norm2();
   break;
  }

  return value;
}

double NOX::TSF::Vector::norm(const NOX::Abstract::Vector& weights) const
{
  return norm(dynamic_cast<const NOX::TSF::Vector&>(weights));
}

double NOX::TSF::Vector::norm(const NOX::TSF::Vector& weights) const
{
  if (weights.length() != this->length()) 
  {
    cerr << "NOX::TSF::Vector::norm - size mismatch for weights vector" << endl;
    throw "NOX::TSF Error";
  }
  return x.norm2(weights.getTSFVector());
  // change here when TSFExtended gets a weighted 2 norm
}

double NOX::TSF::Vector::dot(const NOX::Abstract::Vector& y) const
{
  return dot(dynamic_cast<const NOX::TSF::Vector&>(y));
}

double NOX::TSF::Vector::dot(const NOX::TSF::Vector& y) const
{
  if (y.length() != this->length()) 
  {
    cerr << "NOX::TSF::Vector::dot - size mismatch for y vector" << endl;
    throw "NOX::TSF Error";
  }

  return x.dot(y.x);
}

int NOX::TSF::Vector::length() const
{
  return (x.space()).dim();
}

ostream& NOX::TSF::Vector::leftshift(ostream& stream) const
{
  stream << "[ ";
  for (int i = 0; i < this->length(); i ++) 
    stream << NOX::Utils::sciformat(x.getElement(i),8)<< " ";
  stream << "]";
  return stream;
}

ostream& operator<<(ostream& stream, const NOX::TSF::Vector& v)
{
  return v.leftshift(stream);
}

void NOX::TSF::Vector::print() const
{
  cout << *this << endl;
}

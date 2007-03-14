// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// @HEADER


#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Handle.hpp"
#include "Teuchos_Handleable.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_MPIContainerComm.hpp"
#include "Teuchos_ErrorPolling.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_Version.hpp"


using namespace Teuchos;
using std::string;

//#ifndef DOXYGEN_SHOULD_SKIP_THIS

class VecSpaceBase;

class VecBase : public Handleable<VecBase>
{
public:
  VecBase(){;}
  ~VecBase(){;}

  virtual RefCountPtr<const VecSpaceBase> space() const = 0 ;

  virtual void add(const VecBase* other, 
                   RefCountPtr<VecBase>& result) const = 0 ;


  virtual double dot(const VecBase* other) const = 0 ;

  virtual void scale(const double& a) = 0 ;

  virtual RefCountPtr<VecBase> copy() const = 0 ;

  virtual void print(std::ostream& os) const = 0 ;

  virtual void setElement(int i, const double& x) = 0 ;

  virtual const double& getElement(int i) const = 0 ;

  virtual int dim() const = 0 ;
};


class VecSpaceBase : public ConstHandleable<VecSpaceBase>
{
public:
  VecSpaceBase(){;}
  ~VecSpaceBase(){;}

  virtual RefCountPtr<VecBase> create(const RefCountPtr<const VecSpaceBase>& s) const = 0;

};



class VecSpaceA : public VecSpaceBase
{
public:
  VecSpaceA(int n) : n_(n) {;}
  virtual RefCountPtr<VecBase> create(const RefCountPtr<const VecSpaceBase>& s) const ;
  TEUCHOS_GET_CONST_RCP(VecSpaceBase);
private:
  int n_;
};

class VecA : public VecBase
{
public:
  VecA(int n, const RefCountPtr<const VecSpaceBase>& sp) : x_(n),
                                                           sp_(sp)
  {for (unsigned int i=0; i<x_.size(); i++) x_[i] = 0.0;}

  RefCountPtr<const VecSpaceBase> space() const {return sp_;}
  
  void add(const VecBase* other, 
           RefCountPtr<VecBase>& result) const 
  {
    const VecA* va = dynamic_cast<const VecA*>(other);
    VecA* vr = dynamic_cast<VecA*>(result.get());
    TEST_FOR_EXCEPT(va==0);
    TEST_FOR_EXCEPT(vr==0);
    for (unsigned int i=0; i<x_.size(); i++) 
      {
        vr->x_[i] = x_[i] + va->x_[i];
      }
  }

  double dot(const VecBase* other) const 
  {
    double rtn = 0.0;
    const VecA* va = dynamic_cast<const VecA*>(other);
    TEST_FOR_EXCEPT(va==0);
    for (unsigned int i=0; i<x_.size(); i++) 
      {
        rtn += x_[i] * va->x_[i];
      }
    return rtn;
  }

  void scale(const double& a) 
  {
    for (unsigned int i=0; i<x_.size(); i++) 
      {
        x_[i] *= a;
      }
  }

  RefCountPtr<VecBase> copy() const
  {
    RefCountPtr<VecBase> rtn = space()->create(space());
    VecA* va = dynamic_cast<VecA*>(rtn.get());
    TEST_FOR_EXCEPT(va==0);
    for (unsigned int i=0; i<x_.size(); i++) 
      {
        va->x_[i] = x_[i];
      }
    return rtn;
  }

  void print(std::ostream& os) const 
  {
    for (unsigned int i=0; i<x_.size(); i++) 
      {
        os << i << " " << x_[i] << endl;
      }
  }

  void setElement(int i, const double& x) 
  {x_[i] = x;}

  const double& getElement(int i) const {return x_[i];}

  int dim() const {return x_.size();}

  TEUCHOS_GET_RCP(VecBase);
private:
  Array<double> x_;
  RefCountPtr<const VecSpaceBase> sp_;
};


RefCountPtr<VecBase> VecSpaceA::create(const RefCountPtr<const VecSpaceBase>& s) const
{
  return rcp(new VecA(n_, s));
}

class Vector;



class ConstVector : public virtual ConstHandle<VecBase>
{
public:
  TEUCHOS_CONST_HANDLE_CTORS(ConstVector, VecBase);
  

  RefCountPtr<const VecSpaceBase> space() const {return constPtr()->space();}

  double operator*(const ConstVector& other) const ;  

  void add(const ConstVector& other, Vector& result) const ;

  void copyInto(Vector& x) const ;

  int dim() const {return constPtr()->dim();}

  const double& getElement(int i) const {return constPtr()->getElement(i);}
};

class Vector : public ConstVector,
               public Handle<VecBase>
{
public:
  TEUCHOS_HANDLE_CTORS(Vector, VecBase);

  void scale(const double& a) {ptr()->scale(a);}

  void setElement(int i, const double& x) {ptr()->setElement(i,  x);}

};


Vector copy(const ConstVector& x)
{
  Vector rtn;
  x.copyInto(rtn);
  return rtn;
}

Vector operator+(const ConstVector& a, const ConstVector& b)
{
  Vector result;
  a.add(b, result);
  return result;
}

Vector operator*(const ConstVector& x, const double& a) 
{
  Vector result = copy(x);
  result.scale(a);
  return result;
}

Vector operator*(const double& a, const ConstVector& x) 
{
  return x*a;
}


void ConstVector::add(const ConstVector& other, Vector& result) const 
{
  result = space()->create(space());
  RefCountPtr<VecBase> tmp = result.ptr();
  constPtr()->add(other.constPtr().get(), tmp);
}

void ConstVector::copyInto(Vector& result) const 
{
  result = constPtr()->copy();
}




class VectorSpace : public ConstHandle<VecSpaceBase>
{
public:
  TEUCHOS_CONST_HANDLE_CTORS(VectorSpace, VecSpaceBase);
  
  Vector create() const {return constPtr()->create(constPtr());}
};


std::ostream& operator<<(std::ostream& os, const ConstVector& v)
{
  v.constPtr()->print(os);
  return os;
}

//#endif

/* Test of Teuchos generic handle classes */

int main(int argc, char** argv)
{
  int state = 0;
  cout << Teuchos::Teuchos_Version() << endl << endl;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  try
    {
      VectorSpace space = new VecSpaceA(3);
      Vector x = space.create();
      Vector y = space.create();

      
      for (int i=0; i<x.dim(); i++)
        {
          x.setElement(i, i);
          y.setElement(i, i+2.0);
        }
      Vector z = copy(x);

      cout << "x = " << x << endl;
      cout << "y = " << y << endl;
      cout << "z = " << z << endl;

      cout << "mess = " << 2.0*(x+y+3.0*z) << endl;

      Vector a = 2.0*(x+y+3.0*z);
      cout << "a=" << endl;
      double err = 0.0;
      for (int i=0; i<a.dim(); i++)
        {
          cout << i << " " << a.getElement(i) << endl;
          double x_i = x.getElement(i);
          double y_i = y.getElement(i);
          double z_i = z.getElement(i);
          double t = 2.0*(x_i + y_i + 3.0*z_i);
          err += ::fabs(t - a.getElement(i));
        }
      
      VectorSpace s2 = new VecSpaceA(5);
      VecBase* vb = new VecA(5, s2.constPtr());
      Vector b = vb;

      cout << "b = " << b << endl;

      if (err > 1.0e-12) state = 1;
    }
  catch(std::exception& e)
    {
      cerr << e.what() << endl;
      state = 1;
    }

  if (state != 0)
    {
      cout << "TEST FAILED" << endl;
      return -1;
    }

  
  cout << "TEST PASSED" << endl;
  return state;
}

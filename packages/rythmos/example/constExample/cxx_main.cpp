//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Teuchos_RefCountPtr.hpp"
#include<iostream>

using namespace Teuchos;

// Dummy class to contain a double 
// ------------------------------------------------------------
class Bar 
{
  public:
    // Constructor
    Bar();
    Bar(double x);
    // Destructor
    ~Bar();
    void setx(double x);
    double getx() const; 
  protected:
    double x_;
};
Bar::Bar()
{ 
  x_ = 5.0; 
}
Bar::Bar(double x)
{ 
  x_ = x; 
}
Bar::~Bar() 
{ 
  x_ = 0.0; 
}
void Bar::setx(double x)
{
  x_ = x;
}
double Bar::getx() const
{ 
  return(x_); 
}
// ------------------------------------------------------------

// Class to test out const
class Foo
{
  public:
    Foo();
    ~Foo();
    void setBar(const RefCountPtr<Bar> &Bptr);
    void setConstBar(const RefCountPtr<const Bar> &Bptr);
    const RefCountPtr<Bar> &getBar() const;
    void test1(const RefCountPtr<Bar> &Bptr);
    void test2();
    void test3(RefCountPtr<Bar> &Bptr);
    void test4(const RefCountPtr<Bar> &Bptr);
    RefCountPtr<Bar> &test5();
    RefCountPtr<const Bar> &test6();
    const RefCountPtr<Bar> &test7();
    const RefCountPtr<Bar> test8();
    RefCountPtr<Bar> test9();
    RefCountPtr<Bar> &test10();
    //RefCountPtr<Bar> &test11();
    const RefCountPtr<const Bar> &test12();
    const RefCountPtr<const Bar> &test13() const;
    const RefCountPtr<Bar> &test14() const;
    const RefCountPtr<Bar> test15() const;
    double test16(const Bar &B);

  protected:
    double t_;
    RefCountPtr<Bar> Bptr_;
    RefCountPtr<const Bar> BptrConst_;

};
Foo::Foo() { }
Foo::~Foo() { }
void Foo::setBar(const RefCountPtr<Bar> &Bptr) 
{ 
  Bptr_ = Bptr;  
}
void Foo::setConstBar(const RefCountPtr<const Bar> &Bptr) 
{ 
  BptrConst_ = Bptr; 
  //Bptr_ = Bptr;  // not allowed because Bptr_ is nonconst
  //Bptr->setx(12.0); // not allowed because Bptr is const Bar RefCountPtr
}
const RefCountPtr<Bar> &Foo::getBar() const
{
  return(Bptr_);
}
void Foo::test1(const RefCountPtr<Bar> &Bptr)
{
  Bptr->setx(15.0);
}
void Foo::test2()
{
//  BptrConst_->setx(20.0); // not allowed because BptrConst_ is const
}
void Foo::test3(RefCountPtr<Bar> &Bptr)
{
  RefCountPtr<Bar> Bptr_temp = rcp(new Bar);
  Bptr_temp->setx(25.0);
  Bptr = Bptr_temp;
}
void Foo::test4(const RefCountPtr<Bar> &Bptr)
{
  RefCountPtr<Bar> Bptr_temp = rcp(new Bar);
  Bptr_temp->setx(30.0);
//  Bptr = Bptr_temp; // not allowed because input is const
}
RefCountPtr<Bar> &Foo::test5()
{
  return(Bptr_);
}
RefCountPtr<const Bar> &Foo::test6()
{
  //return(Bptr_); // not allowed because Bptr is nonconst
  return(BptrConst_);
}
const RefCountPtr<Bar> &Foo::test7()
{
  return(Bptr_);
}
const RefCountPtr<Bar> Foo::test8()
{
  return(rcp(new Bar(45.0)));
}
RefCountPtr<Bar> Foo::test9()
{
//  return(rcp(new Bar(55.0))); // allowed also
  return(Bptr_);
}
RefCountPtr<Bar> &Foo::test10()
{
  return(Bptr_);
}
/*
RefCountPtr<Bar> &Foo::test11()
{
//  return(rcp(new Bar(75.0))); // not allowed because I'm passing out a
  // reference to an internal object
}
*/
const RefCountPtr<const Bar> &Foo::test12()
{
  return(BptrConst_);
}
const RefCountPtr<const Bar> &Foo::test13() const
{
  return(BptrConst_);
}
const RefCountPtr<Bar> &Foo::test14() const
{
  return(Bptr_);
}
const RefCountPtr<Bar> Foo::test15() const
{
  return(Bptr_);
}
double Foo::test16(const Bar &B)
{
  return(B.getx());
}
// ------------------------------------------------------------
void somefunction(const Foo &F)
{
  RefCountPtr<Bar> Bptr = F.getBar(); // only allowed if getBar has const after it in definition
}
// ------------------------------------------------------------


int main(int argc, char *argv[])
{
  cout << "main:  This routine tests const conditions." << endl;

  Foo F;
  RefCountPtr<Bar> Bptr = rcp(new Bar);
  Bptr->setx(10.0);
  F.setBar(Bptr);
  cout << "Correct output is x = 10." << endl;
  cout << "x = " << Bptr->getx() << endl;

  F.test1(Bptr);
  cout << "Correct output is x = 15." << endl;
  cout << "x = " << Bptr->getx() << endl;

  // This demonstrates that Bptr was overwritten in this member fcn call with a
  // different RefCountPtr object, which is only possible because the argument
  // list definition was nonconst.
  F.test3(Bptr);
  cout << "Correct output is x = 25." << endl;
  cout << "x = " << Bptr->getx() << endl;

  Bptr = F.test5();
  cout << "Correct output is x = 15." << endl;
  cout << "x = " << Bptr->getx() << endl;

  Bptr->setx(35.0);
//  Bptr = F.test6(); // not allowed because Bptr is nonconst and output of test6 is const
  RefCountPtr<const Bar> BptrConst = rcp(new Bar(40.0));
  F.setConstBar(BptrConst);
  cout << "Correct output is x = 40." << endl;
  cout << "x = " << F.test6()->getx() << endl;  // valid because we put const after Bar::getx

  Bptr = F.test7();
  cout << "Correct output is x = 35." << endl;
  cout << "x = " << Bptr->getx() << endl;  

  Bptr = F.test8();
  cout << "Correct output is x = 45." << endl;
  cout << "x = " << Bptr->getx() << endl;  
  //F.test8() = rcp(new Bar(50.0)); // not allowed because output is const.

  F.test9() = rcp(new Bar(60.0)); // allowed because output is nonconst, but
  // this is crazy.  Especially since the object that is modified is a copied
  // object, so we don't have access to it after this call.

  F.test10() = rcp(new Bar(65.0));  // allowed because output is nonconst, and
  //this does something strange, as it modifies the internal data to the class
  //  through a nonconst reference passed out.
  Bptr = F.test7(); // grab it back out with const output
  cout << "Correct output is x = 65." << endl;
  cout << "x = " << Bptr->getx() << endl;  

//  F.test7() = rcp(new Bar(70.0));  // not allowed because output is const, this
  // is the expected behavior.  You shouldn't be able to do assignments like this.

  BptrConst = F.test12();
  cout << "Correct output is x = 40." << endl;
  cout << "x = " << BptrConst->getx() << endl;  
  //Bptr = F.test12(); // not allowed because Bptr is nonconst

  BptrConst = F.test13(); // this is okay.
  cout << "Correct output is x = 40." << endl;
  cout << "x = " << BptrConst->getx() << endl;  

  Bptr = F.test14(); // this is okay.
  cout << "Correct output is x = 65." << endl;
  cout << "x = " << Bptr->getx() << endl;  

  Bptr = F.test15(); // this is okay.
  cout << "Correct output is x = 65." << endl;
  cout << "x = " << Bptr->getx() << endl;  


  return(0);
}

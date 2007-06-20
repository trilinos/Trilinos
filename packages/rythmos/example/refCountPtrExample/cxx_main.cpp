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

#include "Teuchos_RCP.hpp"
#include<iostream>

class Foo
{
  public:
    Foo() {};
    virtual ~Foo() {};
    virtual double getx() = 0;
};

class Bar : public virtual Foo
{
  public:
    Bar() {};
    ~Bar() {};
    double getx() { return(x_); };
    void setx(double x) { x_ = x; };
  protected:
    double x_;
};

class PrintFoo
{
  public:
    PrintFoo() {};
    PrintFoo(Teuchos::RCP<Foo> &F) { F_ = F; };
    ~PrintFoo() {};
    void setFoo(Teuchos::RCP<Foo> &F) { F = Teuchos::null; };
    void setConstFoo(const Teuchos::RCP<Foo> &F) { F_ = F; };
    void print()
      { std::cout << "x = " << F_->getx() << "!" << std::endl; };
  protected:
    Teuchos::RCP<Foo> F_;
};


// Print foo's x value when passed by RCP
void printFooDirect(Teuchos::RCP<Foo> F)
{
  std::cout << "x = " << F->getx() << "!" << std::endl;
}
// Print foo's x value when passed by reference
void printFooDirect(Foo &F)
{
  std::cout << "x = " << F.getx() << "!" << std::endl;
}

int main(int argc, char *argv[])
{
  std::cout << "Output should be 5: " << std::endl;

  /*
  // This works because this is exactly the type of cast the compiler is expecting.
  Bar B;
  B.setx(5.0);
  printFooDirect(B);
  */

  /*
  // This works because the function printFooDirect is able to make the
  // RCP cast correctly somehow.
  Teuchos::RCP<Bar> B = Teuchos::rcp(new Bar);
  B->setx(5.0);
  printFooDirect(B);
  */
  
  /*
  // This fails because the PrintFoo constructor is not able to make the
  // RCP cast correctly.
  Teuchos::RCP<Bar> B = Teuchos::rcp(new Bar);
  B->setx(5.0);
  PrintFoo PF(B); // fails because PrintFoo takes RCP<Foo> not RCP<Bar>
  PF.print();
  */
  
  /*
  // This fails because PrintFoo.setFoo is not able to make the RCP
  // cast correctly either.
  Teuchos::RCP<Bar> B = Teuchos::rcp(new Bar);
  B->setx(5.0);
  PrintFoo PF;
  PF.setFoo(B);
  PF.print();
  */
  
  /*
  // This does work because the RCP on input to PrintFoo::setConstFoo
  // is a const reference.  This allows the compiler to do the inheritance
  // through the RCP (most likely due to Ross only implementing
  // inheritance for the const RCP constructor).  Contrary to my
  // intuition, the const in front of the RCP, simply says that the
  // RCP itself will not be swapped out with another inside this
  // routine.  It still allows the RCP to modify itself with the count,
  // if copied.  This is the correct way to pass RCP objects, by const
  // reference, provided you're not creating a new RCP and you want to
  // pass it out.
  Teuchos::RCP<Bar> B = Teuchos::rcp(new Bar);
  B->setx(5.0);
  PrintFoo PF;
  PF.setConstFoo(B);
  PF.print();
  */

  
  /*
  // This fails because B is cast (by RCP) as an object of type Foo which
  // doesn't have a member function setx.
  Teuchos::RCP<Foo> B = Teuchos::rcp(new Bar);
  B->setx(5.0); // fails because B is RCP<Foo> not RCP<Bar>
  PrintFoo PF(B);
  PF.print();
  */
  
  
  /*  
  // This does work because the constructor gets the right object and we
  // dynamic cast to get access to the new member function setx.
  Teuchos::RCP<Foo> B = Teuchos::rcp(new Bar);
  (Teuchos::rcp_dynamic_cast<Bar>(B))->setx(5.0);
  PrintFoo PF(B);
  PF.print();
  */  

  /*
  // This does work because we've made the correct RCP from a
  // reference.  In this case, the RCP is not managing the object, so
  // if the function that created the reference goes out of scope, then the
  // RCP will be pointing at garbage.
  Bar B;
  B.setx(5.0);
  Teuchos::RCP<Foo> F = Teuchos::rcp( &B, false );
  PrintFoo PF(F);
  PF.print();
  */

  
  /*
  // This does work because everything is of the right type, and we extracted a
  // reference to the object from the RCP, so when this routine goes
  // out of scope the reference goes away harmlessly.
  Teuchos::RCP<Foo> F = Teuchos::rcp(new Bar);
  Bar &B = *Teuchos::rcp_dynamic_cast<Bar>(F);
  B.setx(5.0);
  PrintFoo PF(F);
  PF.print();
  */

  // This does not work.  It compiles, but has a segfault at runtime because
  // the setFoo member function takes a nonconst reference to the RCP
  // and therefore is allowed to change the variable on exit, in this case
  // making it null, so that when it tries to print, it gets a null pointer to
  // the value.  So using const as in setConstFoo, does not allow the routine
  // to change the variable passed in, which would result in a compile error if
  // we tried to change it.  
  Teuchos::RCP<Foo> F = Teuchos::rcp(new Bar);
  Bar &B = *Teuchos::rcp_dynamic_cast<Bar>(F);
  B.setx(5.0);
  PrintFoo PF;
  PF.setFoo(F);
  PF.print();
  
  return 0;
}

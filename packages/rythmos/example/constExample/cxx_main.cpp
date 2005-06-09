
#include "Teuchos_RefCountPtr.hpp"
#include<iostream>

using namespace Teuchos;

// Dummy class to contain a double 
class Bar 
{
  public:
    // Constructor
    Bar();
    // Destructor
    ~Bar();
    void setx(double x);
    double getx(); 
  protected:
    double x_;
};
Bar::Bar()
{ 
  cout << "Bar::Bar  - address = " << this << endl;
  x_ = 5.0; 
};
Bar::~Bar() 
{ 
  cout << "Bar::~Bar - address = " << this << endl;
  x_ = 0.0; 
};
void Bar::setx(double x)
{
  x_ = x;
};
double Bar::getx() 
{ 
  return(x_); 
};

// Class to test out const
class Foo
{
  public:
    Foo();
    ~Foo();
    void setBar(RefCountPtr<Bar> &Bptr);
//    void test1(RefCountPtr<Bar> &Bptr);
//    void test2(const Bar * Bptr); // also valid:  test2( Bar const * Bptr )
//    void test3(Bar * Bptr);
//    void test4(const Bar * const B) const;

  protected:
    double t_;
    RefCountPtr<Bar> Bptr_;

};
Foo::Foo() 
{ 
  cout << "Foo::Foo  - address = " << this << endl;
  cout << "Foo::Foo  - address of internal Bar pointer = " << &*Bptr_ << endl;
};
Foo::~Foo() 
{
  cout << "Foo::~Foo - address = " << this << endl;
};
void Foo::setBar(RefCountPtr<Bar> &Bptr) 
{ 
  cout << "Foo::setBar" << endl;
  Bptr_ = Bptr; 
  cout << "Foo::setBar - address of internal Bar pointer = " << &*Bptr_ << endl;
};
// test1 shows that if const comes after the *, the function can still modify
// the underlying object.
/*
void Foo::test1(Bar * const Bptr)
{
  cout << "Foo::test1" << endl;
  Bptr->setx(15.0);
};
*/
/*
// test2 shows that if const comes before the *, the function cannot modify the
// underlying object.  This fails to compile.
void Foo::test2(const Bar * Bptr)
{
  cout << "Foo::test2" << endl;
  Bptr->setx(15.0);
};
*/
/*
void Foo::test3(Bar * Bptr)
{
  cout << "Foo::test3" << endl;
  cout << "Foo::test3 - Bptr address = " << Bptr << endl;
  Bar *Bptr_tmp = new Bar;
  Bptr = Bptr_tmp;
  cout << "Foo::test3 - Bptr address = " << Bptr << endl;
};
*/
/*
void Foo::test2(const Bar * Bptr)
{
  cout << "Foo::test2" << endl;
  Bptr->setx(20.0);
};
*/


int main(int argc, char *argv[])
{
  cout << "main:  This routine tests const conditions." << endl;

  Foo F;
  RefCountPtr<Bar> Bptr = rcp(new Bar);
  Bptr->setx(10.0);
  F.setBar(Bptr);
  cout << "Correct output is x = 10." << endl;
  cout << "x = " << Bptr->getx() << endl;

  /*
  F.test1(Bptr);
  cout << "address of Bptr = " << Bptr << endl;
  cout << "Correct output is x = 15." << endl;
  cout << "x = " << Bptr->getx() << endl;

  Bar *Bptr_tmp(NULL);
  cout << "address of Bptr_tmp = " << Bptr_tmp << endl;
  F.test3(Bptr_tmp);
  cout << "address of Bptr_tmp = " << Bptr_tmp << endl;
  */
  /*
  F.test2(B);
  */

  return(0);
};

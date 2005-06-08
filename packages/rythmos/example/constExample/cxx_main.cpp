
//#include "Teuchos_RefCountPtr.hpp"
#include<iostream>

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
  std::cout << "Bar::Bar  - address = " << this << std::endl;
  x_ = 5.0; 
};
Bar::~Bar() 
{ 
  std::cout << "Bar::~Bar - address = " << this << std::endl;
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
    void setBar(Bar * Bptr);
    void test1(Bar * const Bptr);
//    void test2(const Bar * Bptr); // also valid:  test2( Bar const * Bptr )
    void test3(Bar * Bptr);
//    void test4(const Bar * const B) const;

  protected:
    double t_;
    Bar *Bptr_;

};
Foo::Foo() 
{ 
  std::cout << "Foo::Foo  - address = " << this << std::endl;
  std::cout << "Foo::Foo  - address of internal Bar pointer = " << Bptr_ << std::endl;
};
Foo::~Foo() 
{
  std::cout << "Foo::~Foo - address = " << this << std::endl;
};
void Foo::setBar(Bar * Bptr) 
{ 
  std::cout << "Foo::setBar" << std::endl;
  Bptr_ = Bptr; 
  std::cout << "Foo::setBar - address of internal Bar pointer = " << Bptr_ << std::endl;
};
// test1 shows that if const comes after the *, the function can still modify
// the underlying object.
void Foo::test1(Bar * const Bptr)
{
  std::cout << "Foo::test1" << std::endl;
  Bptr->setx(15.0);
};
/*
// test2 shows that if const comes before the *, the function cannot modify the
// underlying object.  This fails to compile.
void Foo::test2(const Bar * Bptr)
{
  std::cout << "Foo::test2" << std::endl;
  Bptr->setx(15.0);
};
*/
void Foo::test3(Bar * Bptr)
{
  std::cout << "Foo::test3" << std::endl;
  std::cout << "Foo::test3 - Bptr address = " << Bptr << std::endl;
  Bar *Bptr_tmp = new Bar;
  Bptr = Bptr_tmp;
  std::cout << "Foo::test3 - Bptr address = " << Bptr << std::endl;
};
/*
void Foo::test2(const Bar * Bptr)
{
  std::cout << "Foo::test2" << std::endl;
  Bptr->setx(20.0);
};
*/


int main(int argc, char *argv[])
{
  std::cout << "main:  This routine tests const conditions." << std::endl;

  Foo F;
  Bar *Bptr = new Bar;
  Bptr->setx(10.0);
  F.setBar(Bptr);
  std::cout << "Correct output is x = 10." << std::endl;
  std::cout << "x = " << Bptr->getx() << std::endl;

  F.test1(Bptr);
  std::cout << "address of Bptr = " << Bptr << std::endl;
  std::cout << "Correct output is x = 15." << std::endl;
  std::cout << "x = " << Bptr->getx() << std::endl;

  Bar *Bptr_tmp(NULL);
  std::cout << "address of Bptr_tmp = " << Bptr_tmp << std::endl;
  F.test3(Bptr_tmp);
  std::cout << "address of Bptr_tmp = " << Bptr_tmp << std::endl;
  /*
  F.test2(B);
  */

  return(0);
};

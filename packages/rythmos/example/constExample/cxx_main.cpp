
//#include "Teuchos_RefCountPtr.hpp"
#include<iostream>

// Dummy class to contain a double 
class Bar 
{
  public:
    // Constructor
    Bar() 
    { 
      std::cout << "Bar::Bar  - address = " << this << std::endl;
      x_ = 5.0; 
    };
    // Destructor
    ~Bar() 
    { 
      std::cout << "Bar::~Bar - address = " << this << std::endl;
      x_ = 0.0; 
    };
    void setx(double x) { x_ = x; };
    double getx() { return(x_); };
  protected:
    double x_;
};

// Class to test out const
class Foo
{
  public:
    Foo() 
    { 
      std::cout << "Foo::Foo  - address = " << this << std::endl;
      std::cout << "Foo::Foo  - address of internal Bar pointer = " << Bptr_ << std::endl;
    };
    ~Foo() 
    {
      std::cout << "Foo::~Foo - address = " << this << std::endl;
    };
    void setBar(Bar * Bptr) 
    { 
      std::cout << "Foo::setBar" << std::endl;
      Bptr_ = Bptr; 
      std::cout << "Foo::setBar - address of internal Bar pointer = " << Bptr_ << std::endl;
    };
    void test1(Bar *Bptr1, Bar *Bptr2);
//    void test2(const Bar *B);
//    void test3(const Bar * const B);
//    void test4(const Bar * const B) const;

  protected:
    double t_;
    Bar *Bptr_;

};
void Foo::test1(Bar *Bptr1, Bar *Bptr2)
{
  std::cout << "Foo::test1" << std::endl;
  Bar *Bptr_tmp;
  Bptr_tmp = Bptr1;
  Bptr1 = Bptr2; 
  Bptr2 = Bptr_tmp;
};
/*
void Foo::test2(const Bar *B)
{
  std::cout << "Foo::test2" << std::endl;
  Bar B2;
  B2.setx(10.0);
  B = B2; // copy of objects
};
*/


int main(int argc, char *argv[])
{
  std::cout << "main:  This routine tests const conditions." << std::endl;

  Foo F;
  Bar *Bptr1 = new Bar; 
  Bar *Bptr2 = new Bar; 
  std::cout << "main:  Address of Bar pointer 1 before test1 = " << Bptr1 << std::endl;
  std::cout << "main:  Address of Bar pointer 2 before test1 = " << Bptr2 << std::endl;

  F.test1(Bptr1,Bptr2);
  std::cout << "main:  Address of Bar pointer 1  after test1  = " << Bptr1 << std::endl;
  std::cout << "main:  Address of Bar pointer 2  after test1  = " << Bptr2 << std::endl;

  Bar *Bptr_tmp;
  Bptr_tmp = Bptr1;
  Bptr1 = Bptr2;
  Bptr2 = Bptr_tmp;
  std::cout << "main:  Address of Bar pointer 1  after swap   = " << Bptr1 << std::endl;
  std::cout << "main:  Address of Bar pointer 2  after swap   = " << Bptr2 << std::endl;

  delete Bptr1;
  delete Bptr2;

  /*
  F.test2(B);
  */

  return(0);
};


//#include "Teuchos_RefCountPtr.hpp"
#include<iostream>

class Foo
{
  public:
    Foo() {};
    ~Foo() {};
    void setBar(const Bar &B) { B_ = B; };
    void test1(Bar &B);
    void test2(const Bar &B);
    void test3(const Bar * const B);
    void test4(const Bar * const B) const;

  protected:
    Bar &B_;

};
void Foo:test1(Bar &B)
{
  Bar B2;
  B2.setx(10.0);
  B = B2;
};

class Bar 
{
  public:
    Bar() { x_ = 5.0; };
    ~Bar() { x_ = 0.0; };
    setx(double x) { x_ = x; };
    double getx() { return(x_); };
  protected:
    double x_;
};

int main(int argc, char *argv[])
{
  std::cout << "This routine tests const conditions." << std::endl;

  Foo F;
  Bar B;
  F.setBar(B);

  F.test1(B);
  std::cout << "x = " << B.getx() << "." << std::endl;

  return 0;
};


#include "Teuchos_RefCountPtr.hpp"
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
    PrintFoo(Teuchos::RefCountPtr<Foo> &F) { F_ = F; };
    ~PrintFoo() {};
    void print()
      { std::cout << "x = " << F_->getx() << "!" << std::endl; };
  protected:
    Teuchos::RefCountPtr<Foo> F_;
};


// Print foo's x value when passed by RefCountPtr
void printFooDirect(Teuchos::RefCountPtr<Foo> F)
{
std::cout << "x = " << F->getx() << "!" << std::endl;
};
// Print foo's x value when passed by reference
void printFooDirect(Foo &F)
{
  std::cout << "x = " << F.getx() << "!" << std::endl;
};

int main(int argc, char *argv[])
{
  std::cout << "Output should be 5: " << std::endl;

  /*
  Bar B;
  B.setx(5.0);
  printFooDirect(B);
  */

  
  /*
  Teuchos::RefCountPtr<Bar> B = Teuchos::rcp(new Bar);
  B->setx(5.0);
  PrintFoo PF(B); // fails because PrintFoo takes RefCountPtr<Foo> not RefCountPtr<Bar>
  PF.print();
  */
  
  
  /*
  Teuchos::RefCountPtr<Foo> B = Teuchos::rcp(new Bar);
  B->setx(5.0); // fails because B is RefCountPtr<Foo> not RefCountPtr<Bar>
  PrintFoo PF(B);
  PF.print();
  */
  
  
  Teuchos::RefCountPtr<Foo> B = Teuchos::rcp(new Bar);
  (Teuchos::rcp_dynamic_cast<Bar>(B))->setx(5.0);
  PrintFoo PF(B);
  PF.print();

  return(0);
};

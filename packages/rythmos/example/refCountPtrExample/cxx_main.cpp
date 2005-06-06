
#include "Teuchos_RefCountPtr.hpp"
#include<iostream>


class Foo
{
  public:
    Foo();
    virtual ~Foo();
    virtual double getx() = 0;
};

class Bar : public virtual Foo
{
  public:
    Bar();
    ~Bar();
    double getx() { return(x_); };
    void setx(double x) { x_ = x; };
  protected:
    double x_;
};


// Print foo's x value when passed by RefCountPtr
void printFoo(Teuchos::RefCountPtr<Foo> F)
{
  cout << "x = " << F->getx() << "!" << endl;
};
// Print foo's x value when passed by reference
void printFoo(Foo &F)
{
  std::cout << "x = " << F.getx() << "!" << std::endl;
};

int main(int argc, char *argv[])
{
  std::cout << "Output should be 5.0: " << endl;

  Bar B;
  B.setx(5.0);
  printFoo(B);

  /*
  Teuchos::RefCountPtr<Bar> B = Teuchos::rcp(new Bar);
  B->setx(5.0);
  printFoo(B); // fails because printFoo takes RefCountPtr<Foo> not RefCountPtr<Bar>
  */
 

  /*
  Teuchos::RefCountPtr<Foo> B = Teuchos::rcp(new Bar);
  B->setx(5.0); // fails because B is RefCountPtr<Foo> not RefCountPtr<Bar>
  printFoo(B);
  */
  
  /*
  Teuchos::RefCountPtr<Foo> B = Teuchos::rcp(new Bar);
  (Teuchos::rcp_dynamic_cast<Bar>(B))->setx(5.0);
  printFoo(B);
  */

  return(0);
};

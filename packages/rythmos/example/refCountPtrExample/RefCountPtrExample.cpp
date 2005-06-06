
#include "Teuchos_RefCountPtr.hpp"

class foo
{
  public:
    foo();
    virtual ~foo();
    virtual double getx() = 0;
}

class bar : public virtual foo
{
  public:
    bar();
    ~bar();
    double getx() { return(x_); };
    void setx(double x) { x_ = x };
  protected:
    double x_;
}

void printfoo(Teuchos::RefCountPtr<foo> foo_)
{
  cout << "x = " << foo_->getx() << "!" << endl;
}

int main(int argc, char *argv[])
{
  Teuchos::RefCountPtr<foo> B = Teuchos::rcp(new bar);
  rcp_dynamic_cast<bar>(B)->setx(5.0);
  printfoo(B);

  return(0);
}

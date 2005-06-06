
#include "Teuchos_RefCountPtr.hpp"

class UtilityFoo
{
  public:
    UtilityFoo();
    virtual ~UtilityFoo();
    virtual double getx() = 0;
};

class UtilityBar : public virtual UtilityFoo
{
  public:
    UtilityBar();
    ~UtilityBar();
    double getx() { return(x_); };
    void setx(double x) { x_ = x; };
  protected:
    double x_;
};

void printFoo(Teuchos::RefCountPtr<UtilityFoo> foo)
{
  cout << "x = " << foo->getx() << "!" << endl;
};

int main(int argc, char *argv[])
{
  /*
  Teuchos::RefCountPtr<UtilityBar> bar = Teuchos::rcp(new UtilityBar);
  bar->setx(5.0);
  printFoo(bar); // fails because printFoo takes RefCountPtr<UtilityFoo> not RefCountPtr<UtilityBar>
  */

  /*
  Teuchos::RefCountPtr<UtilityFoo> bar = Teuchos::rcp(new UtilityBar);
  bar->setx(5.0); // fails because bar is RefCountPtr<UtilityFoo> not RefCountPtr<UtilityBar>
  printFoo(bar);
  */
  
  Teuchos::RefCountPtr<UtilityFoo> bar = Teuchos::rcp(new UtilityBar);
  (Teuchos::rcp_dynamic_cast<UtilityBar>(bar))->setx(5.0);
  printFoo(bar);
  

  return(0);
};

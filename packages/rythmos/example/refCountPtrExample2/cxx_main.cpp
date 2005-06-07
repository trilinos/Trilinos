
#include "Teuchos_RefCountPtr.hpp"
#include<iostream>

class objectA
{
  public:
    objectA() { valueA_ = 5.0; valueB_ = 10.0; };
    ~objectA() { valueA_ = 0.0;  valueB_ = 0.0; };
    double getA() { return(valueA_); };
    double getB() { return(valueB_); };
  protected:
    double valueA_;
    double valueB_;
};

class PrintA
{
  public:
    PrintA() {};
    ~PrintA() {};
    void setA(Teuchos::RefCountPtr<objectA> &A) { A_ = A; };
    void print() { std::cout << "valueA = " << A_->getA() << "." << std::endl; }; 
  protected:
    Teuchos::RefCountPtr<objectA> A_;
};

class PrintB
{
  public:
    PrintB() {};
    ~PrintB() {};
    void setA(Teuchos::RefCountPtr<objectA> &A) { A_ = A; };
    void print() { std::cout << "valueB = " << A_->getB() << "." << std::endl; }; 
  protected:
    Teuchos::RefCountPtr<objectA> A_;
};


void createObjectA(PrintA &PA, PrintB &PB, int TEST)
{
  if (TEST == 1) 
  {
    std::cout <<
      "TEST 1: " << 
      "This creates a RefCountPtr to objectA and gives it to the two print " <<
      "objects and then goes out of scope.  Everything is okay because it is a " <<
      "RefCountPtr instead of a raw pointer." << std::endl;
    Teuchos::RefCountPtr<objectA> OAptr = Teuchos::rcp( new objectA );
    PA.setA(OAptr);
    PB.setA(OAptr);
  } else if (TEST == 2)
  {
    std::cout <<
      "TEST 2: " <<
      "This creates a reference to an objectA, creates an unmanaged RefCountPtr " << 
      "to it and passes it to the two print objects and then goes out of scope.  " << 
      "This is not okay and will result in garbage in the RefCountPtr for the " << 
      "subsequent print operations." << std::endl;
    objectA OA;
    Teuchos::RefCountPtr<objectA> OAptr = Teuchos::rcp( &OA, false );
    PA.setA(OAptr);
    PB.setA(OAptr);
  } else if (TEST == 3)
  {
    std::cout <<
      "This creates a valid RefCountPtr pointing to objectA, and then extracts a " << 
      "reference to it, then gives the RefCountPtr to the two print objects and " << 
      "then goes out of scope.  This works because a reference is just a new name " <<
      "for the underlying object." << std::endl;
    Teuchos::RefCountPtr<objectA> OAptr = Teuchos::rcp( new objectA );
    objectA &OA = *OAptr;
    PA.setA(OAptr);
    PB.setA(OAptr);
  } else if (TEST == 4)
  {
    std::cout <<
      "This creates a valid RefCountPtr pointint to objectA, and then extracts " << 
      "its raw pointer and deletes it.  This is one way to shoot yourself " << 
      "in the foot with RefCountPtr.  " << std::endl;
    Teuchos::RefCountPtr<objectA> OAptr = Teuchos::rcp( new objectA );
    objectA *OA = &*OAptr;
    delete(OA);
    PA.setA(OAptr);
    PB.setA(OAptr);
  } 
};


int main(int argc, char *argv[])
{
  for(int i=1;i<=4;++i)
  {
    PrintA PA;
    PrintB PB;
    createObjectA(PA,PB,i);
    std::cout << "Output should be:  valueA = 5 and valueB = 10." << std::endl;
    PA.print();
    PB.print();
    std::cout << std::endl;
  }

  return 0;
}

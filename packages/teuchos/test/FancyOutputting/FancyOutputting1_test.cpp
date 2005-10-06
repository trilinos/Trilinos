#include "Teuchos_FancyOStream1.hpp"
#include "Teuchos_VerboseObject1.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

class AlgorithmA : public Teuchos::VerboseObject {
public:
  void doAlgorithm()
    {
      Teuchos::RefCountPtr<std::ostream> out = this->getOStream();
      *out << "\nEntering AlgorithmA::doAlgorithm()\n";
      if(1) {
        Teuchos::OSTab tab = this->getOSTab();
        *out << "\nTabbing and doing some stuff ...\n";
      }
      *out << "\nLeaving AlgorithmA::doAlgorithm()\n";
    }

};

int main(int argc, char* argv[])
{

  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObject;
  using Teuchos::OSTab;

  bool success = true;

  try {

    VerboseObject::setDefaultOStream(rcp(new FancyOStream(rcp(&std::cout,false),"  ")));

    std::cout << "\n***\n*** Testing basic FancyOStream and OSTab classes\n***\n";

    RefCountPtr<std::ostream>
      out = VerboseObject::getDefaultOStream();

    *out << "\nThis is very good output\nand I like it a lot!\n";
    *out << "";
    *out << "\n";
    *out << "This should";
    *out << " all be";
    *out << " printed on";
    *out << " the same";
    *out << " line two lines below the above output!\n";
    const int tag2 = OSTab::getNextTag();
    FancyOStream out2(rcp(&std::cout,false),"  ",tag2);
    if(1) {
      OSTab tab;
      *out << "This should be indented one tab!\n";
      if(1) {
        OSTab tab;
        *out << "This should be indented two tabs!\n";
        out2 << "This should be indented zero tabs from out2!\n";
        if(1) {
          OSTab tab(tag2);
          *out << "This should be indented two tabs!\n";
          out2 << "This should be indented one tab from out2!\n";
        }
      }
      *out << "This should be indented one tab!\n";
    }
    *out << "This should be indented zero tabs!\n";
    
    *out << std::endl; // This required overflow() to be overridden!
  
    std::cout << "\n***\n*** Testing VerboseObject base class use\n***\n";
  
    AlgorithmA algoA;
  
    if(1) {
      OSTab tab;
      algoA.doAlgorithm();
    }

    *out << std::endl;

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cout,success)

  return ( success ? 0 : 1 );
  
}

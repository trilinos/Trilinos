#include "Teuchos_VerboseObject2.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_dyn_cast.hpp"

// This is a typical function that would be present in Trilinos right now what
// does not know about FancyOStream and does not derive from VerboseObject.
// However, because of the magic of FancyOStream, this output will be indented
// correctly!
void someDumbFunction( std::ostream &out, const std::string &indentSpacer )
{
  out << "\nEntering someDumbFunction(...)\n";
  if(1) {
    out << std::endl << indentSpacer << "I am \"dumb\" code that knows nothing of FancyOStream and does indenting manually! ...\n";
  }
  out << "\nLeaving someDumbFunction(...)\n";
  // Note that this output will be indented correctly even through it knows nothing of FancyOStream
}

// This is a typical numerical class that derives from VerboseObject and does
// outputting.  Note that the use of the OSTab class requires initialization
// using VerboseObject::getOSTab(...) which takes care of the hassles and is
// easy to use.
class AlgorithmA : public Teuchos::VerboseObject {
public:
  //
  AlgorithmA()
    {
      this->setLinePrefix("ALGO_A"); // I tell me who I am for line prefix outputting
    }
  //
  void doAlgorithm()
    {
      // Here I grab the stream that I will use for outputting.  It is a good
      // idea to grab the RCP to this object just to be safe.
      Teuchos::RefCountPtr<Teuchos::FancyOStream> out = this->getOStream();
      // Here I make sure that my line prefix is set but I do not what to
      // indent since I might assume that the caller will do the indentation.
      // I am not sure who should be responsible for this but we just need to
      // come up with a convention
      Teuchos::OSTab tab = this->getOSTab(0); // This sets the line prefix but does not change the indent!
      *out << "\nEntering AlgorithmA::doAlgorithm()\n";
      if(1) {
        // Here I use a simple macro for the typical case of one tab indent to
        // save typing.  The idea is that this should be as easy to write as
        // OSTab tab; but is more general.
        TEUCHOS_OSTAB;
        *out << "\nI am \"smart\" code that knows about FancyOStream and OSTab ...\n";
        if(1) {
          // Here I am going to be calling a dumb piece of code that does not
          // know about the FancyOStream system and will not use tabs or
          // anything like that.  There is a lot of code in Trilinos that
          // falls in this category.  The first thing I do is manually indent
          // the stream one tab and set a line prefix for the dumb code since
          // it may not do this itself.
          Teuchos::OSTab tab = this->getOSTab(1,"DUMB_ALGO");
          // Now a Pass in the updated FancyOStream object, which is properly
          // indented now, through the std::ostream interface.  I also pass in
          // the string that is being used for creating tabs.  The output from
          // this function will be indented correctly without the dumb code
          // knowing it!
          someDumbFunction(*out,out->getTabIndentStr());
        }
      }
      *out << "\nLeaving AlgorithmA::doAlgorithm()\n";
    }
};

// Here is a simple driver function that I call over and over to show
// different features of FancyOStream
void doAlgorithmStuff()
{

  // Here I just create the algorithm object that derives from VerboseObject.
  // By default, this object will print to *Verbose::getDefaultOStream()
  AlgorithmA algoA;
  // Note that here I could change the stream just this object prints to
  // by calling algoA.setOStream(...).
  
  // Now I call the algorithm which will print to its default output stream
  algoA.doAlgorithm();
  
  *algoA.getOStream() << std::endl;
  
}

int main(int argc, char* argv[])
{

  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObject;
  using Teuchos::OSTab;
  using Teuchos::dyn_cast;

  bool success = true;

  try {

    // Start by setting up a defualt FancyOStream with a new indent string.
    // This output stream object will be used by default for all VerboseObject outputting
    VerboseObject::setDefaultOStream(rcp(new FancyOStream(rcp(&std::cout,false),"  ")));

    // Here I am just grabbing the default output stream
    RefCountPtr<FancyOStream>
      out = VerboseObject::getDefaultOStream();
    // Note that the VerboseObject manages FancyOStream objects and not just
    // std::ostream objects.  This is important to the design and very
    // resonable I think.
    
    //
    // Now I call doAlgorithmStuff() a bunch of times with different setups to
    // show the different kinds of line prefix options
    //
  
    std::cout << "\n***\n*** Testing VerboseObject base class use\n***\n";
  
    std::cout << "\n*** Algorithm output with no line prefixes\n\n";

    doAlgorithmStuff();
  
    std::cout << "\n*** Algorithm output with line prefixes\n\n";

    VerboseObject::getDefaultOStream()->setShowLinePrefix(true);

    doAlgorithmStuff();
  
    std::cout << "\n*** Algorithm output with line prefixes and tab counts\n\n";

    VerboseObject::getDefaultOStream()->setShowTabCount(true);

    doAlgorithmStuff();
  
    std::cout << "\n*** Algorithm output with tab counts\n\n";

    VerboseObject::getDefaultOStream()->setShowLinePrefix(false);

    doAlgorithmStuff();

    //
    // Now I do some other simple tests just to see that FancyOStream is working
    // correctly
    //

    std::cout << "\n***\n*** Testing basic FancyOStream and OSTab classes\n***\n\n";

    *out << "\nThis is very good output\nand I like it a lot!\n";
    *out << "";
    *out << "\n";
    *out << "This should";
    *out << " all be";
    *out << " printed on";
    *out << " the same";
    *out << " line two lines below the above output!\n";
    RefCountPtr<FancyOStream>
      out2 = rcp(new FancyOStream(rcp(new std::ostringstream),"  "));
    if(1) {
      OSTab tab(out);
      *out << "This should be indented one tab!\n";
      if(1) {
        OSTab tab(out);
        *out << "This should be indented two tabs!\n";
        *out2 << "This should be indented zero tabs from out2!\n";
        if(1) {
          OSTab tab(out2);
          *out << "This should be indented two tabs!\n";
          *out2 << "This should be indented one tab from out2!\n";
        }
      }
      *out << "This should be indented one tab!\n";
    }
    *out << "This should be indented zero tabs!\n";
    
    *out << std::endl; // This required overflow() to be overridden!

    *out << "\n***\n*** Now outputting the latent output that was sent to out2\n***\n\n"
         << dyn_cast<std::ostringstream>(*out2->getOStream()).str();
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cout,success)
    
  return ( success ? 0 : 1 );
  
}

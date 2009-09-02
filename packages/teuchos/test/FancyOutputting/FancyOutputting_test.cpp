#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Version.hpp"

// This is a typical function that would be present in Trilinos right now what
// does not know about FancyOStream and does not derive from VerboseObject.
// However, because of the magic of FancyOStream, this output will be indented
// correctly!
void someDumbFunction( std::ostream &out, const std::string &indentSpacer )
{
  out << "\nEntering someDumbFunction(...)\n";
  {
    out << std::endl << indentSpacer << "I am \"dumb\" code that knows nothing of FancyOStream and does indenting manually! ...\n";
  }
  out << "\nLeaving someDumbFunction(...)\n";
  // Note that this output will be indented correctly even through it knows nothing of FancyOStream
}

// This is a function who's interface was written before there was a
// FancyOStream and therefore is written in terms of std::ostream.  However,
// in this case the implementation has been modifed to use FancyOStream as
// shown.
void someLessDumbFunction( std::ostream &out_arg )
{
  using Teuchos::OSTab;
  // Get a FancyOStream from out_arg or create a new one ...
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::getFancyOStream(Teuchos::rcp(&out_arg,false));
  // Do our tab indent and our name.
  OSTab tab(out,1,"LDUMBALGO");
  *out << "\nEntering someLessDumbFunction(...)\n";
  {
    Teuchos::OSTab(out_arg).o()
      << std::endl << "I am less \"dumb\" code that knows about FancyOStream but my interface does not support it directly! ...\n";
    *Teuchos::tab(out)
      << std::endl << "Another print from this less \"dumb\" code ...\n";
  }
  *out << "\nLeaving someLessDumbFunction(...)\n";
}

// This is a typical numerical class that derives from VerboseObject and does
// outputting.  Note that the use of the OSTab class requires initialization
// using VerboseObject::getOSTab(...) which takes care of the hassles and is
// easy to use.
//
// This class also derives from ParameterListAcceptor and uses helper
// functio  ns to read options for VerboseObject from a parameter sublist.
class AlgorithmA
  : public Teuchos::VerboseObject<AlgorithmA>,
    public Teuchos::ParameterListAcceptor
{
public:

  // Constructor(s)

  AlgorithmA();

  // Overridden from ParameterListAccpetor

  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);

  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();

  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();

  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;

  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  // Other functions

  void doAlgorithm();

private:

  enum EAlgoType { ALGO_BOB, ALGO_JOHN, ALGO_HARRY };

  static const std::string toString( AlgorithmA::EAlgoType algoType );

  Teuchos::RCP<Teuchos::ParameterList> paramList_;
  EAlgoType algoType_;
  double algoTol_;
  
};


// Implementations for AlgorithmA

namespace {

const std::string AlgoType_name = "Algo Type";
const std::string AlgoType_default = "Bob";

const std::string AlgoTol_name = "Algo Tol";
const double AlgoTol_default = 1e-5;

} // namespace

const std::string AlgorithmA::toString( AlgorithmA::EAlgoType algoType )
{
  switch(algoType) {
    case ALGO_BOB: return "Bob";
    case ALGO_JOHN: return "John";
    case ALGO_HARRY: return "Harry";
    default: TEST_FOR_EXCEPT("Should never get here!");
  }
  return ""; // never be called!
}


AlgorithmA::AlgorithmA()
  : algoType_(ALGO_BOB), algoTol_(AlgoTol_default)
{
  this->setLinePrefix("ALGO_A"); // I tell me who I am for line prefix outputting
}


void AlgorithmA::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& paramList
  )
{
  TEST_FOR_EXCEPT(is_null(paramList));
  // Validate and set the parameter defaults.  Here, the parameters are
  // validated and the state of *this is not changed unless the parameter
  // validation succeeds.  Also, any validators that are defined for various
  // parameters are passed along so that they can be used in extracting
  // values!
  paramList->validateParametersAndSetDefaults(*this->getValidParameters(),0);
  paramList_ = paramList;
  // Get the enum value for the algorithm type. Here, the actual type stored
  // for the algorithm type in the parameter list is an std::string but this
  // helper function does all the work of extracting the validator set in
  // getValidParameters() and set on *paramList_ through the
  // validateParametersAndSetDefaults(...) function above!
  algoType_ = Teuchos::getIntegralValue<EAlgoType>(*paramList_,AlgoType_name);
  // Get the tolerance for the algorithm.  Here, the actual type of the
  // parameter stored on input could be many different types.  Here, I can
  // just assume that it is a double since it would have been converted to a
  // double above in validateParametersAndSetDefaults(...).
  algoTol_ = Teuchos::getParameter<double>(*paramList_,AlgoTol_name);
  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&*paramList_,this);
#ifdef TEUCHOS_DEBUG
  paramList_->validateParameters(*this->getValidParameters());
#endif
}


Teuchos::RCP<Teuchos::ParameterList>
AlgorithmA::getNonconstParameterList()
{
  return paramList_;
}


Teuchos::RCP<Teuchos::ParameterList>
AlgorithmA::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> paramList = paramList_;
  paramList_ = Teuchos::null;
  return paramList;
}


Teuchos::RCP<const Teuchos::ParameterList>
AlgorithmA::getParameterList() const
{
  return paramList_;
}


Teuchos::RCP<const Teuchos::ParameterList>
AlgorithmA::getValidParameters() const
{
  using Teuchos::RCP; using Teuchos::ParameterList;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;
  static RCP<const ParameterList> validParams;
  if (is_null(validParams)) {
    RCP<ParameterList>
      pl = Teuchos::rcp(new ParameterList("AlgorithmA"));
    setStringToIntegralParameter<EAlgoType>(
      AlgoType_name, AlgoType_default,
      "The algorithm type to use",
      tuple<std::string>("Bob", "John", "Harry"),
      tuple<EAlgoType>(ALGO_BOB, ALGO_JOHN, ALGO_HARRY),
      &*pl
      );
    Teuchos::setDoubleParameter(
      AlgoTol_name, AlgoTol_default,
      "The tolerance for the algorithm.",
      &*pl
      );
    Teuchos::setupVerboseObjectSublist(&*pl);
    validParams = pl;
  }
  return validParams;
}


void AlgorithmA::doAlgorithm()
{
  using Teuchos::OSTab;
  // Get the verbosity that we are going to use
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  // Here I grab the stream that I will use for outputting.  It is a good
  // idea to grab the RCP to this object just to be safe.
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  // Here I set my line prefix and a single indent.  The convention will
  // be that a called function will set its own indent.  This convention makes
  // the most sense.
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab
  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true))
    *out << "\nEntering AlgorithmA::doAlgorithm() with verbLevel="<<Teuchos::toString(verbLevel)<<"\n";
  {
    // Here I use a simple macro for the typical case of one tab indent to
    // save typing.  The idea is that this should be as easy to write as
    // OSTab tab; but is more general.
    TEUCHOS_OSTAB;
    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true))
      *out
        << "\nI am \"smart\" code that knows about FancyOStream and OSTab ...\n"
        << "\nDoing algorithm of type \""<<toString(algoType_)<<"\""
        << "\nUsing tolerance of " << algoTol_ << "\n";
    {
      // Here I temporaraly turn off tabbing so that I can print an imporant warning message.
      OSTab tab2 = this->getOSTab(OSTab::DISABLE_TABBING);
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true))
        *out << "\n***\n*** Warning, I am doing something very dangerous so watch out!!!\n***\n";
    }
    if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true))
      *out << "\nHere I am doing some more stuff and printing with indenting turned back on!\n";
    {
      // Here I am going to be calling a dumb piece of code that does not
      // know about the FancyOStream system and will not use tabs or
      // anything like that.  There is a lot of code in Trilinos that
      // falls in this category.  The first thing I do is manually indent
      // the stream one tab and set a line prefix for the dumb code since
      // it may not do this itself.
      OSTab tab2 = this->getOSTab(1,"DUMBALGO");
      // Now a Pass in the updated FancyOStream object, which is properly
      // indented now, through the std::ostream interface.  I also pass in
      // the std::string that is being used for creating tabs.  The output from
      // this function will be indented correctly without the dumb code
      // knowing it!
      someDumbFunction(*out,out->getTabIndentStr());
    }
    // Here I am calling a less dumb piece of code who's interface does
    // not support FancyOStream but the implementation does.  Note that
    // this function also follows the convention of doing an initial
    // indent.
    someLessDumbFunction(*out);
  }
  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true))
    *out << "\nLeaving AlgorithmA::doAlgorithm()\n";
}


//
// Here is a simple driver function that I call over and over to show
// different features of FancyOStream
//

void doAlgorithmStuff( Teuchos::ParameterList *algoParams = 0 )
{

  // Here I just create the algorithm object that derives from VerboseObject.
  // By default, this object will print to *Verbose::getDefaultOStream()
  AlgorithmA algoA;
  if(algoParams)
    algoA.setParameterList(Teuchos::rcp(algoParams,false));
  // Note that here I could change the stream just this object prints to
  // by calling algoA.setOStream(...).
  
  // Now I call the algorithm which will print to its default output stream
  algoA.doAlgorithm();
  
  *algoA.getOStream() << std::endl;
  
}

//
// Test that static initailziation of VerboseObjectBase and VerboseObject works!
//

class TestVerboseObjectBaseInitialization {
public:
  TestVerboseObjectBaseInitialization()
    {
      // Get the verbosity level for AlgorithmA
      Teuchos::EVerbosityLevel verbLevel = Teuchos::VerboseObject<AlgorithmA>::getDefaultVerbLevel();
      TEST_FOR_EXCEPT_PRINT(verbLevel!=Teuchos::VERB_DEFAULT,&std::cerr);
      // Print to the default default OStream to make sure that the initialization
      // trick worked!
      *Teuchos::VerboseObjectBase::getDefaultOStream()
        << "\n***\n*** Printing to default OStream before main() even starts!\n***\n\n"
        << std::flush;
    }
};

static TestVerboseObjectBaseInitialization testVerboseObjectBaseInitialization;

//
// Main driver program
//

int main(int argc, char* argv[])
{

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::OSTab;
  using Teuchos::dyn_cast;
  using Teuchos::CommandLineProcessor;

  bool success = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  const int numProcs = Teuchos::GlobalMPISession::getNProc();

  try {

    // Get some commandline options
		CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
		CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) return parse_return;

    // Here I am just grabbing the default output stream
    RCP<FancyOStream>
      out = VerboseObjectBase::getDefaultOStream();
    // Note that the VerboseObject manages FancyOStream objects and not just
    // std::ostream objects.  This is important to the design and very
    // resonable I think.

    *out << std::endl << Teuchos::Teuchos_Version() << std::endl << std::endl;

    //
    // Now I call doAlgorithmStuff() a bunch of times with different setups to
    // show the different kinds of line prefix options
    //
  
    *out << "\n***\n*** Testing VerboseObject base class use\n***\n";
  
    *out << "\n*** Algorithm output with default formatting\n\n";
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with no front matter\n\n";
    out->setShowAllFrontMatter(false);
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with processor ranks\n\n";
    out->setShowAllFrontMatter(false).setShowProcRank(true);
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with line prefix names\n\n";
    out->setShowAllFrontMatter(false).setShowLinePrefix(true);
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with tab counts\n\n";
    out->setShowAllFrontMatter(false).setShowTabCount(true);
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with line prefix names and tab counts\n\n";
    out->setShowAllFrontMatter(false).setShowLinePrefix(true).setShowTabCount(true);
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with processor ranks and line prefix names\n\n";
    out->setShowAllFrontMatter(false).setShowProcRank(true).setShowLinePrefix(true);
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with processor ranks and tab counts\n\n";
    out->setShowAllFrontMatter(false).setShowProcRank(true).setShowTabCount(true);
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with processor ranks, line prefix names, and tab counts\n\n";
    out->setShowAllFrontMatter(false).setShowProcRank(true).setShowLinePrefix(true).setShowTabCount(true);
    doAlgorithmStuff();
  
    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Algorithm output with processor ranks, line prefix names, and tab counts but no output for AlgorithmA\n\n";
    Teuchos::VerboseObject<AlgorithmA>::setDefaultVerbLevel(Teuchos::VERB_NONE);
    out->setShowAllFrontMatter(false).setShowProcRank(true).setShowLinePrefix(true).setShowTabCount(true);
    doAlgorithmStuff();
    Teuchos::VerboseObject<AlgorithmA>::setDefaultVerbLevel(Teuchos::VERB_DEFAULT);

    *out << "\n*** Running the algorithm by setting parameters in the parameter list ...\n";

    Teuchos::ParameterList algoParams("AlgorithmA");

    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Set AlgorithmA verbosity level to extreme through a parameter list\n\n";
    algoParams.sublist("VerboseObject").set("Verbosity Level","extreme");
    algoParams.set("Algo Type","Harry");
    algoParams.set("Algo Tol",0.3);
    doAlgorithmStuff(&algoParams);

    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Set AlgorithmA verbosity level to medium and the output file \"AlgorithmA.out\" through a parameter list\n\n";
    algoParams.sublist("VerboseObject").set("Verbosity Level","medium");
    algoParams.sublist("VerboseObject").set("Output File","AlgorithmA.out");
    algoParams.set("Algo Type","John");
    algoParams.set("Algo Tol",10);
    doAlgorithmStuff(&algoParams);

    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n*** Set AlgorithmA verbosity level to low and the output back to default through a parameter list\n\n";
    algoParams.sublist("VerboseObject").set("Verbosity Level","low");
    algoParams.sublist("VerboseObject").set("Output File","none");
    algoParams.set("Algo Tol","20");
    doAlgorithmStuff(&algoParams);

    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1);
    *out << "\n***\n*** Do some more simple tests to make sure things work correctly\n***\n\n";

    //
    // Now I do some other simple tests just to see that FancyOStream is working
    // correctly
    //

    out->setShowAllFrontMatter(false).setShowProcRank(numProcs>1).setShowTabCount(true);
    out->setProcRankAndSize(mpiSession.getRank(),mpiSession.getNProc());
    
    *out << "\n***\n*** Testing basic FancyOStream and OSTab classes\n***\n\n";
    
    *out << "\nThis is very good output\nand I like it a lot!\n";
    *out << "";
    *out << "\n";
    *out << "This should";
    *out << " all be";
    *out << " printed on";
    *out << " the same";
    *out << " line two lines below the above output!\n";
    RCP<FancyOStream>
      out2 = rcp(new FancyOStream(rcp(new std::ostringstream),"  "));
    {
      OSTab tab1(out);
      *out << "This should be indented one tab!\n";
      {
        OSTab tab2(out);
        *out << "This should be indented two tabs!\n";
        *out2 << "This should be indented zero tabs from out2!\n";
        {
          OSTab tab3(out2);
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

    if(success)
      *out << "\nEnd Result: TEST PASSED" << std::endl;
    
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);
    
  return ( success ? 0 : 1 );
  
}

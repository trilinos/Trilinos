/*
// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

#include "AlgorithmA.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"


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

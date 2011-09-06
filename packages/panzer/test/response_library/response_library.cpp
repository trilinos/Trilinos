#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterEntryValidator.hpp"
#include "Teuchos_StrUtils.hpp"

#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_PauseToAttach.hpp"
#include "Panzer_ResponseContainer.hpp"
#include "Panzer_ResponseLibrary.hpp"

#include "TestEvaluators.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

namespace panzer {

TEUCHOS_UNIT_TEST(response_library, test)
{
  typedef Traits::Residual EvalT;

  // build global (or serial communicator)
  #ifdef HAVE_MPI
     Teuchos::RCP<Teuchos::Comm<int> > comm 
           = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
  #else
     Teuchos::RCP<Teuchos::Comm<int> > comm = Teuchos::rcp(new Teuchos::SerialComm<int>);
  #endif
 
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  // panzer::pauseToAttach();

  int worksetSize = 10;
  Teuchos::ParameterList mainParam;
  mainParam.set<int>("Workset Size", worksetSize);

  // build basic field manager
  PHX::FieldManager<Traits> fm;
  {
     RCP<PHX::Evaluator<Traits> > testEval 
        = rcp(new TestEvaluator<EvalT,Traits>(mainParam));
     fm.registerEvaluator<EvalT>(testEval);
  }

  ResponseId dResp  = buildResponse("Dog","Functional");
  ResponseId hResp  = buildResponse("Horse","Functional");
  RCP<ResponseLibrary<Traits> > rLibrary 
        = Teuchos::rcp(new ResponseLibrary<Traits>());
  rLibrary->defineDefaultAggregators();

  rLibrary->reserveVolumeResponse<EvalT>(dResp,"block_0");
  rLibrary->reserveVolumeResponse<EvalT>(hResp,"block_1");

  // test uninitialized access
  TEST_THROW(rLibrary->getVolumeResponse(dResp,"block_0"),std::logic_error);

  std::vector<std::string> eBlocks;
  rLibrary->getRequiredElementBlocks(eBlocks);
  std::sort(eBlocks.begin(),eBlocks.end());
  TEST_EQUALITY(eBlocks.size(),2);
  TEST_EQUALITY(eBlocks[0],"block_0");
  TEST_EQUALITY(eBlocks[1],"block_1");
}

class ResponseEntryValidator : public Teuchos::ParameterEntryValidator {
public:
  ResponseEntryValidator() {}

  ValidStringsList validStringValues() const
  { return Teuchos::null; }

  void validate(const Teuchos::ParameterEntry & entry, 
                const std::string & paramName,
                const std::string & sublistName) const;

  const std::string getXMLTypeName() const
  { return "response"; }

  void printDoc(const std::string & docString, std::ostream &out) const;

  //! Utility function for tokenizing
  static void split(const std::string & str,
                    const std::string & delim,
                    std::vector<std::string> & tokens);
};

void ResponseEntryValidator::split(const std::string & str,const std::string & delim,
                                   std::vector<std::string> & output)
{
   output.clear();

   typedef boost::tokenizer<boost::char_separator<char> > 
           tokenizer;

   boost::char_separator<char> sep(delim.c_str());
   tokenizer tokens(str, sep);
   for(tokenizer::iterator tok_iter = tokens.begin();
       tok_iter != tokens.end(); ++tok_iter) {
      // extract token, remove spaces
      std::string s = *tok_iter;
      boost::trim(s);
      output.push_back(s);
   }
}

void ResponseEntryValidator::validate(Teuchos::ParameterEntry const &entry,  
                                      const std::string & paramName,
                                      const std::string & sublistName) const
{
  const std::string &entryName = entry.getAny(false).typeName();
  Teuchos::any anyValue = entry.getAny(true);
  
  // type passed, validate value
  TEST_FOR_EXCEPTION(!(anyValue.type() == typeid(std::string) ),
    Teuchos::Exceptions::InvalidParameterType,
    "Sorry but it looks like the \"" << paramName << "\"" <<
    " parameter in the \"" << sublistName <<
    "\" sublist does not exist." << std::endl << std::endl <<
    "Error: The value that you entered was the wrong type." << std::endl <<
    "Parameter: " << paramName << std::endl <<
    "Type specified: " << entryName << std::endl <<
    "Type accepted: " << typeid(std::string).name() <<
    std::endl << std::endl);

  const std::string & value = Teuchos::any_cast<std::string>(anyValue);

  std::vector<std::string> tokens;
  split(value,":",tokens);  

  // verify that their is a response type and a evaluation type
  TEST_FOR_EXCEPTION(tokens.size()!=2,
     Teuchos::Exceptions::InvalidParameterValue,
     "The value for \"response\" type parameter named \""+paramName+"\" "
     "is incorrectly formatted. The expected format is\n"
     "   \"<Response type>: <Evaluation type> [, <Evaluation type>]*\""
     "your value is \""+value+"\"");

  std::string respType = tokens[0];  boost::trim(respType);
  std::string evalTypes = tokens[1]; boost::trim(evalTypes);

  TEST_FOR_EXCEPTION(respType.length()==0 || evalTypes.length()==0,
     Teuchos::Exceptions::InvalidParameterValue,
     "The value for \"response\" type parameter named \""+paramName+"\" "
     "in sublist \"" + sublistName + "\" "
     "is incorrectly formatted. The expected format is\n"
     "   \"<Response type>: <Evaluation type> [, <Evaluation type>]*\"\n"
     "your value is \""+value+"\"");

  // must have at least one evaluation type
  split(evalTypes,",",tokens);  
  TEST_FOR_EXCEPTION(tokens.size()==0,
     Teuchos::Exceptions::InvalidParameterValue,
     "The value for \"response\" type parameter named \""+paramName+"\" "
     "in sublist \"" + sublistName + "\" "
     "is incorrectly formatted. The expected format is\n"
     "   \"<Response type>: <Evaluation type> [, <Evaluation type>]*\"\n"
     "your value is \""+value+"\"");
}


void ResponseEntryValidator::printDoc(
  std::string const &docString, std::ostream &out) const
{
  Teuchos::StrUtils::printLines(out,"# ",docString);
  out << "#  Validator Used: " << std::endl;
  out << "#  ResponseEntry Validator" << std::endl;
}


TEUCHOS_UNIT_TEST(response_library, pl_reader)
{
   Teuchos::ParameterList pl;
   pl.sublist("block_0");
   pl.set("Cat","None:Residual,Jacobian","",Teuchos::rcp(new ResponseEntryValidator));
   pl.sublist("block_1");

   std::vector<std::string> tokens;
   ResponseEntryValidator::split(" None :Residual, Jacobian",":",tokens);
    
   TEST_EQUALITY(tokens.size(),2);
   TEST_EQUALITY(tokens[0],"None");
   TEST_EQUALITY(tokens[1],"Residual, Jacobian");

   std::string next = tokens[1];
   ResponseEntryValidator::split(next,",",tokens);

   TEST_EQUALITY(tokens.size(),2);
   TEST_EQUALITY(tokens[0],"Residual");
   TEST_EQUALITY(tokens[1],"Jacobian");

   ResponseEntryValidator rev;

   // standard operating
   {
      Teuchos::ParameterEntry entry;

      entry.setValue<std::string>("Functional : Residual, Jacobian");
      TEST_NOTHROW(rev.validate(entry,"FuncValue","ResponseList"));

      entry.setValue<std::string>("Functional : Residual");
      TEST_NOTHROW(rev.validate(entry,"FuncValue","ResponseList"));
   }

   // failing tests
   {
      Teuchos::ParameterEntry entry;

      entry.setValue<std::string>("Functional : ");
      TEST_THROW(rev.validate(entry,"FuncValue","ResponseList"),Teuchos::Exceptions::InvalidParameterValue);

      entry.setValue<std::string>(" : Residual");
      TEST_THROW(rev.validate(entry,"FuncValue","ResponseList"),Teuchos::Exceptions::InvalidParameterValue);

      entry.setValue<std::string>(" Residual");
      TEST_THROW(rev.validate(entry,"FuncValue","ResponseList"),Teuchos::Exceptions::InvalidParameterValue);
   }
}


}

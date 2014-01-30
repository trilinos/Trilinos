#include "Pike_BlackBoxModelEvaluator_DefaultBase.hpp"

namespace pike {

  BlackBoxModelEvaluatorDefaultBase::~BlackBoxModelEvaluatorDefaultBase() {}

  bool BlackBoxModelEvaluatorDefaultBase::supportsParameter(const std::string& pName) const
  {
    return false;
  }
  
  int BlackBoxModelEvaluatorDefaultBase::getNumberOfParameters() const
  {
    return 0;
  }
  
  std::string BlackBoxModelEvaluatorDefaultBase::getParameterName(const int l) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Error: BlackBoxModelEvaluatorDefaultBase::getParameterName(j) "
			       << "The BlackBoxModelEvaluator named \"" << this->name() 
			       << "\" does not support parameters!");
    return "";
  }
  
  int BlackBoxModelEvaluatorDefaultBase::getParameterIndex(const std::string& pName) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Error: BlackBoxModelEvaluatorDefaultBase::getParameterIndex(name) "
			       << "The BlackBoxModelEvaluator named \"" << this->name() 
			       << "\" does not support parameters!");
    return 0;
  }
  
  
  void BlackBoxModelEvaluatorDefaultBase::setParameter(const int l, const Teuchos::ArrayView<const double>& p)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Error: BlackBoxModelEvaluatorDefaultBase::getParameter(j) "
			       << "The BlackBoxModelEvaluator named \"" << this->name() 
			       << "\" does not support parameters!");
  }
  
  bool BlackBoxModelEvaluatorDefaultBase::supportsResponse(const std::string& rName) const
  {
    return false;
  }
  
  int BlackBoxModelEvaluatorDefaultBase::getNumberOfResponses() const
  {
    return 0;
  }
  
  std::string BlackBoxModelEvaluatorDefaultBase::getResponseName(const int j) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Error: BlackBoxModelEvaluatorDefaultBase::getResponseName(j) "
			       << "The BlackBoxModelEvaluator named \"" << this->name() 
			       << "\" does not support responses!");
    return "";
  }
  
  int BlackBoxModelEvaluatorDefaultBase::getResponseIndex(const std::string& rName) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Error: BlackBoxModelEvaluatorDefaultBase::getResponseIndex(name) "
			       << "The BlackBoxModelEvaluator named \"" << this->name() 
			       << "\" does not support responses!");
    return 0;
  }
  
  Teuchos::ArrayView<const double> BlackBoxModelEvaluatorDefaultBase::getResponse(const int j) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Error: BlackBoxModelEvaluatorDefaultBase::getResponse(j) "
			       << "The BlackBoxModelEvaluator named \"" << this->name() 
			       << "\" does not support responses!");
    return Teuchos::ArrayView<const double>(std::vector<double>(0));
  }

}

#include "Pike_BlackBoxModelEvaluator.hpp"

namespace pike {

  // ***********************
  // Base methods
  // ***********************

  BlackBoxModelEvaluator::~BlackBoxModelEvaluator()
  {}
  
  bool BlackBoxModelEvaluator::isGloballyConverged() const
  { return true; }

  // ***********************
  // Parameter Support
  // ***********************

  bool BlackBoxModelEvaluator::supportsParameter(const std::string& pName) const
  {
    return false;
  }
  
  int BlackBoxModelEvaluator::getNumberOfParameters() const
  {
    return 0;
  }
  
  std::string BlackBoxModelEvaluator::getParameterName(const int l) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Error: pike::BlackBoxModelEvaluator::getParameterName(j) "
			       << "The BlackBoxModelEvaluator named \"" << this->name() 
			       << "\" does not support parameters!");
    return "";
  }
  
  int BlackBoxModelEvaluator::getParameterIndex(const std::string& pName) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Error: pike::BlackBoxModelEvaluator::getParameterIndex(name) "
			       << "The BlackBoxModelEvaluator named \"" << this->name() 
			       << "\" does not support parameters!");
    return 0;
  }
  
  
  void BlackBoxModelEvaluator::setParameter(const int l, const Teuchos::ArrayView<const double>& p)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Error: pike::BlackBoxModelEvaluator::getParameter(j) "
			       << "The BlackBoxModelEvaluator named \"" << this->name() 
			       << "\" does not support parameters!");
  }

  // ***********************
  // Response Support
  // ***********************

  bool BlackBoxModelEvaluator::supportsResponse(const std::string& rName) const
  {
    return false;
  }
  
  int BlackBoxModelEvaluator::getNumberOfResponses() const
  {
    return 0;
  }
  
  std::string BlackBoxModelEvaluator::getResponseName(const int j) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Error: pike::BlackBoxModelEvaluator::getResponseName(j) "
			       << "The BlackBoxModelEvaluator named \"" << this->name() 
			       << "\" does not support responses!");
    return "";
  }
  
  int BlackBoxModelEvaluator::getResponseIndex(const std::string& rName) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Error: pike::BlackBoxModelEvaluator::getResponseIndex(name) "
			       << "The BlackBoxModelEvaluator named \"" << this->name() 
			       << "\" does not support responses!");
    return 0;
  }
  
  Teuchos::ArrayView<const double> BlackBoxModelEvaluator::getResponse(const int j) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Error: pike::BlackBoxModelEvaluator::getResponse(j) "
			       << "The BlackBoxModelEvaluator named \"" << this->name() 
			       << "\" does not support responses!");
    return Teuchos::ArrayView<const double>(std::vector<double>(0));
  }

  // ***********************
  // Transient Support
  // ***********************

  bool BlackBoxModelEvaluator::isTransient() const
  { return false; }
  
  double BlackBoxModelEvaluator::getCurrentTime() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
			       "Error: pike::BlackBoxModelEvaluator::getCurrentTime() is not implemented for "
			       << "the BlackBoxModelEvaluator named \"" << this->name() << "." << std::endl);
    return 0.0;
  }

  double BlackBoxModelEvaluator::getTentativeTime() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
			       "Error: pike::BlackBoxModelEvaluator::getTentativeTime() is not implemented for "
			       << "the BlackBoxModelEvaluator named \"" << this->name() << "." << std::endl);
    return 0.0;
  }
  
  bool BlackBoxModelEvaluator::solvedTentativeStep() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
			       "Error: pike::BlackBoxModelEvaluator::solvedTentativeTime() is not implemented for "
			       << "the BlackBoxModelEvaluator named \"" << this->name() << "." << std::endl);
    return false;
  }
  
  double BlackBoxModelEvaluator::getCurrentTimeStepSize() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
			       "Error: pike::BlackBoxModelEvaluator::solvedTentativeTime() is not implemented for "
			       << "the BlackBoxModelEvaluator named \"" << this->name() << "." << std::endl);
    return 0.0;
  }
  
  double BlackBoxModelEvaluator::getDesiredTimeStepSize() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
			       "Error: pike::BlackBoxModelEvaluator::getDesiredTimeStepSize() is not implemented for "
			       << "the BlackBoxModelEvaluator named \"" << this->name() << "." << std::endl);
    return 0.0;
  }
  
  double BlackBoxModelEvaluator::getMaxTimeStepSize() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
			       "Error: pike::BlackBoxModelEvaluator::getMaxTimeStepSize() is not implemented for "
			       << "the BlackBoxModelEvaluator named \"" << this->name() << "." << std::endl);
    return 0.0;
  }
  
  void BlackBoxModelEvaluator::setNextTimeStepSize(const double& dt)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
			       "Error: pike::BlackBoxModelEvaluator::setNextTimeStepSize() is not implemented for "
			       << "the BlackBoxModelEvaluator named \"" << this->name() << "." << std::endl);
  }
  
  void BlackBoxModelEvaluator::acceptTimeStep()
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
			       "Error: pike::BlackBoxModelEvaluator::acceptTimeStep() is not implemented for "
			       << "the BlackBoxModelEvaluator named \"" << this->name() << "." << std::endl);
  }
  
}

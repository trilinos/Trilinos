#ifndef PANZER_MODEL_EVALUATOR_UTILITIES_HPP
#define PANZER_MODEL_EVALUATOR_UTILITIES_HPP

#include <string>
#include <tuple>

namespace Thyra {
  template<typename ScalarT> class ModelEvaluator;
}

namespace panzer {
  /// Given a parameter name and a ModelEvaluator, returns the index and subindex of the parameter in the ModelEvaluator.
  std::tuple<int,int> findParameterIndex(const std::string& p_name,const Thyra::ModelEvaluator<double>& me, const bool& throw_if_not_found=true);

  /// Given a response name and a ModelEvaluator, returns the index and subindex of the response in the ModelEvaluator.
  std::tuple<int,int> findResponseIndex(const std::string& g_name,const Thyra::ModelEvaluator<double>& me, const bool& throw_if_not_found=true);
}

#endif

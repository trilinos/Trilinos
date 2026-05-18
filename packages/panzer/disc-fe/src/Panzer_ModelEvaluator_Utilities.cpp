#include <Panzer_ModelEvaluator_Utilities.hpp>
#include <Thyra_ModelEvaluator.hpp>
#include <Teuchos_Assert.hpp>

std::tuple<int,int> panzer::findParameterIndex(const std::string& p_name,const Thyra::ModelEvaluator<double>& me, const bool& throw_if_not_found)
{
  bool found = false;
  int index = -1;
  int sub_index = -1;
  for (int p = 0; p < me.Np(); ++p) {
    auto p_names = me.get_p_names(p);
    for (Teuchos::Ordinal p_sub=0; p_sub < p_names->size(); ++p_sub) {
      if ((*p_names)[p_sub] == p_name) {
        found = true;
        index = p;
        sub_index = p_sub;
        break;
      }
    }
    if (found) break;
  }

  TEUCHOS_TEST_FOR_EXCEPTION( (!found and throw_if_not_found), std::runtime_error,
    "ERROR: the parameter \"" << p_name << "\" is not a valid parameter name in the model evaluator!");

  return {index,sub_index};
}

std::tuple<int,int> panzer::findResponseIndex(const std::string& g_name,const Thyra::ModelEvaluator<double>& me, const bool& throw_if_not_found)
{
  bool found = false;
  int index = -1;
  int sub_index = -1;
  for (int g = 0; g < me.Ng(); ++g) {
    auto g_names = me.get_g_names(g);
    for (Teuchos::Ordinal g_sub=0; g_sub < g_names.size(); ++g_sub) {
      std::cout << "g(" << g << "," << g_sub << ")=" << g_names[g_sub] << std::endl;
      if (g_names[g_sub] == g_name) {
        found = true;
        index = g;
        sub_index = g_sub;
        break;
      }
    }
    if (found) break;
  }

  TEUCHOS_TEST_FOR_EXCEPTION( (!found and throw_if_not_found), std::runtime_error,
    "ERROR: the response \"" << g_name << "\" is not a valid response name in the model evaluator!");

  return {index,sub_index};
}

#ifndef KRINO_KRINO_SURFACE_AKRI_ML_MODELS_HPP_
#define KRINO_KRINO_SURFACE_AKRI_ML_MODELS_HPP_

#include <memory>
#include <string>

namespace krino {

bool have_ML_Model_support();

class MLModelInterface
{
public:
  virtual ~MLModelInterface();
};

class ML_Model
{
public:
  ML_Model();
  ~ML_Model();
  void create(const std::string & modelBaseName);
  void evaluate(const size_t dimX, const double * x, const size_t dimY, double * y) const;

private:
  std::unique_ptr<MLModelInterface> myModel;
};

template<size_t NUMVAR>
class ML_Model_with_Gradient
{
public:
  ML_Model_with_Gradient();
  ~ML_Model_with_Gradient();
  void create(const std::string & modelBaseName);
  void evaluate(const size_t dimX, const double * x, const size_t dimY, double * y, double * dydx) const;
  bool is_enabled() const;

private:
  std::unique_ptr<MLModelInterface> myModel;
};

}

#endif /* KRINO_KRINO_SURFACE_AKRI_ML_MODELS_HPP_ */

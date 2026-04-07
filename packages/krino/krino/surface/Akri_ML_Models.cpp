#include <Akri_ML_Models.hpp>

#include <stdexcept>
#include <Akri_config.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include "Sacado.hpp"

#ifdef KRINO_HAVE_POCKET_TENSOR
#include "pt_model.h"
#include "pt_tensor.h"
#include "pt_dispatcher.h"
#define POCKET_TENSOR_SUPPORTS_FAD 0
#endif

namespace krino {

#ifdef KRINO_HAVE_POCKET_TENSOR
static constexpr bool have_pocket_tensor_support() { return true; }
#else
static constexpr bool have_pocket_tensor_support() { return false; }
#endif

bool have_ML_Model_support() { return have_pocket_tensor_support(); }

template<typename DATATYPE>
class UnsupportedPocketTensorInterface : public MLModelInterface
{
public:
  ~UnsupportedPocketTensorInterface() = default;

  static std::unique_ptr<UnsupportedPocketTensorInterface<DATATYPE>> create(const std::string & modelBaseName)
  {
    throw std::runtime_error("ML_Model support not enabled for type " + std::string(typeid(DATATYPE).name()));
    std::unique_ptr<UnsupportedPocketTensorInterface<DATATYPE>> ptr;
    return ptr;
  }
  void evaluate(const size_t dimX, const DATATYPE * x, const size_t dimY, DATATYPE * y) const
  {
    throw std::runtime_error("ML_Model support not enabled for type " + std::string(typeid(DATATYPE).name()));
  }
};

#ifdef KRINO_HAVE_POCKET_TENSOR
template<typename DATATYPE>
class PocketTensorInterface : public MLModelInterface
{
public:
  PocketTensorInterface(std::unique_ptr<pt::Model<DATATYPE>> ptModel) : myPtModel(std::move(ptModel)) {}
  ~PocketTensorInterface() = default;

  static std::unique_ptr<PocketTensorInterface<DATATYPE>> create(const std::string & modelBaseName)
  {
    const std::string & modelFileName = modelBaseName + ".model";
    std::unique_ptr<pt::Model<DATATYPE>> ptModel = pt::Model<DATATYPE>::create(modelFileName);
    if (ptModel == nullptr)
      throw std::runtime_error("Model file '" + modelFileName + "' not found");
    return std::make_unique<PocketTensorInterface<DATATYPE>>(std::move(ptModel));
  }
  void evaluate(const size_t dimX, const DATATYPE * x, const size_t dimY, DATATYPE * y) const
  {
    pt::Tensor<DATATYPE> in{dimX};
    for (size_t i=0; i<dimX; ++i)
      in(i) = x[i];

    pt::Tensor<DATATYPE> out;
    const bool success = myPtModel->predict(mySerialDispatcher, in, out);

    if (!success)
      throw std::runtime_error("Model evaluation failed.");
    if (out.getSize() != dimY)
      throw std::runtime_error("Model prediction size does not match.");

    for (size_t i=0; i<dimY; ++i)
      y[i] = out(i);
  }
private:
  std::unique_ptr<pt::Model<DATATYPE>> myPtModel;
  mutable pt::Dispatcher mySerialDispatcher{1}; // mutable to satisfy predict method that takes non-const Dispatcher&
};
#else
template<typename DATATYPE>
using PocketTensorInterface = UnsupportedPocketTensorInterface<DATATYPE>;
#endif

MLModelInterface::~MLModelInterface() = default;

ML_Model::ML_Model() = default;
ML_Model::~ML_Model() = default;

void ML_Model::create(const std::string & modelBaseName)
{
  myModel = PocketTensorInterface<double>::create(modelBaseName);
}

void ML_Model::evaluate(const size_t dimX, const double * x, const size_t dimY, double * y) const
{
  const auto * model = dynamic_cast<const PocketTensorInterface<double>*>(myModel.get());
  STK_ThrowAssertMsg(model != nullptr, "Model is null, has it been created?");
  model->evaluate(dimX, x, dimY, y);
}

template<size_t NUMVAR>
ML_Model_with_Gradient<NUMVAR>::ML_Model_with_Gradient() {}

template<size_t NUMVAR>
ML_Model_with_Gradient<NUMVAR>::~ML_Model_with_Gradient() {}

template<size_t NUMVAR>
bool ML_Model_with_Gradient<NUMVAR>::is_enabled() const
{
#if POCKET_TENSOR_SUPPORTS_FAD
  return true;
#else
  return false;
#endif
}

template<size_t NUMVAR>
void ML_Model_with_Gradient<NUMVAR>::create(const std::string & modelBaseName)
{
  using FADType = Sacado::Fad::SFad<double,NUMVAR>;
#if POCKET_TENSOR_SUPPORTS_FAD
  myModel = PocketTensorInterface<FADType>::create(modelBaseName);
#else
  myModel = UnsupportedPocketTensorInterface<FADType>::create(modelBaseName);
#endif
}

template<size_t NUMVAR>
void ML_Model_with_Gradient<NUMVAR>::evaluate(const size_t dimX, const double * x, const size_t dimY, double * y, double * dydx) const
{
  using FADType = Sacado::Fad::SFad<double,NUMVAR>;
#if POCKET_TENSOR_SUPPORTS_FAD
  const auto * model = dynamic_cast<const PocketTensorInterface<FADType>*>(myModel.get());
#else
  const auto * model = dynamic_cast<const UnsupportedPocketTensorInterface<FADType>*>(myModel.get());
#endif
  STK_ThrowAssertMsg(model != nullptr, "Model is null, has it been created?");
  STK_ThrowAssert(dimX == NUMVAR);

  std::array<FADType,NUMVAR> xFAD;
  for (size_t i=0; i<NUMVAR; ++i)
    xFAD[i] = FADType(NUMVAR, i, x[i]);

  FADType yFAD;
  model->evaluate(dimX, xFAD.data(), dimY, &yFAD);

  for (size_t j=0; j<dimY; ++j)
  {
    y[j] = yFAD.val();
    for (size_t i=0; i<NUMVAR; ++i)
      dydx[j*NUMVAR+i] = yFAD.dx(i);
  }
}

template class ML_Model_with_Gradient<3>;

}


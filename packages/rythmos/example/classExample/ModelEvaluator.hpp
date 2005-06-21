
#ifndef RYTHMOS_MODELEVALUATOR
#define RYTHMOS_MODELEVALUATOR

namespace Rythmos {
template<class Scalar>
class ModelEvaluator
{
  public:
    ModelEvaluator() {};
    virtual ~ModelEvaluator() {};
    virtual double evalModel(double x, double t) const = 0;
    virtual double get_vector() const = 0;
};
} // namespace Rythmos
#endif // RYTHMOS_MODELEVALUATOR

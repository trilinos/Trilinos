
#ifndef RYTHMOS_STEPPER
#define RYTHMOS_STEPPER

namespace Rythmos {
template<class Scalar>
class Stepper
{
  public:
    Stepper() {};
    virtual ~Stepper() {};
    virtual Scalar TakeStep(Scalar dt) = 0;
    virtual Scalar get_solution() = 0;
};
} // namespace Rythmos
#endif // RYTHMOS_STEPPER

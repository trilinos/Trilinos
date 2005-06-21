
#ifndef RYTHMOS_STEPPER
#define RYTHMOS_STEPPER

namespace Rythmos {
class Stepper
{
  public:
    Stepper() {};
    virtual ~Stepper() {};
    virtual double TakeStep(double dt) = 0;
    virtual double get_solution() = 0;
};
} // namespace Rythmos
#endif // RYTHMOS_STEPPER

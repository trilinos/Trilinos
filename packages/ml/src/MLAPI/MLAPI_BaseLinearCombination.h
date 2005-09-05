#ifndef ML_BASELINEARCOMBINATION_H
#define ML_BASELINEARCOMBINATION_H

namespace MLAPI {

class Space;
class BaseOperator;
class MultiVector;

class BaseLinearCombination 
{
public:
  virtual ~BaseLinearCombination() {};

  virtual const Space GetVectorSpace() const = 0;
  // Computes v += <operations>
  virtual void Update(MultiVector& v) const = 0;
  // Computes v = <operations>
  virtual void Set(MultiVector& v) const = 0;
};

} // namespace MLAPI

#endif // ifdef ML_BASELINEARCOMBINATION_H

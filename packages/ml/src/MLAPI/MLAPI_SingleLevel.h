#ifndef MLAPI_SINGLELEVEL_H
#define MLAPI_SINGLELEVEL_H

#include "ml_config.h"

namespace MLAPI {

class Operator;
class Smoother;

class SingleLevel
{
public:
  SingleLevel(int LevelID = 0) :
    LevelID_(LevelID),
    A_(0),
    P_(0),
    R_(0),
    Spre_(0),
    Spost_(0)
  {
    LevelID_ = LevelID;
  }

  SingleLevel(int LevelID, const Operator* A, const Operator* P, const Operator* R, 
              const Smoother* Spre, const Smoother* Spost)
  {
    A_ = A;
    P_ = P;
    R_ = R;
    Spre_ = Spre;
    Spost_ = Spost;
    LevelID_ = LevelID;
  }

  void SetA(const Operator* A)
  {
    A_ = A;
  }

  void SetP(const Operator* P)
  {
    P_ = P;
  }

  void SetR(const Operator* R)
  {
    R_ = R;
  }

  void SetPreSmoother(const Smoother* Spre)
  {
    Spre_ = Spre;
  }

  void SetPostSmoother(const Smoother* Spost)
  {
    Spost_ = Spost;
  }

  const Operator* A() const
  {
    return(A_);
  }

  const Operator* P() const
  {
    return(P_);
  }

  const Operator* R() const
  {
    return(R_);
  }

  const Smoother* PreSmoother() const
  {
    return(Spre_);
  }

  const Smoother* PostSmoother() const
  {
    return(Spost_);
  }

  int LevelID() const
  {
    return(LevelID_);
  }

private:
  const Operator* A_;
  const Operator* P_;
  const Operator* R_;
  const Smoother* Spre_;
  const Smoother* Spost_;
  int LevelID_;

};
} // namespace MLAPI

#endif

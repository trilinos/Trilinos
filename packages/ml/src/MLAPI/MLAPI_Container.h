#ifndef MLAPI_CONTAINER_H
#define MLAPI_CONTAINER_H


namespace MLAPI {

class Operator;
class InverseOperator;

class Container
{
public:
  Container(int LevelID = 0)
  {
    LevelID_ = LevelID;
  }

  Container(int LevelID, const Operator& A, const Operator& P, 
            const Operator& R, const InverseOperator& S)
  {
    A_ = A;
    P_ = P;
    R_ = R;
    S_ = S;
    LevelID_ = LevelID;
  }

  void SetA(const Operator& A)
  {
    A_ = A;
  }

  void SetP(const Operator& P)
  {
    P_ = P;
  }

  void SetR(const Operator& R)
  {
    R_ = R;
  }

  void SetS(const InverseOperator& S)
  {
    S_ = S;
  }

  const Operator& A() const
  {
    return(A_);
  }

  const Operator& P() const
  {
    return(P_);
  }

  const Operator& R() const
  {
    return(R_);
  }

  const InverseOperator& S() const
  {
    return(S_);
  }

  int LevelID() const
  {
    return(LevelID_);
  }

private:
  Operator A_;
  Operator P_;
  Operator R_;
  InverseOperator S_;
  int LevelID_;

};
} // namespace MLAPI

#endif

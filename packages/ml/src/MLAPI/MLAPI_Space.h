#ifndef ML_SPACE_H
#define ML_SPACE_H
#include "ml_config.h"
#include <iostream>
#include "Epetra_Comm.h"

namespace MLAPI {

class Space {
public:
  Space(int NumMyElements, const Epetra_Comm& Comm) :
    Comm_(&Comm),
    NumGlobalElements_(-1),
    NumMyElements_(NumMyElements)
  {
    Comm.SumAll(&NumMyElements_,&NumGlobalElements_,1);
    Comm.ScanSum(&NumMyElements_,&Offset_,1);
    Offset_ -= NumMyElements_;
  }

  const Epetra_Comm& Comm() const {
    return(*Comm_);
  }

  Space(const Space& RHS)
  {
    Comm_ = &RHS.Comm();
    NumMyElements_ = RHS.NumMyElements();
    NumGlobalElements_ = RHS.NumGlobalElements();
    Offset_ = RHS.Offset();
  }

  ~Space() {};

  int NumMyElements() const
  {
    return(NumMyElements_);
  }

  int NumGlobalElements() const
  {
    return(NumGlobalElements_);
  }

  int GID(int i) const
  {
    return(i + Offset_);
  }

  int Offset() const
  {
    return(Offset_);
  }

  bool operator== (const Space& rhs) const
  {
    bool ok = true;
    if (NumGlobalElements_ != rhs.NumGlobalElements())
      ok = false;
    if (NumMyElements_ != rhs.NumMyElements())
      ok = false;
    return(ok);
  }

  bool operator!= (const Space& rhs) const
  {
    return(!(*this != rhs));
  }

private:
  int NumMyElements_;
  int NumGlobalElements_;
  int Offset_;
  const Epetra_Comm* Comm_;
};

} // namespace MLAPI

#endif

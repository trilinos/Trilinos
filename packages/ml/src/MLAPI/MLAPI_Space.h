#ifndef ML_SPACE_H
#define ML_SPACE_H
#include "ml_include.h"
#include "ml_comm.h"
#include "MLAPI_Workspace.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_IntSerialDenseVector.h"

namespace MLAPI {

class Space {
public:
  Space()
  {
    NumMyElements_ = 0;
    NumGlobalElements_ = 0;
    IsLinear_ = false;
    Offset_ = 0;
  }

  Space(int NumMyElements)
  {
    NumMyElements_ = NumMyElements;
    NumGlobalElements_ = ML_Comm_GsumInt(GetMLComm(),NumMyElements);
    Offset_ = ML_gpartialsum_int(NumMyElements,GetMLComm());
    IsLinear_ = true;
  }

  Space(int NumMyElements, const int* MyGlobalElements)
  {
    NumMyElements_ = NumMyElements;
    NumGlobalElements_ = ML_Comm_GsumInt(GetMLComm(),NumMyElements);
    Offset_ = ML_gpartialsum_int(NumMyElements,GetMLComm());
    MyGlobalElements_ = Teuchos::rcp(new Epetra_IntSerialDenseVector);
    MyGlobalElements_->Resize(NumMyElements);
    for (int i = 0 ; i < NumMyElements ; ++i)
      (*MyGlobalElements_)[i] = MyGlobalElements[i];
    IsLinear_ = false;
  }

  Space(const Space& RHS)
  {
    NumMyElements_ = RHS.NumMyElements();
    NumGlobalElements_ = RHS.NumGlobalElements();
    Offset_ = RHS.Offset();
    IsLinear_ = RHS.IsLinear();
    MyGlobalElements_ = RHS.MyGlobalElements();
  }

  ~Space() {};

  Space& operator=(const Space& RHS)
  {
    NumMyElements_ = RHS.NumMyElements();
    NumGlobalElements_ = RHS.NumGlobalElements();
    Offset_ = RHS.Offset();
    IsLinear_ = RHS.IsLinear();
    MyGlobalElements_ = RHS.MyGlobalElements();
    return(*this);
  }

  void Reshape(int NumMyElements)
  {
    NumMyElements_ = NumMyElements;
    NumGlobalElements_ = ML_Comm_GsumInt(GetMLComm(),NumMyElements);
    Offset_ = ML_gpartialsum_int(NumMyElements,GetMLComm());
    IsLinear_ = false;
    MyGlobalElements_ = Teuchos::null;
  }

  void Reshape(int NumMyElements, const int* MyGlobalElements)
  {
    NumMyElements_ = NumMyElements;
    NumGlobalElements_ = ML_Comm_GsumInt(GetMLComm(),NumMyElements);
    Offset_ = ML_gpartialsum_int(NumMyElements,GetMLComm());
    MyGlobalElements_ = Teuchos::rcp(new Epetra_IntSerialDenseVector);
    MyGlobalElements_->Resize(NumMyElements);
    for (int i = 0 ; i < NumMyElements ; ++i)
      (*MyGlobalElements_)[i] = MyGlobalElements[i];
    IsLinear_ = false;
  }

  int NumMyElements() const
  {
    return(NumMyElements_);
  }

  int NumGlobalElements() const
  {
    return(NumGlobalElements_);
  }

  //! 
  inline int Offset() const
  {
    return(Offset_);
  }

  inline bool operator== (const Space& RHS) const
  {
    if (IsLinear() != RHS.IsLinear())
      return(false);
    if (NumGlobalElements() != RHS.NumGlobalElements())
      return(false);
    if (NumMyElements() != RHS.NumMyElements())
      return(false);

    return(true);
  }

  inline bool operator!= (const Space& RHS) const
  {
    return(!(*this != RHS));
  }

  //! Returns \c true if the decomposition among processors is linear.
  inline bool IsLinear() const
  {
    return(IsLinear_);
  }

  //! Returns a pointer to the list of global nodes.
  inline const Teuchos::RefCountPtr<Epetra_IntSerialDenseVector> 
    MyGlobalElements() const
  {
    return(MyGlobalElements_);
  }

  //! Returns the global ID of local element \c i.
  inline int operator() (int i) const
  {
    if (IsLinear())
      return(i + Offset_);
    else
      return((*MyGlobalElements_)[i]);
  }

private:
  //! Number of elements assigned to the calling processor.
  int NumMyElements_;
  //! Total number of elements.
  int NumGlobalElements_;
  //! If \c true, the decomposition among processors is linear.
  bool IsLinear_;
  //! GID of the first local element (for linear decompositions only).
  int Offset_;
  //! Container of global numbering of local elements.
  Teuchos::RefCountPtr<Epetra_IntSerialDenseVector> MyGlobalElements_;
};

} // namespace MLAPI

#endif

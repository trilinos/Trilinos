#ifndef ML_SPACE_H
#define ML_SPACE_H
#include "ml_include.h"
#include "ml_comm.h"
#include "MLAPI_Workspace.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_IntSerialDenseVector.h"

namespace MLAPI {

/*!
\class Space

\brief Specifies the number of elements in a space and its layout among processors.

Class MLAPI::Space defines the number of elements of linear objects (like
DoubleVector's, Operator's, InverseOperator's and Preconditioner's), and
their layout among the processors. Important methods are NumMyElements(),
NumGlobalElements(), IsLinear(), and the operator().

A simple example of usage follows.
\code
int NumMyElements = 2;
MLAPI::Space MySpace(NumMyElements);
int LID = 2;
int GID = MySpace(LID);
\endcode

\author Marzio Sala, SNL 9214

\date Last updated on 07-Jan-05.
*/

class Space {
public:
  //! Default constructor.
  Space()
  {
    NumMyElements_ = 0;
    NumGlobalElements_ = 0;
    IsLinear_ = false;
    Offset_ = 0;
  }

  //! Constructor with specified number of local element on the calling process.
  Space(int NumMyElements)
  {
    NumMyElements_ = NumMyElements;
    NumGlobalElements_ = ML_Comm_GsumInt(GetMLComm(),NumMyElements);
    Offset_ = ML_gpartialsum_int(NumMyElements,GetMLComm());
    IsLinear_ = true;
  }

  //! Constructor with specified number of local element on the calling process and their global numbering (starting from 0).
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

  //! Copy constructor.
  Space(const Space& RHS)
  {
    NumMyElements_ = RHS.NumMyElements();
    NumGlobalElements_ = RHS.NumGlobalElements();
    Offset_ = RHS.Offset();
    IsLinear_ = RHS.IsLinear();
    MyGlobalElements_ = RHS.MyGlobalElements();
  }

  //! Destructor.
  ~Space() {};

  //! Operator =.
  Space& operator=(const Space& RHS)
  {
    NumMyElements_ = RHS.NumMyElements();
    NumGlobalElements_ = RHS.NumGlobalElements();
    Offset_ = RHS.Offset();
    IsLinear_ = RHS.IsLinear();
    MyGlobalElements_ = RHS.MyGlobalElements();
    return(*this);
  }

  //! Returns \c true if \c this Space is equivalent to \c RHS.
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

  //! Returns \c true if \c this Space is not equivalent to \c RHS.
  inline bool operator!= (const Space& RHS) const
  {
    return(!(*this != RHS));
  }

  //! Reset the dimension of the space by specifying the local number of elements.
  void Reshape(int NumMyElements)
  {
    NumMyElements_ = NumMyElements;
    NumGlobalElements_ = ML_Comm_GsumInt(GetMLComm(),NumMyElements);
    Offset_ = ML_gpartialsum_int(NumMyElements,GetMLComm());
    IsLinear_ = false;
    MyGlobalElements_ = Teuchos::null;
  }

  //! Returns the global ID of local element \c i.
  inline int operator() (int i) const
  {
    if (IsLinear())
      return(i + Offset_);
    else
      return((*MyGlobalElements_)[i]);
  }

  //! Reset the dimension of the space by specifying the local number of elements and their global numbering (starting from 0).
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

  //! Returns the local number of elements on the calling process.
  int NumMyElements() const
  {
    return(NumMyElements_);
  }

  //! Returns the global number of elements.
  int NumGlobalElements() const
  {
    return(NumGlobalElements_);
  }

  //! Returns the global ID of the first element on the calling process (for linear distributions only).
  inline int Offset() const
  {
    return(Offset_);
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

std::ostream& operator<< (std::ostream& os, const Space& v) 
{
  os << std::endl;
  os << "MLAPI::Space" << std::endl;
  os << "NumMyElements() = " << v.NumMyElements() << std::endl;
  os << "NumGlobalElements() = " << v.NumGlobalElements() << std::endl;

  os << "ProcID\t\tLID\t\tGID" << std::endl;
  for (int i = 0 ; i < v.NumMyElements() ; ++i)
    os << 0 << "\t\t" << i << "\t\t" << v(i) << std::endl;
  os << std::endl;
  return(os);
}

} // namespace MLAPI

#endif

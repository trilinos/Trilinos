#ifndef ML_SPACE_H
#define ML_SPACE_H

#include "ml_include.h"
#include "ml_comm.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_BaseObject.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_IntSerialDenseVector.h"
#include <iomanip>

namespace MLAPI {

/*!
\class Space

\brief Specifies the number and distribution among processes of elements.

\author Marzio Sala, SNL 9214

\date Last updated on Feb-05.
*/

class Space : public BaseObject {

public:

  //@{ \name Constructors and Destructors
  
  //! Default constructor, defines an empty space.
  Space()
  {
    NumMyElements_     = 0;
    NumGlobalElements_ = 0;
    IsLinear_          = false;
    Offset_            = 0;
  }

  //! Constructor with specified number of global and local elements.
  /*! 
    Constructs a space with linear distribution.

    \param NumGlobalElements - (In) number of global elements.

    \param NumMyElements - (In) number of local elements. If different
                           from -1, then eithere NumGlobalElements == -1,
                           or the sum across of processors of NumMyElements
                           equals NumGlobalElements
   */
  Space(const int NumGlobalElements, const int NumMyElements = -1)
  {
    Reshape(NumGlobalElements, NumMyElements);
  }

  //! Constructor for non-linear distributions
  /*! 
    \param NumGlobalElements - (In) number of global elements. Set to
                               -1 to compute it automatically.
    \param NumMyElement - (In) number of local elements. Cannot be set
                          to -1.
    \param MyGlobalElements - (In) contains the global ID of each local node.

    \note Global ID always starts from 0.
   */
  Space(const int NumGlobalElements, const int NumMyElements, 
        const int* MyGlobalElements)
  {
    Reshape(NumGlobalElements, NumMyElements, MyGlobalElements);
  }

  //! Copy constructor.
  Space(const Space& RHS)
  {
    NumMyElements_       = RHS.GetNumMyElements();
    NumGlobalElements_   = RHS.GetNumGlobalElements();
    Offset_              = RHS.GetOffset();
    IsLinear_            = RHS.IsLinear();
    RCPMyGlobalElements_ = RHS.GetRCPMyGlobalElements();
  }

  //! Destructor.
  ~Space() {};

  //@}
  //@{ \name Reshape methods

  //! Reset the dimension of the space by specifying the local number of elements.
  void Reshape(const int NumGlobalElements, const int NumMyElements = -1)
  {

    if (NumGlobalElements <= 0 && NumMyElements < 0)
      ML_THROW("NumGlobalElements = " + GetString(NumGlobalElements) +
               " and NumMyElements = " + GetString(NumMyElements), -1);

    if (NumMyElements == -1) {
      NumMyElements_ = NumGlobalElements / GetNumProcs();
      if (GetMyPID() == 0)
        NumMyElements_ += NumGlobalElements % GetNumProcs();
    }
    else
      NumMyElements_ = NumMyElements;

    NumGlobalElements_ = ML_Comm_GsumInt(GetML_Comm(),NumMyElements_);

    if (NumGlobalElements != -1) {
      if (NumGlobalElements != NumGlobalElements_) 
        ML_THROW("Specified # of global elements the sum of local elements (" +
                 GetString(NumGlobalElements) + " vs. " +
                 GetString(NumGlobalElements_), -1);
    }

    Offset_   = ML_gpartialsum_int(NumMyElements_,GetML_Comm());
    IsLinear_ = true;
    
  }

  //! Reset the dimension of the space by specifying the local number of elements and their global numbering (starting from 0).
  void Reshape(const int NumGlobalElements, const int NumMyElements, 
               const int* MyGlobalElements)
  {
    if (NumGlobalElements <= 0 && NumMyElements < 0)
      ML_THROW("NumGlobalElements = " + GetString(NumGlobalElements) +
               " and NumMyElements = " + GetString(NumMyElements), -1);

    if (NumMyElements == -1) {
      NumMyElements_ = NumGlobalElements / GetNumProcs();
      if (GetMyPID() == 0)
        NumMyElements_ += NumGlobalElements % GetNumProcs();
    }
    else
      NumMyElements_ = NumMyElements;

    NumGlobalElements_ = ML_Comm_GsumInt(GetML_Comm(),NumMyElements_);

    if (NumGlobalElements != -1) {
      if (NumGlobalElements != NumGlobalElements_)
        ML_THROW("Specified # of global elements the sum of local elements (" +
                 GetString(NumGlobalElements) + " vs. " +
                 GetString(NumGlobalElements_), -1);
    }

    RCPMyGlobalElements_ = Teuchos::rcp(new Epetra_IntSerialDenseVector);
    RCPMyGlobalElements_->Resize(NumMyElements_);
    for (int i = 0 ; i < NumMyElements_ ; ++i)
      (*RCPMyGlobalElements_)[i] = MyGlobalElements[i];

    Offset_   = -1; // not set
    IsLinear_ = false;
  }

  //@}
  //@{ \name Overloaded operators

  //! Operator =.
  Space& operator=(const Space& RHS)
  {
    NumMyElements_       = RHS.GetNumMyElements();
    NumGlobalElements_   = RHS.GetNumGlobalElements();
    Offset_              = RHS.GetOffset();
    IsLinear_            = RHS.IsLinear();
    RCPMyGlobalElements_ = RHS.GetRCPMyGlobalElements();
    return(*this);
  }

  //! Returns \c true if \c this Space is equivalent to \c RHS.
  /*
    \note this is a cheap check
   */
  inline bool operator== (const Space& RHS) const
  {
    if (IsLinear() != RHS.IsLinear())
      return(false);
    if (GetNumGlobalElements() != RHS.GetNumGlobalElements())
      return(false);
    if (GetNumMyElements() != RHS.GetNumMyElements())
      return(false);

    return(true);
  }

  //! Returns \c true if \c this Space is not equivalent to \c RHS.
  /*
    \note this is a cheap check
   */
  inline bool operator!= (const Space& RHS) const
  {
    return(!(*this == RHS));
  }

  //! Sets the Label of \c this object.
  Space& operator=(const string& Label)
  {
    SetLabel(Label);
    return(*this);
  }

  //! Returns the global ID of local element \c i.
  inline int operator() (int i) const
  {
    if (IsLinear())
      return(i + GetOffset());
    else
      return((*RCPMyGlobalElements_)[i]);
  }

  // @}
  // @{ \name Get and Set methods

  //! Returns the local number of elements on the calling process.
  int GetNumMyElements() const
  {
    return(NumMyElements_);
  }

  //! Returns the global number of elements.
  int GetNumGlobalElements() const
  {
    return(NumGlobalElements_);
  }

  //! Returns the global ID of the first element on the calling process (for linear distributions only).
  inline int GetOffset() const
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
    GetRCPMyGlobalElements() const
  {
    return(RCPMyGlobalElements_);
  }

  // @}
  // @{ \name Miscellanous methods

  //! Prints on ostream basic information about \c this object.
  std::ostream& Print(std::ostream& os,
                      const bool verbose = true) const
  {
    os << std::endl;
    if (GetMyPID() == 0) {
      os << "*** MLAPI::Space ***" << endl;
      os << "Label               = " << GetLabel() << endl;
      os << "NumMyElements()     = " << GetNumMyElements() << endl;
      os << "NumGlobalElements() = " << GetNumGlobalElements() << endl;
      os << "Offset              = " << GetOffset() << endl;
      if (IsLinear())
        os << "Distribution is linear" << endl;
      else
        os << "Distribution is not linear" << endl;
      os << endl;
    }

    if (verbose) {
      if (GetMyPID() == 0) {
        os.width(10);
        os << "ProcID";
        os.width(20);
        os << "LID";
        os.width(20);
        os << "GID" << endl << endl;
      }

      Barrier();

      for (int i = 0 ; i < GetNumMyElements() ; ++i) {
        os.width(10);
        os << GetMyPID();
        os.width(20);
        os << i;
        os.width(20);
        os << (*this)(i) << endl;
      }
    }

    return(os);
  }

  // @}
  
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
  Teuchos::RefCountPtr<Epetra_IntSerialDenseVector> RCPMyGlobalElements_;

};

} // namespace MLAPI

#endif

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

\brief Specifies the number of elements in a space and its layout among processors.

\author Marzio Sala, SNL 9214

\date Last updated on Feb-05.
*/

class Space : public BaseObject {

public:
  //@{ Constructors and Destructors
  //
  //! Default constructor, defines an empty space.
  Space()
  {
    NumMyElements_ = 0;
    NumGlobalElements_ = 0;
    IsLinear_ = false;
    Offset_ = 0;
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
    NumMyElements_     = RHS.NumMyElements();
    NumGlobalElements_ = RHS.NumGlobalElements();
    Offset_            = RHS.Offset();
    IsLinear_          = RHS.IsLinear();
    MyGlobalElements_  = RHS.MyGlobalElements();
  }

  //! Destructor.
  ~Space() {};

  //@}
  //@{ Overloaded operators

  //! Operator =.
  Space& operator=(const Space& RHS)
  {
    NumMyElements_     = RHS.NumMyElements();
    NumGlobalElements_ = RHS.NumGlobalElements();
    Offset_            = RHS.Offset();
    IsLinear_          = RHS.IsLinear();
    MyGlobalElements_  = RHS.MyGlobalElements();
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
    if (NumGlobalElements() != RHS.NumGlobalElements())
      return(false);
    if (NumMyElements() != RHS.NumMyElements())
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
      return(i + Offset_);
    else
      return((*MyGlobalElements_)[i]);
  }

  // @}
  // @{ Reshape methods

  //! Reset the dimension of the space by specifying the local number of elements.
  void Reshape(const int NumGlobalElements, const int NumMyElements = -1)
  {

    if (NumGlobalElements <= 0 && NumMyElements < 0) {
      cerr << "ERROR: In Space::Reshape()" << endl;
      cerr << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
      cerr << "ERROR: wrong input parameters," << endl;
      cerr << "ERROR: NumGlobalElements = " << NumGlobalElements << " and "
           << "ERROR: NumMyElements = " << NumMyElements << endl;
      throw(-1);
    }

    if (NumMyElements == -1) {
      NumMyElements_ = NumGlobalElements / NumProc();
      if (MyPID() == 0)
        NumMyElements_ += NumGlobalElements % NumProc();
    }
    else
      NumMyElements_ = NumMyElements;

    NumGlobalElements_ = ML_Comm_GsumInt(GetML_Comm(),NumMyElements_);

    if (NumGlobalElements != -1) {
      if (NumGlobalElements != NumGlobalElements_) {
        cerr << "ERROR: In Space()::Reshape()" << endl;
        cerr << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
        cerr << "ERROR: Specified number of global elements does not match" << endl;
        cerr << "ERROR: the sum of local elements (" << NumGlobalElements
             << " vs. " << NumGlobalElements_ << ")" << endl;
        throw(-1);
      }
    }

    Offset_   = ML_gpartialsum_int(NumMyElements_,GetML_Comm());
    IsLinear_ = true;
    
  }

  //! Reset the dimension of the space by specifying the local number of elements and their global numbering (starting from 0).
  void Reshape(const int NumGlobalElements, const int NumMyElements, 
               const int* MyGlobalElements)
  {
    if (NumGlobalElements <= 0 && NumMyElements < 0) {
      cerr << "ERROR: In Space::Reshape() (file " << __FILE__
           << ", line " << __LINE__ << ")" << endl;
      cerr << "ERROR: wrong input parameters," << endl;
      cerr << "ERROR: NumGlobalElements = " << NumGlobalElements << " and "
           << "ERROR: NumMyElements = " << NumMyElements << endl;
      throw(-1);
    }

    if (NumMyElements == -1) {
      NumMyElements_ = NumGlobalElements / NumProc();
      if (MyPID() == 0)
        NumMyElements_ += NumGlobalElements % NumProc();
    }
    else
      NumMyElements_ = NumMyElements;

    NumGlobalElements_ = ML_Comm_GsumInt(GetML_Comm(),NumMyElements_);

    if (NumGlobalElements != -1) {
      if (NumGlobalElements != NumGlobalElements_) {
        cerr << "ERROR: In Space() constructor (file " << __FILE__
             << ", line " << __LINE__ << ")" << endl;
        cerr << "ERROR: Specified number of global elements does not match" << endl;
        cerr << "ERROR: the sum of local elements (" << NumGlobalElements
             << " vs. " << NumGlobalElements_ << ")" << endl;
        throw(-1);
      }
    }

    MyGlobalElements_ = Teuchos::rcp(new Epetra_IntSerialDenseVector);
    MyGlobalElements_->Resize(NumMyElements_);
    for (int i = 0 ; i < NumMyElements_ ; ++i)
      (*MyGlobalElements_)[i] = MyGlobalElements[i];

    Offset_ = -1; // not set
    IsLinear_ = false;
  }

  // @}
  // @{ Query methods

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

  //! Prints on ostream basic information about \c this object.
  std::ostream& Print(std::ostream& os,
                      const bool verbose = true) const
  {
    os << std::endl;
    if (MyPID() == 0) {
      os << "*** MLAPI::Space ***" << endl;
      os << "Label               = " << GetLabel() << endl;
      os << "NumMyElements()     = " << NumMyElements() << endl;
      os << "NumGlobalElements() = " << NumGlobalElements() << endl;
      os << "Offset              = " << Offset() << endl;
      if (IsLinear())
        os << "Distribution is linear" << endl;
      else
        os << "Distribution is not linear" << endl;
      os << endl;
    }

    if (verbose) {
      if (MyPID() == 0) {
        os.width(10);
        os << "ProcID";
        os.width(20);
        os << "LID";
        os.width(20);
        os << "GID" << endl << endl;
      }

      GetEpetra_Comm().Barrier();

      for (int i = 0 ; i < NumMyElements() ; ++i) {
        os.width(10);
        os << MyPID();
        os.width(20);
        os << i;
        os.width(20);
        os << (*this)(i) << endl;
      }
    }

    return(os);
  }

private:

  // @}
  // @{ Private data

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

  // @}
};

} // namespace MLAPI

#endif

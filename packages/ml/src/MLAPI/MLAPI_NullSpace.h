#ifndef MLAPI_NULLSPACE_H
#define MLAPI_NULLSPACE_H

#include "Teuchos_RefCountPtr.hpp"
#include "MLAPI_BaseObject.h"
#include "MLAPI_Space.h"
#include "MLAPI_Workspace.h"

using namespace std;

namespace MLAPI {

class NullSpace : public BaseObject {

public:

  // @{ Constructors and destructors

  //! Empty constructor.
  NullSpace()
  {
    Dimension_ = 0;
    RCPValues_    = Teuchos::null;
  }

  //! Builds default null space, gives ownership.
  NullSpace(const Space& VectorSpace, const int Dimension)
  {
    Reshape(VectorSpace, Dimension);
  }

  //! Builds non-default null space, given as one vector.
  NullSpace(const Space& VectorSpace, const int Dimension,
            double* Values, bool ownership = true)
  {
    Reshape(VectorSpace, Dimension, Values, ownership);
  }

  //! Copy constructor.
  NullSpace(const NullSpace& rhs)
  {
    VectorSpace_ = rhs.GetVectorSpace();
    Dimension_   = rhs.GetDimension();
    RCPValues_      = rhs.GetRCPValues();
  }

  // @}
  // @{ Overloaded operators

  //! operator=
  NullSpace& operator=(const NullSpace& rhs) 
  {
    if (this != &rhs) {
      VectorSpace_ = rhs.GetVectorSpace();
      Dimension_   = rhs.GetDimension();
      RCPValues_      = rhs.GetRCPValues();
    }

    return(*this);
  }

  //! Returns the \c row component of the \c null space vector.
  double operator()(const int row, const int col) const
  {
    return(GetValues()[row + col * GetMyLength()]);
  }

  ~NullSpace()
  {
    RCPValues_ = Teuchos::null;
  }

  // @}
  // @{ Reshape methods

  //! Reshapes \c this object with default null space.
  void Reshape(const Space& VectorSpace, const int InputDimension)
  {
    if (InputDimension == 0)
      ML_THROW("Zero null space dimension", -1);

    if (VectorSpace.GetNumGlobalElements() == 0)
      ML_THROW("Empty vector space", -1);

    VectorSpace_ = VectorSpace;
    Dimension_   = InputDimension;
    RCPValues_      = Teuchos::null;
    RCPValues_      = Teuchos::rcp(new double[GetDimension() * GetMyLength()]);

    // populate the null space with default 0's and 1's
    for (int v = 0 ; v < GetDimension() ; ++v) {
      for (int i = 0 ; i < GetMyLength() ; ++i) {
        if ((i + 1) % (v + 1) == 0)
          GetValues()[i + v * GetMyLength()] = 1.0;
        else
          GetValues()[i + v * GetMyLength()] = 0.0;
      }
    }
  }

  //! Reshapes \c this object with non-default null space.
  void Reshape(const Space& VectorSpace, const int Dimension,
               double* Values, bool ownership = true)
  {
    VectorSpace_ = VectorSpace;
    Dimension_   = Dimension;
    RCPValues_      = Teuchos::rcp(Values, ownership);
  }

  // @}
  // @{ Query methods
           
  //! Returns the Space of \c this operator.
  const Space& GetVectorSpace() const 
  {
    return(VectorSpace_);
  }

  //! Returns the local length of the null space vector.
  int GetMyLength() const
  {
    return(VectorSpace_.GetNumMyElements());
  }

  //! Returns the global length of the null space vector.
  int GetGlobalLength() const
  {
    return(VectorSpace_.GetNumGlobalElements());
  }

  //! Returns the dimension of the null space.
  int GetDimension() const 
  {
    return(Dimension_);
  }

  //! Returns a pointer to the internally stores null space vector (non-const).
  double* GetValues()
  {
    return(RCPValues_.get());
  }

  //! Returns a pointer to the internally stores null space vector (const).
  const double* GetValues() const 
  {
    return(RCPValues_.get());
  }

  //! Returns the RefCountPtr for the internally stored null space vector.
  const Teuchos::RefCountPtr<double> GetRCPValues() const 
  {
    return(RCPValues_);
  }

  //! Prints basic information about \c this object.
  virtual std::ostream& Print(std::ostream& os, 
                              const bool verbose = true) const
  {
    os << std::endl;
    if (GetMyPID() == 0) {
      os << "*** MLAPI::NullSpace ***" << endl;
      os << "Label        = " << GetLabel() << endl;
      os << "MyLength     = " << GetMyLength() << endl;
      os << "GlobalLength = " << GetGlobalLength() << endl;
      os << "Dimension    = " << GetDimension() << endl;
      os << endl << endl;
    }

    if (verbose) {

      if (GetMyPID() == 0) {
        os.width(10);
        os << "ProcID";
        os.width(20);
        os << "LID";
        os.width(20);
        os << "GID";
        for (int v = 0 ; v < GetDimension() ; ++v) {
          os.width(20);
          os << "comp " << v;
        }
        os << endl << endl;
      }
      
      for (int iproc = 0 ; iproc < GetNumProcs() ; ++iproc) {

        if (GetMyPID() == iproc) {

          for (int i = 0 ; i < GetMyLength() ; ++i) {
            os.width(10);
            os << GetMyPID();
            os.width(20);
            os << i;
            os.width(20);
            os << GetVectorSpace()(i);
            os.width(20);
            for (int v = 0 ; v < GetDimension() ; ++v) {
              os.width(20);
              os << (*this)(i,v);
            }
            os << endl;
          }
        }

        Barrier();
      }
    }

    return(os);
  }

private:
  // @}
  // @{ Internally used data

  //! Space (defines the local and global lengths).
  Space VectorSpace_;
  //! Dimension of the null space.
  int Dimension_;
  //! Contains the null space pointer (as one vector).
  Teuchos::RefCountPtr<double> RCPValues_;
};

} // namespace MLAPI

#endif // MLAPI_NULLSPACE_H

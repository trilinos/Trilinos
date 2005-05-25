#ifndef AMESOS_TIME_H
#define AMESOS_TIME_H

#include "Epetra_Comm.h"
#include "Epetra_Time.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

using namespace Teuchos;

/*!
 \class Amesos_Time
 
 \brief Amesos_Time: Container for timing information.

 \author Marzio Sala, SNL 9214

 \date Last updated on 24-May-05 (Champions' League Final day)
*/
class Amesos_Time
{
public:

  //! Default constructor.
  Amesos_Time()
  {}

  //! Default destructor.
  ~Amesos_Time()
  {}

  //! Initializes the Time object.
  inline void InitTime(const Epetra_Comm& Comm)
  {
    Time_ = rcp(new Epetra_Time(Comm));
  }

  //! Resets the internally stored time object.
  inline void ResetTime()
  {
    Time_->ResetStartTime();
  }

  //! Adds to field \c what the time elapsed since last call to ResetTime().
  inline void AddTime(const string what)
  {
    data_.set(what, data_.get(what,0.0) + Time_->ElapsedTime());
  }

  //! Gets the comulative time for field \c what.
  inline double GetTime(const string what) const
  {
    return(data_.get(what, 0.0));
  }

private:
  //! Time object.
  RefCountPtr<Epetra_Time> Time_;
  
  //! Container for all fields.
  mutable ParameterList data_;
};

#endif

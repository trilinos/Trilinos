#ifndef AMESOS_TIME_H
#define AMESOS_TIME_H

#include "Epetra_Comm.h"
#include "Epetra_Time.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include <vector>

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

  //! Default constructor to create \c size timers.
  Amesos_Time() :
    size_(1)
  {}

  //! Default destructor.
  ~Amesos_Time()
  {}

  //! Initializes the Time object.
  inline void InitTime(const Epetra_Comm& Comm, int size = 1)
  {
    size_ = size;
    Time_.resize(size_);

    for (int i = 0 ; i < size_ ; ++i)
      Time_[i] = rcp(new Epetra_Time(Comm));
  }

  //! Resets the internally stored time object.
  inline void ResetTime(const int i = 0)
  {
    Time_[i]->ResetStartTime();
  }

  //! Adds to field \c what the time elapsed since last call to ResetTime().
  inline void AddTime(const string what, const int i = 0)
  {
    data_.set(what, data_.get(what,0.0) + Time_[i]->ElapsedTime());
  }

  //! Gets the comulative time for field \c what.
  inline double GetTime(const string what) const
  {
    return(data_.get(what, 0.0));
  }

private:
  //! Number of Epetra_Time objects allocated in \c this object.
  int size_;

  //! Time object.
  std::vector<RefCountPtr<Epetra_Time> > Time_;
  
  //! Container for all fields.
  mutable ParameterList data_;
};

#endif

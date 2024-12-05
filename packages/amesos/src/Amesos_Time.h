#ifndef AMESOS_TIME_H
#define AMESOS_TIME_H

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

#include "Epetra_Comm.h"
#include "Epetra_Time.h"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Array;

/*!
  \struct Amesos_Time_Data
  
  \brief Amesos_Time_Data: Simple struct for storing associated data for Amesos_Time.
  
  \author Heidi Thornquist, SNL 1437
  \date Last updated on 26-March-07
*/
struct Amesos_Time_Data {
  //! Character string identifying this timing data.
  std::string timeName_;
  //! Current timing data.
  double timeVal_;
  
  //! Constructor
  Amesos_Time_Data( std::string timeName, double timeVal ) : 
    timeName_(timeName), 
    timeVal_(timeVal)
  {}
  
  //! Destructor
  virtual ~Amesos_Time_Data() 
  {}
  
};

/*!
  \class Amesos_Time
  
  \brief Amesos_Time: Container for timing information.
  
  \author Marzio Sala, SNL 9214
  
  \date Last updated on 26-March-07 by Heidi Thornquist
*/
class Amesos_Time
{
 public:
  
  //! Default constructor to create \c size timers.
  Amesos_Time() :
    size_(1)
  {}
    
  //! Default destructor.
  virtual ~Amesos_Time()
  {}
    
  //! Initializes the Time object.
  inline void CreateTimer(const Epetra_Comm& Comm, int size = 1)
  {
    size_ = size;
    time_.resize(size_);

    for (int i = 0 ; i < size_ ; ++i)
      time_[i] = rcp(new Epetra_Time(Comm));
  }

  //! Resets the internally stored time object.
  inline void ResetTimer(const int timerID = 0)
  {
    time_[timerID]->ResetStartTime();
  }

  //! Adds to field \c what the time elapsed since last call to ResetTimer().
  inline int AddTime(const std::string what, int dataID, const int timerID = 0)
  {
    // A valid data id is assumed to be > 0, if the id < 0, 
    // then a new entry in the array is created.
    if (dataID < 0) {
      data_.push_back( Amesos_Time_Data( what, time_[timerID]->ElapsedTime() ) );
      return data_.size()-1;
    }
    
    // Check to make sure the data id is valid
    TEUCHOS_TEST_FOR_EXCEPTION(
      timerID >=  (int)(data_.size()), std::logic_error,
      "Amesos_Time::AddTime(...): Error, dataID="<<dataID
      <<" is >= data_.size()="<<data_.size() <<" for dataName=\""<<what<<"\"!"
    );
   
    // The id is valid and the current elapsed time from the indicated timer will be added in. 
    data_[dataID].timeVal_ += time_[timerID]->ElapsedTime();
    return dataID;
  }

  //! Gets the cumulative time using the string.
  inline double GetTime(const std::string what) const
  {
    int dataSize = (int)(data_.size());
    for (int i=0; i<dataSize; ++i) {
      if ( data_[i].timeName_ == what ) {
        return data_[i].timeVal_;
      }
    }
    return 0.0;
  }

  //! Gets the cumulative time using the dataID.
  inline double GetTime(const int dataID) const
  {
    // Return zero if the dataID is not valid
    if ( dataID < 0 || dataID >= (int)(data_.size()) ) {
      return 0.0;
    }
    return data_[dataID].timeVal_;
  }
	
  //! Load up the current timing information into the parameter list.
  inline void GetTiming( Teuchos::ParameterList& list ) const
  {
    int dataSize = (int)(data_.size());
    for (int i=0; i<dataSize; ++i) {
      list.set( data_[i].timeName_, data_[i].timeVal_ );
    }
  }

private:
 
  //! Number of Epetra_Time objects allocated in \c this object.
  int size_;

  //! Time object.
  Array<RCP<Epetra_Time> > time_;
  
  //! Fast accessable container for timing data.
  Array< Amesos_Time_Data > data_;
};

#endif

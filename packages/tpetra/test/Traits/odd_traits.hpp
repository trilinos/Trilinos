/*
Odd traits classes used for testing Traits mechanisms.

Time of Day

*/

#include <iostream>
#include <Teuchos_OrdinalTraits.hpp>
#include "Tpetra_PacketTraits.hpp"

class TimeOfDay {
public:
  
  // constructors and destructors
  TimeOfDay()
    : hour_(0)
    , minute_(0)
  {};
  
  TimeOfDay(int hour, int minute)
    : hour_(hour)
    , minute_(minute)
  {
    if(hourOutOfBounds(hour))
      throw(-1);
    if(minuteOutOfBounds(minute))
      throw(-2);
  };
  
  TimeOfDay(TimeOfDay rhs)
    : hour_(rhs.hour_)
    , minute_(rhs.minute_)
  {};
  
  ~TimeOfDay() {};
  
  // get & set methods
  int getHour() {return(hour_);};
  
  int getMinute() {return(minute_);};
  
  void setHour(int hour) {
    if(hourOutOfBounds(hour))
      throw(-1);
    hour_ = hour;
  };
  
  void setMinute(int minute) {
    if(minuteOutOfBounds(minute))
      throw(-2);
    minute_ = minute;
  };
  
  // overloaded operators
  bool operator<(TimeOfDay const& lhs, TimeOfDay const& rhs) {
    if(lhs.hour_ < rhs.hour_)
      return(true);
    if(lhs.hour_ > rhs.hour_)
      return(false);
    if(lhs.hour_ == rhs.hour_)
      return(lhs.minute_ < rhs.minute_);
  };
  
  bool operator>(TimeOfDay const& lhs, TimeOfDay const& rhs) {
    if(lhs.hour_ > rhs.hour_)
      return(true);
    if(lhs.hour_ < rhs.hour_)
      return(false);
    if(lhs.hour_ == rhs.hour_)
      return(lhs.minute_ > rhs.minute_);
  };

  bool operator==(TimeOfDay const& rhs) {return((hour_ == rhs.hour_) && (minute_ == rhs.minute_));};
  bool operator!=(TimeOfDay const& rhs) {return(!operator==(rhs));};
  bool operator<=(TimeOfDay const& lhs, TimeOfDay const& rhs) {return((lhs < rhs) || (lhs == rhs));};
  bool operator>=(TimeOfDay const& lhs, TimeOfDay const& rhs) {return((lhs > rhs) || (lhs == rhs));};

  TimeOfDay& operator+=(TimeOfDay const& rhs) {
    minute_ += rhs.minute_;
    hour_ += rhs.hour_;
    adjustTime();
    return(this);
  };
  
  TimeOfDay& operator-=(TimeOfDay const& rhs) {
    minute_ -= rhs.minute_;
    hour_ -= rhs.hour_;
    adjustTime();
    return(this);
  };
  
  // NOTE: Very bad compilers may treat this as a const_cast of lhs instead of as an unnamed copy construction.
  // See "More Efficient C++", page 109.
  TimeOfDay& operator+(TimeOfDay const& lhs, TimeOfDay const& rhs) {return(TimeOfDay(lhs) += rhs);};
  TimeOfDay& operator-(TimeOfDay const& lhs, TimeOfDay const& rhs) {return(TimeOfDay(lhs) -= rhs);};

  // overloaded operators that are declared but not defined
  TimeOfDay& operator*=(TimeOfDay const& rhs);
  TimeOfDay& operator/=(TimeOfDay const& rhs);
  TimeOfDay& operator%=(TimeOfDay const& rhs);
  TimeOfDay& operator*(TimeOfDay const& lhs, TimeOfDay const& rhs);
  TimeOfDay& operator/(TimeOfDay const& lhs, TimeOfDay const& rhs);
  TimeOfDay& operator%(TimeOfDay const& lhs, TimeOfDay const& rhs);
  
private:
  int hour_;
  int minute_;
  
  void adjustTime() {
    while(minute_ < 0) {
      minute_ += 60;
      hour_--;
    }
    while(minute_ > 59) {
      minute_ -= 60;
      hour_++;
    }
    while(hour_ < 0)
      hour_ += 24;
    while(hour_ > 23)
      hour_ -= 24;
  };
  
  bool minuteOutOfBounds(int minute) {return(minute < 0 || minute > 59);};
  bool hourOutOfBounds(int hour) {return(hour < 0 || hour > 23);};

};

// overloaded << operator for output
inline std::ostream& operator<<(std::ostream& os, TimeOfDay const& time)
{
  os << time.getHour() << ":" << time.getMinute();
  return(os);
}

// specializations for OrdinalTraits and PacketTraits
namespace Teuchos {
  template<>
  struct OrdinalTraits<TimeOfDay> {
    static inline bool haveMachineParameters() {return(false);};
    static inline TimeOfDay zero()             {return(TimeOfDay(0,0));}; // zero is defined as midnight (0:00)
    static inline TimeOfDay one()              {return(TimeOfDay(0,1));}; // one is defined as one minute past midnight (0:01)
    static inline std::string name()           {return("TimeOfDay");};
  };
}
namespace Tpetra {
  template<>
  struct PacketTraits<TimeOfDay> {
		static inline int packetSize()             {return(sizeof(TimeOfDay));};
    static inline std::string name()           {return("TimeOfDay");};
  };
}

#ifndef TIME_MEMORY_H_
#define TIME_MEMORY_H_
/*
  My hope is that we can measure time and memory usage on multiple
  architectures.  However, some architectures will not support all
  time and memory usage measurements.  Some may support memory 
  measurements only on one processor.  

  In the end there will be two measurement calls (First and Next)
  and a variety of reporting calls.

  Each call to Time_Memory_Object.Time_First() will initialize each of the 
  time and memory fields with a value for the time and memory usage
  on this particular processor.  The time will be from the most recent 
  previous call to .Time_First() from any object.

  Each call to Time_Memory_Object.Time_Next() will update each of the 
  time and memory fields based on the time and memory usage
  on this particular processor since the most recent previous call to 
  .Time_First() from any object.

 */
#include "Epetra_Object.h"
const int UnUsedInt = - 13 ; 
const double UnUsedDbl = -13.0 ; 
class Time_Memory: public Epetra_Object {
    
  public:
  //! Time_Memory Constructor.
  Time_Memory( ) :
    WallTime_(UnUsedDbl) ,
    UserTime_(UnUsedDbl) ,
    SystemTime_(UnUsedDbl) {}  ;

  void Reset_Time();

  void Start_Time();

  void Stop_Time();

  void Time_First();

  void Time_Next();

  double WallTime() const;

  double UserTime() const;

  double SystemTime() const;

  long VirtualMem() const;

  long ResidentMem() const;

  long PageFaults() const;

  virtual ~Time_Memory();

 private:

  double WallStartTime_, UserStartTime_, SystemStartTime_ ;
  double WallTime_, UserTime_, SystemTime_ ;
  long VirtualMem_, ResidentMem_, PageFaults_ ;
  
};

#endif /* TIME_MEMORY_H_ */

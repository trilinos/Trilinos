#include "Time_Memory.h"
#include "TimeMemory.h"

/*
  The Time_Memory class provides two ways of timing execution.  
  1)  Each call to Time_First() or Time_Next() assigns the 
  time since the most recent previous call to one of these
  two functions to the calling object.  This guarantees that
  all time is assigned to exactly one object.
  2)  Reset_Time(), Start_Time() and Stop_Time() allow one to 
  time intervals.  
  
 */

/* 
  Time of the most recent previous call to Time_First (or Time_Next) 
 */
static double WallTime_Last_ = 0.0;
static double UserTime_Last_ = 0.0;
static double SystemTime_Last_ = 0.0;

// Time_Memory::Time_Memory() { }
//  Time_Memory::Time_Memory(const Epetra_Comm* Comm): MyComm(Comm){}
 ; 
Time_Memory::~Time_Memory() { } ; 

void Time_Memory::Time_First( ) {
  double usertime, systemtime, walltime ; 
  long virtualmem, residentmem, pagefaults ;

  getmemusage( &virtualmem, &residentmem, &pagefaults ) ; 
  VirtualMem_ = virtualmem ;
  ResidentMem_ = residentmem ;
  PageFaults_ = pagefaults ; 

  gettimes( &usertime, &systemtime, &walltime ) ; 
  UserTime_ = usertime - UserTime_Last_ ; 
  SystemTime_ = systemtime - SystemTime_Last_ ; 
  WallTime_ = walltime - WallTime_Last_ ; 

  UserTime_Last_ = usertime ; 
  SystemTime_Last_ = systemtime ; 
  WallTime_Last_ = walltime ; 
} 

void Time_Memory::Time_Next(  ) {
  double usertime, systemtime, walltime ; 
  long virtualmem, residentmem, pagefaults ; 

  getmemusage( &virtualmem, &residentmem, &pagefaults ) ; 
  VirtualMem_ = EPETRA_MAX( VirtualMem_, virtualmem) ;
  ResidentMem_ =  EPETRA_MAX( ResidentMem_, residentmem) ;
  PageFaults_ =  EPETRA_MAX( PageFaults_, pagefaults) ; 



  gettimes( &usertime, &systemtime, &walltime ) ; 
  UserTime_ += usertime - UserTime_Last_ ; 
  SystemTime_ += systemtime - SystemTime_Last_ ; 
  WallTime_ += walltime - WallTime_Last_ ; 

  UserTime_Last_ = usertime ; 
  SystemTime_Last_ = systemtime ; 
  WallTime_Last_ = walltime ; 
} 
double Time_Memory::WallTime(  ) const {
  return  WallTime_ ; 
} 
double Time_Memory::UserTime(  ) const {
  return  UserTime_ ; 
}
double Time_Memory::SystemTime(  ) const {
  return  SystemTime_ ; 
}

long Time_Memory::VirtualMem(  ) const {
  return  VirtualMem_ ; 
} 
long Time_Memory::ResidentMem(  ) const {
  return  ResidentMem_ ; 
} 
long Time_Memory::PageFaults(  ) const {
  return  PageFaults_ ; 
} 


void Time_Memory::Reset_Time(  ) {
  WallTime_ = 0.0 ; 
  UserTime_ = 0.0 ; 
  SystemTime_ = 0.0 ; 
  VirtualMem_ = 0 ; 
  PageFaults_ = 0 ;
  ResidentMem_ = 0 ; 

  gettimes( &UserStartTime_, &SystemStartTime_, &WallStartTime_ ) ; 

}
void Time_Memory::Start_Time(  ) {

  gettimes( &UserStartTime_, &SystemStartTime_, &WallStartTime_ ) ; 

}

void Time_Memory::Stop_Time(  ) {
  double usertime, systemtime, walltime ; 
  long virtualmem, residentmem, pagefaults ; 

  getmemusage( &virtualmem, &residentmem, &pagefaults ) ; 
  VirtualMem_ = EPETRA_MAX( VirtualMem_, virtualmem) ;
  ResidentMem_ =  EPETRA_MAX( ResidentMem_, residentmem) ;
  PageFaults_ =  EPETRA_MAX( PageFaults_, pagefaults) ; 



  gettimes( &usertime, &systemtime, &walltime ) ; 
  UserTime_ += usertime - UserStartTime_ ; 
  SystemTime_ += systemtime - SystemStartTime_ ; 
  WallTime_ += walltime - WallStartTime_ ; 

  UserTime_Last_ = usertime ; 
  SystemTime_Last_ = systemtime ; 
  WallTime_Last_ = walltime ; 
} 


#include "Petra_Time.h"


//=============================================================================
Petra_Time::Petra_Time(const Petra_Comm& Comm) 
  : Comm_(&Comm)
{
  StartTime_ = WallTime();
}
//=============================================================================
Petra_Time::Petra_Time(const Petra_Time& Time) 
  : Comm_(Time.Comm_)
{
}
//=============================================================================
Petra_Time::~Petra_Time(void)  
{
}
//=============================================================================
double Petra_Time::WallTime(void) const
{
#ifdef PETRA_MPI

  return(MPI_Wtime());

#elif ICL

   clock_t start;
   double duration;

   start = clock();
  return (double)( start ) / CLOCKS_PER_SEC;

#else

   struct timeval tp;
   static long start=0, startu;
   if (!start)
   {
      gettimeofday(&tp, NULL);
      start = tp.tv_sec;
      startu = tp.tv_usec;
      return(0.0);
   }
   gettimeofday(&tp, NULL);
   return( ((double) (tp.tv_sec - start)) + (tp.tv_usec-startu)/1000000.0 );

#endif

}
//=============================================================================
void Petra_Time::ResetStartTime(void)
{
  StartTime_ = WallTime();
  return;
}
//=============================================================================
double Petra_Time::ElapsedTime(void) const
{
  return(WallTime()-StartTime_);
}

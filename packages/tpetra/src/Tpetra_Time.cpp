namespace Tpetra {
//=============================================================================
Time::Time(const Comm<double, int>& Comm) : Object("Tpetra::Time"), Comm_(&Comm) {
  startTime_ = wallTime();
}

//=============================================================================
Time::Time(const Time& Time) : Object(Time.label()), Comm_(Time.Comm_) {
  startTime_ = Time.startTime_; // We do want to do this in a copy constructor, after all, right?
}

//=============================================================================
double Time::wallTime() const {
  struct timeval tp;
  static long start = 0, startu;
  if (!start)
  {
    gettimeofday(&tp, NULL);
    start = tp.tv_sec;
    startu = tp.tv_usec;
    return(0.0);
  }
  gettimeofday(&tp, NULL);
  return( ((double) (tp.tv_sec - start)) + (tp.tv_usec-startu)/1000000.0 );
}

//=============================================================================
void Time::resetStartTime() {
  startTime_ = wallTime();
}

//=============================================================================
double Time::elapsedTime() const {
  return(wallTime() - startTime_);
}

}  // namespace Tpetra

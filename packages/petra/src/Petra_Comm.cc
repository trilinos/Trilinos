
#include "Petra_Comm.h"


//=============================================================================
#ifdef PETRA_MPI
Petra_Comm::Petra_Comm(MPI_Comm Comm) 
  : Petra_Object("Petra::Comm"), 
    Comm_(Comm),
    ThreadID_(0)
{
    Serial_ = false;
    MPI_Comm_size(Comm, &size_);
    MPI_Comm_rank(Comm, &rank_);
    NodeID_ = rank_;
}
#endif
 
//=============================================================================
Petra_Comm::Petra_Comm(void)
  : Petra_Object("Petra::Comm"), 
    ThreadID_(0)
{
    Serial_ = true;
    size_ = 1;
    rank_ = 0;
    NodeID_ = rank_;
}
 
//=============================================================================
Petra_Comm::Petra_Comm(const Petra_Comm& Comm) : 
  Petra_Object(Comm.Label()),
#ifdef PETRA_MPI
  Comm_(Comm.Comm_),
#endif
  Serial_(Comm.Serial_),
  rank_(Comm.rank_),
  size_(Comm.size_),
  ThreadID_(Comm.ThreadID_),
  NodeID_(Comm.NodeID_)
{
}

//=============================================================================
void Petra_Comm::Barrier(void) const {

#ifdef PETRA_MPI
    if (!Serial_) MPI_Barrier(Comm_);
#endif
}
//=============================================================================
void Petra_Comm::NodeBarrier(void) const {

#ifdef PETRA_MPI
    if (!Serial_) MPI_Barrier(Comm_);
#endif
}
//=============================================================================
int Petra_Comm::Broadcast(double * Values, int Count, int Root) const {

  if (Serial_) return(0);
#ifdef PETRA_MPI
    else PETRA_CHK_ERR(MPI_Bcast(Values, Count, MPI_DOUBLE, Root, Comm_));
#endif
 
  return(0);
}
//=============================================================================
int Petra_Comm::Broadcast(int * Values, int Count, int Root) const {

  if (Serial_) return(0);
#ifdef PETRA_MPI
    else PETRA_CHK_ERR(MPI_Bcast(Values, Count, MPI_INT, Root, Comm_));
#endif
 
  return(0);
}
//=============================================================================
int Petra_Comm::GatherAll(double * MyVals, double * AllVals, int Count) const {

  if (Serial_) 
      for (int i=0; i<Count; i++) AllVals[i] = MyVals[i];
#ifdef PETRA_MPI
    else PETRA_CHK_ERR(MPI_Allgather(MyVals, Count, MPI_DOUBLE, AllVals, Count, MPI_DOUBLE, Comm_));
#endif
 
  return(0);
}
//=============================================================================
int Petra_Comm::GatherAll(int * MyVals, int * AllVals, int Count) const {

  if (Serial_) 
      for (int i=0; i<Count; i++) AllVals[i] = MyVals[i];
#ifdef PETRA_MPI
    else PETRA_CHK_ERR(MPI_Allgather(MyVals, Count, MPI_INT, AllVals, Count, MPI_INT, Comm_));
#endif
 
  return(0);
}
//=============================================================================
int Petra_Comm::SumAll(double * PartialSums, double * GlobalSums, int Count) const {

  if (Serial_) 
      for (int i=0; i<Count; i++) GlobalSums[i] = PartialSums[i];
#ifdef PETRA_MPI
    else PETRA_CHK_ERR(MPI_Allreduce(PartialSums, GlobalSums, Count, MPI_DOUBLE, MPI_SUM, Comm_));
#endif
 
  return(0);
}
//=============================================================================
int Petra_Comm::SumAll(int * PartialSums, int * GlobalSums, int Count) const {

  if (Serial_) 
      for (int i=0; i<Count; i++) GlobalSums[i] = PartialSums[i];
#ifdef PETRA_MPI
    else PETRA_CHK_ERR(MPI_Allreduce(PartialSums, GlobalSums, Count, MPI_INT, MPI_SUM, Comm_));
#endif
 
  return(0);
}
//=============================================================================
int Petra_Comm::MaxAll(double * PartialMaxs, double * GlobalMaxs, int Count) const {

  if (Serial_) 
      for (int i=0; i<Count; i++) GlobalMaxs[i] = PartialMaxs[i];
#ifdef PETRA_MPI
    else PETRA_CHK_ERR(MPI_Allreduce(PartialMaxs, GlobalMaxs, Count, MPI_DOUBLE, MPI_MAX, Comm_));
#endif
 
  return(0);
}
//=============================================================================
int Petra_Comm::MaxAll(int * PartialMaxs, int * GlobalMaxs, int Count) const {

  if (Serial_) 
      for (int i=0; i<Count; i++) GlobalMaxs[i] = PartialMaxs[i];
#ifdef PETRA_MPI
    else PETRA_CHK_ERR(MPI_Allreduce(PartialMaxs, GlobalMaxs, Count, MPI_INT, MPI_MAX, Comm_));
#endif
 
  return(0);
}
//=============================================================================
int Petra_Comm::MinAll(double * PartialMins, double * GlobalMins, int Count) const {

  if (Serial_) 
      for (int i=0; i<Count; i++) GlobalMins[i] = PartialMins[i];
#ifdef PETRA_MPI
    else PETRA_CHK_ERR(MPI_Allreduce(PartialMins, GlobalMins, Count, MPI_DOUBLE, MPI_MIN, Comm_));
#endif
 
  return(0);
}
//=============================================================================
int Petra_Comm::MinAll(int * PartialMins, int * GlobalMins, int Count) const {

  if (Serial_) 
      for (int i=0; i<Count; i++) GlobalMins[i] = PartialMins[i];
#ifdef PETRA_MPI
    else PETRA_CHK_ERR(MPI_Allreduce(PartialMins, GlobalMins, Count, MPI_INT, MPI_MIN, Comm_));
#endif
 
  return(0);
}
//=============================================================================
int Petra_Comm::ScanSum(double * MyVals, double * ScanSums, int Count) const {

  if (Serial_) 
      for (int i=0; i<Count; i++) ScanSums[i] = MyVals[i];
#ifdef PETRA_MPI
    else PETRA_CHK_ERR(MPI_Scan(MyVals, ScanSums, Count, MPI_DOUBLE, MPI_SUM, Comm_));
#endif
 
  return(0);
}
//=============================================================================
int Petra_Comm::ScanSum(int * MyVals, int * ScanSums, int Count) const {

  if (Serial_) 
      for (int i=0; i<Count; i++) ScanSums[i] = MyVals[i];
#ifdef PETRA_MPI
    else PETRA_CHK_ERR(MPI_Scan(MyVals, ScanSums, Count, MPI_INT, MPI_SUM, Comm_));
#endif
 
  return(0);
}
//=============================================================================
Petra_Comm::~Petra_Comm(void)  {

  //delete [] proc_config_;
  //proc_config_ = NULL;
}
//=============================================================================
void Petra_Comm::Print(ostream & os) const
{
  os << "::Processor "<< MyPID()<<" of " << NumProc() << " total processors"; 
  return;
}

// @HEADER BEGIN
// @HEADER END

#include "Teuchos_StackedTimer.hpp"
#include <limits>


namespace Teuchos {

Teuchos::RCP<StackedTimer> StackedTimer::timer;
  
StackedTimer::LevelTimer::LevelTimer() :
    level_(std::numeric_limits<unsigned>::max()),name_("INVALID"),parent_(NULL)
{}

void error_out(const std::string& msg, const bool)
{
  TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,msg);
}
  
#if 0
// Here for MPI stuff, not working/tested
void
StackedTimer::LevelTimer::merge(LevelTimer &partner){

  // Make the timers at this level match
  for (unsigned i=0;i<sub_timers_.size();++i) {
    if ( i == partner.sub_timers_.size())  //ran out of partner ones
      partner.sub_timers_.push_back(LevelTimer(sub_timers_[i].level_,sub_timers_[i].name_.c_str()));
    if (sub_timers_[i].name_.compare( partner.sub_timers_[i].name_)==0)
      continue;
    // So, we need to merge something here.  first see if it is just out of order.
    for (unsigned j = i+1;j<partner.sub_timers_.size();++j)
      if (sub_timers_[i].name_.compare(partner.sub_timers_[j].name_)==0) {
        std::swap(partner.sub_timers_[i],partner.sub_timers_[j]);
        //  partner.sub_timers[i].parent = partner.sub_timers[j].parent = &partner;
        break;
      }
    // If they are not equal still, insert copy
    if (sub_timers_[i].name_.compare(partner.sub_timers_[i].name_) != 0) {
      std::vector< LevelTimer>::iterator iter= partner.sub_timers_.begin();
      for (unsigned j = 0; j < i; ++j) ++iter;
      partner.sub_timers_.insert(iter,LevelTimer(sub_timers_[i].level_,sub_timers_[i].name_.c_str()));
    }
  }
  // Now just tack on the end of the partner list to the main list
  for (unsigned i=sub_timers_.size();i < partner.sub_timers_.size();++i)
  {
    sub_timers_.push_back(LevelTimer(partner.sub_timers_[i].level_,partner.sub_timers_[i].name_.c_str()));
    sub_timers_[i].parent_ = this;
  }

  assert(sub_timers_.size() == partner.sub_timers_.size());

  for (unsigned i=0;i<sub_timers_.size();++i)
    assert ( sub_timers_[i].name_ == partner.sub_timers_[i].name_ ) ;

  for (unsigned i=0;i<sub_timers_.size();++i)
    sub_timers_[i].merge(partner.sub_timers_[i]);

  for (unsigned i=0;i<sub_timers_.size();++i) {
    partner.sub_timers_[i].parent_ = &partner;
  }

}
#endif

void
StackedTimer::LevelTimer::pack() {
  std::vector<char> buff;
  pack_level(buff);

  // wait for rank 0 to send ready to recv
  // MPI_Status stat;
  // int ready_to_recv;
  // MPI_Recv(&ready_to_recv,1,MPI_INT,0, ready_to_recv_tag, mpi_world_,&stat);

  // int size = buff.size();
  // MPI_Send(&size,1,MPI_INT,0,send_size_tag,mpi_world_);
  // MPI_Send(&buff[0],size,MPI_CHAR,0,send_buffer_tag,mpi_world_);
}

StackedTimer::LevelTimer *
StackedTimer::LevelTimer::unpack(unsigned from) {
  // MPI_Status stat;
  // unsigned recv_size;
  // MPI_Recv(&recv_size,1,MPI_INT,from,send_size_tag,mpi_world_,&stat);
  // std::vector<char> buff(recv_size);
  // MPI_Recv(&buff[0], recv_size, MPI_CHAR, from, send_buffer_tag, mpi_world_,&stat);
  // level_ = 0;
  // int position=0;
  // unpack_level(buff, position);
  return this;
}


void
StackedTimer::LevelTimer::unpack_level(std::vector<char> &buff, int &position) {
  // mpi_unpack(name_, buff, position);
  // mpi_unpack(accumulation_, buff, position);
  // mpi_unpack(count_started_, buff, position);
  // mpi_unpack(count_updates_, buff, position);

  // // Sub timer stuff
  // unsigned num_sub_timers;
  // mpi_unpack(num_sub_timers, buff, position);
  // sub_timers_.resize(num_sub_timers);
  // for (unsigned i=0;i<num_sub_timers;++i) {
  //   sub_timers_[i].level_ = level_+1;
  //   sub_timers_[i].unpack_level(buff, position);
  //   sub_timers_[i].parent_ = this;
  // }
}

void
StackedTimer::LevelTimer::pack_level(std::vector<char> &buff) {
  // mpi_pack(name_, buff);
  // mpi_pack(accumulation_, buff);
  // mpi_pack(count_started_, buff);
  // mpi_pack(count_updates_, buff);

  // // Sub timer stuff
  // unsigned num_sub_timers = sub_timers_.size();
  // mpi_pack(num_sub_timers, buff);
  // for (unsigned i=0;i<num_sub_timers;++i)
  //   sub_timers_[i].pack_level(buff);
}

void
StackedTimer::LevelTimer::report(std::ostream &os) {
  for (unsigned i=0; i<level_; ++i)
    os << "   ";
  os << name_<<":"<<accumulatedTime()<< " [" << count_started_<<"] ("<< count_updates_ <<")"<<std::endl;
  double t_total = 0;
  for (size_t i=0; i<sub_timers_.size(); ++i) {
    t_total += sub_timers_[i].accumulatedTime();
    sub_timers_[i].report(os);
  }
  if ( sub_timers_.size() == 0 )
    return;
  for (unsigned i=0; i<=level_; ++i)
    os << "   ";
  os << "Remainder: " << accumulatedTime() - t_total<<std::endl;

}

} //namespace Teuchos

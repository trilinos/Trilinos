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

BaseTimer::TimeInfo
StackedTimer::LevelTimer::findTimer(std::string &name) {
  BaseTimer::TimeInfo t;
  if (get_full_name() == name)
    t = BaseTimer::TimeInfo(this);
  else {
    for (unsigned i=0;i<sub_timers_.size(); ++i){
      t = sub_timers_[i].findTimer(name);
      if (t.time > 0.)
        return t;
    }
  }
  return t;
}
void
StackedTimer::flatten() {
  int num_timers = top_->countTimers();
  flat_names_.resize(num_timers);
  unsigned pos=0;
  top_->addTimerNames(flat_names_, pos);
}

void
StackedTimer::merge(Teuchos::RCP<Teuchos::Comm<int> > comm){
  Array<std::string> all_names;
  mergeCounterNames(*comm, flat_names_, all_names, Union);
  flat_names_ = all_names;
}

void
StackedTimer::collectRemoteData(Teuchos::RCP<Teuchos::Comm<int> > comm) {
  int comm_size = size(*comm);
  int comm_rank = rank(*comm);
  Array<BaseTimer::TimeInfo> local_time_data;
  if (comm_rank == 0) {
    time_data_.resize(comm_size);
    for (int r=0; r<comm_size; ++r)
      time_data_[r].resize(flat_names_.size());
  } else {
    local_time_data.resize(flat_names_.size());
  }

  for (unsigned i=0; i<flat_names_.size(); ++i) {
    auto t = top_->findTimer(flat_names_[i]);
    if ( comm_rank == 0 )
      time_data_[0][i] = t;
    else
      local_time_data[i] = t;
  }

  int buffer_size = sizeof(BaseTimer::TimeInfo);
  if ( comm_rank == 0 )
    buffer_size += ((std::size_t)&time_data_[0][flat_names_.size()-1] - (std::size_t)&time_data_[0][0]);
  else
    buffer_size += ((std::size_t)&local_time_data[flat_names_.size()-1] - (std::size_t)&local_time_data[0]);

  if (comm_rank == 0 )  {
    for (int r=1; r<comm_size; ++r)
      receive(*comm, r, buffer_size, (char*)&time_data_[r][0]);
  } else
    send(*comm, buffer_size,(char*)&local_time_data[0] ,0);
}

std::pair<std::string, std::string> getPrefix(std::string &name) {
  for (std::size_t i=name.size()-1; i>0; --i)
    if (name[i] == ':') {
        return std::pair<std::string, std::string>(name.substr(0,i), name.substr(i+1, name.size()));
    }
  return std::pair<std::string, std::string>(std::string(""), name);
}


double
StackedTimer::printLevel (std::string prefix, int print_level, std::ostream &os, std::vector<bool> &printed) {
  double total_time = 0.0;
  std::size_t num_entries = time_data_.size();
  for (std::size_t i=0; i<flat_names_.size(); ++i ) {
    if (printed[i])
      continue;
    int level = std::count(flat_names_[i].begin(), flat_names_[i].end(), ':');
    if (level != print_level)
      continue;
    auto split_names = getPrefix(flat_names_[i]);
    if ( prefix == split_names.first) {
      double average_t=0;
      long average_count=0;
      long long average_updates=0;
      for (std::size_t r=0; r<num_entries; ++r) {
        average_t += time_data_[r][i].time;
        average_updates += time_data_[r][i].updates;
        average_count += time_data_[r][i].count;
      }
      for (int l=0; l<level; ++l)
        os << "    ";
      os << split_names.second << ": ";
      os << average_t/num_entries <<
          " ["<<average_count/num_entries<<"] ("<<average_updates/num_entries<<")\n";
      printed[i] = true;
      double sub_time = printLevel(flat_names_[i], level+1, os, printed);
      if (sub_time > 0 ) {
        for (int l=0; l<=level; ++l)
          os << "    ";
        os << "Remainder: " << average_t/num_entries - sub_time<<std::endl;
      }
      total_time += average_t/num_entries;
    }
  }
  return total_time;
}

void
StackedTimer::report(std::ostream &os, Teuchos::RCP<Teuchos::Comm<int> > comm) {
  flatten();
  merge(comm);
  collectRemoteData(comm);
  if (rank(*comm) == 0 ) {
    std::vector<bool> printed(flat_names_.size(), false);
    printLevel("", 0, os, printed);
  }
}

} //namespace Teuchos

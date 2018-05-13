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
StackedTimer::merge(Teuchos::RCP<const Teuchos::Comm<int> > comm){
  Array<std::string> all_names;
  mergeCounterNames(*comm, flat_names_, all_names, Union);
  flat_names_ = all_names;
}

void
StackedTimer::collectRemoteData(Teuchos::RCP<const Teuchos::Comm<int> > comm) {
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

  int buffer_size = sizeof(BaseTimer::TimeInfo[flat_names_.size()]);

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
StackedTimer::printLevel (std::string prefix, int print_level, std::ostream &os, std::vector<bool> &printed, double parent_time, const OutputOptions &options) {
  double total_time = 0.0;
  std::size_t num_entries = time_data_.size();
  for (std::size_t i=0; i<flat_names_.size(); ++i ) {
    if (printed[i])
      continue;
    int level = std::count(flat_names_[i].begin(), flat_names_[i].end(), ':');
    if (level != print_level)
      continue;
    auto split_names = getPrefix(flat_names_[i]);
    if ( prefix != split_names.first)
      continue;

    int num_ranks=0;
    double sum_t=0., sum_sq=0., min_t=std::numeric_limits<double>::max(), max_t=0.;
    long sum_count=0;
    long long sum_updates=0;
    for (std::size_t r=0; r<num_entries; ++r) {
      if (time_data_[r][i].count) {
        num_ranks++;
        min_t = std::min(min_t, time_data_[r][i].time);
        max_t = std::max(max_t, time_data_[r][i].time);
        sum_t += time_data_[r][i].time;
        sum_sq += time_data_[r][i].time*time_data_[r][i].time;
        sum_updates += time_data_[r][i].updates;
        sum_count += time_data_[r][i].count;
      }
    }
    for (int l=0; l<level; ++l)
      os << "    ";
    os << split_names.second << ": ";
    // output averge time
    os << sum_t/num_ranks;
    // output percentage
    if ( options.output_fraction && parent_time>0)
      os << " - "<<sum_t/num_ranks/parent_time*100<<"%";
    // output count
    os << " ["<<sum_count/num_ranks<<"]";
    // output total counts
    if ( options.output_total_updates )
      os << " ("<<sum_updates/num_ranks<<")";
    // Output min and maxs
    if ( options.output_minmax ) {
      os << " {min="<<min_t<<", max="<<max_t;
      if (num_ranks>1)
        os<<", std dev="<<sqrt((sum_sq-sum_t*sum_t/num_ranks)/(num_ranks-1));
      os << "}";
    }
    // Output histogram
    if ( options.output_histogram && num_entries >1 ) {
      std::vector<int> hist(options.num_histogram, 0);
      double dh = (max_t-min_t)/options.num_histogram;
      if (dh==0) // Put everything into bin 1
        dh=1;
      for (std::size_t r=0; r<num_entries; ++r) {
        if (time_data_[r][i].count) {
          int bin=(time_data_[r][i].time - min_t)/dh;
          bin = std::max(std::min(bin,options.num_histogram-1) , 0);
          hist[bin]++;
        }
      }
      // dump the histogram
      os << " <";
      for (int h=0;h<options.num_histogram; ++h) {
        if (h)
          os <<", "<<hist[h];
        else
          os << hist[h];
      }
      os << ">";
    }
    os << std::endl;
    printed[i] = true;
    double sub_time = printLevel(flat_names_[i], level+1, os, printed, sum_t/num_ranks, options);
    if (sub_time > 0 ) {
      for (int l=0; l<=level; ++l)
        os << "    ";
      os << "Remainder: " << sum_t/num_ranks - sub_time<<std::endl;
    }
    total_time += sum_t/num_ranks;
  }
  return total_time;
}

void
StackedTimer::report(std::ostream &os, Teuchos::RCP<const Teuchos::Comm<int> > comm, OutputOptions options) {
  flatten();
  merge(comm);
  collectRemoteData(comm);
  if (rank(*comm) == 0 ) {
    std::vector<bool> printed(flat_names_.size(), false);
    printLevel("", 0, os, printed, 0., options);
  }
}

} //namespace Teuchos

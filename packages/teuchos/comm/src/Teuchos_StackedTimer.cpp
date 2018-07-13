// @HEADER BEGIN
// @HEADER END

#include "Teuchos_StackedTimer.hpp"
#include <limits>


namespace Teuchos {


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
StackedTimer::LevelTimer::findTimer(const std::string &name, bool& found) {
  BaseTimer::TimeInfo t;
  if (get_full_name() == name) {
    t = BaseTimer::TimeInfo(this);
    found = true;
  }
  else {
    for (unsigned i=0;i<sub_timers_.size(); ++i){
      t = sub_timers_[i].findTimer(name,found);
      if (found)
        return t;
    }
  }
  return t;
}

void
StackedTimer::flatten() {
  int num_timers = timer_.countTimers();
  flat_names_.resize(num_timers);
  unsigned pos=0;
  timer_.addTimerNames(flat_names_, pos);
}

void
StackedTimer::merge(Teuchos::RCP<const Teuchos::Comm<int> > comm){
  Array<std::string> all_names;
  mergeCounterNames(*comm, flat_names_, all_names, Union);
  flat_names_ = all_names;
}

void
StackedTimer::collectRemoteData(Teuchos::RCP<const Teuchos::Comm<int> > comm, const OutputOptions &options) {
  // allocate everything
  int num_names = flat_names_.size();
  sum_.resize(num_names);
  count_.resize(num_names);
  active_.resize(num_names);

  if (options.output_minmax || options.output_histogram) {
    min_.resize(num_names);
    max_.resize(num_names);
    if ( options.output_minmax )
      sum_sq_.resize(num_names);
  }


  if (options.output_histogram ) {
    hist_.resize(options.num_histogram);
    for (int i=0;i<options.num_histogram ; ++i)
      hist_[i].resize(num_names);
  }

  // Temp data
  Array<double> time(num_names);
  Array<unsigned long> count(num_names);
  Array<unsigned long long> updates;
  if (options.output_total_updates)
    updates.resize(num_names);
  Array<int> used(num_names);
  Array<int> bins;

  if (options.output_histogram)
    bins.resize(num_names);

  // set initial values
  for (int i=0;i<num_names; ++i) {
    bool found = false; // ignore result here
    auto t = timer_.findTimer(flat_names_[i],found);
    time[i] = t.time;
    count[i] = t.count;
    used[i] = t.count==0? 0:1;
    if (options.output_total_updates)
      updates[i] = t.updates;
  }

  // Now reduce the data
  reduce<int, double>(time.getRawPtr(), sum_.getRawPtr(), num_names, REDUCE_SUM, 0, *comm);
  reduce(count.getRawPtr(), count_.getRawPtr(), num_names, REDUCE_SUM, 0, *comm);
  reduce(used.getRawPtr(), active_.getRawPtr(), num_names, REDUCE_SUM, 0, *comm);

  if (min_.size()) {
    reduceAll(*comm, REDUCE_MIN, num_names, time.getRawPtr(), min_.getRawPtr());
    reduceAll(*comm, REDUCE_MAX, num_names, time.getRawPtr(), max_.getRawPtr());
  }

  if (options.output_histogram) {
    for (int i=0;i<num_names; ++i) {

      double dh = (max_[i]-min_[i])/options.num_histogram;
      if (dh==0) // Put everything into bin 1
        dh=1;
      if (used[i]) {
        int bin=(time[i]- min_[i])/dh;
        bins[i] = std::max(std::min(bin,options.num_histogram-1) , 0);
      }
    }
    // Recycle the used array for the temp bin array
    for (int j=0; j<options.num_histogram; ++j){
      for (int i=0;i<num_names; ++i) {
        if (bins[i] == j )
          used[i]=1;
        else
          used[i]=0;
      }
      reduce(used.getRawPtr(), hist_[j].getRawPtr(), num_names, REDUCE_SUM, 0, *comm);
    }
  }

  if (sum_sq_.size()) {
    for (int i=0;i<num_names; ++i)
      time[i] *= time[i];
    reduce(time.getRawPtr(), sum_sq_.getRawPtr(), num_names, REDUCE_SUM, 0, *comm);
  }

}

std::pair<std::string, std::string> getPrefix(const std::string &name) {
  for (std::size_t i=name.size()-1; i>0; --i)
    if (name[i] == '@') {
      return std::pair<std::string, std::string>(name.substr(0,i), name.substr(i+1, name.size()));
    }
  return std::pair<std::string, std::string>(std::string(""), name);
}


double
StackedTimer::printLevel (std::string prefix, int print_level, std::ostream &os, std::vector<bool> &printed, double parent_time, const OutputOptions &options) {
  double total_time = 0.0;

  for (int i=0; i<flat_names_.size(); ++i ) {
    if (printed[i])
      continue;
    int level = std::count(flat_names_[i].begin(), flat_names_[i].end(), '@');
    if ( (level != print_level) || (level >= options.max_levels) )
      continue;
    auto split_names = getPrefix(flat_names_[i]);
    if ( prefix != split_names.first)
      continue;

    // Output the data
    for (int l=0; l<level; ++l)
      os << "|   ";
    os << split_names.second << ": ";
    // output averge time
    os << sum_[i]/active_[i];
    // output percentage
    if ( options.output_fraction && parent_time>0)
      os << " - "<<sum_[i]/active_[i]/parent_time*100<<"%";
    // output count
    os << " ["<<count_[i]/active_[i]<<"]";
    // output total counts
    if ( options.output_total_updates )
      os << " ("<<updates_[i]/active_[i]<<")";
    // Output min and maxs
    if ( options.output_minmax && active_[i]>1) {
      os << " {min="<<min_[i]<<", max="<<max_[i];
      if (active_[i]>1)
        os<<", std dev="<<sqrt((sum_sq_[i]-sum_[i]*sum_[i]/active_[i])/(active_[i]-1));
      os << "}";
    }
    // Output histogram
    if ( options.output_histogram && active_[i] >1 ) {
      // dump the histogram
      os << " <";
      for (int h=0;h<options.num_histogram; ++h) {
        if (h)
          os <<", "<<hist_[h][i];
        else
          os << hist_[h][i];
      }
      os << ">";
    }
    os << std::endl;
    printed[i] = true;
    double sub_time = printLevel(flat_names_[i], level+1, os, printed, sum_[i]/active_[i], options);
    if (sub_time > 0 ) {
      for (int l=0; l<=level; ++l)
        os << "|   ";
      os << "Remainder: " <<  sum_[i]/active_[i]- sub_time;
      if ( options.output_fraction && (sum_[i]/active_[i] > 0.) )
        os << " - "<< (sum_[i]/active_[i]- sub_time)/(sum_[i]/active_[i])*100 << "%";
      os <<std::endl;
    }
    total_time += sum_[i]/active_[i];
  }
  return total_time;
}

void
StackedTimer::report(std::ostream &os, Teuchos::RCP<const Teuchos::Comm<int> > comm, OutputOptions options) {
  flatten();
  merge(comm);
  collectRemoteData(comm, options);
  if (rank(*comm) == 0 ) {
    if (options.print_warnings) {
      os << "*** Teuchos::StackedTimer::report() - Remainder for a block will be ***"
         << "\n*** incorrect if a timer in the block does not exist on every rank  ***"
         << "\n*** of the MPI Communicator.                                        ***"
         << std::endl;
    }
    if ( (options.max_levels != INT_MAX) && options.print_warnings) {
      os << "Teuchos::StackedTimer::report() - max_levels set to " << options.max_levels
         << ", to print more levels, increase value of OutputOptions::max_levels." << std::endl;
    }
    std::vector<bool> printed(flat_names_.size(), false);
    printLevel("", 0, os, printed, 0., options);
  }
}

} //namespace Teuchos

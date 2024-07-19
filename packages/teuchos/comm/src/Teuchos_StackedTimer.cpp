// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_StackedTimer.hpp"
#include <limits>
#include <ctime>
#include <cctype>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <sstream>

namespace Teuchos {


StackedTimer::LevelTimer::LevelTimer() :
    level_(std::numeric_limits<unsigned>::max()),name_("INVALID"),parent_(nullptr)
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

const BaseTimer*
StackedTimer::LevelTimer::findBaseTimer(const std::string &name) const {
  const BaseTimer* t = nullptr;
  if (get_full_name() == name) {
    return this;
  }
  else {
    for (unsigned i=0;i<sub_timers_.size(); ++i){
      t = sub_timers_[i].findBaseTimer(name);
      if (t != nullptr)
        return t;
    }
  }
  return t;
}
  
BaseTimer::TimeInfo
StackedTimer::LevelTimer::findTimer(const std::string &name, bool& found) {
  BaseTimer::TimeInfo t;
  auto full_name = get_full_name();
  if (full_name.size() > name.size())
    return t;
  if ( strncmp(full_name.c_str(), name.c_str(), full_name.size()))
    return t;
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
  updates_.resize(num_names);
  active_.resize(num_names);

  if (options.output_minmax || options.output_histogram || options.output_proc_minmax) {
    min_.resize(num_names);
    max_.resize(num_names);
    if ( options.output_minmax )
      sum_sq_.resize(num_names);
    else
      sum_sq_.resize(0);
  } else {
    min_.resize(0);
    max_.resize(0);
    sum_sq_.resize(0);
  }

  if (options.output_proc_minmax) {
    procmin_.resize(num_names);
    procmax_.resize(num_names);
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
    reduceAll(*comm, REDUCE_MAX, num_names, time.getRawPtr(), max_.getRawPtr());
    for (int i=0;i<num_names;++i)
      if (!used[i])
        time[i] = max_[i];
    reduceAll(*comm, REDUCE_MIN, num_names, time.getRawPtr(), min_.getRawPtr());
    for (int i=0;i<num_names;++i)
      if (!used[i])
        time[i] = 0.;
    if (procmin_.size()) {
      Array<int> procmin(num_names);
      Array<int> procmax(num_names);
      int commRank = comm->getRank();
      for (int i=0;i<num_names; ++i) {
        if (used[i] && (min_[i]==time[i]))
          procmin[i] = commRank;
        else
          procmin[i] = -1;
        if (used[i] && (max_[i]==time[i]))
          procmax[i] = commRank;
        else
          procmax[i] = -1;
      }
      reduceAll(*comm, REDUCE_MAX, num_names, procmin.getRawPtr(), procmin_.getRawPtr());
      reduceAll(*comm, REDUCE_MAX, num_names, procmax.getRawPtr(), procmax_.getRawPtr());
    }
  }

  if (options.output_histogram) {
    for (int i=0;i<num_names; ++i) {

      double dh = (max_[i]-min_[i])/options.num_histogram;
      if (dh==0) // Put everything into bin 1
        dh=1;
      if (used[i]) {
        int bin=(time[i]- min_[i])/dh;
        bins[i] = std::max(std::min(bin,options.num_histogram-1) , 0);
      } else
        bins[i] = -1;
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
      return std::pair<std::string, std::string>(name.substr(0,i), name.substr(i+1));
    }
  return std::pair<std::string, std::string>(std::string(""), name);
}

double
StackedTimer::computeColumnWidthsForAligment(std::string prefix,
                                             int print_level,
                                             std::vector<bool> &printed,
                                             double parent_time,
                                             const OutputOptions &options)
{
  // This replicates printLevel but counts column width instead of
  // printing to ostream. This must be kept in sync with printLevel()
  double total_time = 0.0;

  for (int i=0; i<flat_names_.size(); ++i ) {
    if (sum_[i]/active_[i] <= options.drop_time)
      continue;
    if (printed[i])
      continue;
    int level = std::count(flat_names_[i].begin(), flat_names_[i].end(), '@');
    if ( (level != print_level) || (level >= options.max_levels) )
      continue;
    auto split_names = getPrefix(flat_names_[i]);
    if ( prefix != split_names.first)
      continue;

    // Output the indentation level and timer name
    {
      std::ostringstream os;
      for (int l=0; l<level; ++l)
        os << "|   ";
      // Output the timer name
      os << split_names.second << ": ";
      alignments_.timer_names_= std::max(alignments_.timer_names_,os.str().size());
    }

    // output averge time
    {
      std::ostringstream os;
      os << sum_[i]/active_[i];
      alignments_.average_time_ = std::max(alignments_.average_time_,os.str().size());
    }

    // output percentage
    if ( options.output_fraction && parent_time>0) {
      std::ostringstream os;
      os << " - "<<sum_[i]/active_[i]/parent_time*100<<"%";
      alignments_.fraction_ = std::max(alignments_.fraction_,os.str().size());
    }

    // output count
    {
      std::ostringstream os;
      os << " ["<<count_[i]/active_[i]<<"]";
      alignments_.count_ = std::max(alignments_.count_,os.str().size());
    }

    // output total counts
    if ( options.output_total_updates) {
      std::ostringstream os;
      os << " ("<<updates_[i]/active_[i]<<")";
      alignments_.total_updates_ = std::max(alignments_.total_updates_,os.str().size());
    }

    // Output min and maxs
    if ( options.output_minmax && active_[i]>1) {
      {
        std::ostringstream os;
        os << " {min=" << min_[i];
        alignments_.min_ = std::max(alignments_.min_,os.str().size());
      }
      {
        std::ostringstream os;
        os << ", max=" << max_[i];
        if (active_[i] <= 1)
          os << "}";
        alignments_.max_ = std::max(alignments_.max_,os.str().size());
      }
      if (procmin_.size()) {
        std::ostringstream os;
        os << ", proc min=" << procmin_[i];
        if (active_[i] <= 1)
          os << "}";
        alignments_.procmin_ = std::min(alignments_.procmin_,os.str().size());
      }
      if (procmax_.size()) {
        std::ostringstream os;
        os << ", proc max=" << procmax_[i];
        if (active_[i] <= 1)
          os << "}";
        alignments_.procmax_ = std::max(alignments_.procmax_,os.str().size());
      }
      if (active_[i]>1) {
        std::ostringstream os;
        os << ", std dev=" << sqrt(std::max<double>(sum_sq_[i]-sum_[i]*sum_[i]/active_[i],0.0)/(active_[i]-1));
        os << "}";
        alignments_.stddev_ = std::max(alignments_.stddev_,os.str().size());
      }
    }
    // Output histogram
    if ( options.output_histogram && active_[i] >1 ) {
      std::ostringstream os;
      os << " <";
      for (int h=0;h<options.num_histogram; ++h) {
        if (h)
          os <<", "<<hist_[h][i];
        else
          os << hist_[h][i];
      }
      os << ">";
      alignments_.histogram_ = std::max(alignments_.histogram_,os.str().size());
    }

    printed[i] = true;
    double sub_time = computeColumnWidthsForAligment(flat_names_[i], level+1, printed, sum_[i]/active_[i], options);

    // Print Remainder
    if (sub_time > 0 ) {
      if (options.print_names_before_values) {
        std::ostringstream tmp;
        for (int l=0; l<=level; ++l)
          tmp << "|   ";
        tmp << "Remainder: ";
        alignments_.timer_names_ = std::max(alignments_.timer_names_,tmp.str().size());
      }
      {
        std::ostringstream tmp;
        tmp << sum_[i]/active_[i]- sub_time;
        alignments_.average_time_ = std::max(alignments_.average_time_,tmp.str().size());
      }
      if ( options.output_fraction && (sum_[i]/active_[i] > 0.) ) {
        std::ostringstream tmp;
        tmp << " - "<< (sum_[i]/active_[i]- sub_time)/(sum_[i]/active_[i])*100 << "%";
        alignments_.fraction_ = std::max(alignments_.fraction_,tmp.str().size());
      }
    }

    total_time += sum_[i]/active_[i];
  }
  return total_time;
}

double
StackedTimer::printLevel (std::string prefix, int print_level, std::ostream &os, std::vector<bool> &printed, double parent_time, const OutputOptions &options)
{
  // NOTE: If you change the outputting format or logic in this
  // function, you must make a corresponding change to the function
  // computeColumnWidthsForAlignment() or the alignments will be
  // incorrect if the user requests aligned output!

  double total_time = 0.0;

  for (int i=0; i<flat_names_.size(); ++i ) {
    if (sum_[i]/active_[i] <= options.drop_time) {
      continue;
    }
    if (printed[i])
      continue;
    int level = std::count(flat_names_[i].begin(), flat_names_[i].end(), '@');
    if ( (level != print_level) || (level >= options.max_levels) )
      continue;
    auto split_names = getPrefix(flat_names_[i]);
    if ( prefix != split_names.first)
      continue;

    // Output the indentation level
    if (options.print_names_before_values) {
      std::ostringstream tmp;
      for (int l=0; l<level; ++l) {
        tmp << "|   ";
      }
      // Output the timer name
      tmp << split_names.second << ": ";
      if (options.align_columns)
        os << std::left << std::setw(alignments_.timer_names_);
      os << tmp.str();
    }
    // output averge time
    {
      std::ostringstream tmp;
      tmp << sum_[i]/active_[i];
      if (options.align_columns)
        os << std::left << std::setw(alignments_.average_time_);
      os << tmp.str();
    }
    // output percentage
    if ( options.output_fraction && parent_time>0) {
      std::ostringstream tmp;
      tmp << " - "<<sum_[i]/active_[i]/parent_time*100<<"%";
      if (options.align_columns)
        os << std::left << std::setw(alignments_.fraction_);
      os << tmp.str();
    }
    // to keep alignment for later columns if requested
    else if (options.output_fraction) {
      if (options.align_columns)
        os << std::setw(alignments_.fraction_) << " ";
    }
    // output count
    {
      std::ostringstream tmp;
      tmp << " ["<<count_[i]/active_[i]<<"]";
      if (options.align_columns)
        os << std::left << std::setw(alignments_.count_);
      os << tmp.str();
    }
    // output total counts
    if ( options.output_total_updates ) {
      std::ostringstream tmp;
      tmp << " ("<<updates_[i]/active_[i]<<")";
      if (options.align_columns)
        os << std::left << std::setw(alignments_.total_updates_);
      os << tmp.str();
    }
    // Output min and maxs
    if ( options.output_minmax && active_[i]>1) {
      {
        std::ostringstream tmp;
        tmp << " {min="<<min_[i];
        if (options.align_columns)
          os << std::left << std::setw(alignments_.min_);
        os << tmp.str();
      }
      {
        std::ostringstream tmp;
        tmp <<", max="<<max_[i];
        if (active_[i] <= 1)
          tmp << "}";
        if (options.align_columns)
          os << std::left << std::setw(alignments_.max_);
        os << tmp.str();
      }
      if (procmin_.size()) {
        std::ostringstream tmp;
        tmp <<", proc min="<<procmin_[i];
        if (active_[i] <= 1)
          tmp << "}";
        if (options.align_columns)
          os << std::left << std::setw(alignments_.procmin_);
        os << tmp.str();
      }
      if (procmax_.size()) {
        std::ostringstream tmp;
        tmp <<", proc max="<<procmax_[i];
        if (active_[i] <= 1)
          tmp << "}";
        if (options.align_columns)
          os << std::left << std::setw(alignments_.procmax_);
        os << tmp.str();
      }
      if (active_[i]>1) {
        std::ostringstream tmp;
        tmp << ", std dev="<<sqrt(std::max<double>(sum_sq_[i]-sum_[i]*sum_[i]/active_[i],0.0)/(active_[i]-1));
        tmp << "}";
        if (options.align_columns)
          os << std::left << std::setw(alignments_.stddev_);
        os << tmp.str();
      }
    }
    else if ( options.output_minmax) {
      // this block keeps alignment for single rank timers
      size_t offset = alignments_.min_ + alignments_.max_ + alignments_.stddev_;
      for (size_t j=0; j < offset; ++j)
        os << " ";
    }

    // Output histogram
    if ( options.output_histogram && active_[i] >1 ) {
      std::ostringstream tmp;
      tmp << " <";
      for (int h=0;h<options.num_histogram; ++h) {
        if (h)
          tmp <<", "<<hist_[h][i];
        else
          tmp << hist_[h][i];
      }
      tmp << ">";
      if (options.align_columns)
        os << std::left << std::setw(alignments_.histogram_);
      os << tmp.str();
    }
    else if ( options.output_histogram) {
      // this block keeps alignment for single rank timers
      for (size_t j=0; j < alignments_.histogram_; ++j)
        os << " ";
    }

    if (! options.print_names_before_values) {
      std::ostringstream tmp;
      tmp << " ";
      for (int l=0; l<level; ++l) {
        tmp << "|   ";
      }
      // Output the timer name
      tmp << split_names.second << ": ";
      os << tmp.str();
    }

    os << std::endl;
    printed[i] = true;
    double sub_time = printLevel(flat_names_[i], level+1, os, printed, sum_[i]/active_[i], options);

    // Print Remainder
    if (sub_time > 0 ) {
      if (options.print_names_before_values) {
        std::ostringstream tmp;
        for (int l=0; l<=level; ++l)
          tmp << "|   ";
        tmp << "Remainder: ";
        if (options.align_columns)
          os << std::left << std::setw(alignments_.timer_names_);
        os << tmp.str();
      }
      {
        std::ostringstream tmp;
        tmp << sum_[i]/active_[i]- sub_time;
        if (options.align_columns)
          os << std::left << std::setw(alignments_.average_time_);
        os << tmp.str();
      }
      if ( options.output_fraction && (sum_[i]/active_[i] > 0.) ) {
        if (options.align_columns)
          os << std::left << std::setw(alignments_.fraction_);
        std::ostringstream tmp;
        tmp << " - "<< (sum_[i]/active_[i]- sub_time)/(sum_[i]/active_[i])*100 << "%";
        os << tmp.str();
      }
      if (! options.print_names_before_values) {
        {
          size_t offset = 0;
          offset += alignments_.count_;
          if (options.output_total_updates)
            offset += alignments_.total_updates_;
          if (options.output_minmax)
            offset += alignments_.min_ + alignments_.max_ + alignments_.stddev_;
          if (options.output_histogram)
            offset += alignments_.histogram_;
          for (size_t j=0; j < offset; ++j)
            os << " ";
        }
        std::ostringstream tmp;
        tmp << " ";
        for (int l=0; l<=level; ++l)
          tmp << "|   ";
        tmp << "Remainder: ";
        if (options.align_columns)
          os << std::left << std::setw(alignments_.timer_names_);
        os << tmp.str();
      }
      os << std::endl;
    }
    total_time += sum_[i]/active_[i];
  }
  return total_time;
}

static void printXMLEscapedString(std::ostream& os, const std::string& str)
{
  for(char c : str)
  {
    switch(c)
    {
      case '<':
        os << "&lt;";
        break;
      case '>':
        os << "&gt;";
        break;
      case '\'':
        os << "&apos;";
        break;
      case '"':
        os << "&quot;";
        break;
      case '&':
        os << "&amp;";
        break;
      //NOTE: unescaped curly braces {} are valid in XML,
      //however Watchr has a bug with parsing them
      case '{':
        os << '(';
        break;
      case '}':
        os << ')';
        break;
      default:
        os << c;
    }
  }
}

double
StackedTimer::printLevelXML (std::string prefix, int print_level, std::ostream& os, std::vector<bool> &printed, double parent_time, const std::string& rootName)
{
  constexpr int indSpaces = 2;
  int indent = indSpaces * print_level;

  double total_time = 0.0;

  for (int i=0; i<flat_names_.size(); ++i) {
    if (printed[i])
      continue;
    int level = std::count(flat_names_[i].begin(), flat_names_[i].end(), '@');
    if ( level != print_level)
      continue;
    auto split_names = getPrefix(flat_names_[i]);
    if ( prefix != split_names.first)
      continue;
    // Output the indentation level
    for (int j = 0; j < indent; j++)
      os << " ";
    os << "<timing name=\"";
    if(level == 0 && rootName.length())
      printXMLEscapedString(os, rootName);
    else
      printXMLEscapedString(os, split_names.second);
    os << "\" value=\"" << sum_[i]/active_[i] << "\"";
    printed[i] = true;
    //note: don't need to pass in prependRoot, since the recursive calls don't apply to the root level
    //Print the children to a temporary string. If it's empty, can close the current XML element on the same line.
    std::ostringstream osInner;
    double sub_time = printLevelXML(flat_names_[i], print_level+1, osInner, printed, sum_[i]/active_[i]);
    std::string innerContents = osInner.str();
    if(innerContents.length())
    {
      os << ">\n";
      os << innerContents;
      // Print Remainder
      if (sub_time > 0 ) {
        for (int j = 0; j < indent + indSpaces; j++)
          os << " ";
        os << "<timing name=\"Remainder\" value=\"" << (sum_[i]/active_[i] - sub_time) << "\"/>\n";
      }
      //having printed child nodes, close the XML element on its own line
      for (int j = 0; j < indent; j++)
        os << " ";
      os << "</timing>\n";
    }
    else
    {
      //Just a leaf node.
      os << "/>\n";
    }
    total_time += sum_[i]/active_[i];
  }
  return total_time;
}

void
StackedTimer::report(std::ostream &os, Teuchos::RCP<const Teuchos::Comm<int> > comm, OutputOptions options) {
  aggregateMpiData(comm, options);
  if (rank(*comm) == 0 ) {
    if (options.print_warnings) {
      os << "*** Teuchos::StackedTimer::report() - Remainder for a level will be ***"
         << "\n*** incorrect if a timer in the level does not exist on every rank  ***"
         << "\n*** of the MPI Communicator.                                        ***"
         << std::endl;
    }
    if ( (options.max_levels != INT_MAX) && options.print_warnings) {
      os << "Teuchos::StackedTimer::report() - max_levels manually set to " << options.max_levels
         << ". \nTo print more levels, increase value of OutputOptions::max_levels." << std::endl;
    }
    if ( (! options.print_names_before_values) && (! options.align_columns)) {
      options.align_columns = true;
      if (options.print_warnings)
        os << "Teuchos::StackedTimer::report() - option print_names_before_values=false "
           << "\nrequires that the option align_columns=true too. Setting the value for "
           << "\nalign_column to true."
           << std::endl;
    }
    if (options.align_columns) {
      std::vector<bool> printed(flat_names_.size(), false);
      computeColumnWidthsForAligment("", 0, printed, 0., options);
    }

    std::vector<bool> printed(flat_names_.size(), false);
    printLevel("", 0, os, printed, 0., options);
  }
}

void
StackedTimer::reportXML(std::ostream &os, const std::string& datestamp, const std::string& timestamp, Teuchos::RCP<const Teuchos::Comm<int> > comm)
{
  OutputOptions defaultOptions;
  aggregateMpiData(comm, defaultOptions);
  if (rank(*comm) == 0 ) {
    std::vector<bool> printed(flat_names_.size(), false);
    os << "<?xml version=\"1.0\"?>\n";
    os << "<performance-report date=\"" << timestamp << "\" name=\"nightly_run_" << datestamp << "\" time-units=\"seconds\">\n";
    printLevelXML("", 0, os, printed, 0.0);
    os << "</performance-report>\n";
  }
}

std::string
StackedTimer::reportWatchrXML(const std::string& name, Teuchos::RCP<const Teuchos::Comm<int> > comm) {
  const char* rawWatchrDir = getenv("WATCHR_PERF_DIR");
  const char* rawBuildName = getenv("WATCHR_BUILD_NAME");
  const char* rawGitSHA = getenv("TRILINOS_GIT_SHA");
  const char* rawBuildDateOverride = getenv("WATCHR_BUILD_DATE");
  //WATCHR_PERF_DIR is required (will also check nonempty below)
  if(!rawWatchrDir)
    return "";
  std::string watchrDir = rawWatchrDir;
  if(!watchrDir.length())
  {
    //Output directory has not been set, so don't produce output.
    return "";
  }
  //But the build name is optional (may be empty)
  std::string buildName = rawBuildName ? rawBuildName : "";
  std::string datestamp;
  std::string timestamp;
  {
    char buf[256];
    time_t t;
    struct tm* tstruct;
    time(&t);
    tstruct = gmtime(&t);
    if(rawBuildDateOverride)
    {
      //Parse the year, month, day
      int year = 0, month = 0, day = 0;
      sscanf(rawBuildDateOverride, "%d_%d_%d", &year, &month, &day);
      //Sanity check the values
      if(year <= 2000 || year > 2100)
        throw std::invalid_argument("$WATCHR_BUILD_DATE has invalid year or is not in YYYY_MM_DD format.");
      if(month < 1 || month > 12)
        throw std::invalid_argument("$WATCHR_BUILD_DATE has invalid month or is not in YYYY_MM_DD format.");
      if(day < 1 || day > 31)
        throw std::invalid_argument("$WATCHR_BUILD_DATE has invalid day or is not in YYYY_MM_DD format.");
      snprintf(buf, 256, "%04d_%02d_%02d", year, month, day);
      datestamp = buf;
      strftime(buf, 256, "T%H:%M:%S", tstruct);
      std::string justTime = buf;
      snprintf(buf, 256, "%04d-%02d-%02d", year, month, day);
      timestamp = std::string(buf) + justTime;
    }
    else
    {
      strftime(buf, 256, "%Y_%m_%d", tstruct);
      datestamp = buf;
      strftime(buf, 256, "%FT%H:%M:%S", tstruct);
      timestamp = buf;
    }
  }
  OutputOptions defaultOptions;
  aggregateMpiData(comm,defaultOptions);
  std::string fullFile;
  //only open the file on rank 0
  if(rank(*comm) == 0) {
    std::string nameNoSpaces = name;
    for(char& c : nameNoSpaces)
    {
      if(isspace(c))
        c = '_';
    }
    if(buildName.length())
    {
      //In filename, replace all whitespace with underscores
      std::string buildNameNoSpaces = buildName;
      for(char& c : buildNameNoSpaces)
      {
        if(isspace(c))
          c = '_';
      }
      fullFile = watchrDir + '/' + buildNameNoSpaces + "-" + nameNoSpaces + '_' + datestamp + ".xml";
    }
    else
      fullFile = watchrDir + '/' + nameNoSpaces + '_' + datestamp + ".xml";
    std::ofstream os(fullFile);
    std::vector<bool> printed(flat_names_.size(), false);
    os << "<?xml version=\"1.0\"?>\n";
    os << "<performance-report date=\"" << timestamp << "\" name=\"nightly_run_" << datestamp << "\" time-units=\"seconds\">\n";
    if(rawGitSHA)
    {
      std::string gitSHA(rawGitSHA);
      //Output the first 10 (hex) characters
      if(gitSHA.length() > 10)
        gitSHA = gitSHA.substr(0, 10);
      os << "  <metadata key=\"Trilinos Version\" value=\"" << gitSHA << "\"/>\n";
    }
    printLevelXML("", 0, os, printed, 0.0, buildName + ": " + name);
    os << "</performance-report>\n";
  }
  return fullFile;
}

void StackedTimer::enableVerbose(const bool enable_verbose)
{enable_verbose_ = enable_verbose;}

void StackedTimer::enableVerboseTimestamps(const unsigned levels)
{verbose_timestamp_levels_ = levels;}

void StackedTimer::setVerboseOstream(const Teuchos::RCP<std::ostream>& os)
{verbose_ostream_ = os;}

void StackedTimer::disableTimers()
{enable_timers_ = false;}

void StackedTimer::enableTimers()
{enable_timers_ = true;}

void StackedTimer::aggregateMpiData(Teuchos::RCP<const Teuchos::Comm<int> > comm, OutputOptions options)
{
  flatten();
  merge(comm);
  collectRemoteData(comm, options);
  global_mpi_aggregation_called_ = true;
}

double StackedTimer::getMpiAverageTime(const std::string& flat_timer_name)
{
  auto i = getFlatNameIndex(flat_timer_name);
  return sum_[i] / active_[i];
}

double StackedTimer::getMpiAverageCount(const std::string& flat_timer_name)
{
  auto i = getFlatNameIndex(flat_timer_name);
  return static_cast<double>(count_[i]) / static_cast<double>(active_[i]);
}

int StackedTimer::getFlatNameIndex(const std::string& flat_timer_name)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!global_mpi_aggregation_called_, std::runtime_error,
                             "ERROR: StackedTimer::getAverageMpiTime() - must call aggregateMpiData() first!");

  auto search = std::find(flat_names_.begin(),flat_names_.end(),flat_timer_name);

  TEUCHOS_TEST_FOR_EXCEPTION(search == flat_names_.end(),std::runtime_error,
                             "ERROR: StackedTimer::getAverageMpiTime() - the timer named \""
                             << flat_timer_name << "\" does not exist!");

  auto i = std::distance(flat_names_.begin(),search);
  return static_cast<int>(i);
}

bool StackedTimer::isTimer(const std::string& flat_timer_name)
{
  TEUCHOS_TEST_FOR_EXCEPTION(!global_mpi_aggregation_called_,std::runtime_error,
                             "ERROR: StackedTimer::isTimer() - must call aggregateMpiData() before using this query!");

  auto search = std::find(flat_names_.begin(),flat_names_.end(),flat_timer_name);
  return (search == flat_names_.end()) ? false : true;
}

} //namespace Teuchos

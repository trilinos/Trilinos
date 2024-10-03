// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_TESTINGXMLUTILITIES_HPP_
#define TPETRA_TESTINGXMLUTILITIES_HPP_

/// \file Tpetra_TestingXMLUtilities.hpp
/// \brief Internal utilities for testing Tpetra.
///
/// \warning This header file and its contents are implementation
///   details of Tpetra.  Users must not rely on this file existing,
///   or on any contents of this file.

#include <ctime>
#include "Teuchos_DefaultSerialComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_StackedTimer.hpp"
#include "Teuchos_Array.hpp"

namespace Tpetra {
  template<class T>
  class TestingXMLUtilities {
  public:
    /// \brief Reporting of data in a Watchr-compatible XML format
    ///
    /// \param label_class [in] - XML object type to output
    ///
    /// \param units [in] - Units to output in reporting
    ///
    /// \param base_name [in] - base name to prepend to the XML "name" in the individual line output
    ///
    /// \param names [in] - List of individual outputs
    ///
    /// \param values [in] - Values for individual output (these will be gathered across all MPI ranks)
    ///
    /// \param comm [in] - Teuchos::Comm for reductions
    ///
    /// \warning This code is unlikely to work correctly if the "names" don't all appear on each rank in the same order"
    ///
    ///
    /// Sample output (items in parentheses are input parameters, items in braces are enviornment varaibels
    ///
    /// <?xml version="1.0"?>
    /// <performance-report date="2023-08-07T20:38:11" name="nightly_run_2023_08_07" time-units="(units)">
    /// <metadata key="Trilinos Version" value="{TRILINOS_GIT_SHA}"/>
    /// <memory name="{WATCHR_BUILD_NAME}: (base_name) (names[0])" value="values[0]"/>
    /// <memory name="{WATCHR_BUILD_NAME}: (base_name) (names[1])" value="values[1]"/>
    /// <memory name="{WATCHR_BUILD_NAME}: (base_name) (names[2])" value="values[2]"/>
    /// </performance-report>
    std::string
    reportWatchrXML(const std::string label_class,const std::string units, const std::string base_name, const Teuchos::Array<std::string> & names, const Teuchos::Array<T> & values, Teuchos::RCP<const Teuchos::Comm<int> > comm) {
      using Teuchos::Array;
      
      TEUCHOS_TEST_FOR_EXCEPTION(names.size()!=values.size(), std::runtime_error, "reportWatchrXML: names and values are not the same size");
      flat_names_=names;
      flat_times_=values;

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

      Teuchos::StackedTimer::OutputOptions defaultOptions;
      collectRemoteData(comm, defaultOptions);

      std::string fullFile;
      //only open the file on rank 0
      if(rank(*comm) == 0) {
        std::string nameNoSpaces = base_name;
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
        Teuchos::Array<bool> printed(flat_names_.size(), false);
        os << "<?xml version=\"1.0\"?>\n";
        os << "<performance-report date=\"" << timestamp << "\" name=\"nightly_run_" << datestamp << "\" time-units=\""<<units<<"\">\n";
        if(rawGitSHA)
          {
            std::string gitSHA(rawGitSHA);
            //Output the first 10 (hex) characters
            if(gitSHA.length() > 10)
              gitSHA = gitSHA.substr(0, 10);
            os << "  <metadata key=\"Trilinos Version\" value=\"" << gitSHA << "\"/>\n";
          }
        printLevelXML(label_class,"", 0, os, printed, 0.0, buildName + ": " + base_name);
        os << "</performance-report>\n";
      }
      return fullFile;
     
    }


  private:

    void
    collectRemoteData(Teuchos::RCP<const Teuchos::Comm<int> > comm, const Teuchos::StackedTimer::OutputOptions &options) {
      using Teuchos::Array;
      using Teuchos::reduce;
      using Teuchos::reduceAll;
      using Teuchos::REDUCE_MAX;
      using Teuchos::REDUCE_MIN;
      using Teuchos::REDUCE_SUM;

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
      /*
        for (int i=0;i<num_names; ++i) {
        bool found = false; // ignore result here
        auto t = timer_.findTimer(flat_names_[i],found);
        time[i] = t.time;
        count[i] = t.count;
        used[i] = t.count==0? 0:1;
        if (options.output_total_updates)
        updates[i] = t.updates;
        }
      */
      for (int i=0;i<num_names; ++i) {
        time[i] = flat_times_[i];
        count[i] = 1;
        used[i] = 1;
        if (options.output_total_updates)
          updates[i] = 1;
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


    double
    printLevelXML (std::string label_class, std::string prefix, int print_level, std::ostream& os, Teuchos::Array<bool> &printed, double parent_time, const std::string& rootName="")
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
        os << "<"<<label_class<<" name=\"";
        if(level == 0 && rootName.length())
          printXMLEscapedString(os, rootName);
        else
          printXMLEscapedString(os, split_names.second);
        os << " "<<flat_names_[i];
        os << "\" value=\"" << sum_[i]/active_[i] << "\"";
        printed[i] = true;
        //note: don't need to pass in prependRoot, since the recursive calls don't apply to the root level
        //Print the children to a temporary string. If it's empty, can close the current XML element on the same line.
        std::ostringstream osInner;
        double sub_time = printLevelXML(label_class,flat_names_[i], print_level+1, osInner, printed, sum_[i]/active_[i]);
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


    std::pair<std::string, std::string> getPrefix(const std::string &name) {
      for (std::size_t i=name.size()-1; i>0; --i)
        if (name[i] == '@') {
          return std::pair<std::string, std::string>(name.substr(0,i), name.substr(i+1));
        }
      return std::pair<std::string, std::string>(std::string(""), name);
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



    Teuchos::Array<std::string> flat_names_;
    Teuchos::Array<T> flat_times_;
    Teuchos::Array<double> min_;
    Teuchos::Array<double> max_;
    Teuchos::Array<int> procmin_;
    Teuchos::Array<int> procmax_;
    Teuchos::Array<double> sum_;
    Teuchos::Array<double> sum_sq_;
    Teuchos::Array<Teuchos::Array<int>> hist_;
    Teuchos::Array<unsigned long> count_;
    Teuchos::Array<unsigned long long> updates_;
    Teuchos::Array<int> active_;
    
  }; // class TestingXMLUtilities
} // namespace Tpetra

#endif // TPETRA_TESTINGXMLUTILITIES_HPP_

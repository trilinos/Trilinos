/*
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan2 Directory for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */

#ifndef ZOLTAN2_DIRECTORY_CLOCK_H_
#define ZOLTAN2_DIRECTORY_CLOCK_H_

#include <Teuchos_CommHelpers.hpp>
#include <time.h>

// a temporary class to get some timing information out - to be deleted
class Zoltan2_Directory_Clock {
  public:
    struct time_data {
      time_data() : count(0), sum(0.0) {}
      int count;
      double sum;
      int sublevel;
    };

    Zoltan2_Directory_Clock(std::string name_, int sublevel_ = 0)
      : name(name_), sublevel(sublevel_), bCompleted(false) {
      startTime = getTime();
    }

    ~Zoltan2_Directory_Clock() {
      if(!bCompleted) {
        complete();
      }
    }

    void complete() {
      if(get_time_data().find(name) == get_time_data().end()) {
        get_time_data()[name] = time_data(); // create with 0 sum and 0 count
      }
      get_time_data()[name].sublevel = sublevel;
      get_time_data()[name].count += 1;
      get_time_data()[name].sum += (getTime() - startTime);
      bCompleted = true;
    }

    static void clearTimes() {
      get_time_data().clear();
    }

    static void logAndClearTimes(const std::string& message,
      Teuchos::RCP<const Teuchos::Comm<int> > comm) {

      // convert maps to vectors
      std::vector<std::string> map_names(get_time_data().size());
      std::vector<double> localSums(get_time_data().size()); // for reduce
      std::vector<time_data> map_info(get_time_data().size());

      // each proc has a summed time for each clock name
      // reduce to get the average across all procs
      size_t index = 0;
      for(auto itr = get_time_data().begin(); itr != get_time_data().end(); ++itr) {
        map_names[index] = itr->first;
        localSums[index] = itr->second.sum; // for reduce
        map_info[index] = itr->second;
        ++index;
      }

      std::vector<double> avSums(get_time_data().size(), 0.0);
      Teuchos::reduceAll<int, double>(*comm, Teuchos::REDUCE_SUM,
        get_time_data().size(), &localSums[0], &avSums[0]);
      for(size_t n = 0; n < avSums.size(); ++n) {
        avSums[n] /= static_cast<double>(comm->getSize());
      }

      // Now we can plot some messaging
      // First plot the average fo all ranks, then plot the average just for
      // rank 0 so we have some idea about fluctuations. We expect the values to
      // be similar though for edges cases (very low total IDs) it's possible
      // the ranks will be managing very different loads
      if(comm->getRank() == 0) {
        std::cout << message << std::endl;
        for(size_t n = 0; n < map_names.size(); ++n) {
          for(int indent = 0; indent < map_info[n].sublevel; ++indent) {
            std::cout << "  ";
          }
          std::cout << "  " << std::setw(20-map_info[n].sublevel*2)
            << std::left << map_names[n] << " "
            << "(" << map_info[n].count << ")  "
            << std::fixed << std::setprecision(2)
            << avSums[n] << " (" << localSums[n] << ")" << std::endl;
        }
      }

      clearTimes();
      comm->barrier();
    }

  private:
    typedef std::map<std::string, time_data> time_data_t;

    static time_data_t & get_time_data() {
      static time_data_t static_time_data;
      return static_time_data;
    }

    double getTime() { return static_cast<double>(clock()) /
      static_cast<double>(CLOCKS_PER_SEC); }

    std::string name;
    int sublevel;
    double startTime;
    bool bCompleted;
};

#endif

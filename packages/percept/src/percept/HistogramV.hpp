// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef HistogramV_hpp
#define HistogramV_hpp

#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <stk_util/diag/PrintTable.hpp>

namespace percept
{

#if defined(STK_HAS_MPI)

    inline void my_all_reduce_max( stk::ParallelMachine comm ,
                                           const double  * local , double * global , unsigned count )
    {
      double * tmp = const_cast<double*>( local );
      MPI_Allreduce( tmp , global , count , MPI_DOUBLE , MPI_MAX , comm );
    }

    inline void my_all_reduce_max( stk::ParallelMachine comm ,
                                        const int  * local , int * global , unsigned count )
    {
      int * tmp = const_cast<int*>( local );
      MPI_Allreduce( tmp , global , count , MPI_INT , MPI_MAX , comm );
    }

    // min
    inline void my_all_reduce_min( stk::ParallelMachine comm ,
                                   const double  * local , double * global , unsigned count )
    {
      double * tmp = const_cast<double*>( local );
      MPI_Allreduce( tmp , global , count , MPI_DOUBLE , MPI_MIN , comm );
    }

    inline void my_all_reduce_min( stk::ParallelMachine comm ,
                                   const int  * local , int * global , unsigned count )
    {
      int * tmp = const_cast<int*>( local );
      MPI_Allreduce( tmp , global , count , MPI_INT , MPI_MIN , comm );
    }

    // sum
    inline void my_all_reduce_sum( stk::ParallelMachine comm ,
                                   const unsigned  * local , unsigned * global , unsigned count )
    {
      unsigned * tmp = const_cast<unsigned*>( local );
      MPI_Allreduce( tmp , global , count , MPI_UNSIGNED , MPI_SUM , comm );
    }



#endif


template<typename T>
class HistogramV
{
public:
  typedef std::vector<T> RangesType;
  typedef std::pair<unsigned,unsigned> CountStatsType;
  typedef unsigned CountsType;

protected:
  CountStatsType m_count_stats;
  std::vector<T> m_ranges;
  std::vector<CountsType> m_counts;
  const std::vector<T>& m_data;

public:
  // create histogram using passed in data
  HistogramV(const std::vector<T>& data)
    : m_data(data)
  {
    reset();
  }

  virtual ~HistogramV() {}

  // compute ranges for bins in the histogram
  void compute_uniform_bins(unsigned num_bins, bool use_log_scale=false)
  {
    compute_ranges();
    unsigned num_ranges = num_bins+1;
    T min = m_ranges[0];
    T max = m_ranges[1];
    if (use_log_scale && min <= 0 )
      throw std::runtime_error("Histogram::compute_uniform_bins request for log scale for data with nonpositive values");

    m_ranges.resize(num_ranges);
    m_ranges[0] = min;
    m_ranges[num_bins] = max;
    for (unsigned i=1; i < num_ranges-1; ++i)
      {
        double di = double(i)/double(num_bins);
        //  di = std::log(di);
        if (use_log_scale)
          {
            m_ranges[i] = std::exp(T(std::log(double(min)) + di*(std::log(double(max))-std::log(double(min)))));
          }
        else
          {
            m_ranges[i] = T(double(min) + di*(double(max)-double(min)));
          }
      }
    m_count_stats = compute_count();
  }

  const std::vector<T>& get_data() const
  {
    return m_data;
  }


  // get the ranges of the histogram bins
  const RangesType& get_ranges() const
  {
    return m_ranges;
  }

  const CountStatsType& get_count_stats() const
  {
    return m_count_stats;
  }

  const std::vector<CountsType>& get_counts() const
  {
    return m_counts;
  }

  // set the ranges of the histogram bins
  void set_ranges(const RangesType& ranges)
  {
    for (unsigned i=1; i < ranges.size(); ++i)
      {
        //std::cout << "ranges[" << i-1 << "] = " << ranges[i-1] << " ranges[" << i << "]= " << ranges[i] << std::endl;
        if (ranges[i] <= ranges[i-1])
          {
            throw std::runtime_error("Histogram::set_ranges ranges must be monotonically increasing");
          }
      }
    m_ranges = ranges;
    m_count_stats = compute_count();
  }


protected:
  /** resets values computed from the data */
  void reset() {
    m_ranges.resize(2);
    m_counts.resize(0);
  }

  /** gets min/max of the data */
  void compute_ranges()
  {
    reset();
    m_ranges[0] = std::numeric_limits<T>::max();
    m_ranges[1] = -(0.1*std::numeric_limits<T>::max());
    for (size_t i=0; i < m_data.size(); ++i)
      {
        m_ranges[0] = std::min(m_ranges[0], m_data[i]);
        m_ranges[1] = std::max(m_ranges[1], m_data[i]);
      }
#if defined( STK_HAS_MPI )
    T minv=m_ranges[0];
    T maxv=m_ranges[1];
    my_all_reduce_min( MPI_COMM_WORLD , &minv , &m_ranges[0] , 1);
    my_all_reduce_max( MPI_COMM_WORLD , &maxv , &m_ranges[1] , 1);
#endif
  }

  std::pair<unsigned,unsigned> compute_count()
  {
    m_counts.resize(m_ranges.size() - 1);
    m_counts.assign(m_ranges.size() - 1, 0);
    for (unsigned i=0; i < m_data.size(); ++i)
      {
        T data = m_data[i];
        for (size_t j=1; j < m_ranges.size(); ++j)
          {
            if ((j == m_ranges.size() -1 && m_ranges[j-1] <= data && data <= m_ranges[j])
                || ( m_ranges[j-1] <= data && data < m_ranges[j]))
              {
                ++m_counts[j-1];
              }
          }
      }
#if defined( STK_HAS_MPI )
    std::vector<unsigned> cc(m_counts);
    my_all_reduce_sum( MPI_COMM_WORLD , &cc[0] , &m_counts[0] , m_counts.size());
#endif

    unsigned max_count = m_counts[0];
    unsigned total_count = m_counts[0];
    for (unsigned i=1; i < m_counts.size(); i++)
      {
        max_count = std::max(max_count, m_counts[i]);
        total_count += m_counts[i];
      }
    return std::pair<unsigned,unsigned>(max_count,total_count);
  }

public:
  static void print_histogram(const HistogramV<T>& h, std::ostream& stream)
  {
    const std::vector<T>& data= h.get_data();
    const std::vector<typename HistogramV<T>::CountsType>& counts = h.get_counts();
    const typename HistogramV<T>::RangesType& ranges = h.get_ranges();

    stream << "Histogram::print:: data= " << std::endl;
    for (unsigned i=0; i < data.size(); ++i)
      {
        stream << " " << data[i];
      }
    stream << std::endl;
    stream << "Histogram::print:: ranges= " << std::endl;
    for (unsigned i=0; i < ranges.size(); ++i)
      {
        stream << " " << ranges[i];
      }
    stream << std::endl;

    stream << "Histogram::print:: counts= " << std::endl;
    for (unsigned i=0; i < counts.size(); ++i)
      {
        stream << " " << counts[i];
      }
    stream << std::endl;
  }


  static void print_simple_table(const HistogramV<T>& h, const std::string titles[4], std::ostream& stream)
  {
    const std::pair<unsigned,unsigned>& max_tot_count = h.get_count_stats();
    size_t width1 = titles[1].length()+2; // for "# "
    size_t width2 = titles[2].length();
    size_t width3 = titles[3].length();
    size_t width10 = width1;
    size_t width20 = width2;
    size_t width30 = width3;
    size_t w0=0, w1=0, w11=0;

    const std::vector<typename HistogramV<T>::CountsType>& counts = h.get_counts();
    const typename HistogramV<T>::RangesType& ranges = h.get_ranges();

    for (unsigned i=0; i < counts.size(); ++i)
      {
        std::ostringstream ostr0, ostr1, ostr11, ostr2, ostr3;
        ostr0 << ranges[i];
        ostr1 << ranges[i+1];
        ostr11 << (ranges[i]+ranges[i+1])/2.;
        w0 = std::max(w0, ostr0.str().length());
        w1 = std::max(w1, ostr1.str().length());
        w11 = std::max(w11, ostr11.str().length());

        unsigned count = counts[i];
        double percent = 100.0*double(count)/double(max_tot_count.second);
        //stream << " " << count << " " << percent << std::endl;
        ostr2 << count;
        ostr3 << percent;
        width2 = std::max(width2, ostr2.str().length());
        width3 = std::max(width3, ostr3.str().length());
      }
    size_t w012=w0+w1+w11+3;
    if (w012 < width1)
      {
        if (width1 % 3) ++width1;
        if (width1 % 3) ++width1;
        if (width1 % 3) ++width1;
        w0=width1/3;
        w1=width1/3;
        w11=width1/3;
      }
    else
      {
        width1 = w012;
      }

    std::string title1_filled = "# "+titles[1]+ (width10<width1?std::string(width1-width10,' '):"");
    std::string title2_filled = titles[2] + (width20<width2?std::string(width2-width20,' '):"");
    std::string title3_filled = titles[3] + (width30<width3?std::string(width3-width30,' '):"");
    stream << "# Histogram: " << titles[0] << std::endl;
    stream << std::setw(width1) << std::right << title1_filled << " " << std::setw(width2) <<  title2_filled << " " << std::setw(width3) << title3_filled << std::endl;
    std::string line(width1+width2+width3+2,'-');
    line = "# "+line;
    stream << line << std::endl;
    for (unsigned i=0; i < counts.size(); ++i)
      {
        stream << std::setw(w0) << ranges[i] << " ";
        stream << std::setw(w1) << ranges[i+1] << " ";
        stream << std::setw(w11) << (ranges[i]+ranges[i+1])/2.;
        unsigned count = counts[i];
        double percent = 100.0*double(count)/double(max_tot_count.second);
        stream << " " << std::setw(width2) << count << " " << std::setw(width3) << percent << std::endl;
      }
  }

  static void print_table(const HistogramV<T>& h, const std::string titles[4], std::ostream& os, unsigned max_column_width = 80, const char bar_symbol = '=', bool use_percentage=false)
  {
    std::string title=titles[0];
    const std::pair<unsigned,unsigned>& max_tot_count = h.get_count_stats();
    if (!use_percentage && max_tot_count.first > max_column_width)
      {
        os << "Histogram: using percentages since count is greater than screen width";
        print_table(h, titles, os, max_column_width, bar_symbol, true);
        return;
      }
    stk::PrintTable table;
    table.setTitle(title+"\n");
    //table.setAutoEndCol(false);

    std::ostringstream ostr0;
    double npsym = 100.0/double(max_column_width);
    if (use_percentage)
      {
        ostr0 << "[0:100]" << "  each '" << bar_symbol << "' is " << npsym << "%";
      }
    else
      ostr0 << "[0:" << max_tot_count.first << "]";

    table << "|" << stk::justify(stk::PrintTable::Cell::CENTER) << "Ranges" <<  "|" << stk::justify(stk::PrintTable::Cell::CENTER)  << (use_percentage?"Percentage In Range":"Count In Bucket") << "|" << stk::end_header;
    table << "|" << stk::justify(stk::PrintTable::Cell::CENTER) << "      " <<  "|" << stk::justify(stk::PrintTable::Cell::CENTER)  << ostr0.str() << "|" << stk::end_header;
    if (use_percentage)
      {
        std::string line_equal(max_column_width, bar_symbol);
        table << "|" << "       " <<  "|" << stk::justify(stk::PrintTable::Cell::CENTER)  << line_equal << "|" << stk::end_header;
      }

    const std::vector<typename HistogramV<T>::CountsType>& counts = h.get_counts();
    const typename HistogramV<T>::RangesType& ranges = h.get_ranges();

    //os << "debug: m_ranges.size() = " << m_ranges.size() << " m_counts.size() = " << m_counts.size();
    for (unsigned i=0; i < counts.size(); ++i)
      {
        std::ostringstream ostr;
        ostr << std::setw(8) << ranges[i];
        ostr << " - ";
        ostr << std::setw(8) << ranges[i+1];
        //table << "|" << m_ranges[i-1] << " - " << m_ranges[i] << "|" ;
        table << "|" << ostr.str() << "|" ;
        unsigned end_loop = counts[i];
        if (use_percentage)
          end_loop = unsigned(double(max_column_width)*double(end_loop)/double(max_tot_count.second));
        std::string count(end_loop, bar_symbol);
        table << stk::justify(stk::PrintTable::Cell::LEFT) <<  count << "|" << stk::end_row;
      }

    os << "\n" << table;

  }

};


}

#endif // HistogramV_hpp

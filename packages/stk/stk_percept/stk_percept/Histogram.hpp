#ifndef stk_percept_Histogram_h
#define stk_percept_Histogram_h

#include <stk_util/diag/PrintTable.hpp>

#if defined( STK_HAS_MPI )
#include <stk_util/parallel/ParallelReduce.hpp>
#endif

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <utility>
#include <string>
#include <cmath>
#include <sstream>
#include <vector>
#include <limits>

namespace stk {
  namespace percept {


    /**
     * Usage:
     *
     *   Histogram<double> h;
     *   double data[] = {1,2,3,5,8,13};
     *   // create the data - any std::vector methods usable
     *   h.assign(data,data+6);
     *   h.compute_uniform_bins(2);
     *   h.print(std::cout);
     *   h.compute_uniform_bins(2, true); // use log scale
     *   h.print(std::cout);
     *   double ranges[] = {1,4,9,13};
     *   std::vector<double> vranges(ranges,ranges+4);
     *   // user-supplied bins
     *   h.set_ranges(vranges);
     *   h.print(std::cout);
     *   // print in table format
     *   h.print_table(std::cout);
     *   // print in table format with percentages
     *   h.print_table(std::cout, true);
     */
#if defined(STK_HAS_MPI)

    inline void my_all_reduce_max( ParallelMachine comm ,
                                           const double  * local , double * global , unsigned count )
    {
      double * tmp = const_cast<double*>( local );
      MPI_Allreduce( tmp , global , count , MPI_DOUBLE , MPI_MAX , comm );
    }

    inline void my_all_reduce_max( ParallelMachine comm ,
                                        const int  * local , int * global , unsigned count )
    {
      int * tmp = const_cast<int*>( local );
      MPI_Allreduce( tmp , global , count , MPI_INT , MPI_MAX , comm );
    }

    // min
    inline void my_all_reduce_min( ParallelMachine comm ,
                                   const double  * local , double * global , unsigned count )
    {
      double * tmp = const_cast<double*>( local );
      MPI_Allreduce( tmp , global , count , MPI_DOUBLE , MPI_MIN , comm );
    }

    inline void my_all_reduce_min( ParallelMachine comm ,
                                   const int  * local , int * global , unsigned count )
    {
      int * tmp = const_cast<int*>( local );
      MPI_Allreduce( tmp , global , count , MPI_INT , MPI_MIN , comm );
    }

    // sum
    inline void my_all_reduce_sum( ParallelMachine comm ,
                                   const unsigned  * local , unsigned * global , unsigned count )
    {
      unsigned * tmp = const_cast<unsigned*>( local );
      MPI_Allreduce( tmp , global , count , MPI_UNSIGNED , MPI_SUM , comm );
    }



#endif

    template<typename T>
    class Histogram : public std::vector<T>
    {
      typedef unsigned CountsType;
      typedef std::pair<unsigned,unsigned> CountStatsType;
      CountStatsType m_count_stats;
      std::vector<T> m_ranges;
      std::vector<CountsType> m_counts;
      unsigned m_max_column_width;
      const char m_bar_symbol;
      std::string m_titles[4];
    public:
      typedef std::vector<T> RangesType;

      Histogram(unsigned max_column_width=40, const char bar_symbol = '=') :
        m_max_column_width(max_column_width), m_bar_symbol(bar_symbol)
      {
        reset();
        set_titles();
      }

      //virtual ~Histogram() {}
      void set_titles(std::string title="Histogram",
                      std::string title1="Ranges / Mid-Point", std::string title2="Counts", std::string title3="Percentage")
      {
        m_titles[0]=title;
        m_titles[1]=title1;
        m_titles[2]=title2;
        m_titles[3]=title3;
      }

    private:
      /** resets values computed from the data */
      void reset() {
        m_ranges.resize(2);
        m_counts.resize(0);
      }

      /** gets min/max of the data */
      void compute_ranges()
      {
        if (!this->size()) return;
        reset();
        m_ranges[0] = m_ranges[1] = (*this)[0];
        for (unsigned i=1; i < this->size(); ++i)
          {
            m_ranges[0] = std::min(m_ranges[0], (*this)[i]);
            m_ranges[1] = std::max(m_ranges[1], (*this)[i]);
          }
#if defined( STK_HAS_MPI )
        T minv=m_ranges[0];
        T maxv=m_ranges[1];
        my_all_reduce_min( MPI_COMM_WORLD , &minv , &m_ranges[0] , 1);
        my_all_reduce_max( MPI_COMM_WORLD , &maxv , &m_ranges[1] , 1);
#endif
      }
    public:

      /** returns a copy of current ranges */
      RangesType get_ranges() { return m_ranges; }
      void set_ranges(RangesType& ranges) {
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

    private:
      std::pair<unsigned,unsigned> compute_count()
      {
        m_counts.resize(m_ranges.size() - 1);
        m_counts.assign(m_ranges.size() - 1, 0);
        for (unsigned i=0; i < this->size(); ++i)
          {
            T data = (*this)[i];
            for (unsigned j=1; j < m_ranges.size(); ++j)
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

      void print(std::ostream& stream)
      {
        stream << "Histogram::print:: data= " << std::endl;
        for (unsigned i=0; i < this->size(); ++i)
          {
            stream << " " << (*this)[i];
          }
        stream << std::endl;
        stream << "Histogram::print:: ranges= " << std::endl;
        for (unsigned i=0; i < m_ranges.size(); ++i)
          {
            stream << " " << m_ranges[i];
          }
        stream << std::endl;

        stream << "Histogram::print:: counts= " << std::endl;
        for (unsigned i=0; i < m_counts.size(); ++i)
          {
            stream << " " << m_counts[i];
          }
        stream << std::endl;
      }

      void print_simple_table(std::ostream& stream)
      {
        std::pair<unsigned,unsigned> max_tot_count = m_count_stats;
        size_t width1 = m_titles[1].length()+2; // for "# "
        size_t width2 = m_titles[2].length();
        size_t width3 = m_titles[3].length();
        size_t width10 = width1;
        size_t width20 = width2;
        size_t width30 = width3;
        size_t w0=0, w1=0, w11=0;

        for (unsigned i=0; i < m_counts.size(); ++i)
          {
            std::ostringstream ostr0, ostr1, ostr11, ostr2, ostr3;
            ostr0 << m_ranges[i];
            ostr1 << m_ranges[i+1];
            ostr11 << (m_ranges[i]+m_ranges[i+1])/2.;
            w0 = std::max(w0, ostr0.str().length());
            w1 = std::max(w1, ostr1.str().length());
            w11 = std::max(w11, ostr11.str().length());

            unsigned count = m_counts[i];
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

        std::string title1_filled = "# "+m_titles[1]+ (width10<width1?std::string(width1-width10,' '):"");
        std::string title2_filled = m_titles[2] + (width20<width2?std::string(width2-width20,' '):"");
        std::string title3_filled = m_titles[3] + (width30<width3?std::string(width3-width30,' '):"");
        stream << "# Histogram: " << m_titles[0] << std::endl;
        stream << std::setw(width1) << std::right << title1_filled << " " << std::setw(width2) <<  title2_filled << " " << std::setw(width3) << title3_filled << std::endl;
        std::string line(width1+width2+width3+2,'-');
        line = "# "+line;
        stream << line << std::endl;
        for (unsigned i=0; i < m_counts.size(); ++i)
          {
            stream << std::setw(w0) << m_ranges[i] << " ";
            stream << std::setw(w1) << m_ranges[i+1] << " ";
            stream << std::setw(w11) << (m_ranges[i]+m_ranges[i+1])/2.;
            unsigned count = m_counts[i];
            double percent = 100.0*double(count)/double(max_tot_count.second);
            stream << " " << std::setw(width2) << count << " " << std::setw(width3) << percent << std::endl;
          }
      }

      void print_table(std::ostream& os, bool use_percentage=false)
      {
        std::string title=m_titles[0];
        std::pair<unsigned,unsigned> max_tot_count = m_count_stats;
        if (!use_percentage && max_tot_count.first > m_max_column_width)
          {
            os << "Histogram: using percentages since count is greater than screen width";
            print_table(os, true);
          }
        stk::PrintTable table;
        table.setTitle(title+"\n");
        //table.setAutoEndCol(false);

        std::ostringstream ostr0;
        double npsym = 100.0/double(m_max_column_width);
        if (use_percentage)
          {
            ostr0 << "[0:100]" << "  each '" << m_bar_symbol << "' is " << npsym << "%";
          }
        else
          ostr0 << "[0:" << max_tot_count.first << "]";

        table << "|" << justify(PrintTable::Cell::CENTER) << "Ranges" <<  "|" << justify(PrintTable::Cell::CENTER)  << (use_percentage?"Percentage In Range":"Count In Bucket") << "|" << stk::end_header;
        table << "|" << justify(PrintTable::Cell::CENTER) << "      " <<  "|" << justify(PrintTable::Cell::CENTER)  << ostr0.str() << "|" << stk::end_header;
        if (use_percentage)
          {
            std::string line_equal(m_max_column_width, m_bar_symbol);
            table << "|" << "       " <<  "|" << justify(PrintTable::Cell::CENTER)  << line_equal << "|" << stk::end_header;
          }

        //os << "debug: m_ranges.size() = " << m_ranges.size() << " m_counts.size() = " << m_counts.size();
        for (unsigned i=0; i < m_counts.size(); ++i)
          {
            std::ostringstream ostr;
            ostr << m_ranges[i] << " - " << m_ranges[i+1];
            //table << "|" << m_ranges[i-1] << " - " << m_ranges[i] << "|" ;
            table << "|" << ostr.str() << "|" ;
            unsigned end_loop = m_counts[i];
            if (use_percentage)
              end_loop = unsigned(double(m_max_column_width)*double(end_loop)/double(max_tot_count.second));
            std::string count(end_loop, m_bar_symbol);
            table << justify(PrintTable::Cell::LEFT) <<  count << "|" << stk::end_row;
          }

        os << "\n" << table;

      }

    };

  }
}

#endif

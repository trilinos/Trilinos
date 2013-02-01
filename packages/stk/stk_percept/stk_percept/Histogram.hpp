#ifndef stk_percept_Histogram_h
#define stk_percept_Histogram_h

#include <stk_util/diag/PrintTable.hpp>

#include <iostream>
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
     *   h.compute_ranges();
     *   h.compute_uniform_buckets(2);
     *   h.print(std::cout);
     *   h.compute_uniform_buckets(2, true); // use log scale
     *   h.print(std::cout);
     *   double ranges[] = {1,4,9,13};
     *   std::vector<double> vranges(ranges,ranges+4);
     *   // user-supplied buckets
     *   h.set_ranges(vranges);
     *   h.print(std::cout);
     *   // print in table format
     *   h.print_table(std::cout);
     *   // print in table format with percentages
     *   h.print_table(std::cout, true);
     */
    template<typename T>
    class Histogram : public std::vector<T>
    {
      std::vector<T> m_ranges;
      std::vector<unsigned> m_counts;
      unsigned m_max_column_width;
    public:
      typedef std::vector<T> RangesType;

      Histogram(unsigned max_column_width=40) : m_max_column_width(max_column_width)
      {
        reset();
      }
      ~Histogram() {}

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
      }

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
      }

      void compute_uniform_buckets(unsigned num_buckets, bool use_log_scale=false)
      {
        compute_ranges();
        unsigned num_ranges = num_buckets+1;
        T min = m_ranges[0];
        T max = m_ranges[1];
        if (use_log_scale && min <= 0 )
          throw std::runtime_error("Histogram::compute_uniform_buckets request for log scale for data with nonpositive values");

        m_ranges.resize(num_ranges);
        m_ranges[0] = min;
        m_ranges[num_buckets] = max;
        for (unsigned i=1; i < num_ranges-1; ++i)
          {
            double di = double(i)/double(num_buckets);
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
      }

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
                  ++m_counts[j-1];
              }
          }
        unsigned max_count = m_counts[0];
        unsigned total_count = m_counts[0];
        for (unsigned i=1; i < m_counts.size(); i++)
          {
            max_count = std::max(max_count, m_counts[i]);
            total_count += m_counts[i];
          }
        return std::pair<unsigned,unsigned>(max_count,total_count);
      }

      void print(std::ostream& stream)
      {
        compute_count();
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

      void print_table(std::ostream& os, bool use_percentage=false)
      {
        std::pair<unsigned,unsigned> max_tot_count = compute_count();
        if (!use_percentage && max_tot_count.first > m_max_column_width)
          {
            os << "Histogram: using percentages since count is greater than screen width";
            print_table(os, true);
          }
        stk::PrintTable table;
        table.setTitle("Histogram\n");
        //table.setAutoEndCol(false);

        std::ostringstream ostr0;
        if (use_percentage)
          ostr0 << "[0:100]";
        else
          ostr0 << "[0:" << max_tot_count.first << "]";

        table << "|" << justify(PrintTable::Cell::CENTER) << "Buckets" <<  "|" << justify(PrintTable::Cell::CENTER)  << (use_percentage?"Percentage In Bucket":"Count In Bucket") << "|" << stk::end_header;
        table << "|" << justify(PrintTable::Cell::CENTER) << "       " <<  "|" << justify(PrintTable::Cell::CENTER)  << ostr0.str() << "|" << stk::end_header;
        if (use_percentage)
          {
            std::string line_equal(m_max_column_width, '=');
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
            std::string count(end_loop, '=');
            table << justify(PrintTable::Cell::LEFT) <<  count << "|" << stk::end_row;
          }

        os << "\n" << table;

      }

    };

  }
}

#endif

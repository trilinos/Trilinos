// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#pragma once
#ifndef ROL2_TABLE_HPP
#define ROL2_TABLE_HPP

namespace ROL2 {

class Table {
public:

  Table( std::string                        name_in,
         std::initializer_list<std::string> column_names_in,
         bool use_hline_in = true,
         int precision = 6 ) :
    name_(name_in),
    column_names_(column_names_in),
    vecstring_(precision),
    use_hline_(use_hline_in) {
    for(auto e : column_names_) column_widths_.emplace_back(e.size());
  }

  template<typename...Columns>
  void add_row( Columns&&...cols ) {
    auto row = vecstring_( std::forward<Columns>(cols)... );
    for( std::size_t k=0u; k<column_names_.size(); ++k )
      if( static_cast<int>(row[k].size()) > column_widths_[k] )
        column_widths_[k] = static_cast<int>(row[k].size());
    values_.emplace_back(row);
  }

  template<typename T>
  void add_row( const std::vector<T>& vec ) {
     for( auto& elem : vec ) values_.emplace_back(vecstring_(elem));
  }

  void clear_values() {
    values_.clear();
    for(auto e : column_names_) column_widths_.emplace_back(e.size());
  }

  friend std::ostream& operator << ( std::ostream&, const Table& );

private:

  class VectorString {
  public:
    VectorString( int p=4 ) : precision_(p) {}

    template<typename T>
    std::vector<std::string> operator() ( const std::vector<T>& vec ) const {
      std::vector<std::string> result;
      for( auto& e : vec ) result.emplace_back(as_string(e));
      return result;
    }

    template<typename...Args>
    std::vector<std::string> operator() ( Args&&...args ) const {
      std::vector<std::string> result;
      int unpack[]{0, (result.emplace_back(as_string(std::forward<Args>(args))), 0)...};
      static_cast<void>(unpack);
      return result;
    }

  private:
    template<typename T>
    std::string as_string( T&& t, std::true_type ) const {
      std::stringstream s;
      std::string buf;
      auto value = std::forward<T>(t);
      if( value < 0 ) {
        value *= -1; buf = "-";
      }
      else buf = " ";
      s << std::setprecision(precision_) << std::scientific << value;
      return buf+s.str();
    }

    template<typename T>
    std::string as_string( T&& t, std::false_type ) const {
      std::stringstream s;
      s << std::forward<T>(t);
      return s.str();
    }

    template<typename T>
    std::string as_string( T&& t ) const {
      return as_string(std::forward<T>(t), std::is_signed<T>{} );
    }

    inline std::string as_string( bool b ) { return b ? std::string("true") :
      std::string("false"); }

    int precision_ = 4;

  }; // class VectorString

  VectorString                          vecstring_;
  std::string                           name_;
  std::string                           hline_;
  std::vector<std::string>              column_names_;
  std::vector<std::vector<std::string>> values_;
  std::vector<int>                      column_widths_;
  int gap_width_  = 3;
  bool use_hline_ = true;

}; // Table


inline std::ostream& operator << ( std::ostream& os, const Table& table ) {
  os << table.name_ << std::endl;
  os << std::string(table.gap_width_,' ');
  for( std::size_t k=0u; k<table.column_names_.size(); ++k ) {
    int cw = table.column_widths_[k];
    int hw = (cw-table.column_names_[k].size())/2;
    os << std::setw(cw+table.gap_width_) << std::left
       << (std::string(hw,' ')+table.column_names_[k]);
  }
  os << std::endl;

  if( table.use_hline_ ) {
    for( auto k : table.column_widths_ ) {
      os << std::setw(table.gap_width_+k) << std::left << std::string(table.gap_width_,' ') + std::string(k,'-');
    }
    os << std::endl;
  }
  for( const auto& row : table.values_ ) {
    os << std::string(table.gap_width_,' ');
    for( std::size_t k=0u; k<row.size(); ++k )
      os << std::setw(table.gap_width_+table.column_widths_[k]) << std::left << row[k];
    os << std::endl;
  }
  return os;
}

} // namespace ROL2


#endif //ROL2_TABLE_HPP


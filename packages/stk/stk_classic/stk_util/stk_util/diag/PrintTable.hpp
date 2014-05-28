/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_UTIL_DIAG_PrintTable_hpp
#define STK_UTIL_DIAG_PrintTable_hpp

#include <vector>
#include <string>
#include <sstream>

#include <stk_util/diag/Writer_fwd.hpp>

namespace stk {

class PrintTable
{
  template<typename T>
  friend PrintTable &operator<<(PrintTable &table, const T &t);

public:
  typedef std::string::size_type ColumnWidth;
  typedef std::vector<ColumnWidth> ColumnWidthVector;

  struct Cell
  {
    enum Flags {
      SPAN = 0x01
    };

    enum Justification {
      LEFT		= 1,
      RIGHT		= 2,
      CENTER		= 3,
      JUSTIFY_MASK	= 0x0F,
      TRUNC		= 0x10,
      ENDS		= 0x20
    };

    Cell()
      : m_string(),
	m_flags(0),
	m_justification(RIGHT | TRUNC),
	m_indent(0),
	m_width(0)
    {}

    Cell(const Cell &cell)
      : m_string(cell.m_string),
 	m_flags(cell.m_flags),
 	m_justification(cell.m_justification),
	m_indent(cell.m_indent),
	m_width(cell.m_width)
    {}

    Cell &operator=(const Cell &cell) {
      m_string = cell.m_string;
      m_flags = cell.m_flags;
      m_justification = cell.m_justification;
      m_indent = cell.m_indent;
      m_width = cell.m_width;

      return *this;
    }
        
    std::string			m_string;
    int				m_flags;
    int				m_justification;
    ColumnWidth	                m_indent;
    ColumnWidth	                m_width;
  };

  typedef std::vector<Cell> Row;
  typedef std::vector<Row> Table;

  enum Flags {
    AUTO_END_COL                = 0x01,
    COMMA_SEPARATED_VALUES      = 0x02,
    PRINT_TRANSPOSED            = 0x04
  };

  PrintTable()
    : m_ostream(0),
      m_flags(AUTO_END_COL),
      m_tableWidth(0)
  {
    m_table.push_back(Row());
  }

  explicit PrintTable(std::ostream &os)
    : m_ostream(&os),
      m_flags(AUTO_END_COL),
      m_commentPrefix(),
      m_tableWidth(0)
  {
    m_table.push_back(Row());
  }

private:
  PrintTable(const PrintTable &);
  PrintTable &operator=(const PrintTable &);

public:
  ~PrintTable()
  {}
  
  Row::size_type headerSize() const {
    return m_header.empty() ? 0 : m_header.begin()->size();
  }

  inline Table::size_type size() const {
    return m_table.size();
  }

  inline std::ostringstream &getCurrentString() {
    return m_currentString;
  }

  inline bool autoEndCol() const {
    return m_flags & AUTO_END_COL;
  }

  inline PrintTable &setAutoEndCol(bool auto_end_col = true) {
    if (auto_end_col)
      m_flags |= AUTO_END_COL;
    else
      m_flags &= ~AUTO_END_COL;
    return *this;
  }

  inline bool commaSeparatedValues() const {
    return m_flags & COMMA_SEPARATED_VALUES;
  }

  inline PrintTable &setCommaSeparatedValues(bool comma_separated_values = true) {
    if (comma_separated_values)
      m_flags |= COMMA_SEPARATED_VALUES;
    else
      m_flags &= ~COMMA_SEPARATED_VALUES;
    return *this;
  }

  inline PrintTable &setCommentPrefix(const std::string &comment_prefix) {
    m_commentPrefix = comment_prefix;
    return *this;
  }

  inline const std::string &getCommentPrefix() const {
    return m_commentPrefix;
  }

  /**
   *  This function sets the title to the given new title. The title can also
   *  be accessed directly.
   */
  inline PrintTable &setTitle(const std::string &title){
    m_title = title;
    return *this;
  }

  inline const std::string &getTitle() const {
    return m_title;
  }

  /**
   * Member function <b>operator&lt;&lt;</b> is the manipulator instantiation
   * function
   *
   * @return			a <b>PrintTable</b> reference to this object
   */
  inline PrintTable& operator<<(PrintTable& (*f)(PrintTable&)) {
    f(*this);
    return *this;

  }

  /**
   * Member function <b>operator&lt;&lt;</b> passes the ios_base manipulator to the
   * output stream.
   *
   * @return			a <b>PrintTable</b> reference to this object
   */
  inline PrintTable& operator<<(std::ios_base& (*f)(std::ios_base&)) {
    f(m_currentString);
    return *this;
  }

  inline PrintTable &push() {
    ++m_currentCell.m_indent;
    return *this;
  }

  inline PrintTable &span() {
    m_currentCell.m_flags |= Cell::SPAN;
    return *this;
  }

  inline PrintTable &pop() {
    if (m_currentCell.m_indent != 0)
      --m_currentCell.m_indent;

    return *this;
  }

  inline PrintTable &cell_width(ColumnWidth my_width) {
    m_currentCell.m_width = my_width;

    return *this;
  }

  inline PrintTable &indent(ColumnWidth my_indent) {
    m_currentCell.m_indent = my_indent;

    return *this;
  }

  inline PrintTable &justify(int justification) {
    m_currentCell.m_justification = justification;

    return *this;
  }

  PrintTable &end_col();

  PrintTable &end_row();

  PrintTable &at(size_t row, size_t col);

  PrintTable &end_header() {
    m_header.push_back(m_table.back());
    m_table.pop_back();
    m_table.push_back(Row());
    return *this;
  }

  PrintTable &end_format() {
    m_format = m_table.back();
    m_table.pop_back();
    m_table.push_back(Row());
    return *this;
  }

  void calculate_column_widths() const;

  void transpose_table() const;

  /**
   *  This function prints out the table to it's PrintTable
   */
  std::ostream &print(std::ostream &os) const;

  /**
   *  This function prints out the table to it's PrintTable
   */
  std::ostream &printRow(std::ostream &os, const Row &row) const;

  /**
   *  This function prints out the table to it's PrintTable
   */
  std::ostream &printHeaderBar(std::ostream &os) const;

  /**
   *  This function prints out the table to it's PrintTable
   */
  std::ostream &csvPrint(std::ostream &os) const;

  /**
   *  This function prints out the table to it's PrintTable
   */
  diag::Writer & verbose_print(diag::Writer &dout) const;

private:
  /**
   * Member function <b>normalize_table</b> makes sure that the table has a field at
   * <b>row</b> and <b>column</b>.
   *
   * @param row		an <b>int</b> value of the row to ensure existence.
   * @param col		an <b>int</b> value of the column to ensure existence.
   */
  void normalize_table(int row, int col);

private:
  std::ostream *		m_ostream;
  std::string			m_title;
  Table				m_header;
  Row				m_format;
  Table				m_table;
  Cell				m_currentCell;
  std::ostringstream		m_currentString;
  int				m_flags;
  std::string                   m_commentPrefix;
  mutable ColumnWidthVector	m_columnWidth;
  mutable ColumnWidth	        m_tableWidth;
};


template<typename T>
inline PrintTable &operator<<(PrintTable &table, const T &t) {
  table.m_currentString << t;
  if (table.autoEndCol())
    table.end_col();

  return table;
}


struct cell_width
{
  cell_width(PrintTable::ColumnWidth width)
    : m_width(width)
  {}

  PrintTable::ColumnWidth	m_width;
};


struct at
{
  at(size_t row, size_t col)
    : m_row(row),
      m_col(col)
  {}

  size_t                        m_row;
  size_t                        m_col; 
};


struct indent
{
  indent(PrintTable::ColumnWidth my_indent)
    : m_indent(my_indent)
  {}

  PrintTable::ColumnWidth	m_indent;
};


struct justify
{
  justify(int my_justify)
    : m_justify(my_justify)
  {}

  int		m_justify;
};


inline PrintTable &operator<<(PrintTable &tout, const at &m) {
  tout.at(m.m_row, m.m_col);
  return tout;
}

inline PrintTable &operator<<(PrintTable &tout, const cell_width &m) {
  tout.cell_width(m.m_width);
  return tout;
}

inline PrintTable &operator<<(PrintTable &tout, const indent &m) {
  tout.indent(m.m_indent);
  return tout;
}

inline PrintTable &operator<<(PrintTable &tout, const justify &m) {
  tout.justify(m.m_justify);
  return tout;
}

inline PrintTable &end_col(PrintTable &tout) {
  return tout.end_col();
}

inline PrintTable &end_row(PrintTable &tout) {
  return tout.end_row();
}

inline PrintTable &end_header(PrintTable &tout) {
  return tout.end_header();
}

inline PrintTable &end_format(PrintTable &tout) {
  return tout.end_format();
}

inline PrintTable &push(PrintTable &tout) {
  return tout.push();
}

inline PrintTable &pop(PrintTable &tout) {
  return tout.pop();
}

inline PrintTable &span(PrintTable &tout) {
  return tout.span();
}

inline std::ostream &operator<<(std::ostream &os, const PrintTable &table){
  return table.print(os);
}

inline diag::Writer &operator<<(diag::Writer &dout, const PrintTable &table){
  return table.verbose_print(dout);
}

} // namespace stk

//namespace sierra {
//
////typedef stk::PrintTable PrintTable;
//using stk::PrintTable;
//typedef stk::cell_width cell_width;
//typedef stk::at at;
//typedef stk::justify justify;
//
//} // namespace sierra

#endif // STK_UTIL_DIAG_PrintTable_hpp

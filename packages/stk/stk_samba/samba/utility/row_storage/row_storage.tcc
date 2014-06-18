#ifndef SAMBA_SAMBA_UTILITY_ROW_STORAGE_ROW_STORAGE_TCC
#define SAMBA_SAMBA_UTILITY_ROW_STORAGE_ROW_STORAGE_TCC

#include <samba/utility/debug_message.hpp>

#include <boost/assert.hpp>

namespace samba { namespace utility {

template <typename T>
inline
row_storage<T>::row_storage( size_t arg_num_rows
                             ,size_t arg_colunm_size_hint
                             ,value_type const& arg_default
                           )
  : m_column_size_hint(arg_colunm_size_hint)
  , m_num_rows(arg_num_rows)
  , m_compressed(false)
  , m_default(arg_default)
  , m_row_begin(arg_num_rows)
  , m_row_end(arg_num_rows)
  , m_column_data(arg_num_rows)
  , m_compressed_column_data()
{}

template <typename T>
inline
row_storage<T>::row_storage( row_storage<T> const& b)
  : m_column_size_hint(b.m_colunm_size_hint)
  , m_num_rows(b.m_num_rows)
  , m_compressed(b.m_compressed)
  , m_default(b.m_default)
  , m_row_begin(m_num_rows)
  , m_row_end(m_num_rows)
  , m_column_data()
  , m_compressed_column_data(b.m_compressed_column_data)
{
  //b is compressed, setup begin/end
  if ( m_compressed)  {
    size_t offset = 0;
    for (size_t i=0; i<m_num_rows; ++i) {
      m_row_begin[i] = &m_compressed_column_data[i+offset];
      m_row_end[i] = m_row_begin[i] + b.num_columns(i);
      offset += b.num_columns(i);
    }
  }
  else {
    compress_helper(b.m_column_data);
  }
}

template <typename T>
inline
row_storage<T> & row_storage<T>::operator=(row_storage<T> const& b)
{
  if (this != &b)
  {
    m_column_size_hint = b.m_column_size_hint;
    m_num_rows = b.m_num_rows;
    m_compressed = b.m_compressed;
    m_default =b.m_default;

    if (b.m_compressed)
    { //clear out m_column_data
      column_data_type temp;
      m_column_data.swap(temp);

      m_compressed_column_data = b.m_compressed_column_data;
      size_t offset = 0;
      for (size_t i=0; i<m_num_rows; ++i) {
        m_row_begin[i] = &m_compressed_column_data[i+offset];
        m_row_end[i] = m_row_begin[i] + b.num_columns(i);
        offset += b.num_columns(i);
      }
    }
    else {
      compress_helper(b.m_column_data);
    }
  }
  return this;
}

template <typename T>
inline
typename row_storage<T>::iterator row_storage<T>::begin(size_t row)
{
  BOOST_ASSERT_MSG(row < m_row_begin.size(),
                   (debug_message() << "Row " << row << " is out-of-bounds, max is " << m_row_begin.size()));
  return m_row_begin[row];
}

template <typename T>
inline
typename row_storage<T>::const_iterator row_storage<T>::begin(size_t row) const
{
  BOOST_ASSERT_MSG(row < m_row_begin.size(),
                   (debug_message() << "Row " << row << " is out-of-bounds, max is " << m_row_begin.size()));
  return m_row_begin[row];
}

template <typename T>
inline
typename row_storage<T>::iterator row_storage<T>::end(size_t row)
{
  BOOST_ASSERT_MSG(row < m_row_end.size(),
                   (debug_message() << "Row " << row << " is out-of-bounds, max is " << m_row_end.size()));
  return m_row_end[row];
}

template <typename T>
inline
typename row_storage<T>::const_iterator row_storage<T>::end(size_t row) const
{
  BOOST_ASSERT_MSG(row < m_row_end.size(),
                   (debug_message() << "Row " << row << " is out-of-bounds, max is " << m_row_end.size()));
  return m_row_end[row];
}

template <typename T>
inline
typename row_storage<T>::iterator row_storage<T>::operator[](size_t row)
{ return begin(row); }

template <typename T>
inline
typename row_storage<T>::const_iterator row_storage<T>::operator[](size_t row) const
{ return begin(row); }

template <typename T>
inline
typename row_storage<T>::iterator row_storage<T>::absolute_col(size_t col)
{
  BOOST_ASSERT_MSG(m_compressed, "This query only valid in compressed state");
  return &m_compressed_column_data[0] + col;
}

template <typename T>
inline
typename row_storage<T>::const_iterator row_storage<T>::absolute_col(size_t col) const
{
  BOOST_ASSERT_MSG(m_compressed, "This query only valid in compressed state");
  return &m_compressed_column_data[0] + col;
}

//resize row
template <typename T>
inline void row_storage<T>::resize_rows(size_t arg_num_rows)
{
  //no-op if resizing to same number of rows
  if (arg_num_rows == m_num_rows) return;

  inflate();

  size_t orig_num_rows = m_num_rows;
  m_num_rows = arg_num_rows;

  m_column_data.resize(m_num_rows);

  if (m_num_rows > orig_num_rows) {
    // growing
    m_row_begin.resize(m_num_rows, NULL);
    m_row_end.resize(m_num_rows,   NULL);

    for (size_t i=0; i<m_num_rows; ++i)
    {
      m_column_data[i].reserve(m_column_size_hint);
      m_row_begin[i] = &m_column_data[i][0];
      m_row_end[i]   = m_row_begin[i] + m_column_data[i].size();
    }
  }
  else {
    // shrinking
    m_row_begin.resize(m_num_rows, NULL);
    m_row_end.resize(m_num_rows,   NULL);
  }
}

//resize columns for row
template <typename T>
inline void row_storage<T>::resize_columns( size_t row
                                            ,size_t arg_num_columns
                                           )
{
  BOOST_ASSERT_MSG(row < num_rows(),
                   (debug_message() << "Row " << row << " is out-of-bounds, max is " << num_rows()));

  //no-op if resizing to same number of columns
  if ( num_columns(row) == arg_num_columns) return;

  inflate();

  m_column_data[row].resize(arg_num_columns,m_default);

  m_row_begin[row] = &m_column_data[row][0];
  m_row_end[row]   = m_row_begin[row] + m_column_data[row].size();
}

//insert column
template <typename T>
inline void row_storage<T>::insert_column( size_t row
                                          ,size_t column
                                         )
{ insert_column(row,column,m_default); }

//insert column
template <typename T>
inline void row_storage<T>::insert_column( size_t row
                                          ,size_t column
                                          ,value_type const& value
                                         )
{
  if (num_columns(row) <= column) {
    resize_columns(row,column+1);
  }

  begin(row)[column] = value;
}

//insert columns
template <typename T>
template <typename ColumnIterator>
inline void row_storage<T>::insert_columns( size_t row
                                           ,ColumnIterator first
                                           ,ColumnIterator last
                                          )
{
  size_t num_column = last -first;

  if (num_columns(row) != num_column) {
    resize_columns(row,num_column);
  }

  std::copy(first,last,begin(row));
}

template <typename T>
inline void row_storage<T>::insert_new_column( size_t row
                                              ,size_t column
                                              ,value_type const& value
                                             )
{
  insert_new_columns(row, column, &value, &value + 1);
}

template <typename T>
inline void row_storage<T>::insert_new_columns( size_t row
                                               ,size_t starting_at_column
                                               ,size_t num_to_insert
                                               ,value_type const& value
                                              )
{
  // kind of inefficient...
  std::vector<value_type> temp_vec(num_to_insert, value);
  insert_new_columns(row, starting_at_column, temp_vec.begin(), temp_vec.end());
}

template <typename T>
template <typename ColumnIterator>
inline void row_storage<T>::insert_new_columns( size_t row
                                               ,size_t starting_at_column
                                               ,ColumnIterator first
                                               ,ColumnIterator last
                                              )
{
  BOOST_ASSERT_MSG(row < num_rows(),
                   (debug_message() << "Row " << row << " is out-of-bounds, max is " << num_rows()));
  BOOST_ASSERT_MSG(starting_at_column <= num_columns(row),
                   (debug_message() << "Column " << starting_at_column << " is out-of-bounds, max is " << num_columns(row)));

  const size_t num_column = last -first;
  if (num_column == 0) return;

  inflate();

  // Control the growth of the vector. This should be a no-op if the col-hint
  // has already caused us to reserve the necessary space
  m_column_data[row].reserve(m_column_data[row].size() + num_column);

  m_column_data[row].insert(m_column_data[row].begin() + starting_at_column, first, last);

  m_row_begin[row] = &m_column_data[row][0];
  m_row_end[row]   = m_row_begin[row] + m_column_data[row].size();
}

template <typename T>
inline
void row_storage<T>::push_back_column( size_t row
                                       ,value_type const& value
                                       )
{
  insert_column(row, num_columns(row), value);
}

template <typename T>
inline
void row_storage<T>::erase_column( size_t row
                                  ,size_t column
                                  )
{
  BOOST_ASSERT_MSG(row < num_rows(),
                   (debug_message() << "Row " << row << " is out-of-bounds, max is " << num_rows()));

  if (column < num_columns(row)) {
    inflate();
    m_column_data[row].erase(m_column_data[row].begin() + column);
    --m_row_end[row];
  }
}

template <typename T>
inline
void row_storage<T>::swap(row_storage & b)
{
  row_storage & a = *this;
  std::swap(a.m_column_size_hint, b.m_column_size_hint);
  std::swap(a.m_num_rows, b.m_num_rows);
  std::swap(a.m_compressed, b.m_compressed);
  std::swap(a.m_default, b.m_default);
  a.m_row_begin.swap(b.m_row_begin);
  a.m_row_end.swap(b.m_row_end);
  a.m_column_data.swap(b.m_column_data);
  a.m_compressed_column_data.swap(b.m_compressed_column_data);
}

//inflate
template <typename T>
inline void row_storage<T>::inflate()
{
  if (!m_compressed) return;
  m_compressed = false;

  m_column_data.resize(m_num_rows);

  for (size_t i=0; i<m_num_rows; ++i) {
    m_column_data[i].insert(m_column_data[i].end(), begin(i), end(i));
    m_row_begin[i] = &m_column_data[i][0];
    m_row_end[i]   = m_row_begin[i] + m_column_data[i].size();
  }

  //clear out m_compressed_column_data
  compressed_column_data_type temp;
  m_compressed_column_data.swap(temp);
}

//compress
template <typename T>
inline void row_storage<T>::compress_helper(column_data_type const& column_data)
{
  if (m_compressed) return;
  m_compressed = true;

  m_num_rows = column_data.size();

  //resize compressed column_data
  {
    size_t resize_size = 0;
    for(size_t i=0; i<m_num_rows; ++i)
      resize_size += column_data[i].size();

    m_compressed_column_data.resize(resize_size);
  }

  //copy row's column to new vector
  size_t offset = 0;
  for(size_t i=0; i<m_num_rows; ++i)
  {
    m_row_begin[i] = &m_compressed_column_data[offset];
    m_row_end[i]   = m_row_begin[i] + column_data[i].size();

    offset += column_data[i].size();

    std::copy( column_data[i].begin()
              ,column_data[i].end()
              ,m_row_begin[i]
             );
  }

  //clear out m_column_data
  column_data_type temp;
  m_column_data.swap(temp);
}

}} //namespace samba::utility

#endif //SAMBA_SAMBA_UTILITY_ROW_STORAGE_ROW_STORAGE_TCC

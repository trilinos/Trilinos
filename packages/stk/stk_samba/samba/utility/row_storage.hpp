#ifndef SAMBA_SAMBA_UTILITY_ROW_STORAGE_HPP
#define SAMBA_SAMBA_UTILITY_ROW_STORAGE_HPP

#include <vector>
#include <algorithm>

namespace samba { namespace utility {

// Container to flatten storage when a data structure
// similar to a vector< vector<T> > is needed
//
// eg given the following vectors
//
//  -----
// |x|x|x|
//  -----
//  -----------
// |y|y|y|y|y|y|
//  -----------
//  -
// |z|
//  -
//
//  which will be compressed to
//
//  --------------------
// |x|x|x|y|y|y|y|y|y||z|
//  --------------------
//
// row_storage can either be in a compressed state or inflated state, not both.
//
// When making modifications, storage will automatically be switched to inflated and
// will remain inflated until user requests a compression.
template <typename T>
class row_storage
{
  typedef std::vector<std::vector<T> > column_data_type;
  typedef std::vector<T>  compressed_column_data_type;

  public:
    typedef T value_type;
    typedef size_t size_type;

    typedef T* iterator;
    typedef const T* const_iterator;

    row_storage( size_t arg_num_rows = 0
                ,size_t arg_colunm_size_hint = 0
                ,value_type const& arg_default = value_type()
               );

    row_storage( row_storage const& b);

    row_storage & operator=(row_storage const& b);

    ////////////////////////////
    // Accessors
    ////////////////////////////

    //begin(row)
    iterator begin(size_t row);

    //begin(row)
    const_iterator begin(size_t row) const;

    //end(row)
    iterator end(size_t row);

    //end(row)
    const_iterator end(size_t row) const;

    //operator[row]
    iterator operator[](size_t row);

    //operator[row]
    const_iterator operator[](size_t row) const;

    iterator absolute_col(size_t col);

    const_iterator absolute_col(size_t col) const;

    ////////////////////////////
    // Queries
    ////////////////////////////

    //num_rows
    size_t num_rows() const
    { return m_num_rows; }

    size_t total_column_count() const
    { return m_compressed_column_data.size(); }

    //num_columns(row)
    size_t num_columns(size_t row) const
    { return end(row)-begin(row); }

    bool is_compressed() const
    { return m_compressed; }

    ////////////////////////////
    // Modification API
    ////////////////////////////

    void resize_rows(size_t arg_num_rows);

    void resize_columns(size_t row, size_t num_columns);

    // WARNING: none of the insert_column(...) methods have typical insert behavior

    //insert the column
    void insert_column( size_t row
                       ,size_t column
                      );

    void insert_column( size_t row
                       ,size_t column
                       ,value_type const& value
                      );

    //this can be used to clear all the columns of
    template <typename ColumnIterator>
    void insert_columns( size_t row
                        ,ColumnIterator first
                        ,ColumnIterator last
                       );

    // Implements true insert behavior
    void insert_new_column( size_t row
                           ,size_t column
                           ,value_type const& value
                          );

    // Implements true insert behavior
    void insert_new_columns( size_t row
                            ,size_t starting_at_column
                            ,size_t num_to_insert
                            ,value_type const& value
                           );

    template <typename ColumnIterator>
    void insert_new_columns( size_t row
                            ,size_t starting_at_column
                            ,ColumnIterator first
                            ,ColumnIterator last
                           );

    void push_back_column( size_t row
                           ,value_type const& value
                           );

    void erase_column( size_t row
                      ,size_t column
                     );

    void compress()
    { compress_helper(m_column_data); }

    void inflate();

    void set_column_size_hint(size_t arg_colunm_size_hint)
    { m_column_size_hint = arg_colunm_size_hint; }

    size_t column_size_hint() const
    { return m_column_size_hint; }

    void swap(row_storage & b);

    friend void swap_rows( row_storage &a, size_t a_row, row_storage& b, size_t b_row )
    {
      a.inflate();
      b.inflate();
      a.m_column_data[a_row].swap(b.m_column_data[b_row]);
      std::swap(a.m_row_begin[a_row], b.m_row_begin[b_row]);
      std::swap(a.m_row_end[a_row], b.m_row_end[b_row]);
    }

    std::vector<value_type> & get_column_vector(size_t row)
    {
      inflate();
      return m_column_data[row];
    }

#ifndef PRINT_SAMBA_PARTITION_MEMORY_FOOTPRINT
  private:
#endif
    void compress_helper(column_data_type const& column_data);

    size_t  m_column_size_hint;
    size_t  m_num_rows;
    bool    m_compressed;

    value_type m_default;

    std::vector<iterator> m_row_begin;
    std::vector<iterator> m_row_end;

    column_data_type            m_column_data;
    compressed_column_data_type m_compressed_column_data;
};

}}

namespace std {

template <typename T>
inline void
swap(samba::utility::row_storage<T> & a, samba::utility::row_storage<T> & b)
{ a.swap(b); }

};

#include <samba/utility/row_storage/row_storage.tcc>

#endif //SAMBA_SAMBA_UTILITY_ROW_STORAGE_HPP

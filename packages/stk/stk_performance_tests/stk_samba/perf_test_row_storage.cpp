#define ENABLE_ROW_STORAGE_VS_VECTOR_VECTOR false

#if ENABLE_ROW_STORAGE_VS_VECTOR_VECTOR

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <samba/utility/row_storage.hpp>

#include <stk_performance_test_includes/cpu_time.hpp>

#include <vector>
#include <algorithm>
#include <numeric>
#include <cstdlib>
#include <iomanip>
#include <iostream>

namespace {

double run_vector( size_t num_rows, size_t num_columns, size_t progress_freq, size_t show_progress)
{
  //timing variables
  double vector_fill = 0;
  double vector_iterate = 0;

  std::ostream_iterator<size_t> out_size_t(std::cout, ", ");

  std::vector< std::vector<size_t> > vec(num_rows);

  std::vector<size_t> columns(num_rows);
  for (size_t i=0; i<num_rows; ++i) columns[i] = i;

  //randomize the column iteration to simulate adding
  //back connectivity in random order
  std::srand(0);
  std::random_shuffle(columns.begin(), columns.end());

  for (size_t i=0; i<num_rows; ++i) {
    double start_time = cpu_time();
    for (size_t j=0; j < num_columns; ++j) {
      vec[columns[i]].push_back(j*j+10);
    }
    double end_time = cpu_time() - start_time;
    vector_fill += end_time;
    if (show_progress && i%progress_freq == 0) {
      std::cout << "\b\b\b\b" << 100.0*i/num_rows  << "%" << std::flush;
    }
  }

  if(show_progress) {
    std::cout << "\b\b\b\b";
  }

  //vector iteration
  {
    for (size_t iter = 0; iter < 100; ++iter)
    {
      double start_time = cpu_time();
      for (size_t i=0; i<num_rows; ++i) {
        size_t s = std::accumulate(vec[i].begin(),vec[i].end(), 0);
        EXPECT_TRUE(s>0);
      }
      double end_time = cpu_time() - start_time;
      vector_iterate += end_time;
      if (show_progress) {
        std::cout << "\b\b\b" << iter  << "%" << std::flush;
      }
    }

    if (show_progress)
      std::cout << "\b\b\b";
  }

  double total_time = vector_fill + vector_iterate;

  if (show_progress) {
  std::cout << std::endl
    << "     Vector:" << std::endl
    << "       Fill: " << vector_fill << " seconds" << std::endl
    << "  Iteration: " << vector_iterate << " seconds" << std::endl
    << "      TOTOL: " << total_time << " seconds" << std::endl;
  }

  return total_time;
}

double run_row_store( size_t num_rows, size_t num_columns, size_t progress_freq, size_t show_progress)
{
  double row_fill = 0;
  double row_compress = 0;
  double row_iterate = 0;
  double row_inflate = 0;

  std::ostream_iterator<size_t> out_size_t(std::cout, ", ");

  samba::utility::row_storage<size_t> row_store(num_rows);

  std::vector<size_t> columns(num_rows);
  for (size_t i=0; i<num_rows; ++i) columns[i] = i;

  //randomize the column iteration to simulate adding
  //back connectivity in random order
  std::srand(0);
  std::random_shuffle(columns.begin(), columns.end());

  //bouce back and forth to simulate so that the vector memory is not alreay linearized

  for (size_t i=0; i<num_rows; ++i) {
    double start_time = cpu_time();
    for (size_t j=0; j < num_columns; ++j) {
      row_store.insert_column(columns[i],j,j*j+10);
    }
    double end_time = cpu_time() - start_time;
    row_fill += end_time;

    if (show_progress && i%progress_freq == 0) {
      std::cout << "\b\b\b\b" << 100.0*i/num_rows  << "%" << std::flush;
    }
  }

  if(show_progress) {
    std::cout << "\b\b\b\b";
  }

  {
    double start_time = cpu_time();

    row_store.compress();

    double end_time = cpu_time() - start_time;
    row_compress = end_time;
  }

  //row_store iteration
  {
    for (size_t iter = 0; iter < 100; ++iter)
    {
      double start_time = cpu_time();
      for (size_t i=0; i<num_rows; ++i) {
        size_t s = std::accumulate(row_store.begin(i),row_store.end(i), 0);
        EXPECT_TRUE(s>0);
      }
      double end_time = cpu_time() - start_time;
      row_iterate += end_time;
      if (show_progress) {
        std::cout << "\b\b\b" << iter  << "%" << std::flush;
      }
    }

    if (show_progress)
      std::cout << "\b\b\b    ";
  }

  {
    double start_time = cpu_time();

    row_store.inflate();

    double end_time = cpu_time() - start_time;
    row_inflate = end_time;
  }

  double total_time = row_fill + row_compress + row_iterate + row_inflate;

  if (show_progress) {
    std::cout << std::endl
      << "  Row Store:" << std::endl
      << "       Fill: " << row_fill << " seconds" << std::endl
      << "   Compress: " << row_compress << " seconds" << std::endl
      << "  Iteration: " << row_iterate << " seconds" << std::endl
      << "    Inflate: " << row_inflate << " seconds" << std::endl
      << "      TOTAL: " << total_time << " seconds" << std::endl;
  }

  return total_time;
}

struct times{
  size_t row;
  size_t column;
  double vector_time;
  double row_store_time;
};

} // namespace

STKUNIT_UNIT_TEST(samba, row_storage_vs_vector_vector)
{
  std::cout << std::setprecision(3);

  size_t min_row =   50000; //50 thousand
  size_t max_row = 5000000; //5 million

  size_t min_col = 5;
  size_t max_col = 500;

  std::vector<times> timings;

  bool show_progress = true;

  for (size_t row = min_row; row <= max_row; row *= 10) {
    for (size_t column = min_col; column <= max_col; column *= 10) {
      //don't run largest case since it swaps to disc
      if (row == max_row && column == max_col) break;

      times t;

      t.row = row;
      t.column = column;

      size_t progress_freq = row/100;

      std::cout << std::endl
                << "Rows: " << row << std::endl
                << "Columns: " << column << std::endl
                << std::endl;

      t.vector_time =
        run_vector(
           row
          ,column
          ,progress_freq
          ,show_progress
          );

      std::cout << std::endl;

      t.row_store_time  =
        run_row_store(
           row
          ,column
          ,progress_freq
          ,show_progress
          );

      timings.push_back(t);

    }
  }

  std::cout << std::endl;
  std::cout << std::endl;

  size_t width = 10;

  std::cout << std::setw(width) << std::right << "Rows";
  std::cout << std::setw(width) << std::right << "Columns";
  std::cout << std::setw(width) << std::right << "Size";
  std::cout << std::setw(width) << std::right << "Speedup";
  std::cout << std::endl;

  for (std::vector<times>::const_iterator itr=timings.begin(),end_itr=timings.end();
       itr != end_itr;
       ++itr
       )
  {
    times t = *itr;
    std::cout << std::setw(width) << std::right << t.row;
    std::cout << std::setw(width) << std::right << t.column;
    std::cout << std::setw(width) << std::right << t.row*t.column;
    std::cout << std::setw(width) << std::right << t.vector_time/t.row_store_time;
    std::cout << std::endl;
  }

  std::cout << std::endl;
}

#endif

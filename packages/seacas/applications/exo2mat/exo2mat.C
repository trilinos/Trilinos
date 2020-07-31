/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/* exodus II to matlab m file, copy of
   exo2mat.
   exo2mat was written by mrtabba
   exo2m modifications by gmreese to permit usage on machines without
     the matlab libraries.

   modified by D. Todd Griffith, 01/12/2006
      to include changes made in versions 1.4 through 1.6 on the
      SEACAS tools repository as the previous modifications dated
      12/08 and 12/15/2005 were made starting with Version 1.3.
      In particular, the only changes are those made here are those
      made in Version 1.6 which include a special provision
      for exodus files which contain no distribution factors
   modified by D. Todd Griffith, 12/15/2005
      to write distribution factors as double precision type
   modified by D. Todd Griffith  12/08/2005
      to include complete writing of side set and node set information
      so it will be available for the mat2exo conversion stage

*/

#include <algorithm>
#include <cstring> // for strlen, etc
#include <iostream>
#include <numeric>
#include <vector>

#include "add_to_log.h" // for add_to_log
#include "exodusII.h"   // for ex_get_variable_param, etc
#include "fmt/chrono.h"
#include "fmt/printf.h"
#include "matio.h" // for Mat_VarCreate, Mat_VarFree, etc
#include <cassert> // for assert
#include <cstddef> // for size_t
#include <cstdio>  // for fprintf, printf, sprintf, etc
#include <cstdlib> // for free, calloc, exit
#include <ctime>
#if MATIO_VERSION < 151
#error "MatIO Version 1.5.1 or greater is required"
#endif

#define EXT ".mat"
static int textfile = 0;

static FILE * m_file   = nullptr; /* file for m file output */
static mat_t *mat_file = nullptr; /* file for binary .mat output */
static bool   debug    = false;

static const char *qainfo[] = {
    "exo2mat",
    "2019/05/18",
    "4.07",
};

std::string time_stamp(const std::string &format)
{
  if (format == "") {
    return std::string("");
  }
  auto  calendar_time = std::time(nullptr);
  auto *local_time    = std::localtime(&calendar_time);

  std::string time_string = fmt::format(format, *local_time);
  return time_string;
}

void logger(const char *message)
{
  const std::string tsFormat = "[%H:%M:%S] ";
  fmt::print(std::clog, "{}: {}\n", time_stamp(tsFormat), message);
}

void usage()
{
  fmt::print("exo2mat [options] exodus_file_name.\n"
             "   the exodus_file_name is required (exodus only).\n"
             "   Options:\n"
             "   -t    write a text (.m) file rather than a binary .mat\n"
             "   -o    output file name (rather than auto generate)\n"
             "   -c    use cell arrays for transient variables.\n"
             "   -v5   output version 5 mat file\n"
             "   -v73  output version 7.3 mat file (hdf5-based) [default]\n"
             "   -v7.3 output version 7.3 mat file (hdf5-based)\n"
             " ** note **\n"
             "Binary files are written by default on all platforms.\n");
}

/* put a string into an m file. If the string has
   line feeds, we put it as ints, and use 'char()' to convert it */
void mPutStr(const std::string &name, const char *str)
{
  assert(m_file != nullptr);
  if (strchr(str, '\n') == nullptr) {
    fmt::fprintf(m_file, "%s='%s';\n", name, str);
  }
  else {
    fmt::fprintf(m_file, "%s=[", name);
    size_t i;
    size_t j;
    for (j = i = 0; i < std::strlen(str); i++, j++) {
      if (j >= 20) {
        j = 0;
        fmt::fprintf(m_file, "...\n");
      }
      fmt::fprintf(m_file, "%d ", str[i]);
    }
    fmt::fprintf(m_file, "];\n");
    fmt::fprintf(m_file, "%s=char(%s);\n", name, name);
  }
}

/* put double array in m file */
void mPutDbl(const std::string &name, int n1, int n2, double *pd)
{
  assert(m_file != nullptr);
  if (n1 == 1 && n2 == 1) {
    fmt::fprintf(m_file, "%s=%15.8e;\n", name, *pd);
    return;
  }
  fmt::fprintf(m_file, "%s=zeros(%d,%d);\n", name, n1, n2);
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      fmt::fprintf(m_file, "%s(%d,%d)=%15.8e;\n", name, i + 1, j + 1, pd[i * n2 + j]);
    }
  }
}

/* put integer array in m file */
void mPutInt(const std::string &name, int pd)
{
  assert(m_file != nullptr);
  fmt::fprintf(m_file, "%s=%d;\n", name, pd);
}

/* put integer array in m file */
void mPutInt(const std::string &name, int n1, int n2, int *pd)
{
  assert(m_file != nullptr);
  if (n1 == 1 && n2 == 1) {
    fmt::fprintf(m_file, "%s=%d;\n", name, *pd);
    return;
  }
  fmt::fprintf(m_file, "%s=zeros(%d,%d);\n", name, n1, n2);
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      fmt::fprintf(m_file, "%s(%d,%d)=%d;\n", name, i + 1, j + 1, pd[i * n2 + j]);
    }
  }
}

/* put string in mat file*/
int matPutStr(const std::string &name, char *str)
{
  int       error  = 0;
  matvar_t *matvar = nullptr;
  size_t    dims[2];

  dims[0] = 1;
  dims[1] = std::strlen(str);

  matvar = Mat_VarCreate(name.c_str(), MAT_C_CHAR, MAT_T_UINT8, 2, dims, str, MAT_F_DONT_COPY_DATA);
  if (matvar != nullptr) {
    error = Mat_VarWrite(mat_file, matvar, MAT_COMPRESSION_NONE);
    Mat_VarFree(matvar);
  }
  else {
    error = 1;
  }
  return error;
}

/* put double in mat file*/
int matPutDbl(const std::string &name, int n1, int n2, double *pd)
{
  int       error  = 0;
  matvar_t *matvar = nullptr;

  size_t dims[2];
  dims[0] = n1;
  dims[1] = n2;

  matvar =
      Mat_VarCreate(name.c_str(), MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, pd, MAT_F_DONT_COPY_DATA);
  if (matvar != nullptr) {
    error = Mat_VarWrite(mat_file, matvar, MAT_COMPRESSION_ZLIB);
    Mat_VarFree(matvar);
  }
  else {
    error = 1;
  }
  return error;
}

/* put integer in mat file*/
int matPutInt(const std::string &name, int n1, int n2, int *pd)
{
  int       error  = 0;
  matvar_t *matvar = nullptr;

  size_t dims[2];
  dims[0] = n1;
  dims[1] = n2;

  matvar = Mat_VarCreate(name.c_str(), MAT_C_INT32, MAT_T_INT32, 2, dims, pd, MAT_F_DONT_COPY_DATA);
  if (matvar != nullptr) {
    error = Mat_VarWrite(mat_file, matvar, MAT_COMPRESSION_ZLIB);
    Mat_VarFree(matvar);
  }
  else {
    error = 1;
  }
  return error;
}

/* wrappers for the output routine types */
void PutStr(const std::string &name, const char *str)
{
  if (textfile != 0) {
    mPutStr(name, str);
  }
  else {
    matPutStr(name, const_cast<char *>(str));
  }
}

int PutInt(const std::string &name, int pd)
{
  int error = 0;
  if (textfile != 0) {
    mPutInt(name, pd);
  }
  else {
    error = matPutInt(name, 1, 1, &pd);
  }
  return error;
}

int PutInt(const std::string &name, int n1, int n2, int *pd)
{
  int error = 0;
  if (textfile != 0) {
    mPutInt(name, n1, n2, pd);
  }
  else {
    error = matPutInt(name, n1, n2, pd);
  }
  return error;
}

int PutDbl(const std::string &name, int n1, int n2, double *pd)
{
  int error = 0;
  if (textfile != 0) {
    mPutDbl(name, n1, n2, pd);
  }
  else {
    error = matPutDbl(name, n1, n2, pd);
  }
  return error;
}

char **get_exodus_names(size_t count, int size)
{
  auto names = new char *[count];
  for (size_t i = 0; i < count; i++) {
    names[i] = new char[size + 1];
    std::memset(names[i], '\0', size + 1);
  }
  return names;
}

void delete_exodus_names(char **names, int count)
{
  for (int i = 0; i < count; i++) {
    delete[] names[i];
  }
  delete[] names;
}

void get_put_user_names(int exo_file, ex_entity_type type, int num_blocks, const char *mname)
{
  int max_name_length = ex_inquire_int(exo_file, EX_INQ_DB_MAX_USED_NAME_LENGTH);
  max_name_length     = max_name_length < 32 ? 32 : max_name_length;
  ex_set_max_name_length(exo_file, max_name_length);
  char **names = get_exodus_names(num_blocks, max_name_length + 1);
  ex_get_names(exo_file, type, names);

  std::string user_names;
  for (int j = 0; j < num_blocks; j++) {
    user_names += names[j];
    user_names += "\n";
  }
  PutStr(mname, user_names.c_str());
  delete_exodus_names(names, num_blocks);
}

void get_put_names(int exo_file, ex_entity_type type, int num_vars, const char *mname)
{
  int max_name_length = ex_inquire_int(exo_file, EX_INQ_DB_MAX_USED_NAME_LENGTH);
  max_name_length     = max_name_length < 32 ? 32 : max_name_length;
  char **names        = get_exodus_names(num_vars, max_name_length + 1);

  if (debug) {
    logger("\tReading variable names");
  }
  ex_get_variable_names(exo_file, type, num_vars, names);

  std::string mat;
  for (int i = 0; i < num_vars; i++) {
    mat += names[i];
    mat += "\n";
  }
  if (debug) {
    logger("\tWriting variable names");
  }
  PutStr(mname, mat.c_str());

  delete_exodus_names(names, num_vars);
}

std::vector<std::string> get_names(int exo_file, ex_entity_type type, int num_vars)
{
  int max_name_length = ex_inquire_int(exo_file, EX_INQ_DB_MAX_USED_NAME_LENGTH);
  max_name_length     = max_name_length < 32 ? 32 : max_name_length;
  char **names        = get_exodus_names(num_vars, max_name_length + 1);

  if (debug) {
    logger("\tReading variable names");
  }
  ex_get_variable_names(exo_file, type, num_vars, names);

  std::vector<std::string> mat(num_vars);
  for (int i = 0; i < num_vars; i++) {
    mat[i] = names[i];
  }
  delete_exodus_names(names, num_vars);
  return mat;
}

void get_put_vars(int exo_file, ex_entity_type type, int num_blocks, int num_vars,
                  int num_time_steps, const std::vector<int> &num_per_block,
                  const std::string &prefix, bool use_cell_arrays)

{
  /* truth table */
  if (debug) {
    logger("\tTruth Table");
  }
  std::vector<int> truth_table(num_vars * num_blocks);
  ex_get_truth_table(exo_file, type, num_blocks, num_vars, truth_table.data());

  std::vector<int> ids(num_blocks);
  ex_get_ids(exo_file, type, ids.data());

  size_t num_entity = std::accumulate(num_per_block.begin(), num_per_block.end(), 0);

  if (use_cell_arrays) {
    std::string var_name = prefix + "var";

    size_t dims[2];
    dims[0] = 2;
    dims[1] = num_vars;
    matvar_t *cell_array =
        Mat_VarCreate(var_name.c_str(), MAT_C_CELL, MAT_T_CELL, 2, dims, nullptr, 0);
    assert(cell_array);

    std::vector<double> scr(num_vars * num_time_steps * num_entity);
    dims[0]       = num_entity;
    dims[1]       = num_time_steps;
    size_t offset = 0;

    // Get vector of variable names...
    auto names = get_names(exo_file, type, num_vars);

    std::vector<matvar_t *> cell_element(num_vars * 2);

    int j = 0;
    for (int i = 0; i < num_vars; i++) {
      size_t sdims[2];
      sdims[0]        = 1;
      sdims[1]        = names[i].length();
      cell_element[j] = Mat_VarCreate(nullptr, MAT_C_CHAR, MAT_T_UINT8, 2, sdims,
                                      (void *)names[i].c_str(), MAT_F_DONT_COPY_DATA);
      Mat_VarSetCell(cell_array, j, cell_element[j]);
      j++;

      cell_element[j] = Mat_VarCreate(nullptr, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &scr[offset],
                                      MAT_F_DONT_COPY_DATA);
      assert(cell_element[j]);
      Mat_VarSetCell(cell_array, j, cell_element[j]);
      size_t n = 0;
      for (int jj = 0; jj < num_time_steps; jj++) {
        for (int k = 0; k < num_blocks; k++) {
          if (truth_table[num_vars * k + i] == 1) {
            ex_get_var(exo_file, jj + 1, type, i + 1, ids[k], num_per_block[k], &scr[n + offset]);
          }
          n += num_per_block[k];
        }
      }
      offset += num_time_steps * num_entity;
      j++;
    }
    Mat_VarWrite(mat_file, cell_array, MAT_COMPRESSION_NONE);
    Mat_VarFree(cell_array);
  }
  else {
    std::string var_name = prefix + "names";
    get_put_names(exo_file, type, num_vars, var_name.c_str());

    std::vector<double> scr(num_entity * num_time_steps);

    std::string format = prefix + "var%02d";
    for (int i = 0; i < num_vars; i++) {
      if (debug) {
        logger("\tReading");
      }
      std::fill(scr.begin(), scr.end(), 0.0);
      size_t      n   = 0;
      std::string str = fmt::sprintf(format.c_str(), i + 1);
      for (int j = 0; j < num_time_steps; j++) {
        for (int k = 0; k < num_blocks; k++) {
          if (truth_table[num_vars * k + i] == 1) {
            ex_get_var(exo_file, j + 1, type, i + 1, ids[k], num_per_block[k], &scr[n]);
          }
          n = n + num_per_block[k];
        }
      }
      if (debug) {
        logger("\tWriting");
      }
      PutDbl(str, num_entity, num_time_steps, scr.data());
    }
  }
}

std::vector<int> handle_element_blocks(int exo_file, int num_blocks, bool use_cell_arrays)
{
  std::vector<int> ids(num_blocks);
  ex_get_ids(exo_file, EX_ELEM_BLOCK, ids.data());

  std::vector<int> num_elem_in_block(num_blocks);

  // Storing:
  // 1) name
  // 2) id
  // 3) block topology type
  // 4) connectivity
  int max_name_length = ex_inquire_int(exo_file, EX_INQ_DB_MAX_USED_NAME_LENGTH);
  max_name_length     = max_name_length < 32 ? 32 : max_name_length;
  if (use_cell_arrays) {
    int    num_field = 4;
    size_t dims[2];
    dims[0] = num_field;
    dims[1] = num_blocks;
    matvar_t *cell_array =
        Mat_VarCreate("element_blocks", MAT_C_CELL, MAT_T_CELL, 2, dims, nullptr, 0);
    assert(cell_array);

    std::vector<matvar_t *> cell_element(num_blocks * num_field);

    std::vector<int> num_node_per_elem(num_blocks);

    size_t                   conn_size = 0;
    std::vector<std::string> types(num_blocks);

    for (int i = 0; i < num_blocks; i++) {
      char type[33];
      int  num_elem = 0;
      int  num_node = 0;
      int  num_attr = 0;
      ex_get_block(exo_file, EX_ELEM_BLOCK, ids[i], type, &num_elem, &num_node, nullptr, nullptr,
                   &num_attr);
      types[i]             = std::string(type);
      num_elem_in_block[i] = num_elem;
      num_node_per_elem[i] = num_node;
      conn_size += num_elem * num_node;
    }

    std::vector<int> connect(conn_size);
    size_t           conn_off = 0;
    for (int i = 0; i < num_blocks; i++) {
      std::vector<char> name(max_name_length + 1);
      ex_get_name(exo_file, EX_ELEM_BLOCK, ids[i], name.data());
      dims[0]      = 1;
      dims[1]      = std::strlen(name.data());
      size_t index = num_field * i + 0;
      cell_element[index] =
          Mat_VarCreate(nullptr, MAT_C_CHAR, MAT_T_UINT8, 2, dims, (void *)name.data(), 0);
      Mat_VarSetCell(cell_array, index, cell_element[index]);

      dims[0] = 1;
      dims[1] = 1;
      index   = num_field * i + 1;
      cell_element[index] =
          Mat_VarCreate(nullptr, MAT_C_INT32, MAT_T_INT32, 2, dims, &ids[i], MAT_F_DONT_COPY_DATA);
      Mat_VarSetCell(cell_array, index, cell_element[index]);

      dims[0] = 1;
      dims[1] = types[i].length();
      index   = num_field * i + 2;
      cell_element[index] =
          Mat_VarCreate(nullptr, MAT_C_CHAR, MAT_T_UINT8, 2, dims, (void *)types[i].c_str(), 0);
      Mat_VarSetCell(cell_array, index, cell_element[index]);

      dims[0] = num_node_per_elem[i];
      dims[1] = num_elem_in_block[i];
      index   = num_field * i + 3;
      ex_get_conn(exo_file, EX_ELEM_BLOCK, ids[i], &connect[conn_off], nullptr, nullptr);
      cell_element[index] = Mat_VarCreate(nullptr, MAT_C_INT32, MAT_T_INT32, 2, dims,
                                          &connect[conn_off], MAT_F_DONT_COPY_DATA);
      assert(cell_element[index]);
      Mat_VarSetCell(cell_array, index, cell_element[index]);

      conn_off += num_node_per_elem[i] * num_elem_in_block[i];
    }
    Mat_VarWrite(mat_file, cell_array, MAT_COMPRESSION_NONE);
    Mat_VarFree(cell_array);
  }
  else {
    std::vector<int>    connect;
    std::vector<double> attr;

    PutInt("blkids", num_blocks, 1, ids.data());
    std::vector<char> type(max_name_length + 1);
    std::string       types;
    for (int i = 0; i < num_blocks; i++) {
      int num_elem = 0;
      int num_node = 0;
      int num_attr = 0;
      ex_get_block(exo_file, EX_ELEM_BLOCK, ids[i], type.data(), &num_elem, &num_node, nullptr,
                   nullptr, &num_attr);
      types += type.data();
      types += "\n";
      num_elem_in_block[i] = num_elem;
      connect.resize(num_elem * num_node);
      ex_get_conn(exo_file, EX_ELEM_BLOCK, ids[i], connect.data(), nullptr, nullptr);
      std::string str = fmt::sprintf("blk%02d", i + 1);
      PutInt(str, num_node, num_elem, connect.data());

      // Handle block attributes (if any...)
      attr.resize(num_elem);
      str = fmt::sprintf("blk%02d_nattr", i + 1);
      PutInt(str, num_attr);
      if (num_attr > 0) {
        std::string attr_names;
        char **     names = get_exodus_names(num_attr, max_name_length + 1);
        ex_get_attr_names(exo_file, EX_ELEM_BLOCK, ids[i], names);
        for (int j = 0; j < num_attr; j++) {
          attr_names += names[j];
          attr_names += "\n";
        }
        str = fmt::sprintf("blk%02d_attrnames", i + 1);
        PutStr(str, attr_names.c_str());
        delete_exodus_names(names, num_attr);

        for (int j = 0; j < num_attr; j++) {
          str = fmt::sprintf("blk%02d_attr%02d", i + 1, j + 1);
          ex_get_one_attr(exo_file, EX_ELEM_BLOCK, ids[i], j + 1, attr.data());
          PutDbl(str, num_elem, 1, attr.data());
        }
      }
    }

    get_put_user_names(exo_file, EX_ELEM_BLOCK, num_blocks, "blkusernames");
    PutStr("blknames", types.c_str());
  }
  return num_elem_in_block;
}

std::vector<int> handle_node_sets(int exo_file, int num_sets, bool use_cell_arrays)
{
  std::vector<int> num_nodes(num_sets);
  if (num_sets > 0) {
    if (debug) {
      logger("Node Sets");
    }
    std::vector<int> ids(num_sets);
    ex_get_ids(exo_file, EX_NODE_SET, ids.data());

    size_t           tot_nodes = 0;
    size_t           tot_dfac  = 0;
    std::vector<int> num_df(num_sets);
    for (int i = 0; i < num_sets; i++) {
      int n1;
      int n2;
      ex_get_set_param(exo_file, EX_NODE_SET, ids[i], &n1, &n2);
      num_nodes[i] = n1;
      num_df[i]    = n2;
      tot_nodes += n1;
      tot_dfac += n2;
    }

    // Storing:
    // 1) name
    // 2) id
    // 3) node list
    // 4) distribution factors
    if (use_cell_arrays) {
      size_t dims[2];
      dims[0] = 4;
      dims[1] = num_sets;
      matvar_t *cell_array =
          Mat_VarCreate("node_sets", MAT_C_CELL, MAT_T_CELL, 2, dims, nullptr, 0);
      assert(cell_array);

      std::vector<matvar_t *> cell_element(num_sets * 4);

      std::vector<int>    node_list(tot_nodes);
      std::vector<double> dist_fac(tot_dfac);
      size_t              nl_off = 0;
      size_t              df_off = 0;

      int max_name_length = ex_inquire_int(exo_file, EX_INQ_DB_MAX_USED_NAME_LENGTH);
      max_name_length     = max_name_length < 32 ? 32 : max_name_length;

      for (int i = 0; i < num_sets; i++) {
        std::vector<char> name(max_name_length + 1);
        ex_get_name(exo_file, EX_NODE_SET, ids[i], name.data());
        dims[0]      = 1;
        dims[1]      = std::strlen(name.data());
        size_t index = 4 * i + 0;
        cell_element[index] =
            Mat_VarCreate(nullptr, MAT_C_CHAR, MAT_T_UINT8, 2, dims, (void *)name.data(), 0);
        Mat_VarSetCell(cell_array, index, cell_element[index]);

        dims[0]             = 1;
        dims[1]             = 1;
        index               = 4 * i + 1;
        cell_element[index] = Mat_VarCreate(nullptr, MAT_C_INT32, MAT_T_INT32, 2, dims, &ids[i],
                                            MAT_F_DONT_COPY_DATA);
        Mat_VarSetCell(cell_array, index, cell_element[index]);

        dims[0] = num_nodes[i];
        dims[1] = 1;
        index   = 4 * i + 2;
        ex_get_set(exo_file, EX_NODE_SET, ids[i], &node_list[nl_off], nullptr);
        /* nodes list */
        cell_element[index] = Mat_VarCreate(nullptr, MAT_C_INT32, MAT_T_INT32, 2, dims,
                                            &node_list[nl_off], MAT_F_DONT_COPY_DATA);
        assert(cell_element[index]);
        Mat_VarSetCell(cell_array, index, cell_element[index]);

        /* distribution-factors list */
        ex_get_set_dist_fact(exo_file, EX_NODE_SET, ids[i], &dist_fac[df_off]);
        index               = 4 * i + 3;
        cell_element[index] = Mat_VarCreate(nullptr, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims,
                                            &dist_fac[df_off], MAT_F_DONT_COPY_DATA);
        assert(cell_element[index]);
        Mat_VarSetCell(cell_array, index, cell_element[index]);

        nl_off += num_nodes[i];
        df_off += num_df[i];
      }
      Mat_VarWrite(mat_file, cell_array, MAT_COMPRESSION_NONE);
      Mat_VarFree(cell_array);
    }
    else {
      PutInt("nsids", num_sets, 1, ids.data());

      for (int i = 0; i < num_sets; i++) {
        std::vector<int> node_list(num_nodes[i]);
        ex_get_set(exo_file, EX_NODE_SET, ids[i], node_list.data(), nullptr);
        /* nodes list */
        std::string str;
        str = fmt::sprintf("nsnod%02d", i + 1);
        PutInt(str, node_list.size(), 1, node_list.data());

        /* distribution-factors list */
        if (num_df[i] > 0) {
          std::vector<double> dist_fac(num_df[i]);
          ex_get_set_dist_fact(exo_file, EX_NODE_SET, ids[i], dist_fac.data());
          str = fmt::sprintf("nsfac%02d", i + 1);
          PutDbl(str, dist_fac.size(), 1, dist_fac.data());
        }
      }
    }

    get_put_user_names(exo_file, EX_NODE_SET, num_sets, "nsusernames");

    /* Store # nodes and # dis. factors per node set */
    PutInt("nnsnodes", num_sets, 1, num_nodes.data());
    PutInt("nnsdfac", num_sets, 1, num_df.data());
  }
  return num_nodes;
}

std::vector<int> handle_side_sets(int exo_file, int num_sets, bool use_cell_arrays)
{
  std::vector<int> num_sideset_sides(num_sets);
  std::vector<int> num_sideset_dfac(num_sets);
  std::vector<int> num_sideset_nodes(num_sets);
  if (num_sets > 0) {
    std::vector<int> ids(num_sets);
    ex_get_ids(exo_file, EX_SIDE_SET, ids.data());

    // Storing:
    // 1) name
    // 2) id
    // 3) element list
    // 4) side list
    // 5) node count per face
    // 6) face node list
    // 7) distribution factors

    if (use_cell_arrays) {
      size_t dims[2];
      dims[0] = 7;
      dims[1] = num_sets;
      matvar_t *cell_array =
          Mat_VarCreate("side_sets", MAT_C_CELL, MAT_T_CELL, 2, dims, nullptr, 0);
      assert(cell_array);

      std::vector<matvar_t *> cell_element(num_sets * 7);

      size_t           num_sides = ex_inquire_int(exo_file, EX_INQ_SS_ELEM_LEN);
      size_t           num_nodes = ex_inquire_int(exo_file, EX_INQ_SS_NODE_LEN);
      std::vector<int> elem_list(num_sides);
      std::vector<int> side_list(num_sides);
      std::vector<int> num_nodes_per_side(num_sides);
      std::vector<int> side_nodes(num_nodes);

      // size_t num_df    = ex_inquire_int(exo_file, EX_INQ_SS_DF_LEN);
      // If `num_df == 0` or if it isn't equal to `num_nodes`, then
      // all df will be set to 1.0, but in any case, the size of
      // `ssdfac` should be num_nodes and not num_df.
      std::vector<double> ssdfac(num_nodes);

      size_t side_off = 0;
      size_t node_off = 0;
      size_t df_off   = 0;

      int max_name_length = ex_inquire_int(exo_file, EX_INQ_DB_MAX_USED_NAME_LENGTH);
      max_name_length     = max_name_length < 32 ? 32 : max_name_length;

      for (int i = 0; i < num_sets; i++) {
        std::vector<char> name(max_name_length + 1);
        ex_get_name(exo_file, EX_SIDE_SET, ids[i], name.data());
        dims[0]      = 1;
        dims[1]      = std::strlen(name.data());
        size_t index = 7 * i + 0;
        cell_element[index] =
            Mat_VarCreate(nullptr, MAT_C_CHAR, MAT_T_UINT8, 2, dims, (void *)name.data(), 0);
        Mat_VarSetCell(cell_array, index, cell_element[index]);

        dims[0]             = 1;
        dims[1]             = 1;
        index               = 7 * i + 1;
        cell_element[index] = Mat_VarCreate(nullptr, MAT_C_INT32, MAT_T_INT32, 2, dims, &ids[i],
                                            MAT_F_DONT_COPY_DATA);
        Mat_VarSetCell(cell_array, index, cell_element[index]);

        int n1;
        int n2;
        ex_get_set_param(exo_file, EX_SIDE_SET, ids[i], &n1, &n2);
        num_sideset_sides[i] = n1;
        num_sideset_dfac[i]  = n2;

        bool has_ss_dfac = !(n2 == 0 || n1 == n2);
        if (!has_ss_dfac) {
          num_sideset_dfac[i] = num_sideset_nodes[i];
        }

        ex_get_side_set_node_list_len(exo_file, ids[i], &num_sideset_nodes[i]);

        /* element and side list for side sets (dgriffi) */
        ex_get_set(exo_file, EX_SIDE_SET, ids[i], &elem_list[side_off], &side_list[side_off]);
        dims[0]             = num_sideset_sides[i];
        dims[1]             = 1;
        index               = 7 * i + 2;
        cell_element[index] = Mat_VarCreate(nullptr, MAT_C_INT32, MAT_T_INT32, 2, dims,
                                            &elem_list[side_off], MAT_F_DONT_COPY_DATA);
        Mat_VarSetCell(cell_array, index, cell_element[index]);

        index               = 7 * i + 3;
        cell_element[index] = Mat_VarCreate(nullptr, MAT_C_INT32, MAT_T_INT32, 2, dims,
                                            &side_list[side_off], MAT_F_DONT_COPY_DATA);
        Mat_VarSetCell(cell_array, index, cell_element[index]);

        ex_get_side_set_node_list(exo_file, ids[i], &num_nodes_per_side[side_off],
                                  &side_nodes[node_off]);

        /* number-of-nodes-per-side list */
        index               = 7 * i + 4;
        cell_element[index] = Mat_VarCreate(nullptr, MAT_C_INT32, MAT_T_INT32, 2, dims,
                                            &num_nodes_per_side[side_off], MAT_F_DONT_COPY_DATA);
        Mat_VarSetCell(cell_array, index, cell_element[index]);

        dims[0]             = num_sideset_nodes[i];
        dims[1]             = 1;
        index               = 7 * i + 5;
        cell_element[index] = Mat_VarCreate(nullptr, MAT_C_INT32, MAT_T_INT32, 2, dims,
                                            &side_nodes[node_off], MAT_F_DONT_COPY_DATA);
        Mat_VarSetCell(cell_array, index, cell_element[index]);

        /* distribution-factors list */
        if (has_ss_dfac) {
          ex_get_set_dist_fact(exo_file, EX_SIDE_SET, ids[i], &ssdfac[df_off]);
        }
        else {
          n2 = num_sideset_dfac[i];
          for (int j = 0; j < n2; j++) {
            ssdfac[j] = 1.0;
          }
        }
        dims[0]             = num_sideset_dfac[i];
        dims[1]             = 1;
        index               = 7 * i + 6;
        cell_element[index] = Mat_VarCreate(nullptr, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims,
                                            &ssdfac[df_off], MAT_F_DONT_COPY_DATA);
        Mat_VarSetCell(cell_array, index, cell_element[index]);

        side_off += num_sideset_sides[i];
        node_off += num_sideset_nodes[i];
        df_off += num_sideset_dfac[i];
      }
      Mat_VarWrite(mat_file, cell_array, MAT_COMPRESSION_NONE);
      Mat_VarFree(cell_array);
    }
    else {
      PutInt("ssids", num_sets, 1, ids.data());
      std::vector<int>    elem_list;
      std::vector<int>    side_list;
      std::vector<int>    num_nodes_per_side;
      std::vector<int>    side_nodes;
      std::vector<double> ssdfac;
      for (int i = 0; i < num_sets; i++) {
        int n1;
        int n2;
        ex_get_set_param(exo_file, EX_SIDE_SET, ids[i], &n1, &n2);
        num_sideset_sides[i] = n1;
        num_sideset_dfac[i]  = n2;

        bool has_ss_dfac = (n2 != 0);
        if (n2 == 0 || n1 == n2) {
          ex_get_side_set_node_list_len(exo_file, ids[i], &n2);
        }

        num_nodes_per_side.resize(n1);
        side_nodes.resize(n2);
        ex_get_side_set_node_list(exo_file, ids[i], num_nodes_per_side.data(), side_nodes.data());

        /* number-of-nodes-per-side list */
        std::string str;
        str = fmt::sprintf("ssnum%02d", i + 1);
        PutInt(str, n1, 1, num_nodes_per_side.data());
        /* nodes list */
        str = fmt::sprintf("ssnod%02d", i + 1);
        PutInt(str, n2, 1, side_nodes.data());

        /* distribution-factors list */
        if (has_ss_dfac) {
          ssdfac.resize(n2);
          ex_get_set_dist_fact(exo_file, EX_SIDE_SET, ids[i], ssdfac.data());
          str = fmt::sprintf("ssfac%02d", i + 1);
          PutDbl(str, n2, 1, ssdfac.data());
        }

        /* element and side list for side sets (dgriffi) */
        elem_list.resize(n1);
        side_list.resize(n1);
        ex_get_set(exo_file, EX_SIDE_SET, ids[i], elem_list.data(), side_list.data());
        str = fmt::sprintf("ssside%02d", i + 1);
        PutInt(str, n1, 1, side_list.data());
        str = fmt::sprintf("sselem%02d", i + 1);
        PutInt(str, n1, 1, elem_list.data());
      }
    }
    get_put_user_names(exo_file, EX_SIDE_SET, num_sets, "ssusernames");

    /* Store # sides and # dis. factors per side set (dgriffi) */
    PutInt("nsssides", num_sets, 1, num_sideset_sides.data());
    PutInt("nssdfac", num_sets, 1, num_sideset_dfac.data());
  }
  return num_sideset_sides;
}

void handle_coordinates(int exo_file, size_t num_nodes, int num_axes)
{
  if (debug) {
    logger("Coordinates");
  }
  std::vector<double> x, y, z;
  x.resize(num_nodes);
  if (num_axes >= 2) {
    y.resize(num_nodes);
  }
  if (num_axes == 3) {
    z.resize(num_nodes);
  }
  ex_get_coord(exo_file, x.data(), y.data(), z.data());
  PutDbl("x0", num_nodes, 1, x.data());
  if (num_axes >= 2) {
    PutDbl("y0", num_nodes, 1, y.data());
  }
  if (num_axes == 3) {
    PutDbl("z0", num_nodes, 1, z.data());
  }
}

/**********************************************************************/
/* remove an argument from the list */
void del_arg(int *argc, char *argv[], int j)
{
  for (int jj = j + 1; jj < *argc; jj++) {
    argv[jj - 1] = argv[jj];
  }
  (*argc)--;
  argv[*argc] = nullptr;
}
/**********************************************************************/
int main(int argc, char *argv[])
{
  std::string oname{};
  std::string filename{};

  std::string str;

  const char *ext = EXT;

  int err;
  int num_axes;
  int num_blocks;
  int num_side_sets;
  int num_node_sets;
  int num_time_steps;
  int num_info_lines;
  int num_global_vars;
  int num_nodal_vars;
  int num_element_vars;
  int num_nodeset_vars;
  int num_sideset_vars;

  size_t num_nodes    = 0;
  size_t num_elements = 0;

  int  mat_version     = 73;
  bool use_cell_arrays = false;

  /* process arguments */
  for (int j = 1; j < argc && argv[j][0] == '-'; j++) {
    if (strcmp(argv[j], "-t") == 0) { /* write text file (*.m) */
      del_arg(&argc, argv, j);
      textfile = 1;
      j--;
      continue;
    }
    if (strcmp(argv[j], "-h") == 0 || strcmp(argv[j], "--help") == 0) { /* write help info */
      del_arg(&argc, argv, j);
      usage();
      exit(1);
    }
    if (strcmp(argv[j], "-d") == 0) { /* write help info */
      del_arg(&argc, argv, j);
      j--;
      debug = true;
      continue;
    }
    if (strcmp(argv[j], "-c") == 0) { /* use cell arrays */
      del_arg(&argc, argv, j);
      j--;
      use_cell_arrays = true;
      continue;
    }
    if (strcmp(argv[j], "-v73") == 0) { /* Version 7.3 */
      del_arg(&argc, argv, j);
      mat_version = 73;
      j--;
      continue;
    }
    // This matches the option used in matlab
    if ((strcmp(argv[j], "-v7.3") == 0) || (strcmp(argv[j], "-V7.3") == 0)) { /* Version 7.3 */
      del_arg(&argc, argv, j);
      mat_version = 73;
      j--;
      continue;
    }
    if (strcmp(argv[j], "-v5") == 0) { /* Version 5 (default) */
      del_arg(&argc, argv, j);
      mat_version = 50;
      j--;
      continue;
    }
    if (strcmp(argv[j], "-o") == 0) { /* specify output file name */
      del_arg(&argc, argv, j);
      if (argv[j] != nullptr) {
        oname = argv[j];
        del_arg(&argc, argv, j);
        fmt::print("output file: {}\n", oname);
      }
      else {
        fmt::print(stderr, "ERROR: Invalid output file specification.\n");
        return 2;
      }
      j--;
      continue;
    }
    // Unrecognized option...
    fmt::print(stderr, "ERROR: Unrecognized option: '{}'\n", argv[j]);
    usage();
    exit(1);
  }

  /* QA Info */
  printf("%s: %s, %s\n", qainfo[0], qainfo[2], qainfo[1]);

  /* usage message*/
  if (argc != 2) {
    usage();
    exit(1);
  }

  /* open output file */
  if (textfile != 0) {
    ext = ".m";
  }
  if (oname.empty()) {
    filename   = argv[1];
    size_t pos = filename.find_last_of('.');
    if (pos != std::string::npos) {
      filename = filename.substr(0, pos);
    }
    filename += std::string(ext);
  }
  else {
    filename = oname;
  }

  if (textfile != 0) {
    m_file = fopen(filename.c_str(), "w");
    if (m_file == nullptr) {
      fmt::print(stderr, "ERROR: Unable to open '{}'\n", filename);
      exit(1);
    }
  }
  else {
    if (mat_version == 50) {
      mat_file = Mat_CreateVer(filename.c_str(), nullptr, MAT_FT_MAT5);
    }
    else if (mat_version == 73) {
      mat_file = Mat_CreateVer(filename.c_str(), nullptr, MAT_FT_MAT73);
    }

    if (mat_file == nullptr) {
      fmt::print(stderr, "ERROR: Unable to create matlab file '{}'\n", filename);
      exit(1);
    }
  }

  ex_opts(EX_VERBOSE);

  /* word sizes */
  int cpu_word_size = sizeof(double);
  int io_word_size  = 0;

  /* open exodus file */
  float exo_version;
  int   exo_file = ex_open(argv[1], EX_READ, &cpu_word_size, &io_word_size, &exo_version);
  if (exo_file < 0) {
    fmt::print(stderr, "ERROR: Cannot open '{}'\n", argv[1]);
    exit(1);
  }

  /* print */
  fmt::print("\ttranslating {} to {}...\n", argv[1], filename);

  /* read database parameters */
  char *line = reinterpret_cast<char *>(calloc((MAX_LINE_LENGTH + 1), sizeof(char)));
  ex_get_init(exo_file, line, &num_axes, &num_nodes, &num_elements, &num_blocks, &num_node_sets,
              &num_side_sets);
  num_info_lines = ex_inquire_int(exo_file, EX_INQ_INFO);
  num_time_steps = ex_inquire_int(exo_file, EX_INQ_TIME);
  ex_get_variable_param(exo_file, EX_GLOBAL, &num_global_vars);
  ex_get_variable_param(exo_file, EX_NODAL, &num_nodal_vars);
  ex_get_variable_param(exo_file, EX_ELEM_BLOCK, &num_element_vars);
  ex_get_variable_param(exo_file, EX_NODE_SET, &num_nodeset_vars);
  ex_get_variable_param(exo_file, EX_SIDE_SET, &num_sideset_vars);

  /* export parameters */
  PutInt("naxes", num_axes);
  PutInt("nnodes", num_nodes);
  PutInt("nelems", num_elements);
  PutInt("nblks", num_blocks);
  PutInt("nnsets", num_node_sets);
  PutInt("nssets", num_side_sets);
  PutInt("nsteps", num_time_steps);
  PutInt("ngvars", num_global_vars);
  PutInt("nnvars", num_nodal_vars);
  PutInt("nevars", num_element_vars);
  PutInt("nnsvars", num_nodeset_vars);
  PutInt("nssvars", num_sideset_vars);

  /* allocate -char- scratch space*/
  int nstr2   = num_info_lines;
  nstr2       = std::max(nstr2, num_blocks);
  nstr2       = std::max(nstr2, num_node_sets);
  nstr2       = std::max(nstr2, num_side_sets);
  char **str2 = get_exodus_names(nstr2, 512);

  /* title */
  PutStr("Title", line);

  /* information records */
  if (num_info_lines > 0) {
    ex_get_info(exo_file, str2);
    std::string ostr;
    for (int i = 0; i < num_info_lines; i++) {
      if (std::strlen(str2[i]) > 0) {
        ostr += str2[i];
        ostr += "\n";
      }
    }
    PutStr("info", ostr.c_str());
    ostr = "";
    for (int i = 0; i < num_info_lines; i++) {
      if (std::strlen(str2[i]) > 0 && strncmp(str2[i], "cavi", 4) == 0) {
        ostr += str2[i];
        ostr += "\n";
      }
    }
    PutStr("cvxp", ostr.c_str());
  }

  /* nodal coordinates */
  handle_coordinates(exo_file, num_nodes, num_axes);

  /* side sets */
  if (debug) {
    logger("Side Sets");
  }
  auto num_sideset_sides = handle_side_sets(exo_file, num_side_sets, use_cell_arrays);

  /* node sets (section by dgriffi) */
  auto num_nodeset_nodes = handle_node_sets(exo_file, num_node_sets, use_cell_arrays);

  /* element blocks */
  if (debug) {
    logger("Element Blocks");
  }
  auto num_elem_in_block = handle_element_blocks(exo_file, num_blocks, use_cell_arrays);

  /* time values */
  if (num_time_steps > 0) {
    if (debug) {
      logger("Time Steps");
    }
    std::vector<double> scr(num_time_steps);
    ex_get_all_times(exo_file, scr.data());
    PutDbl("time", num_time_steps, 1, scr.data());
  }

  /* global variables */
  if (num_global_vars > 0) {
    if (debug) {
      logger("Global Variables");
    }

    if (use_cell_arrays) {
      size_t dims[2];
      dims[0]              = 2;
      dims[1]              = num_global_vars;
      matvar_t *cell_array = Mat_VarCreate("gvar", MAT_C_CELL, MAT_T_CELL, 2, dims, nullptr, 0);
      assert(cell_array);
      std::vector<double> scr(num_time_steps * num_global_vars);
      dims[0]       = num_time_steps;
      dims[1]       = 1;
      size_t offset = 0;
      // Get vector of variable names...
      auto gnames = get_names(exo_file, EX_GLOBAL, num_global_vars);

      std::vector<matvar_t *> cell_element(num_global_vars * 2);
      int                     j = 0;
      for (int i = 0; i < num_global_vars; i++) {
        size_t sdims[2];
        sdims[0]        = 1;
        sdims[1]        = gnames[i].length();
        cell_element[j] = Mat_VarCreate(nullptr, MAT_C_CHAR, MAT_T_UINT8, 2, sdims,
                                        (void *)gnames[i].c_str(), MAT_F_DONT_COPY_DATA);
        Mat_VarSetCell(cell_array, j, cell_element[j]);
        j++;

        ex_get_var_time(exo_file, EX_GLOBAL, i + 1, 1, 1, num_time_steps, &scr[offset]);
        cell_element[j] = Mat_VarCreate(nullptr, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &scr[offset],
                                        MAT_F_DONT_COPY_DATA);
        assert(cell_element[j]);
        Mat_VarSetCell(cell_array, j, cell_element[j]);
        offset += num_time_steps;
        j++;
      }
      Mat_VarWrite(mat_file, cell_array, MAT_COMPRESSION_NONE);
      Mat_VarFree(cell_array);
    }
    else {
      get_put_names(exo_file, EX_GLOBAL, num_global_vars, "gnames");
      std::vector<double> scr(num_time_steps);
      for (int i = 0; i < num_global_vars; i++) {
        str = fmt::sprintf("gvar%02d", i + 1);
        ex_get_var_time(exo_file, EX_GLOBAL, i + 1, 1, 1, num_time_steps, scr.data());
        PutDbl(str, num_time_steps, 1, scr.data());
      }
    }
  }

  /* nodal variables */
  if (num_nodal_vars > 0) {
    if (debug) {
      logger("Nodal Variables");
    }
    if (debug) {
      logger("\tNames");
    }
    if (use_cell_arrays) {
      size_t dims[2];
      dims[0]              = 2;
      dims[1]              = num_nodal_vars;
      matvar_t *cell_array = Mat_VarCreate("nvar", MAT_C_CELL, MAT_T_CELL, 2, dims, nullptr, 0);
      assert(cell_array);
      std::vector<double> scr(num_nodal_vars * num_time_steps * num_nodes);
      dims[0]       = num_nodes;
      dims[1]       = num_time_steps;
      size_t offset = 0;
      // Get vector of variable names...
      auto nnames = get_names(exo_file, EX_NODAL, num_nodal_vars);

      std::vector<matvar_t *> cell_element(num_nodal_vars * 2);
      int                     j = 0;
      for (int i = 0; i < num_nodal_vars; i++) {
        size_t sdims[2];
        sdims[0]        = 1;
        sdims[1]        = nnames[i].length();
        cell_element[j] = Mat_VarCreate(nullptr, MAT_C_CHAR, MAT_T_UINT8, 2, sdims,
                                        (void *)nnames[i].c_str(), MAT_F_DONT_COPY_DATA);
        Mat_VarSetCell(cell_array, j, cell_element[j]);
        j++;

        cell_element[j] = Mat_VarCreate(nullptr, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &scr[offset],
                                        MAT_F_DONT_COPY_DATA);
        assert(cell_element[j]);
        Mat_VarSetCell(cell_array, j, cell_element[j]);
        for (int k = 0; k < num_time_steps; k++) {
          ex_get_var(exo_file, k + 1, EX_NODAL, i + 1, 1, num_nodes, &scr[num_nodes * k + offset]);
        }
        offset += num_time_steps * num_nodes;
        j++;
      }
      Mat_VarWrite(mat_file, cell_array, MAT_COMPRESSION_NONE);
      Mat_VarFree(cell_array);
    }
    else {
      get_put_names(exo_file, EX_NODAL, num_nodal_vars, "nnames");

      std::vector<double> scr(num_nodes * num_time_steps);
      for (int i = 0; i < num_nodal_vars; i++) {
        str = fmt::sprintf("nvar%02d", i + 1);
        if (debug) {
          logger("\tReading");
        }
        for (int j = 0; j < num_time_steps; j++) {
          ex_get_var(exo_file, j + 1, EX_NODAL, i + 1, 1, num_nodes, &scr[num_nodes * j]);
        }
        if (debug) {
          logger("\tWriting");
        }
        PutDbl(str, num_nodes, num_time_steps, scr.data());
      }
    }
  }

  /* element variables */
  if (num_element_vars > 0) {
    if (debug) {
      logger("Element Variables");
    }
    get_put_vars(exo_file, EX_ELEM_BLOCK, num_blocks, num_element_vars, num_time_steps,
                 num_elem_in_block, "e", use_cell_arrays);
  }

  /* nodeset variables */
  if (num_nodeset_vars > 0) {
    if (debug) {
      logger("Nodeset Variables");
    }
    get_put_vars(exo_file, EX_NODE_SET, num_node_sets, num_nodeset_vars, num_time_steps,
                 num_nodeset_nodes, "ns", use_cell_arrays);
  }

  /* sideset variables */
  if (num_sideset_vars > 0) {
    if (debug) {
      logger("Sideset Variables");
    }
    get_put_vars(exo_file, EX_SIDE_SET, num_side_sets, num_sideset_vars, num_time_steps,
                 num_sideset_sides, "ss", use_cell_arrays);
  }

  /* node and element number maps */
  if (debug) {
    logger("Node and Element Number Maps");
  }
  ex_opts(0); /* turn off error reporting. It is not an error to have no map*/
  std::vector<int> ids(num_nodes);
  err = ex_get_id_map(exo_file, EX_NODE_MAP, ids.data());
  if (err == 0) {
    PutInt("node_num_map", num_nodes, 1, ids.data());
  }

  ids.resize(num_elements);
  err = ex_get_id_map(exo_file, EX_ELEM_MAP, ids.data());
  if (err == 0) {
    PutInt("elem_num_map", num_elements, 1, ids.data());
  }

  if (debug) {
    logger("Closing file");
  }
  ex_close(exo_file);

  if (textfile != 0) {
    fclose(m_file);
  }
  else {
    Mat_Close(mat_file);
  }
  free(line);

  delete_exodus_names(str2, nstr2);

  fmt::print("done...\n");

  /* exit status */
  add_to_log("exo2mat", 0);
  return (0);
}

/*
 * Copyright(C) 2011 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following
 *   disclaimer in the documentation and/or other materials provided
 *   with the distribution.
 *
 * * Neither the name of Sandia Corporation nor the names of its
 *   contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
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
#include <cstring> // for strcat, strlen, strcpy, etc
#include <iostream>
#include <numeric>
#include <vector>

#include "add_to_log.h" // for add_to_log
#include "exodusII.h"   // for ex_get_variable_param, etc
#include "matio.h"      // for Mat_VarCreate, Mat_VarFree, etc
#include <assert.h>     // for assert
#include <stddef.h>     // for size_t
#include <stdio.h>      // for fprintf, printf, sprintf, etc
#include <stdlib.h>     // for free, calloc, exit, malloc
#if MATIO_VERSION < 151
#error "MatIO Version 1.5.1 or greater is required"
#endif

#if __cplusplus > 199711L
#define TOPTR(x) x.data()
#else
#define TOPTR(x) (x.empty() ? nullptr : &x[0])
#endif

#define EXT ".mat"
int textfile = 0;

FILE * m_file   = nullptr; /* file for m file output */
mat_t *mat_file = nullptr; /* file for binary .mat output */
int    debug    = 0;

static const char *qainfo[] = {
    "exo2mat", "2015/10/28", "3.02",
};

std::string time_stamp(const std::string &format)
{
  if (format == "") {
    return std::string("");
  }
  else {
    const int   length = 256;
    static char time_string[length];

    time_t     calendar_time = time(nullptr);
    struct tm *local_time    = localtime(&calendar_time);

    int error = strftime(time_string, length, format.c_str(), local_time);
    if (error != 0) {
      time_string[length - 1] = '\0';
      return std::string(time_string);
    }
    else {
      return std::string("[ERROR]");
    }
  }
}

void logger(const char *message)
{
  const std::string tsFormat = "[%H:%M:%S] ";
  std::clog << time_stamp(tsFormat) << ": " << message << "\n";
}

void usage()
{
  std::cout << "exo2mat [options] exodus_file_name.\n"
            << "   the exodus_file_name is required (exodus only).\n"
            << "   Options:\n"
            << "   -t    write a text (.m) file rather than a binary .mat\n"
            << "   -o    output file name (rather than auto generate)\n"
            << "   -v5   output version 5 mat file\n"
            << "   -v73  output version 7.3 mat file (hdf5-based) [default]\n"
            << "   -v7.3 output version 7.3 mat file (hdf5-based)\n"
            << " ** note **\n"
            << "Binary files are written by default on all platforms.\n";
}

/* put a string into an m file. If the string has
   line feeds, we put it as ints, and use 'char()' to convert it */
void mPutStr(const char *name, const char *str)
{
  assert(m_file != nullptr);
  if (strchr(str, '\n') == nullptr)
    fprintf(m_file, "%s='%s';\n", name, str);
  else {
    fprintf(m_file, "%s=[", name);
    size_t i;
    size_t j;
    for (j = i = 0; i < strlen(str); i++, j++) {
      if (j >= 20) {
        j = 0;
        fprintf(m_file, "...\n");
      }
      fprintf(m_file, "%d ", str[i]);
    }
    fprintf(m_file, "];\n");
    fprintf(m_file, "%s=char(%s);\n", name, name);
  }
}

/* put double array in m file */
void mPutDbl(const char *name, int n1, int n2, double *pd)
{
  assert(m_file != nullptr);
  if (n1 == 1 && n2 == 1) {
    fprintf(m_file, "%s=%15.8e;\n", name, *pd);
    return;
  }
  fprintf(m_file, "%s=zeros(%d,%d);\n", name, n1, n2);
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++)
      fprintf(m_file, "%s(%d,%d)=%15.8e;\n", name, i + 1, j + 1, pd[i * n2 + j]);
}

/* put integer array in m file */
void mPutInt(const char *name, int pd)
{
  assert(m_file != nullptr);
  fprintf(m_file, "%s=%d;\n", name, pd);
  return;
}

/* put integer array in m file */
void mPutInt(const char *name, int n1, int n2, int *pd)
{
  assert(m_file != nullptr);
  if (n1 == 1 && n2 == 1) {
    fprintf(m_file, "%s=%d;\n", name, *pd);
    return;
  }
  fprintf(m_file, "%s=zeros(%d,%d);\n", name, n1, n2);
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++)
      fprintf(m_file, "%s(%d,%d)=%d;\n", name, i + 1, j + 1, pd[i * n2 + j]);
}

/* put string in mat file*/
int matPutStr(const char *name, char *str)
{
  int       error  = 0;
  matvar_t *matvar = nullptr;
  size_t    dims[2];

  dims[0] = 1;
  dims[1] = strlen(str);

  matvar = Mat_VarCreate(name, MAT_C_CHAR, MAT_T_UINT8, 2, dims, str, MAT_F_DONT_COPY_DATA);
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
int matPutDbl(const char *name, int n1, int n2, double *pd)
{
  int       error  = 0;
  matvar_t *matvar = nullptr;

  size_t dims[2];
  dims[0] = n1;
  dims[1] = n2;

  matvar = Mat_VarCreate(name, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, pd, MAT_F_DONT_COPY_DATA);
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
int matPutInt(const char *name, int n1, int n2, int *pd)
{
  int       error  = 0;
  matvar_t *matvar = nullptr;

  size_t dims[2];
  dims[0] = n1;
  dims[1] = n2;

  matvar = Mat_VarCreate(name, MAT_C_INT32, MAT_T_INT32, 2, dims, pd, MAT_F_DONT_COPY_DATA);
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
void PutStr(const char *name, const char *str)
{
  if (textfile)
    mPutStr(name, str);
  else
    matPutStr(name, (char *)str);
}

int PutInt(const char *name, int pd)
{
  int error = 0;
  if (textfile)
    mPutInt(name, pd);
  else
    error = matPutInt(name, 1, 1, &pd);
  return error;
}

int PutInt(const char *name, int n1, int n2, int *pd)
{
  int error = 0;
  if (textfile)
    mPutInt(name, n1, n2, pd);
  else
    error = matPutInt(name, n1, n2, pd);
  return error;
}

int PutDbl(const char *name, int n1, int n2, double *pd)
{
  int error = 0;
  if (textfile)
    mPutDbl(name, n1, n2, pd);
  else
    error = matPutDbl(name, n1, n2, pd);
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

void get_put_names(int exo_file, ex_entity_type type, int num_vars, const char *mname)
{
  int max_name_length = ex_inquire_int(exo_file, EX_INQ_DB_MAX_USED_NAME_LENGTH);
  max_name_length     = max_name_length < 32 ? 32 : max_name_length;
  char **names        = get_exodus_names(num_vars, max_name_length + 1);

  if (debug)
    logger("\tReading variable names");
  ex_get_variable_names(exo_file, type, num_vars, names);

  std::string mat;
  for (int i = 0; i < num_vars; i++) {
    mat += names[i];
    mat += "\n";
  }
  if (debug)
    logger("\tWriting variable names");
  PutStr(mname, mat.c_str());

  delete_exodus_names(names, num_vars);
}

void get_put_vars(int exo_file, ex_entity_type type, int num_blocks, int num_vars,
                  int num_time_steps, const std::vector<int> &num_per_block, const char *mname)

{
  /* truth table */
  if (debug)
    logger("\tTruth Table");
  std::vector<int> truth_table(num_vars * num_blocks);
  ex_get_truth_table(exo_file, type, num_blocks, num_vars, TOPTR(truth_table));

  size_t              num_entity = std::accumulate(num_per_block.begin(), num_per_block.end(), 0);
  std::vector<double> scr(num_entity * num_time_steps);

  std::vector<int> ids(num_blocks);
  ex_get_ids(exo_file, type, TOPTR(ids));

  char str[32];
  for (int i = 0; i < num_vars; i++) {
    if (debug)
      logger("\tReading");
    std::fill(scr.begin(), scr.end(), 0.0);
    size_t n = 0;
    sprintf(str, mname, i + 1);
    for (int j = 0; j < num_time_steps; j++) {
      for (int k = 0; k < num_blocks; k++) {
        if (truth_table[num_vars * k + i] == 1) {
          ex_get_var(exo_file, j + 1, type, i + 1, ids[k], num_per_block[k], &scr[n]);
        }
        n = n + num_per_block[k];
      }
    }
    if (debug)
      logger("\tWriting");
    PutDbl(str, num_entity, num_time_steps, TOPTR(scr));
  }
}

/**********************************************************************/
/* remove an argument from the list */
void del_arg(int *argc, char *argv[], int j)
{
  for (int jj    = j + 1; jj < *argc; jj++)
    argv[jj - 1] = argv[jj];
  (*argc)--;
  argv[*argc] = nullptr;
}
/**********************************************************************/
int main(int argc, char *argv[])
{
  char *oname = nullptr, *dot = nullptr, *filename = nullptr;
  char  str[32];

  const char *ext = EXT;

  int n, n1, n2, err, num_axes, num_blocks, num_side_sets, num_node_sets, num_time_steps,
      num_info_lines, num_global_vars, num_nodal_vars, num_element_vars, num_nodeset_vars,
      num_sideset_vars;

  size_t num_nodes    = 0;
  size_t num_elements = 0;

  int mat_version = 73;

  /* process arguments */
  for (int j = 1; j < argc; j++) {
    if (strcmp(argv[j], "-t") == 0) { /* write text file (*.m) */
      del_arg(&argc, argv, j);
      textfile = 1;
      j--;
      continue;
    }
    if (strcmp(argv[j], "-h") == 0) { /* write help info */
      del_arg(&argc, argv, j);
      usage();
      exit(1);
    }
    if (strcmp(argv[j], "-d") == 0) { /* write help info */
      del_arg(&argc, argv, j);
      j--;
      debug = 1;
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
      if (argv[j]) {
        oname = (char *)calloc(strlen(argv[j]) + 10, sizeof(char));
        strcpy(oname, argv[j]);
        del_arg(&argc, argv, j);
        std::cout << "output file: " << oname << "\n";
      }
      else {
        std::cerr << "ERROR: Invalid output file specification.\n";
        return 2;
      }
      j--;

      continue;
    }
  }

  /* QA Info */
  printf("%s: %s, %s\n", qainfo[0], qainfo[2], qainfo[1]);

  /* usage message*/
  if (argc != 2) {
    usage();
    exit(1);
  }

  /* open output file */
  if (textfile)
    ext = ".m";

  if (!oname) {
    filename = (char *)malloc(strlen(argv[1]) + 10);
    strcpy(filename, argv[1]);
    dot = strrchr(filename, '.');
    if (dot)
      *dot = '\0';
    strcat(filename, ext);
  }
  else {
    filename = oname;
  }

  if (textfile) {
    m_file = fopen(filename, "w");
    if (!m_file) {
      std::cerr << "ERROR: Unable to open " << filename << "\n";
      exit(1);
    }
  }
  else {
    if (mat_version == 50) {
      mat_file = Mat_CreateVer(filename, nullptr, MAT_FT_MAT5);
    }
    else if (mat_version == 73) {
      mat_file = Mat_CreateVer(filename, nullptr, MAT_FT_MAT73);
    }

    if (mat_file == nullptr) {
      std::cerr << "ERROR: Unable to create matlab file " << filename << "\n";
      exit(1);
    }
  }

  /* word sizes */
  int cpu_word_size = sizeof(double);
  int io_word_size  = 0;

  /* open exodus file */
  float exo_version;
  int   exo_file = ex_open(argv[1], EX_READ, &cpu_word_size, &io_word_size, &exo_version);
  if (exo_file < 0) {
    std::cerr << "ERROR: Cannot open " << argv[1] << "\n";
    exit(1);
  }

  /* print */
  std::cout << "\ttranslating " << argv[1] << " to " << filename << "...\n";

  /* read database paramters */
  char *line = (char *)calloc((MAX_LINE_LENGTH + 1), sizeof(char));
  ex_get_init(exo_file, line, &num_axes, &num_nodes, &num_elements, &num_blocks, &num_node_sets,
              &num_side_sets);
  num_info_lines = ex_inquire_int(exo_file, EX_INQ_INFO);
  num_time_steps = ex_inquire_int(exo_file, EX_INQ_TIME);
  ex_get_variable_param(exo_file, EX_GLOBAL, &num_global_vars);
  ex_get_variable_param(exo_file, EX_NODAL, &num_nodal_vars);
  ex_get_variable_param(exo_file, EX_ELEM_BLOCK, &num_element_vars);
  ex_get_variable_param(exo_file, EX_NODE_SET, &num_nodeset_vars);
  ex_get_variable_param(exo_file, EX_SIDE_SET, &num_sideset_vars);

  /* export paramters */
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
      if (strlen(str2[i]) > 0) {
        ostr += str2[i];
        ostr += "\n";
      }
    }
    PutStr("info", ostr.c_str());
    ostr = "";
    for (int i = 0; i < num_info_lines; i++) {
      if (strlen(str2[i]) > 0 && strncmp(str2[i], "cavi", 4) == 0) {
        ostr += str2[i];
        ostr += "\n";
      }
    }
    PutStr("cvxp", ostr.c_str());
  }

  /* nodal coordinates */
  {
    if (debug) {
      logger("Coordinates");
    }
    std::vector<double> x, y, z;
    x.resize(num_nodes);
    if (num_axes >= 2)
      y.resize(num_nodes);
    if (num_axes == 3)
      z.resize(num_nodes);
    ex_get_coord(exo_file, TOPTR(x), TOPTR(y), TOPTR(z));
    PutDbl("x0", num_nodes, 1, TOPTR(x));
    if (num_axes >= 2) {
      PutDbl("y0", num_nodes, 1, TOPTR(y));
    }
    if (num_axes == 3) {
      PutDbl("z0", num_nodes, 1, TOPTR(z));
    }
  }

  /* side sets */
  std::vector<int> num_sideset_sides(num_side_sets);
  std::vector<int> ids;
  if (num_side_sets > 0) {
    if (debug) {
      logger("Side Sets");
    }
    ids.resize(num_side_sets);
    ex_get_ids(exo_file, EX_SIDE_SET, TOPTR(ids));
    PutInt("ssids", num_side_sets, 1, TOPTR(ids));
    std::vector<int>    nssdfac(num_side_sets);
    std::vector<int>    iscr;
    std::vector<int>    jscr;
    std::vector<double> scr;
    std::vector<int>    elem_list;
    std::vector<int>    side_list;
    std::vector<int>    junk;
    for (int i = 0; i < num_side_sets; i++) {
      ex_get_set_param(exo_file, EX_SIDE_SET, ids[i], &n1, &n2);
      num_sideset_sides[i] = n1;
      nssdfac[i]           = n2;
      /*
       * the following provision is from Version 1.6 when there are no
       * distribution factors in exodus file
       */
      bool has_ss_dfac = (n2 != 0);
      if (n2 == 0 || n1 == n2) {

        std::cerr << "WARNING: Exodus II file does not contain distribution factors.\n";

        /* n1=number of faces, n2=number of df */
        /* using distribution factors to determine number of nodes in the sideset
           causes a lot grief since some codes do not output distribution factors
           if they are all equal to 1. mkbhard: I am using the function call below
           to figure out the total number of nodes in this sideset. Some redundancy
           exists, but it works for now */

        junk.resize(n1);
        ex_get_side_set_node_count(exo_file, ids[i], TOPTR(junk));
        n2 = 0; /* n2 will be equal to the total number of nodes in the sideset */
        for (int j = 0; j < n1; j++)
          n2 += junk[j];
      }

      iscr.resize(n1);
      jscr.resize(n2);
      ex_get_side_set_node_list(exo_file, ids[i], TOPTR(iscr), TOPTR(jscr));
      /* number-of-nodes-per-side list */
      sprintf(str, "ssnum%02d", i + 1);
      PutInt(str, n1, 1, TOPTR(iscr));
      /* nodes list */
      sprintf(str, "ssnod%02d", i + 1);
      PutInt(str, n2, 1, TOPTR(jscr));

      /* distribution-factors list */
      scr.resize(n2);
      if (has_ss_dfac) {
        ex_get_side_set_dist_fact(exo_file, ids[i], TOPTR(scr));
      }
      else {
        for (int j = 0; j < n2; j++) {
          scr[j] = 1.0;
        }
      }
      sprintf(str, "ssfac%02d", i + 1);
      PutDbl(str, n2, 1, TOPTR(scr));

      /* element and side list for side sets (dgriffi) */
      elem_list.resize(n1);
      side_list.resize(n1);
      ex_get_set(exo_file, EX_SIDE_SET, ids[i], TOPTR(elem_list), TOPTR(side_list));
      sprintf(str, "ssside%02d", i + 1);
      PutInt(str, n1, 1, TOPTR(side_list));
      sprintf(str, "sselem%02d", i + 1);
      PutInt(str, n1, 1, TOPTR(elem_list));
    }
    /* Store # sides and # dis. factors per side set (dgriffi) */
    PutInt("nsssides", num_side_sets, 1, TOPTR(num_sideset_sides));
    PutInt("nssdfac", num_side_sets, 1, TOPTR(nssdfac));
  }

  /* node sets (section by dgriffi) */
  std::vector<int> num_nodeset_nodes(num_node_sets);
  if (num_node_sets > 0) {
    if (debug) {
      logger("Node Sets");
    }
    std::vector<int>    iscr;
    std::vector<double> scr;
    ids.resize(num_node_sets);
    ex_get_ids(exo_file, EX_NODE_SET, TOPTR(ids));
    PutInt("nsids", num_node_sets, 1, TOPTR(ids));

    std::vector<int> num_nodeset_df(num_node_sets);
    for (int i = 0; i < num_node_sets; i++) {
      ex_get_set_param(exo_file, EX_NODE_SET, ids[i], &n1, &n2);
      iscr.resize(n1);
      ex_get_node_set(exo_file, ids[i], TOPTR(iscr));
      /* nodes list */
      sprintf(str, "nsnod%02d", i + 1);
      PutInt(str, n1, 1, TOPTR(iscr));
      {
        /* distribution-factors list */
        scr.resize(n2);
        ex_get_node_set_dist_fact(exo_file, ids[i], TOPTR(scr));
        sprintf(str, "nsfac%02d", i + 1);
        PutDbl(str, n2, 1, TOPTR(scr));
      }
      num_nodeset_nodes[i] = n1;
      num_nodeset_df[i]    = n2;
    }

    /* Store # nodes and # dis. factors per node set */
    PutInt("nnsnodes", num_node_sets, 1, TOPTR(num_nodeset_nodes));
    PutInt("nnsdfac", num_node_sets, 1, TOPTR(num_nodeset_df));
  }

  /* element blocks */
  if (debug) {
    logger("Element Blocks");
  }
  std::vector<int> num_elem_in_block(num_blocks);
  {
    ids.resize(num_blocks);
    std::vector<int> iscr;
    ex_get_ids(exo_file, EX_ELEM_BLOCK, TOPTR(ids));
    PutInt("blkids", num_blocks, 1, TOPTR(ids));
    for (int i = 0; i < num_blocks; i++) {
      ex_get_elem_block(exo_file, ids[i], str2[i], &n, &n1, &n2);
      num_elem_in_block[i] = n;
      iscr.resize(n * n1);
      ex_get_conn(exo_file, EX_ELEM_BLOCK, ids[i], TOPTR(iscr), nullptr, nullptr);
      sprintf(str, "blk%02d", i + 1);
      PutInt(str, n1, n, TOPTR(iscr));
    }
    str[0] = '\0';
    for (int i = 0; i < num_blocks; i++) {
      strcat(str, str2[i]);
      strcat(str, "\n");
    }
    PutStr("blknames", str);
  }

  /* time values */
  if (num_time_steps > 0) {
    if (debug) {
      logger("Time Steps");
    }
    std::vector<double> scr(num_time_steps);
    ex_get_all_times(exo_file, TOPTR(scr));
    PutDbl("time", num_time_steps, 1, TOPTR(scr));
  }

  /* global variables */
  if (num_global_vars > 0) {
    if (debug) {
      logger("Global Variables");
    }
    get_put_names(exo_file, EX_GLOBAL, num_global_vars, "gnames");

    std::vector<double> scr(num_time_steps);
    for (int i = 0; i < num_global_vars; i++) {
      sprintf(str, "gvar%02d", i + 1);
      ex_get_glob_var_time(exo_file, i + 1, 1, num_time_steps, TOPTR(scr));
      PutDbl(str, num_time_steps, 1, TOPTR(scr));
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
    get_put_names(exo_file, EX_NODAL, num_nodal_vars, "nnames");

    std::vector<double> scr(num_nodes * num_time_steps);
    for (int i = 0; i < num_nodal_vars; i++) {
      sprintf(str, "nvar%02d", i + 1);
      if (debug) {
        logger("\tReading");
      }
      for (int j = 0; j < num_time_steps; j++) {
        ex_get_nodal_var(exo_file, j + 1, i + 1, num_nodes, &scr[num_nodes * j]);
      }
      if (debug) {
        logger("\tWriting");
      }
      PutDbl(str, num_nodes, num_time_steps, TOPTR(scr));
    }
  }

  /* element variables */
  if (num_element_vars > 0) {
    if (debug) {
      logger("Element Variables");
    }
    get_put_names(exo_file, EX_ELEM_BLOCK, num_element_vars, "enames");

    get_put_vars(exo_file, EX_ELEM_BLOCK, num_blocks, num_element_vars, num_time_steps,
                 num_elem_in_block, "evar%02d");
  }

  /* nodeset variables */
  if (num_nodeset_vars > 0) {
    if (debug) {
      logger("Nodeset Variables");
    }
    get_put_names(exo_file, EX_NODE_SET, num_nodeset_vars, "nsnames");

    get_put_vars(exo_file, EX_NODE_SET, num_node_sets, num_nodeset_vars, num_time_steps,
                 num_nodeset_nodes, "nsvar%02d");
  }

  /* sideset variables */
  if (num_sideset_vars > 0) {
    if (debug) {
      logger("Sideset Variables");
    }
    get_put_names(exo_file, EX_SIDE_SET, num_sideset_vars, "ssnames");

    get_put_vars(exo_file, EX_SIDE_SET, num_side_sets, num_sideset_vars, num_time_steps,
                 num_sideset_sides, "ssvar%02d");
  }

  /* node and element number maps */
  if (debug) {
    logger("Node and Element Number Maps");
  }
  ex_opts(0); /* turn off error reporting. It is not an error to have no map*/
  ids.resize(num_nodes);
  err = ex_get_node_num_map(exo_file, TOPTR(ids));
  if (err == 0) {
    PutInt("node_num_map", num_nodes, 1, TOPTR(ids));
  }

  ids.resize(num_elements);
  err = ex_get_elem_num_map(exo_file, TOPTR(ids));
  if (err == 0) {
    PutInt("elem_num_map", num_elements, 1, TOPTR(ids));
  }

  if (debug) {
    logger("Closing file");
  }
  ex_close(exo_file);

  if (textfile)
    fclose(m_file);
  else
    Mat_Close(mat_file);

  std::cout << "done...\n";

  free(filename);
  free(line);

  delete_exodus_names(str2, nstr2);

  /* exit status */
  add_to_log("exo2mat", 0);
  return (0);
}

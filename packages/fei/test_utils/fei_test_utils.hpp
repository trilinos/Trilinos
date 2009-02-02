/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_test_utils_hpp_
#define _fei_test_utils_hpp_

#include <fei_macros.hpp>

#include <fei_mpi.h>

#include <fei_fwd.hpp>

#include <fei_FillableMat.hpp>

#include <vector>
#include <string>

/** The fei_test_utils namespace contains general test-utility functions.
*/
namespace fei_test_utils {

  /** Look through argv for a '-i' argument which names an input file,
      and a '-d' argument which names the directory in which the
      input file is to be found. If the '-i' argument is present but
      the '-d' argument is not present, then the path to the input
      file will be assumed to be ".".
  */
  std::string construct_filename(int argc, char** argv);

  /** If the macro FEI_SER is not defined, call MPI_Init and then call
      MPI_Comm_rank to set localProc and MPI_Comm_size to set numProcs.

      If the macro FEI_SER is defined, then simply set localProc = 0 and
      set numProcs = 1.
  */
  int initialize_mpi(int argc, char** argv,
                     int& localProc, int& numProcs);

  /** If flag is in the specified command-line arguments, and
      it is accompanied by either "yes" or "true" (case insensitive), then
      return true. If it is accompanied by "no" or "false", return false.
      If flag is not in the specified command-line arguments, return false
      unless the optional parameter 'default_result' is set to true.
      
      (e.g., return true if some argv[i] == flag and argv[i+1] == "yes".)

  */
  bool bool_arg(const char* flag, int argc, char** argv,
                bool default_result=false);

  /** Return argument string corresponding to flag from the command-line
      arguments.
      Example: if command-line arguments include some argv[i]="-d"
      and argv[i+1] = "path", and this function is called with flag=="-d",
      then the returned string value will be "path".
  */
  std::string get_arg_value(const char* flag, int argc, char** argv);

  /** Broadcast string from 'root' processor so that on exit, the string
      has the same contents on all processors. Does nothing if num-procs==1
      or if the compile-time macro FEI_SER is defined.
  */
  void broadcast_string(MPI_Comm comm, int root, std::string& strg);

  /** Check command-line arguments for an input-file specified by '-i', and
      optionally a path specified by '-d', and then read the contents of
      the input-file into the user-provided parameter-strings.

      Filename is constructed on proc 0, file is read on proc 0, and contents
      are returned on all processors.
  */
  int get_filename_and_read_input(int argc, char** argv,
				  MPI_Comm comm, int localProc,
				  std::vector<std::string>& stdstrings);

  /** read contents of file line by line into a vector of strings. This is a
      purely serial operation.
  */
  void read_file_lines_into_strings(const char* filename,
				    std::vector<std::string>& file_contents);

  /** This function reads the contents of filename on processor 0, and
      broadcasts those contents (strings) to all processors.
      All processors return the file contents in the file_contents argument.

      If the file is not found, or can't be opened, an std::runtime_error will be
      thrown.
  */
  void read_input_file(const char* filename, MPI_Comm comm,
		       std::vector<std::string>& file_contents);

  /** Given a file-name and test-name, return the named benchmark value.

      If anything goes wrong, such as the file can't be read, or the specified
      testname doesn't appear in the file, throw an std::runtime_error.
  */
  double get_file_benchmark(const char* filename,
			    const char* testname);

  /** Given two values determine whether the values are within 'margin' percentage
      points of each other.
      i.e., return true if
      100*abs(value1 - value2)/max(abs(value1),abs(value2)) <= margin

      Note: if max(abs(value1),abs(value2)) < 1.e-14, return true
  */
  bool within_percentage_margin(double value1,
				double value2,
				unsigned margin);

  bool check_and_cout_test_result(std::string testname,
				  double value,
				  double file_value,
				  unsigned margin);

  std::string check_test_result(double value,
				double goldvalue,
				unsigned margin);

  int compare_with_file_benchmark(const char* name,
				  double benchmark,
				  const char* filename);

  int whichArg(int argc, const char*const* argv, const char* findarg);

  int dirname(const char* name, char*& dir);

  void print_args(int argc, char** argv);

  int compareMatrices(fei::FillableMat& mat1, fei::FillableMat& mat2, double tol=1.e-15);

  int readMatrix(const char* baseName, int np, fei::FillableMat& matrix);

  int readMatrix(const char* fileName, fei::FillableMat& matrix);

  int writeMatrix(const char* fileName, fei::FillableMat& matrix);

  int copy_feiMatrix_to_FillableMat(fei::Matrix& feimat, fei::FillableMat& ssmat);

}//namespace fei_test_utils

#endif


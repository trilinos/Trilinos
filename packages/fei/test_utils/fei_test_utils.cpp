/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <fei_sstream.hpp>
#include <fei_fstream.hpp>

#include <fei_test_utils.hpp>

#include <snl_fei_Utils.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_MatrixGraph.hpp>
#include <fei_Matrix.hpp>

#include <cmath>

#undef fei_file
#define fei_file "fei_test_utils.cpp"
#include <fei_ErrMacros.hpp>

namespace fei_test_utils {

std::string construct_filename(int argc, char** argv)
{
  std::string path(".");
  std::string filename;
  std::string result;

  int dashDarg = whichArg(argc, argv, "-d");

  int dashIarg = whichArg(argc, argv, "-i");

  if (dashIarg < 0) {
    fei::console_out() << "fei_test_utils::construct_filename: argument '-i' not found."
             << FEI_ENDL;
    return(result);
  }

  if (dashDarg > -1 && argc > dashDarg+1) {
    if (argv[dashDarg+1] != 0) {
      path = argv[dashDarg+1];
    }
  }

  if (argc > dashIarg+1) {
    if (argv[dashIarg+1] != 0) {
      filename = argv[dashIarg+1];
    }
  }

  FEI_OSTRINGSTREAM osstr;
  osstr << path;

  std::size_t string_length = path.size();
  if (path[string_length-1] != '/') osstr << "/";
  
  osstr << filename;

  return( osstr.str() );
}

int initialize_mpi(int argc, char** argv, int& localProc, int& numProcs)
{
#ifndef FEI_SER
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) ERReturn(-1);
  if (MPI_Comm_rank(MPI_COMM_WORLD, &localProc) != MPI_SUCCESS) ERReturn(-1);
  if (MPI_Comm_size(MPI_COMM_WORLD, &numProcs) != MPI_SUCCESS) ERReturn(-1);
#else
  localProc = 0;
  numProcs = 1;
#endif
  return(0);
}

bool bool_arg(const char* flag, int argc, char** argv,
                       bool default_result)
{
  int arg_index = whichArg(argc, argv, flag);
  bool result = default_result;
  if (arg_index > -1) {
    if (argc > arg_index+1) {
      if (argv[arg_index+1] != 0) {
        std::string strarg(argv[arg_index+1]);
        if (strarg == "yes" || strarg == "Yes" || strarg == "YES" ||
            strarg == "true" || strarg == "True" || strarg == "TRUE") {
          result = true;
        }

        if (strarg == "no" || strarg == "No" || strarg == "NO" ||
            strarg == "false" || strarg == "False" || strarg == "FALSE") {
          result = false;
        }
      }
    }
  }

  return(result);
}

std::string get_arg_value(const char* flag, int argc, char** argv)
{
  std::string result;
  int arg_index = whichArg(argc, argv, flag);
  if (arg_index > -1) {
    if (argc > arg_index+1) {
      if (argv[arg_index+1] != 0) {
        result = argv[arg_index+1];
      }
    }
  }

  return(result);
}

void broadcast_string(MPI_Comm comm, int root, std::string& strg)
{
#ifdef FEI_SER
  return;
#else
  int numprocs = 1, localproc = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &localproc);
  if (numprocs==1) {
    return;
  }

  if (localproc == root) {
    int length = strg.size();
    MPI_Bcast(&length, 1, MPI_INT, root, comm);
    char* cstr = const_cast<char*>(strg.c_str());
    MPI_Bcast(cstr, strg.size(), MPI_CHAR, root, comm);
  }
  else {
    int length;
    MPI_Bcast(&length, 1, MPI_INT, root, comm);
    char* charstr = new char[length+1];
    MPI_Bcast(charstr, length, MPI_CHAR, root, comm);
    charstr[length] = '\0';
    strg = charstr;
    delete [] charstr;
  }

  return;
#endif
}

void read_file_lines_into_strings(const char* filename,
				   std::vector<std::string>& file_contents)
{
  FEI_IFSTREAM infile(filename);
  if (!infile) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei_test_utils::read_file_lines_into_strings ERROR, couldn't open file '"
         << filename << "'";
    throw std::runtime_error(osstr.str());
  }

  file_contents.clear();

  std::string line;

  getline(infile, line);
  while(!infile.eof()) {
    file_contents.push_back(line);

    line.resize(0);
    getline(infile, line);
  }
}

int get_filename_and_read_input(int argc, char** argv,
				 MPI_Comm comm, int localProc,
				 std::vector<std::string>& stdstrings)
{
  std::string filename;
  if (localProc == 0) {
    filename = construct_filename(argc, argv);
  }

  try {
    read_input_file(filename.c_str(), comm, stdstrings);
  }
  catch(std::runtime_error& exc) {
    fei::console_out() << exc.what() << FEI_ENDL;
    ERReturn(-1);
  }

  return(0);
}

void read_input_file(const char* filename,
		      MPI_Comm comm,
		      std::vector<std::string>& file_contents)
{
  int localProc =0;
#ifndef FEI_SER
  int numProcs = 1;
  MPI_Comm_rank(comm, &localProc);
  MPI_Comm_size(comm, &numProcs);
#endif

  if (localProc == 0) {
    read_file_lines_into_strings(filename, file_contents);

#ifndef FEI_SER
    if (numProcs > 1) {
      int num = file_contents.size();

      MPI_Bcast(&num, 1, MPI_INT, 0, comm);

      for(int j=0; j<num; ++j) {
        const char* cptr = file_contents[j].c_str();
        int length = strlen(cptr) + 1;
        MPI_Bcast(&length, 1, MPI_INT, 0, comm);
        MPI_Bcast((void*)cptr, length, MPI_CHAR, 0, comm);
      }
    }
#endif
  }
#ifndef FEI_SER
  else {//localProc != 0
    int num = -1;
    MPI_Bcast(&num, 1, MPI_INT, 0, comm);

    file_contents.resize(0);

    for(int j=0; j<num; ++j) {
      int length = 0;
      MPI_Bcast(&length, 1, MPI_INT, 0, comm);
      char* newstring = new char[length];
      MPI_Bcast(newstring, length, MPI_CHAR, 0, comm);
      std::string strg(newstring);
      file_contents.push_back(strg);
      delete [] newstring;
    }
  }
#endif
}

double get_file_benchmark(const char* filename,
			   const char* testname)
{
  std::vector<std::string> file_contents;
  read_file_lines_into_strings(filename, file_contents);

  double file_benchmark = 0.0;
  int err2 = snl_fei::getDoubleParamValue(testname,
					  file_contents,
					  file_benchmark);
  if (err2 != 0) {
    throw std::runtime_error("fei_test_utils::get_file_benchmark failed to find named benchmark");
  }

  return(file_benchmark);
}

bool within_percentage_margin(double value1,
			       double value2,
			       unsigned margin)
{
  double maxval = std::abs(value1) > std::abs(value2) ? std::abs(value1) : std::abs(value2);
  if (maxval < 1.e-14) {
    //They're both close enough to zero, return true without dividing.
    return(true);
  }

  double difference = std::abs(value2 - value1);
  double lhs = 100.0*(difference/maxval);
  double rhs = 1.0*margin;
  return( lhs <= rhs );
}

int whichArg(int argc, const char*const* argv, const char* findarg)
{
  for(int i=0; i<argc; i++) {
    if (argv[i] != NULL) {
      if (strcmp(findarg, argv[i]) == 0) return(i);
    }
  }
  return(-1);
}

bool check_and_cout_test_result(std::string testname,
				 double value,
				 double file_value,
				 unsigned margin)
{
  bool result = within_percentage_margin(value, file_value, margin);
  FEI_COUT << testname << " " << value << " ";
  if (!result) {
    FEI_COUT << "NOT ";
  }
  FEI_COUT << "within "<<margin<<"% of file value " << file_value << ". ";

  if (result) {
    FEI_COUT << "passed."<<FEI_ENDL;
  }
  else {
    FEI_COUT << "FAILED."<<FEI_ENDL;
  }
  return(result);
}

std::string check_test_result(double value,
				double goldvalue,
				unsigned margin)
{
  bool result = within_percentage_margin(value, goldvalue, margin);
  return( result ? "passed" : "FAILED");
}

int compare_with_file_benchmark(const char* name,
				 double benchmark,
				 const char* filename)
{
  if (name == NULL) return(-1);
  int returnValue = 0;

#if defined(FEI_PLATFORM) && defined(FEI_OPT_LEVEL)

  FEI_OSTRINGSTREAM testname;
  testname << name<<"_"<<FEI_PLATFORM<<"_"<<FEI_OPT_LEVEL;

  double file_benchmark = 0.0;
  bool file_benchmark_available = true;
  try {
    file_benchmark = get_file_benchmark(filename, testname.str().c_str());
  }
  catch (std::runtime_error& exc) {
    file_benchmark_available = false;
  }

  if (file_benchmark_available) {
    bool test_passed = check_and_cout_test_result(testname.str(),
							   benchmark,
							   file_benchmark,
							   10);
    returnValue = test_passed ? 0 : 1;
    if (returnValue != 0) {
      return(returnValue);
    }
  }

#else
  FEI_COUT << "(compile with -DFEI_PLATFORM=<platform> -DFEI_OPT_LEVEL=<opt|dbg>)"<<FEI_ENDL;
#endif

  return(returnValue);
}

int dirname(const char* name, char*& dir)
{
  //given a string which is a file-name, create a new string 'dir' and
  //populate it with the directory-name, or file-name minus the last '/file'
  //section.

  int i, len = strlen(name);
  dir = new char[len];
  for(i=0; i<len; ++i) {
    dir[i] = name[i];
  }

  i = len-1;
  if (i > 1) {
    while(i>0) {
      if (dir[i] == '/') {
	dir[i] = '\0';
	break;
      }
      else {
	dir[i] = '\0';
	--i;
      }
    }
  }

  if (i==0) dir[i] = '.';

  return(0);
}

void print_args(int argc, char** argv){
    FEI_COUT << "argc: " << argc << FEI_ENDL;

    for(int i=0; i<argc; i++){
      if (argv[i] != NULL) {
        FEI_COUT << "argv["<<i<<"]: " << argv[i] << FEI_ENDL;
      }
    }
    FEI_COUT << FEI_ENDL;
}

int compareMatrices(fei::FillableMat& mat1, fei::FillableMat& mat2, double tol)
{
  if (mat1 == mat2) return 0;
  FEI_COUT << "compareMatrices returned not-equal, tol=="<<tol << FEI_ENDL;
  return 1;
}

int readMatrix(const char* baseName, int np, fei::FillableMat& matrix)
{
  for(int i=0; i<np; i++) {
    FEI_OSTRINGSTREAM fileName;
    fileName <<baseName<<"."<<np<<"."<<i;
    FEI_IFSTREAM infile(fileName.str().c_str());

    infile.setf(IOS_SCIENTIFIC, IOS_FLOATFIELD);

    int row, col, tmp, numRows, numCols;
    double value;
    infile >> numRows;
    infile >> numCols;
    infile >> tmp;

    infile >> row;
    while(!infile.eof()) {
      infile >> col;
      infile >> value;

      matrix.putCoef(row, col, value);

      infile >> row;
    }
  }

  return(0);
}

int readMatrix(const char* fileName, fei::FillableMat& matrix)
{
  matrix.clear();

  FEI_IFSTREAM infile(fileName);
  if (!infile) {
    FEI_COUT << "ERROR opening file " << fileName << FEI_ENDL;
    return(1);
  }

  infile.setf(IOS_SCIENTIFIC, IOS_FLOATFIELD);

  int row, col, tmp, numRows, numCols;
  double value;
  infile >> numRows;
  infile >> numCols;
  infile >> tmp;

  infile >> row;
  while(!infile.eof()) {
    infile >> col;
    infile >> value;

    matrix.putCoef(row, col, value);

    infile >> row;
  }

  return(0);
}

int writeMatrix(const char* fileName, fei::FillableMat& matrix)
{
  FEI_OFSTREAM outfile(fileName);
  if (!outfile) {
    FEI_COUT << "ERROR opening file " << fileName << FEI_ENDL;
    return(1);
  }

  outfile.setf(IOS_SCIENTIFIC, IOS_FLOATFIELD);

  outfile << matrix;

  return(0);
}

int copy_feiMatrix_to_FillableMat(fei::Matrix& feimat, fei::FillableMat& fmat)
{
  fmat.clear();

  fei::SharedPtr<fei::VectorSpace> rowspace =
    feimat.getMatrixGraph()->getRowSpace();

  MPI_Comm comm = rowspace->getCommunicator();
  int localProc = fei::localProc(comm);

  std::vector<int> globalOffsets;
  rowspace->getGlobalIndexOffsets(globalOffsets);
  int numLocalRows = globalOffsets[localProc+1]-globalOffsets[localProc];
  int firstLocalRow = globalOffsets[localProc];

  for(int i=0; i<numLocalRows; ++i) {
    int row = firstLocalRow+i;
    int rowLen;
    feimat.getRowLength(row, rowLen);

    std::vector<int> colindices(rowLen);
    std::vector<double> coefs(rowLen);

    feimat.copyOutRow(row, rowLen, &coefs[0], &colindices[0]);

    fmat.putRow(row, &colindices[0], &coefs[0], rowLen);
  }

  return(0);
}

} //namespace fei_test_utils

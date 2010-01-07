#include <Tpetra_MatrixIO.hpp>
#include <Kokkos_ConfigDefs.hpp>
#include <Kokkos_SerialNode.hpp>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <algorithm>
#include <iterator>
#include <exception>
#include <string>
#ifdef HAVE_KOKKOS_TBB
#  include <Kokkos_TBBNode.hpp>
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
#  include <Kokkos_TPINode.hpp>
#endif
#ifdef HAVE_KOKKOS_THRUST
#  include <Kokkos_ThrustGPUNode.hpp>
#endif

#define EXP_INST_NODE_SCALAR_ORDINAL(NODE,SCALAR,ORDINAL) \
  template                                                \
  void \
  Tpetra::Utils::readHBMatrix<SCALAR,ORDINAL,ORDINAL,NODE,Kokkos::DefaultSparseMultiply<SCALAR,ORDINAL,NODE >,Kokkos::DefaultSparseSolve<SCALAR,ORDINAL,NODE > >(     \
      const std::string &, const Teuchos::RCP<const Teuchos::Comm<int> > &, const Teuchos::RCP<NODE > &, \
      Teuchos::RCP< Tpetra::CrsMatrix<SCALAR,ORDINAL,ORDINAL,NODE,Kokkos::DefaultSparseMultiply<SCALAR,ORDINAL,NODE >,Kokkos::DefaultSparseSolve<SCALAR,ORDINAL,NODE > > > &);

#define EXP_INST_NODE_ORDINAL(NODE,ORDINAL) \
  EXP_INST_NODE_SCALAR_ORDINAL(NODE,double,ORDINAL) \
  EXP_INST_NODE_SCALAR_ORDINAL(NODE,float ,ORDINAL)

EXP_INST_NODE_ORDINAL(Kokkos::SerialNode,int)

#ifdef HAVE_KOKKOS_TBB
EXP_INST_NODE_ORDINAL(Kokkos::TBBNode,int)
#endif

#ifdef HAVE_KOKKOS_THREADPOOL
EXP_INST_NODE_ORDINAL(Kokkos::TPINode,int)
#endif

#ifdef HAVE_KOKKOS_THRUST
#ifdef HAVE_KOKKOS_CUDA_DOUBLE
EXP_INST_NODE_SCALAR_ORDINAL(Kokkos::ThrustGPUNode,double,ORDINAL)
#endif
#ifdef HAVE_KOKKOS_CUDA_FLOAT
EXP_INST_NODE_SCALAR_ORDINAL(Kokkos::ThrustGPUNode,float ,ORDINAL)
#endif
#endif

bool Tpetra::Utils::parseIfmt(Teuchos::ArrayRCP<char> fmt, int &perline, int &width) {
  TEST_FOR_EXCEPT(fmt.size() != 0 && fmt[fmt.size()-1] != '\0');
  // parses integers n and d out of (nId)
  bool error = true;
  std::transform(fmt.begin(), fmt.end(), fmt, static_cast < int(*)(int) > (std::toupper));
  if (std::sscanf(fmt.getRawPtr(),"(%dI%d)",&perline,&width) == 2) {
    error = false;
  }
  return error;
}

bool Tpetra::Utils::parseRfmt(Teuchos::ArrayRCP<char> fmt, int &perline, int &width, int &prec, char &valformat) {
  TEST_FOR_EXCEPT(fmt.size() != 0 && fmt[fmt.size()-1] != '\0');
  std::transform(fmt.begin(), fmt.end(), fmt, static_cast < int(*)(int) > (std::toupper));
  // find the first left paren '(' and the last right paren ')'
  Teuchos::ArrayRCP<char>::iterator firstLeftParen = std::find( fmt.begin(),  fmt.end(), '(');
  Teuchos::ArrayRCP<char>::iterator lastRightParen = std::find(std::reverse_iterator<Teuchos::ArrayRCP<char>::iterator>(fmt.end()), 
                                                               std::reverse_iterator<Teuchos::ArrayRCP<char>::iterator>(fmt.begin()), 
                                                               ')').base()-1;
  // select the substring between the parens, including them
  // if neither was found, set the string to empty
  if (firstLeftParen == fmt.end() || lastRightParen == fmt.begin()) {
    fmt.resize(0 + 1);
    fmt[0] = '\0';
  }
  else {
    fmt += (firstLeftParen - fmt.begin());
    size_t newLen = lastRightParen - firstLeftParen + 1;
    fmt.resize(newLen + 1);
    fmt[newLen] = '\0';
  }
  if (std::find(fmt.begin(),fmt.end(),'P') != fmt.end()) {
    // not supported
    return true;
  }
  bool error = true;
  if (std::sscanf(fmt.getRawPtr(),"(%d%c%d.%d)",&perline,&valformat,&width,&prec) == 4) {
    if (valformat == 'E' || valformat == 'D' || valformat == 'F') {
      error = false;
    }
  } 
  return error;
}


void Tpetra::Utils::readHBHeader(std::ifstream &fin, Teuchos::ArrayRCP<char> &Title, Teuchos::ArrayRCP<char> &Key, Teuchos::ArrayRCP<char> &Type, 
                           int &Nrow, int &Ncol, int &Nnzero, int &Nrhs,
                           Teuchos::ArrayRCP<char> &Ptrfmt, Teuchos::ArrayRCP<char> &Indfmt, Teuchos::ArrayRCP<char> &Valfmt, Teuchos::ArrayRCP<char> &Rhsfmt, 
                           int &Ptrcrd, int &Indcrd, int &Valcrd, int &Rhscrd, Teuchos::ArrayRCP<char> &Rhstype) {
  int Totcrd, Neltvl, Nrhsix;
  const int MAXLINE = 81;
  char line[MAXLINE];
  //
  Title.resize(72 + 1);  std::fill(Title.begin(),  Title.end(),  '\0');
  Key.resize(8 + 1);     std::fill(Key.begin(),    Key.end(),    '\0');
  Type.resize(3 + 1);    std::fill(Type.begin(),   Type.end(),   '\0');
  Ptrfmt.resize(16 + 1); std::fill(Ptrfmt.begin(), Ptrfmt.end(), '\0');
  Indfmt.resize(16 + 1); std::fill(Indfmt.begin(), Indfmt.end(), '\0');
  Valfmt.resize(20 + 1); std::fill(Valfmt.begin(), Valfmt.end(), '\0');
  Rhsfmt.resize(20 + 1); std::fill(Rhsfmt.begin(), Rhsfmt.end(), '\0');
  //
  const std::string errStr("Tpetra::Utils::readHBHeader(): Improperly formatted H/B file.");
  /*  First line:   (A72,A8) */
  fin.getline(line,MAXLINE);
  TEST_FOR_EXCEPTION( std::sscanf(line,"%*s") < 0, std::runtime_error, errStr);
  (void)std::sscanf(line, "%72c%8[^\n]", Title.getRawPtr(), Key.getRawPtr());
  /*  Second line:  (5I14) */
  fin.getline(line,MAXLINE);
  TEST_FOR_EXCEPTION(std::sscanf(line,"%*s") < 0, std::runtime_error, errStr);
  TEST_FOR_EXCEPTION(std::sscanf(line,"%14d%14d%14d%14d%14d",&Totcrd,&Ptrcrd,&Indcrd,&Valcrd,&Rhscrd) != 5, std::runtime_error, errStr);
  /*  Third line:   (A3, 11X, 4I14) */
  fin.getline(line,MAXLINE);
  TEST_FOR_EXCEPTION(std::sscanf(line,"%*s") < 0, std::runtime_error, errStr);
  TEST_FOR_EXCEPTION(std::sscanf(line, "%3c%14i%14i%14i%14i", Type.getRawPtr(),&Nrow,&Ncol,&Nnzero,&Neltvl) != 5 , std::runtime_error, errStr);
  std::transform(Type.begin(), Type.end(), Type.begin(), static_cast < int(*)(int) > (std::toupper));
  /*  Fourth line:  */
  fin.getline(line,MAXLINE);
  TEST_FOR_EXCEPTION(std::sscanf(line,"%*s") < 0, std::runtime_error, errStr);
  TEST_FOR_EXCEPTION(std::sscanf(line,"%16c%16c%20c%20c",Ptrfmt.getRawPtr(),Indfmt.getRawPtr(),Valfmt.getRawPtr(),Rhsfmt.getRawPtr()) != 4, std::runtime_error, errStr);
  /*  (Optional) Fifth line: */
  if (Rhscrd != 0 ) { 
    Rhstype.resize(3 + 1,'\0');
    fin.getline(line,MAXLINE);
    TEST_FOR_EXCEPTION(std::sscanf(line,"%*s") < 0, std::runtime_error, errStr);
    TEST_FOR_EXCEPTION(std::sscanf(line,"%3c%14d%14d", Rhstype.getRawPtr(), &Nrhs, &Nrhsix) != 3, std::runtime_error, errStr);
  }
}


void Tpetra::Utils::readHBInfo(const std::string &filename, int &M, int &N, int &nz, Teuchos::ArrayRCP<char> &Type, int &Nrhs) {
  std::ifstream fin;
  int Ptrcrd, Indcrd, Valcrd, Rhscrd; 
  int Nrow, Ncol, Nnzero;
  Teuchos::ArrayRCP<char> Title, Key, Rhstype, Ptrfmt, Indfmt, Valfmt, Rhsfmt;
  try {
    fin.open(filename.c_str(),std::ifstream::in);
    Tpetra::Utils::readHBHeader(fin, Title, Key, Type, Nrow, Ncol, Nnzero, Nrhs,
                          Ptrfmt, Indfmt, Valfmt, Rhsfmt, 
                          Ptrcrd, Indcrd, Valcrd, Rhscrd, Rhstype);
    fin.close();
  }
  catch (std::exception &e) {
    TEST_FOR_EXCEPTION(true, std::runtime_error, 
        "Tpetra::Utils::readHBInfo() of filename \"" << filename << "\" caught exception: " << std::endl
        << e.what() << std::endl);
  }
}


void Tpetra::Utils::readHBMatDouble(const std::string &filename, int &numRows, int &numCols, int &numNZ, std::string &type, Teuchos::ArrayRCP<int> &colPtrs, Teuchos::ArrayRCP<int> &rowInds, Teuchos::ArrayRCP<double> &vals) {
  // NOTE: if complex, allocate enough space for 2*NNZ and interleave real and imaginary parts (real,imag)
  //       if pattern, don't touch parameter vals; do not allocate space space for values
  try {
    std::ifstream fin;
    int ptrCrd, indCrd, valCrd;
    Teuchos::ArrayRCP<char> Title, Key, Rhstype, Ptrfmt, Indfmt, Valfmt, Rhsfmt;
    const int MAXSIZE = 81;
    char lineBuf[MAXSIZE];
    // nitty gritty
    int ptrsPerLine, ptrWidth, indsPerLine, indWidth, valsPerLine, valWidth, valPrec;
    char valFlag;
    // 
    fin.open(filename.c_str(),std::ifstream::in);
    {
      // we don't care about RHS-related stuff, so declare those vars in an expiring scope
      int Nrhs, rhsCrd;
      Teuchos::ArrayRCP<char> Rhstype, Rhsfmt;
      Teuchos::ArrayRCP<char> TypeArray;
      Tpetra::Utils::readHBHeader(fin, Title, Key, TypeArray, numRows, numCols, numNZ, Nrhs,
                                  Ptrfmt, Indfmt, Valfmt, Rhsfmt, 
                                  ptrCrd, indCrd, valCrd, rhsCrd, Rhstype);
      if (TypeArray.size() > 0) {
        type.resize(TypeArray.size()-1);
        std::copy(TypeArray.begin(), TypeArray.end(), type.begin());
      }
    }
    const std::string errStr("Tpetra::Utils::readHBHeader(): Improperly formatted H/B file.");
    const bool readPatternOnly = (type[0] == 'P' || type[0] == 'p');
    const bool readComplex     = (type[0] == 'C' || type[0] == 'c');
    /*  Parse the array input formats from Line 3 of HB file  */
    TEST_FOR_EXCEPTION( Tpetra::Utils::parseIfmt(Ptrfmt,ptrsPerLine,ptrWidth) == true, std::runtime_error,
        "Tpetra::Utils::readHBMatDouble(): error parsing. Invalid/unsupported file format.");
    TEST_FOR_EXCEPTION( Tpetra::Utils::parseIfmt(Indfmt,indsPerLine,indWidth) == true, std::runtime_error,
        "Tpetra::Utils::readHBMatDouble(): error parsing. Invalid/unsupported file format.");
    if (readPatternOnly == false) {
      TEST_FOR_EXCEPTION( Tpetra::Utils::parseRfmt(Valfmt,valsPerLine,valWidth,valPrec,valFlag) == true, std::runtime_error,
          "Tpetra::Utils::readHBMatDouble(): error parsing. Invalid/unsupported file format.");
    }
    // special case this: the reason is that the number of colPtrs read is numCols+1, which is non-zero even if numCols == 0
    // furthermore, if numCols == 0, there is nothing of interest to read
    if (numCols == 0) return;
    // allocate space for column pointers, row indices, and values
    // if the file is empty, do not touch these ARCPs
    colPtrs = Teuchos::arcp<int>(numCols+1);
    if (numNZ > 0) {
      rowInds = Teuchos::arcp<int>(numNZ);  
      if (readPatternOnly == false) {
        if (readComplex) {
          vals = Teuchos::arcp<double>(2*numNZ); 
        }
        else {
          vals = Teuchos::arcp<double>(numNZ); 
        }
      }
    }
    /* Read column pointer array:
       Specifically, read ptrCrd number of lines, and on each line, read ptrsPerLine number of integers, each of width ptrWidth
       Store these in colPtrs */
    {
      int colPtrsRead = 0;
      char NullSub = '\0';
      for (int lno=0; lno < ptrCrd; ++lno) {
        fin.getline(lineBuf, MAXSIZE);
        TEST_FOR_EXCEPTION(std::sscanf(lineBuf,"%*s") < 0, std::runtime_error, errStr);
        char *linePtr = lineBuf;
        for (int ptr=0; ptr < ptrsPerLine; ++ptr) {
          if (colPtrsRead == numCols + 1) break;
          int cptr;
          // terminate the string at the end of the current ptr block, saving the character in that location
          std::swap(NullSub,linePtr[ptrWidth]);
          // read the ptr
          std::sscanf(linePtr, "%d", &cptr);
          // put the saved character back, and put the '\0' back into NullSub for use again
          std::swap(NullSub,linePtr[ptrWidth]);
          linePtr += ptrWidth;
          colPtrs[colPtrsRead++] = cptr;
        }
      }
      TEST_FOR_EXCEPT(colPtrsRead != numCols + 1);
    }
    /* Read row index array:
       Specifically, read indCrd number of lines, and on each line, read indsPerLine number of integers, each of width indWidth
       Store these in rowInds */
    {
      char NullSub = '\0';
      int indicesRead = 0;
      for (int lno=0; lno < indCrd; ++lno) {
        fin.getline(lineBuf, MAXSIZE);
        TEST_FOR_EXCEPTION(std::sscanf(lineBuf,"%*s") < 0, std::runtime_error, errStr);
        char *linePtr = lineBuf;
        for (int ind=0; ind < indsPerLine; ++ind) {
          if (indicesRead == numNZ) break;
          int ind;
          // terminate the string at the end of the current ind block, saving the character in that location
          std::swap(NullSub,linePtr[indWidth]);
          // read the ind
          std::sscanf(linePtr, "%d", &ind);
          // put the saved character back, and put the '\0' back into NullSub for use again
          std::swap(NullSub,linePtr[indWidth]);
          linePtr += indWidth;
          rowInds[indicesRead++] = ind;
        }
      }
      TEST_FOR_EXCEPT(indicesRead != numNZ);
    }
    /* Read array of values:
       Specifically, read valCrd number of lines, and on each line, read valsPerLine number of real values, each of width/precision valWidth/valPrec
       Store these in vals
       If readComplex, then read twice as many non-zeros, and interleave real,imag into vals */
    if (readPatternOnly == false) {
      int totalNumVals;
      if (readComplex) {
        totalNumVals = 2*numNZ;
      }
      else {
        totalNumVals = numNZ;
      }
      char NullSub = '\0';
      int valsRead = 0;
      for (int lno=0; lno < valCrd; ++lno) {
        fin.getline(lineBuf, MAXSIZE);
        TEST_FOR_EXCEPTION(std::sscanf(lineBuf,"%*s") < 0, std::runtime_error, errStr);
        // if valFlag == 'D', then we need to convert [dD] in fp vals into [eE] that scanf can parse
        if (valFlag == 'D') std::replace_if(lineBuf, lineBuf+MAXSIZE, std::bind2nd(std::equal_to<char>(), 'D'), 'E'); 
        char *linePtr = lineBuf;
        for (int val=0; val < valsPerLine; ++val) {
          if (valsRead == totalNumVals) break;
          double val;
          // terminate the string at the end of the current val block, saving the character in that location
          std::swap(NullSub,linePtr[valWidth]);
          // read the val
          std::sscanf(linePtr, "%le", &val);
          // put the saved character back, and put the '\0' back into NullSub for use again
          std::swap(NullSub,linePtr[valWidth]);
          linePtr += valWidth;
          vals[valsRead++] = val;
        }
      }
      TEST_FOR_EXCEPT(valsRead != totalNumVals);
    }
    fin.close();
  }
  catch (std::exception &e) {
    TEST_FOR_EXCEPTION(true, std::runtime_error, 
        "Tpetra::Utils::readHBInfo() of filename \"" << filename << "\" caught exception: " << std::endl
        << e.what() << std::endl);
  }
}

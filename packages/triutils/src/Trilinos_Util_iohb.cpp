/*
Fri Aug 15 16:29:47 EDT 1997

                       Harwell-Boeing File I/O in C
                                V. 1.0

           National Institute of Standards and Technology, MD.
                             K.A. Remington

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                NOTICE

 Permission to use, copy, modify, and distribute this software and
 its documentation for any purpose and without fee is hereby granted
 provided that the above copyright notice appear in all copies and
 that both the copyright notice and this permission notice appear in
 supporting documentation.

 Neither the Author nor the Institution (National Institute of Standards
 and Technology) make any representations about the suitability of this
 software for any purpose. This software is provided "as is" without
 expressed or implied warranty.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                         ---------------------
                         INTERFACE DESCRIPTION
                         ---------------------
  ---------------
  QUERY FUNCTIONS
  ---------------

  FUNCTION:

  int readHB_info(const char *filename, int *M, int *N, int *nz,
  char **Type, int *Nrhs)

  DESCRIPTION:

  The readHB_info function opens and reads the header information from
  the specified Harwell-Boeing file, and reports back the number of rows
  and columns in the stored matrix (M and N), the number of nonzeros in
  the matrix (nz), the 3-character matrix type(Type), and the number of
  right-hand-sides stored along with the matrix (Nrhs).  This function
  is designed to retrieve basic size information which can be used to
  allocate arrays.

  FUNCTION:

  int  readHB_header(std::FILE* in_file, char* Title, char* Key, char* Type,
                    int* Nrow, int* Ncol, int* Nnzero, int* Nrhs,
                    char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
                    int* Ptrcrd, int* Indcrd, int* Valcrd, int* Rhscrd,
                    char *Rhstype)

  DESCRIPTION:

  More detailed than the readHB_info function, readHB_header() reads from
  the specified Harwell-Boeing file all of the header information.


  ------------------------------
  DOUBLE PRECISION I/O FUNCTIONS
  ------------------------------
  FUNCTION:

  int readHB_newmat_double(const char *filename, int *M, int *N, *int nz,
  int **colptr, int **rowind,  double**val)

  int readHB_mat_double(const char *filename, int *colptr, int *rowind,
  double*val)


  DESCRIPTION:

  This function opens and reads the specified file, interpreting its
  contents as a sparse matrix stored in the Harwell/Boeing standard
  format.  (See readHB_aux_double to read auxillary vectors.)
        -- Values are interpreted as double precision numbers. --

  The "mat" function uses _pre-allocated_ vectors to hold the index and
  nonzero value information.

  The "newmat" function allocates vectors to hold the index and nonzero
  value information, and returns pointers to these vectors along with
  matrix dimension and number of nonzeros.

  FUNCTION:

  int readHB_aux_double(const char* filename, const char AuxType, double b[])

  int readHB_newaux_double(const char* filename, const char AuxType, double** b)

  DESCRIPTION:

  This function opens and reads from the specified file auxillary vector(s).
  The char argument Auxtype determines which type of auxillary vector(s)
  will be read (if present in the file).

                  AuxType = 'F'   right-hand-side
                  AuxType = 'G'   initial estimate (Guess)
                  AuxType = 'X'   eXact solution

  If Nrhs > 1, all of the Nrhs vectors of the given type are read and
  stored in column-major order in the vector b.

  The "newaux" function allocates a vector to hold the values retrieved.
  The "mat" function uses a _pre-allocated_ vector to hold the values.

  FUNCTION:

  int writeHB_mat_double(const char* filename, int M, int N,
                        int nz, const int colptr[], const int rowind[],
                        const double val[], int Nrhs, const double rhs[],
                        const double guess[], const double exact[],
                        const char* Title, const char* Key, const char* Type,
                        char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
                        const char* Rhstype)

  DESCRIPTION:

  The writeHB_mat_double function opens the named file and writes the specified
  matrix and optional auxillary vector(s) to that file in Harwell-Boeing
  format.  The format arguments (Ptrfmt,Indfmt,Valfmt, and Rhsfmt) are
  character strings specifying "Fortran-style" output formats -- as they
  would appear in a Harwell-Boeing file.  They are used to produce output
  which is as close as possible to what would be produced by Fortran code,
  but note that "D" and "P" edit descriptors are not supported.
  If NULL, the following defaults will be used:
                    Ptrfmt = Indfmt = "(8I10)"
                    Valfmt = Rhsfmt = "(4E20.13)"

  -----------------------
  CHARACTER I/O FUNCTIONS
  -----------------------
  FUNCTION:

  int readHB_mat_char(const char* filename, int colptr[], int rowind[],
                                           char val[], char* Valfmt)
  int readHB_newmat_char(const char* filename, int* M, int* N, int* nonzeros,
                          int** colptr, int** rowind, char** val, char** Valfmt)

  DESCRIPTION:

  This function opens and reads the specified file, interpreting its
  contents as a sparse matrix stored in the Harwell/Boeing standard
  format.  (See readHB_aux_char to read auxillary vectors.)
              -- Values are interpreted as char strings.     --
  (Used to translate exact values from the file into a new storage format.)

  The "mat" function uses _pre-allocated_ arrays to hold the index and
  nonzero value information.

  The "newmat" function allocates char arrays to hold the index
  and nonzero value information, and returns pointers to these arrays
  along with matrix dimension and number of nonzeros.

  FUNCTION:

  int readHB_aux_char(const char* filename, const char AuxType, char b[])
  int readHB_newaux_char(const char* filename, const char AuxType, char** b,
                         char** Rhsfmt)

  DESCRIPTION:

  This function opens and reads from the specified file auxillary vector(s).
  The char argument Auxtype determines which type of auxillary vector(s)
  will be read (if present in the file).

                  AuxType = 'F'   right-hand-side
                  AuxType = 'G'   initial estimate (Guess)
                  AuxType = 'X'   eXact solution

  If Nrhs > 1, all of the Nrhs vectors of the given type are read and
  stored in column-major order in the vector b.

  The "newaux" function allocates a character array to hold the values
                retrieved.
  The "mat" function uses a _pre-allocated_ array to hold the values.

  FUNCTION:

  int writeHB_mat_char(const char* filename, int M, int N,
                        int nz, const int colptr[], const int rowind[],
                        const char val[], int Nrhs, const char rhs[],
                        const char guess[], const char exact[],
                        const char* Title, const char* Key, const char* Type,
                        char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
                        const char* Rhstype)

  DESCRIPTION:

  The writeHB_mat_char function opens the named file and writes the specified
  matrix and optional auxillary vector(s) to that file in Harwell-Boeing
  format.  The format arguments (Ptrfmt,Indfmt,Valfmt, and Rhsfmt) are
  character strings specifying "Fortran-style" output formats -- as they
  would appear in a Harwell-Boeing file.  Valfmt and Rhsfmt must accurately
  represent the character representation of the values stored in val[]
  and rhs[].

  If NULL, the following defaults will be used for the integer vectors:
                    Ptrfmt = Indfmt = "(8I10)"
                    Valfmt = Rhsfmt = "(4E20.13)"


*/

/*---------------------------------------------------------------------*/
/* If zero-based indexing is desired, _SP_base should be set to 0      */
/* This will cause indices read from H-B files to be decremented by 1  */
/*             and indices written to H-B files to be incremented by 1 */
/*            <<<  Standard usage is _SP_base = 1  >>>                 */
#ifndef _SP_base
#define _SP_base 1
#endif
/*---------------------------------------------------------------------*/

#include "Trilinos_Util_iohb.h"

#include<cstring>
#include<cmath>
#include <cstdlib>
using std::malloc;
using std::free;
using std::size_t;

char* substr(const char* S, const int pos, const int len);
void upcase(char* S);
void IOHBTerminate(const char* message);

int readHB_info(const char* filename, int* M, int* N, int* nz, char** Type,
                                                      int* Nrhs)
{
/****************************************************************************/
/*  The readHB_info function opens and reads the header information from    */
/*  the specified Harwell-Boeing file, and reports back the number of rows  */
/*  and columns in the stored matrix (M and N), the number of nonzeros in   */
/*  the matrix (nz), and the number of right-hand-sides stored along with   */
/*  the matrix (Nrhs).                                                      */
/*                                                                          */
/*  For a description of the Harwell Boeing standard, see:                  */
/*            Duff, et al.,  ACM TOMS Vol.15, No.1, March 1989              */
/*                                                                          */
/*    ----------                                                            */
/*    **CAVEAT**                                                            */
/*    ----------                                                            */
/*  **  If the input file does not adhere to the H/B format, the  **        */
/*  **             results will be unpredictable.                 **        */
/*                                                                          */
/****************************************************************************/
    std::FILE *in_file;
    int Ptrcrd, Indcrd, Valcrd, Rhscrd;
    int Nrow, Ncol, Nnzero;
    char *mat_type;
    char Title[73], Key[9], Rhstype[4];
    char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];

    mat_type = (char *) malloc(4);
    if ( mat_type == NULL ) IOHBTerminate("Insufficient memory for mat_type\n");

    if ( (in_file = std::fopen( filename, "r")) == NULL ) {
       std::fprintf(stderr,"Error: Cannot open file: %s\n",filename);
       return 0;
    }

    readHB_header(in_file, Title, Key, mat_type, &Nrow, &Ncol, &Nnzero, Nrhs,
                  Ptrfmt, Indfmt, Valfmt, Rhsfmt,
                  &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);
    std::fclose(in_file);
    *Type = mat_type;
    *(*Type+3) = '\0';
    *M    = Nrow;
    *N    = Ncol;
    *nz   = Nnzero;
    if (Rhscrd == 0) {*Nrhs = 0;}

/*  In verbose mode, print some of the header information:   */
/*
    if (verbose == 1)
    {
        printf("Reading from Harwell-Boeing file %s (verbose on)...\n",filename);
        printf("  Title: %s\n",Title);
        printf("  Key:   %s\n",Key);
        printf("  The stored matrix is %i by %i with %i nonzeros.\n",
                *M, *N, *nz );
        printf("  %i right-hand--side(s) stored.\n",*Nrhs);
    }
*/

    return 1;

}



int readHB_header(std::FILE* in_file, char* Title, char* Key, char* Type,
                    int* Nrow, int* Ncol, int* Nnzero, int* Nrhs,
                    char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
                    int* Ptrcrd, int* Indcrd, int* Valcrd, int* Rhscrd,
                    char *Rhstype)
{
/*************************************************************************/
/*  Read header information from the named H/B file...                   */
/*************************************************************************/
    int Totcrd,Neltvl,Nrhsix;
    char line[BUFSIZ];

/*  First line:   */
    if (std::fgets(line, BUFSIZ, in_file) == NULL) {
      std::fprintf(stderr,"Error: Failed to read from file.\n");
      return 0;
    }
    if ( std::sscanf(line,"%*s") < 0 )
        IOHBTerminate("Trilinos_Util_iohb.cpp: Null (or blank) first line of HB file.\n");
    (void) std::sscanf(line, "%72c%8[^\n]", Title, Key);
    *(Key+8) = '\0';
    *(Title+72) = '\0';

/*  Second line:  */
    if (std::fgets(line, BUFSIZ, in_file) == NULL) {
      std::fprintf(stderr,"Error: Failed to read from file.\n");
      return 0;
    }
    if ( std::sscanf(line,"%*s") < 0 )
        IOHBTerminate("Trilinos_Util_iohb.cpp: Null (or blank) second line of HB file.\n");
    if ( std::sscanf(line,"%i",&Totcrd) != 1) Totcrd = 0;
    if ( std::sscanf(line,"%*i%i",Ptrcrd) != 1) *Ptrcrd = 0;
    if ( std::sscanf(line,"%*i%*i%i",Indcrd) != 1) *Indcrd = 0;
    if ( std::sscanf(line,"%*i%*i%*i%i",Valcrd) != 1) *Valcrd = 0;
    if ( std::sscanf(line,"%*i%*i%*i%*i%i",Rhscrd) != 1) *Rhscrd = 0;

/*  Third line:   */
    if (std::fgets(line, BUFSIZ, in_file) == NULL) {
      std::fprintf(stderr,"Error: Failed to read from file.\n");
      return 0;
    }
    if ( std::sscanf(line,"%*s") < 0 )
        IOHBTerminate("Trilinos_Util_iohb.cpp: Null (or blank) third line of HB file.\n");
    if ( std::sscanf(line, "%3c", Type) != 1)
        IOHBTerminate("Trilinos_Util_iohb.cpp: Invalid Type info, line 3 of Harwell-Boeing file.\n");
    upcase(Type);
    if ( std::sscanf(line,"%*3c%i",Nrow) != 1) *Nrow = 0 ;
    if ( std::sscanf(line,"%*3c%*i%i",Ncol) != 1) *Ncol = 0 ;
    if ( std::sscanf(line,"%*3c%*i%*i%i",Nnzero) != 1) *Nnzero = 0 ;
    if ( std::sscanf(line,"%*3c%*i%*i%*i%i",&Neltvl) != 1) Neltvl = 0 ;

/*  Fourth line:  */
    if (std::fgets(line, BUFSIZ, in_file) == NULL) {
      std::fprintf(stderr,"Error: Failed to read from file.\n");
      return 0;
    }
    if ( std::sscanf(line,"%*s") < 0 )
        IOHBTerminate("Trilinos_Util_iohb.cpp: Null (or blank) fourth line of HB file.\n");
    if ( std::sscanf(line, "%16c",Ptrfmt) != 1)
        IOHBTerminate("Trilinos_Util_iohb.cpp: Invalid format info, line 4 of Harwell-Boeing file.\n");
    if ( std::sscanf(line, "%*16c%16c",Indfmt) != 1)
        IOHBTerminate("Trilinos_Util_iohb.cpp: Invalid format info, line 4 of Harwell-Boeing file.\n");
    if ( std::sscanf(line, "%*16c%*16c%20c",Valfmt) != 1)
        IOHBTerminate("Trilinos_Util_iohb.cpp: Invalid format info, line 4 of Harwell-Boeing file.\n");
    std::sscanf(line, "%*16c%*16c%*20c%20c",Rhsfmt);
    *(Ptrfmt+16) = '\0';
    *(Indfmt+16) = '\0';
    *(Valfmt+20) = '\0';
    *(Rhsfmt+20) = '\0';

/*  (Optional) Fifth line: */
    if (*Rhscrd != 0 )
    {
       if (std::fgets(line, BUFSIZ, in_file) == NULL) {
         std::fprintf(stderr,"Error: Failed to read from file.\n");
         return 0;
       }
       if ( std::sscanf(line,"%*s") < 0 )
           IOHBTerminate("Trilinos_Util_iohb.cpp: Null (or blank) fifth line of HB file.\n");
       if ( std::sscanf(line, "%3c", Rhstype) != 1)
         IOHBTerminate("Trilinos_Util_iohb.cpp: Invalid RHS type information, line 5 of Harwell-Boeing file.\n");
       if ( std::sscanf(line, "%*3c%i", Nrhs) != 1) *Nrhs = 0;
       if ( std::sscanf(line, "%*3c%*i%i", &Nrhsix) != 1) Nrhsix = 0;
    }
    return 1;
}


int readHB_mat_double(const char* filename, int colptr[], int rowind[],
                                                                 double val[])
{
/****************************************************************************/
/*  This function opens and reads the specified file, interpreting its      */
/*  contents as a sparse matrix stored in the Harwell/Boeing standard       */
/*  format and creating compressed column storage scheme vectors to hold    */
/*  the index and nonzero value information.                                */
/*                                                                          */
/*    ----------                                                            */
/*    **CAVEAT**                                                            */
/*    ----------                                                            */
/*  Parsing real formats from Fortran is tricky, and this file reader       */
/*  does not claim to be foolproof.   It has been tested for cases when     */
/*  the real values are printed consistently and evenly spaced on each      */
/*  line, with Fixed (F), and Exponential (E or D) formats.                 */
/*                                                                          */
/*  **  If the input file does not adhere to the H/B format, the  **        */
/*  **             results will be unpredictable.                 **        */
/*                                                                          */
/****************************************************************************/
    std::FILE *in_file;
    int i,j,ind,col,offset,count,last,Nrhs;
    int Ptrcrd, Indcrd, Valcrd, Rhscrd;
    int Nrow, Ncol, Nnzero, Nentries;
    int Ptrperline, Ptrwidth, Indperline, Indwidth;
    int Valperline, Valwidth, Valprec;
    int Valflag;           /* Indicates 'E','D', or 'F' float format */
    char* ThisElement;
    char Title[73], Key[9], Type[4] = "XXX", Rhstype[4];
    char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];
    char line[BUFSIZ];

    if ( (in_file = std::fopen( filename, "r")) == NULL ) {
       std::fprintf(stderr,"Error: Cannot open file: %s\n",filename);
       return 0;
    }

    readHB_header(in_file, Title, Key, Type, &Nrow, &Ncol, &Nnzero, &Nrhs,
                  Ptrfmt, Indfmt, Valfmt, Rhsfmt,
                  &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);

/*  Parse the array input formats from Line 3 of HB file  */
    ParseIfmt(Ptrfmt,&Ptrperline,&Ptrwidth);
    ParseIfmt(Indfmt,&Indperline,&Indwidth);
    if ( Type[0] != 'P' ) {          /* Skip if pattern only  */
    ParseRfmt(Valfmt,&Valperline,&Valwidth,&Valprec,&Valflag);
    }

/*  Read column pointer array:   */

    offset = 1-_SP_base;  /* if base 0 storage is declared (via macro definition), */
                          /* then storage entries are offset by 1                  */

    ThisElement = (char *) malloc(Ptrwidth+1);
    if ( ThisElement == NULL ) IOHBTerminate("Insufficient memory for ThisElement.");
    *(ThisElement+Ptrwidth) = '\0';
    count=0;
    for (i=0;i<Ptrcrd;i++)
    {
       if (std::fgets(line, BUFSIZ, in_file) == NULL) {
         std::fprintf(stderr,"Error: Failed to read from file.\n");
         return 0;
       }
       if ( std::sscanf(line,"%*s") < 0 )
         IOHBTerminate("Trilinos_Util_iohb.cpp: Null (or blank) line in pointer data region of HB file.\n");
       col =  0;
       for (ind = 0;ind<Ptrperline;ind++)
       {
          if (count > Ncol) break;
          std::strncpy(ThisElement,line+col,Ptrwidth);
  /* ThisElement = substr(line,col,Ptrwidth); */
          colptr[count] = std::atoi(ThisElement)-offset;
          count++; col += Ptrwidth;
       }
    }
    free(ThisElement);

/*  Read row index array:  */

    ThisElement = (char *) malloc(Indwidth+1);
    if ( ThisElement == NULL ) IOHBTerminate("Insufficient memory for ThisElement.");
    *(ThisElement+Indwidth) = '\0';
    count = 0;
    for (i=0;i<Indcrd;i++)
    {
       if (std::fgets(line, BUFSIZ, in_file) == NULL) {
         std::fprintf(stderr,"Error: Failed to read from file.\n");
         return 0;
       }
       if ( std::sscanf(line,"%*s") < 0 )
         IOHBTerminate("Trilinos_Util_iohb.cpp: Null (or blank) line in index data region of HB file.\n");
       col =  0;
       for (ind = 0;ind<Indperline;ind++)
       {
          if (count == Nnzero) break;
          std::strncpy(ThisElement,line+col,Indwidth);
/*        ThisElement = substr(line,col,Indwidth); */
          rowind[count] = std::atoi(ThisElement)-offset;
          count++; col += Indwidth;
       }
    }
    free(ThisElement);

/*  Read array of values:  */

    if ( Type[0] != 'P' ) {          /* Skip if pattern only  */

       if ( Type[0] == 'C' ) Nentries = 2*Nnzero;
           else Nentries = Nnzero;

    ThisElement = (char *) malloc(Valwidth+2);
    if ( ThisElement == NULL ) IOHBTerminate("Insufficient memory for ThisElement.");
    *(ThisElement+Valwidth) = '\0';
    *(ThisElement+Valwidth+1) = '\0';
    count = 0;
    for (i=0;i<Valcrd;i++)
    {
       if (std::fgets(line, BUFSIZ, in_file) == NULL) {
         std::fprintf(stderr,"Error: Failed to read from file.\n");
         return 0;
       }
       if ( std::sscanf(line,"%*s") < 0 )
         IOHBTerminate("Trilinos_Util_iohb.cpp: Null (or blank) line in value data region of HB file.\n");
       if (Valflag == 'D')  {
          while( std::strchr(line,'D') ) *std::strchr(line,'D') = 'E';
/*           *std::strchr(Valfmt,'D') = 'E'; */
       }
       col =  0;
       for (ind = 0;ind<Valperline;ind++)
       {
          if (count == Nentries) break;
          std::strncpy(ThisElement,line+col,Valwidth);
          /*ThisElement = substr(line,col,Valwidth);*/
          if ( Valflag != 'F' && std::strchr(ThisElement,'E') == NULL ) {
             /* insert a char prefix for exp */
             last = std::strlen(ThisElement);
             for (j=last+1;j>=0;j--) {
                ThisElement[j] = ThisElement[j-1];
                if ( ThisElement[j] == '+' || ThisElement[j] == '-' ) {
                   ThisElement[j-1] = Valflag;
                   break;
                }
             }
          }
          val[count] = std::atof(ThisElement);
          count++; col += Valwidth;
          *(ThisElement+Valwidth) = '\0';
          *(ThisElement+Valwidth+1) = '\0';
       }
    }
    free(ThisElement);
    }

    std::fclose(in_file);
    return 1;
}

int readHB_newmat_double(const char* filename, int* M, int* N, int* nonzeros,
                         int** colptr, int** rowind, double** val)
{
  int Nrhs;
  char *Type;

  if (readHB_info(filename, M, N, nonzeros, &Type, &Nrhs) == 0) {
    return 0;
  }

  *colptr = (int *)malloc((*N+1)*sizeof(int));
  if ( *colptr == NULL ) IOHBTerminate("Insufficient memory for colptr.\n");
  *rowind = (int *)malloc(*nonzeros*sizeof(int));
  if ( *rowind == NULL ) IOHBTerminate("Insufficient memory for rowind.\n");
  if ( Type[0] == 'C' ) {
/*
   std::fprintf(stderr, "Warning: Reading complex data from HB file %s.\n",filename);
   std::fprintf(stderr, "         Real and imaginary parts will be interlaced in val[].\n");
*/
           /* Malloc enough space for real AND imaginary parts of val[] */
  *val = (double *)malloc(*nonzeros*sizeof(double)*2);
  if ( *val == NULL ) IOHBTerminate("Insufficient memory for val.\n");
  } else {
    if ( Type[0] != 'P' ) {
      /* Malloc enough space for real array val[] */
      *val = (double *)malloc(*nonzeros*sizeof(double));
      if ( *val == NULL ) IOHBTerminate("Insufficient memory for val.\n");
    }
  }  /* No val[] space needed if pattern only */
  return readHB_mat_double(filename, *colptr, *rowind, *val);

}

int readHB_aux_double(const char* filename, const char AuxType, double b[])
{
/****************************************************************************/
/*  This function opens and reads the specified file, placing auxillary     */
/*  vector(s) of the given type (if available) in b.                        */
/*  Return value is the number of vectors successfully read.                */
/*                                                                          */
/*                AuxType = 'F'   full right-hand-side vector(s)            */
/*                AuxType = 'G'   initial Guess vector(s)                   */
/*                AuxType = 'X'   eXact solution vector(s)                  */
/*                                                                          */
/*    ----------                                                            */
/*    **CAVEAT**                                                            */
/*    ----------                                                            */
/*  Parsing real formats from Fortran is tricky, and this file reader       */
/*  does not claim to be foolproof.   It has been tested for cases when     */
/*  the real values are printed consistently and evenly spaced on each      */
/*  line, with Fixed (F), and Exponential (E or D) formats.                 */
/*                                                                          */
/*  **  If the input file does not adhere to the H/B format, the  **        */
/*  **             results will be unpredictable.                 **        */
/*                                                                          */
/****************************************************************************/
    std::FILE *in_file;
    int i,j,n,maxcol,start,stride,col,last,linel;
    int Ptrcrd, Indcrd, Valcrd, Rhscrd;
    int Nrow, Ncol, Nnzero, Nentries;
    int Nrhs, nvecs, rhsi;
    int Rhsperline, Rhswidth, Rhsprec;
    int Rhsflag;
    char *ThisElement;
    char Title[73], Key[9], Type[4] = "XXX", Rhstype[4];
    char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];
    char line[BUFSIZ];

    if ((in_file = std::fopen( filename, "r")) == NULL) {
      std::fprintf(stderr,"Error: Cannot open file: %s\n",filename);
      return 0;
     }

    readHB_header(in_file, Title, Key, Type, &Nrow, &Ncol, &Nnzero, &Nrhs,
                  Ptrfmt, Indfmt, Valfmt, Rhsfmt,
                  &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);

    if (Nrhs <= 0)
    {
      std::fprintf(stderr, "Warn: Attempt to read auxillary vector(s) when none are present.\n");
      return 0;
    }
    if (Rhstype[0] != 'F' )
    {
      std::fprintf(stderr,"Warn: Attempt to read auxillary vector(s) which are not stored in Full form.\n");
      std::fprintf(stderr,"       Rhs must be specified as full. \n");
      return 0;
    }

/* If reading complex data, allow for interleaved real and imaginary values. */
    if ( Type[0] == 'C' ) {
       Nentries = 2*Nrow;
     } else {
       Nentries = Nrow;
    }

    nvecs = 1;

    if ( Rhstype[1] == 'G' ) nvecs++;
    if ( Rhstype[2] == 'X' ) nvecs++;

    if ( AuxType == 'G' && Rhstype[1] != 'G' ) {
      std::fprintf(stderr, "Warn: Attempt to read auxillary Guess vector(s) when none are present.\n");
      return 0;
    }
    if ( AuxType == 'X' && Rhstype[2] != 'X' ) {
      std::fprintf(stderr, "Warn: Attempt to read auxillary eXact solution vector(s) when none are present.\n");
      return 0;
    }

    ParseRfmt(Rhsfmt, &Rhsperline, &Rhswidth, &Rhsprec,&Rhsflag);
    maxcol = Rhsperline*Rhswidth;

/*  Lines to skip before starting to read RHS values... */
    n = Ptrcrd + Indcrd + Valcrd;

    for (i = 0; i < n; i++) {
      if (std::fgets(line, BUFSIZ, in_file) == NULL) {
        std::fprintf(stderr,"Error: Failed to read from file.\n");
        return 0;
      }
    }

/*  start  - number of initial aux vector entries to skip   */
/*           to reach first  vector requested               */
/*  stride - number of aux vector entries to skip between   */
/*           requested vectors                              */
    if ( AuxType == 'F' ) start = 0;
    else if ( AuxType == 'G' ) start = Nentries;
    else start = (nvecs-1)*Nentries;
    stride = (nvecs-1)*Nentries;

    if (std::fgets(line, BUFSIZ, in_file) == NULL) {
      std::fprintf(stderr,"Error: Failed to read from file.\n");
      return 0;
    }
    linel= std::strchr(line,'\n')-line;
    col = 0;
/*  Skip to initial offset */

    for (i=0;i<start;i++) {
       if ( col >=  ( maxcol<linel?maxcol:linel ) ) {
           if (std::fgets(line, BUFSIZ, in_file) == NULL) {
             std::fprintf(stderr,"Error: Failed to read from file.\n");
             return 0;
           }
           linel= std::strchr(line,'\n')-line;
           col = 0;
       }
       col += Rhswidth;
    }
    if (Rhsflag == 'D')  {
       while( std::strchr(line,'D') ) *std::strchr(line,'D') = 'E';
    }

/*  Read a vector of desired type, then skip to next */
/*  repeating to fill Nrhs vectors                   */

  ThisElement = (char *) malloc(Rhswidth+1);
  if ( ThisElement == NULL ) IOHBTerminate("Insufficient memory for ThisElement.");
  *(ThisElement+Rhswidth) = '\0';
  for (rhsi=0;rhsi<Nrhs;rhsi++) {

    for (i=0;i<Nentries;i++) {
       if ( col >= ( maxcol<linel?maxcol:linel ) ) {
           if (std::fgets(line, BUFSIZ, in_file) == NULL) {
             std::fprintf(stderr,"Error: Failed to read from file.\n");
             return 0;
           }
           linel= std::strchr(line,'\n')-line;
           if (Rhsflag == 'D')  {
              while( std::strchr(line,'D') ) *std::strchr(line,'D') = 'E';
           }
           col = 0;
       }
       std::strncpy(ThisElement,line+col,Rhswidth);
       /*ThisElement = substr(line, col, Rhswidth);*/
          if ( Rhsflag != 'F' && std::strchr(ThisElement,'E') == NULL ) {
             /* insert a char prefix for exp */
             last = std::strlen(ThisElement);
             for (j=last+1;j>=0;j--) {
                ThisElement[j] = ThisElement[j-1];
                if ( ThisElement[j] == '+' || ThisElement[j] == '-' ) {
                   ThisElement[j-1] = Rhsflag;
                   break;
                }
             }
          }
       b[i] = std::atof(ThisElement);
       col += Rhswidth;
    }

/*  Skip any interleaved Guess/eXact vectors */

    for (i=0;i<stride;i++) {
       if ( col >= ( maxcol<linel?maxcol:linel ) ) {
           if (std::fgets(line, BUFSIZ, in_file) == NULL) {
             std::fprintf(stderr,"Error: Failed to read from file.\n");
             return 0;
           }
           linel= std::strchr(line,'\n')-line;
           col = 0;
       }
       col += Rhswidth;
    }

  }
  free(ThisElement);


    std::fclose(in_file);
    return Nrhs;
}

int readHB_newaux_double(const char* filename, const char AuxType, double** b)
{
        int Nrhs = 0;
        int M = 0;
        int N = 0;
        int nonzeros = 0;
        char *Type = NULL;

        readHB_info(filename, &M, &N, &nonzeros, &Type, &Nrhs);
        if ( Nrhs <= 0 ) {
          std::fprintf(stderr,"Warn: Requested read of aux vector(s) when none are present.\n");
          return 0;
        } else {
          if ( Type[0] == 'C' ) {
            std::fprintf(stderr, "Warning: Reading complex aux vector(s) from HB file %s.",filename);
            std::fprintf(stderr, "         Real and imaginary parts will be interlaced in b[].");
            *b = (double *)malloc(M*Nrhs*sizeof(double)*2);
            if ( *b == NULL ) IOHBTerminate("Insufficient memory for rhs.\n");
            return readHB_aux_double(filename, AuxType, *b);
          } else {
            *b = (double *)malloc(M*Nrhs*sizeof(double));
            if ( *b == NULL ) IOHBTerminate("Insufficient memory for rhs.\n");
            return readHB_aux_double(filename, AuxType, *b);
          }
        }
}

int writeHB_mat_double(const char* filename, int M, int N,
                        int nz, const int colptr[], const int rowind[],
                        const double val[], int Nrhs, const double rhs[],
                        const double guess[], const double exact[],
                        const char* Title, const char* Key, const char* Type,
                        char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
                        const char* Rhstype)
{
/****************************************************************************/
/*  The writeHB function opens the named file and writes the specified      */
/*  matrix and optional right-hand-side(s) to that file in Harwell-Boeing   */
/*  format.                                                                 */
/*                                                                          */
/*  For a description of the Harwell Boeing standard, see:                  */
/*            Duff, et al.,  ACM TOMS Vol.15, No.1, March 1989              */
/*                                                                          */
/****************************************************************************/
    std::FILE *out_file;
    int i,j,entry,offset,acount,linemod;
    int totcrd, ptrcrd, indcrd, valcrd, rhscrd;
    int nvalentries, nrhsentries;
    int Ptrperline, Ptrwidth, Indperline, Indwidth;
    int Rhsperline, Rhswidth, Rhsprec;
    int Rhsflag;
    int Valperline, Valwidth, Valprec;
    int Valflag;           /* Indicates 'E','D', or 'F' float format */
    char pformat[16],iformat[16],vformat[19],rformat[19];

    if ( Type[0] == 'C' ) {
         nvalentries = 2*nz;
         nrhsentries = 2*M;
    } else {
         nvalentries = nz;
         nrhsentries = M;
    }

    if ( filename != NULL ) {
       if ( (out_file = std::fopen( filename, "w")) == NULL ) {
         std::fprintf(stderr,"Error: Cannot open file: %s\n",filename);
         return 0;
       }
    } else out_file = stdout;

    if ( Ptrfmt == NULL ) strcpy(Ptrfmt, "(8I10)");
    ParseIfmt(Ptrfmt,&Ptrperline,&Ptrwidth);
    std::sprintf(pformat,"%%%dd",Ptrwidth);
    ptrcrd = (N+1)/Ptrperline;
    if ( (N+1)%Ptrperline != 0) ptrcrd++;

    if ( Indfmt == NULL ) Indfmt =  Ptrfmt;
    ParseIfmt(Indfmt,&Indperline,&Indwidth);
    std::sprintf(iformat,"%%%dd",Indwidth);
    indcrd = nz/Indperline;
    if ( nz%Indperline != 0) indcrd++;

    if ( Type[0] != 'P' ) {          /* Skip if pattern only  */
      if ( Valfmt == NULL ) strcpy(Valfmt, "(4E20.13)");
      ParseRfmt(Valfmt,&Valperline,&Valwidth,&Valprec,&Valflag);
      if (Valflag == 'D') *std::strchr(Valfmt,'D') = 'E';
      if (Valflag == 'F')
         std::sprintf(vformat,"%% %d.%df",Valwidth,Valprec);
      else
         std::sprintf(vformat,"%% %d.%dE",Valwidth,Valprec);
      valcrd = nvalentries/Valperline;
      if ( nvalentries%Valperline != 0) valcrd++;
    } else valcrd = 0;

    if ( Nrhs > 0 ) {
       if ( Rhsfmt == NULL ) Rhsfmt = Valfmt;
       ParseRfmt(Rhsfmt,&Rhsperline,&Rhswidth,&Rhsprec, &Rhsflag);
       if (Rhsflag == 'F')
          std::sprintf(rformat,"%% %d.%df",Rhswidth,Rhsprec);
       else
          std::sprintf(rformat,"%% %d.%dE",Rhswidth,Rhsprec);
       if (Rhsflag == 'D') *std::strchr(Rhsfmt,'D') = 'E';
       rhscrd = nrhsentries/Rhsperline;
       if ( nrhsentries%Rhsperline != 0) rhscrd++;
       if ( Rhstype[1] == 'G' ) rhscrd+=rhscrd;
       if ( Rhstype[2] == 'X' ) rhscrd+=rhscrd;
       rhscrd*=Nrhs;
    } else rhscrd = 0;

    totcrd = 4+ptrcrd+indcrd+valcrd+rhscrd;


/*  Print header information:  */

    std::fprintf(out_file,"%-72s%-8s\n%14d%14d%14d%14d%14d\n",Title, Key, totcrd,
            ptrcrd, indcrd, valcrd, rhscrd);
    std::fprintf(out_file,"%3s%11s%14d%14d%14d\n",Type,"          ", M, N, nz);
    std::fprintf(out_file,"%-16s%-16s%-20s", Ptrfmt, Indfmt, Valfmt);
    if ( Nrhs != 0 ) {
/*     Print Rhsfmt on fourth line and                                    */
/*           optional fifth header line for auxillary vector information: */
       std::fprintf(out_file,"%-20s\n%-14s%d\n",Rhsfmt,Rhstype,Nrhs);
    } else std::fprintf(out_file,"\n");

    offset = 1-_SP_base;  /* if base 0 storage is declared (via macro definition), */
                          /* then storage entries are offset by 1                  */

/*  Print column pointers:   */
    for (i=0;i<N+1;i++)
    {
       entry = colptr[i]+offset;
       std::fprintf(out_file,pformat,entry);
       if ( (i+1)%Ptrperline == 0 ) std::fprintf(out_file,"\n");
    }

   if ( (N+1) % Ptrperline != 0 ) std::fprintf(out_file,"\n");

/*  Print row indices:       */
    for (i=0;i<nz;i++)
    {
       entry = rowind[i]+offset;
       std::fprintf(out_file,iformat,entry);
       if ( (i+1)%Indperline == 0 ) std::fprintf(out_file,"\n");
    }

   if ( nz % Indperline != 0 ) std::fprintf(out_file,"\n");

/*  Print values:            */

    if ( Type[0] != 'P' ) {          /* Skip if pattern only  */

    for (i=0;i<nvalentries;i++)
    {
       std::fprintf(out_file,vformat,val[i]);
       if ( (i+1)%Valperline == 0 ) std::fprintf(out_file,"\n");
    }

    if ( nvalentries % Valperline != 0 ) std::fprintf(out_file,"\n");

/*  If available,  print right hand sides,
           guess vectors and exact solution vectors:  */
    acount = 1;
    linemod = 0;
    if ( Nrhs > 0 ) {
       for (i=0;i<Nrhs;i++)
       {
          for ( j=0;j<nrhsentries;j++ ) {
            std::fprintf(out_file,rformat,rhs[j]);
            if ( acount++%Rhsperline == linemod ) std::fprintf(out_file,"\n");
          }
          if ( (acount-1)%Rhsperline != linemod ) {
            std::fprintf(out_file,"\n");
            linemod = (acount-1)%Rhsperline;
          }
          rhs += nrhsentries;
          if ( Rhstype[1] == 'G' ) {
            for ( j=0;j<nrhsentries;j++ ) {
              std::fprintf(out_file,rformat,guess[j]);
              if ( acount++%Rhsperline == linemod ) std::fprintf(out_file,"\n");
            }
            if ( (acount-1)%Rhsperline != linemod ) {
              std::fprintf(out_file,"\n");
              linemod = (acount-1)%Rhsperline;
            }
            guess += nrhsentries;
          }
          if ( Rhstype[2] == 'X' ) {
            for ( j=0;j<nrhsentries;j++ ) {
              std::fprintf(out_file,rformat,exact[j]);
              if ( acount++%Rhsperline == linemod ) std::fprintf(out_file,"\n");
            }
            if ( (acount-1)%Rhsperline != linemod ) {
              std::fprintf(out_file,"\n");
              linemod = (acount-1)%Rhsperline;
            }
            exact += nrhsentries;
          }
       }
    }

    }

    if ( std::fclose(out_file) != 0){
      std::fprintf(stderr,"Error closing file in writeHB_mat_double().\n");
      return 0;
    } else return 1;

}

int readHB_mat_char(const char* filename, int colptr[], int rowind[],
                                           char val[], char* Valfmt)
{
/****************************************************************************/
/*  This function opens and reads the specified file, interpreting its      */
/*  contents as a sparse matrix stored in the Harwell/Boeing standard       */
/*  format and creating compressed column storage scheme vectors to hold    */
/*  the index and nonzero value information.                                */
/*                                                                          */
/*    ----------                                                            */
/*    **CAVEAT**                                                            */
/*    ----------                                                            */
/*  Parsing real formats from Fortran is tricky, and this file reader       */
/*  does not claim to be foolproof.   It has been tested for cases when     */
/*  the real values are printed consistently and evenly spaced on each      */
/*  line, with Fixed (F), and Exponential (E or D) formats.                 */
/*                                                                          */
/*  **  If the input file does not adhere to the H/B format, the  **        */
/*  **             results will be unpredictable.                 **        */
/*                                                                          */
/****************************************************************************/
    std::FILE *in_file;
    int i,j,ind,col,offset,count,last;
    int Nrow,Ncol,Nnzero,Nentries,Nrhs;
    int Ptrcrd, Indcrd, Valcrd, Rhscrd;
    int Ptrperline, Ptrwidth, Indperline, Indwidth;
    int Valperline, Valwidth, Valprec;
    int Valflag;           /* Indicates 'E','D', or 'F' float format */
    char* ThisElement;
    char line[BUFSIZ];
    char Title[73], Key[9], Type[4] = "XXX", Rhstype[4];
    char Ptrfmt[17], Indfmt[17], Rhsfmt[21];

    if ( (in_file = std::fopen( filename, "r")) == NULL ) {
       std::fprintf(stderr,"Error: Cannot open file: %s\n",filename);
       return 0;
    }

    readHB_header(in_file, Title, Key, Type, &Nrow, &Ncol, &Nnzero, &Nrhs,
                  Ptrfmt, Indfmt, Valfmt, Rhsfmt,
                  &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);

/*  Parse the array input formats from Line 3 of HB file  */
    ParseIfmt(Ptrfmt,&Ptrperline,&Ptrwidth);
    ParseIfmt(Indfmt,&Indperline,&Indwidth);
    if ( Type[0] != 'P' ) {          /* Skip if pattern only  */
       ParseRfmt(Valfmt,&Valperline,&Valwidth,&Valprec,&Valflag);
       if (Valflag == 'D') {
          *std::strchr(Valfmt,'D') = 'E';
       }
    }

/*  Read column pointer array:   */

    offset = 1-_SP_base;  /* if base 0 storage is declared (via macro definition), */
                          /* then storage entries are offset by 1                  */

    ThisElement = (char *) malloc(Ptrwidth+1);
    if ( ThisElement == NULL ) IOHBTerminate("Insufficient memory for ThisElement.");
    *(ThisElement+Ptrwidth) = '\0';
    count=0;
    for (i=0;i<Ptrcrd;i++)
    {
       if (std::fgets(line, BUFSIZ, in_file) == NULL) {
         std::fprintf(stderr,"Error: Failed to read from file.\n");
         return 0;
       }
       if ( std::sscanf(line,"%*s") < 0 )
         IOHBTerminate("Trilinos_Util_iohb.cpp: Null (or blank) line in pointer data region of HB file.\n");
       col =  0;
       for (ind = 0;ind<Ptrperline;ind++)
       {
          if (count > Ncol) break;
          std::strncpy(ThisElement,line+col,Ptrwidth);
          /*ThisElement = substr(line,col,Ptrwidth);*/
          colptr[count] = std::atoi(ThisElement)-offset;
          count++; col += Ptrwidth;
       }
    }
    free(ThisElement);

/*  Read row index array:  */

    ThisElement = (char *) malloc(Indwidth+1);
    if ( ThisElement == NULL ) IOHBTerminate("Insufficient memory for ThisElement.");
    *(ThisElement+Indwidth) = '\0';
    count = 0;
    for (i=0;i<Indcrd;i++)
    {
       if (std::fgets(line, BUFSIZ, in_file) == NULL) {
         std::fprintf(stderr,"Error: Failed to read from file.\n");
         return 0;
       }
       if ( std::sscanf(line,"%*s") < 0 )
         IOHBTerminate("Trilinos_Util_iohb.cpp: Null (or blank) line in index data region of HB file.\n");
       col =  0;
       for (ind = 0;ind<Indperline;ind++)
       {
          if (count == Nnzero) break;
          std::strncpy(ThisElement,line+col,Indwidth);
          /*ThisElement = substr(line,col,Indwidth);*/
          rowind[count] = std::atoi(ThisElement)-offset;
          count++; col += Indwidth;
       }
    }
    free(ThisElement);

/*  Read array of values:  AS CHARACTERS*/

    if ( Type[0] != 'P' ) {          /* Skip if pattern only  */

       if ( Type[0] == 'C' ) Nentries = 2*Nnzero;
           else Nentries = Nnzero;

    ThisElement = (char *) malloc(Valwidth+1);
    if ( ThisElement == NULL ) IOHBTerminate("Insufficient memory for ThisElement.");
    *(ThisElement+Valwidth) = '\0';
    count = 0;
    for (i=0;i<Valcrd;i++)
    {
       if (std::fgets(line, BUFSIZ, in_file) == NULL) {
         std::fprintf(stderr,"Error: Failed to read from file.\n");
         return 0;
       }
       if ( std::sscanf(line,"%*s") < 0 )
         IOHBTerminate("Trilinos_Util_iohb.cpp: Null (or blank) line in value data region of HB file.\n");
       if (Valflag == 'D') {
          while( std::strchr(line,'D') ) *std::strchr(line,'D') = 'E';
       }
       col =  0;
       for (ind = 0;ind<Valperline;ind++)
       {
          if (count == Nentries) break;
          ThisElement = &val[count*Valwidth];
          std::strncpy(ThisElement,line+col,Valwidth);
          /*std::strncpy(ThisElement,substr(line,col,Valwidth),Valwidth);*/
          if ( Valflag != 'F' && std::strchr(ThisElement,'E') == NULL ) {
             /* insert a char prefix for exp */
             last = std::strlen(ThisElement);
             for (j=last+1;j>=0;j--) {
                ThisElement[j] = ThisElement[j-1];
                if ( ThisElement[j] == '+' || ThisElement[j] == '-' ) {
                   ThisElement[j-1] = Valflag;
                   break;
                }
             }
          }
          count++; col += Valwidth;
       }
    }
    }

    return 1;
}

int readHB_newmat_char(const char* filename, int* M, int* N, int* nonzeros, int** colptr,
                          int** rowind, char** val, char** Valfmt)
{
    std::FILE *in_file;
    int Nrhs;
    int Ptrcrd, Indcrd, Valcrd, Rhscrd;
    int Valperline, Valwidth, Valprec;
    int Valflag;           /* Indicates 'E','D', or 'F' float format */
    char Title[73], Key[9], Type[4] = "XXX", Rhstype[4];
    char Ptrfmt[17], Indfmt[17], Rhsfmt[21];

    if ((in_file = std::fopen( filename, "r")) == NULL) {
      std::fprintf(stderr,"Error: Cannot open file: %s\n",filename);
      return 0;
     }

    *Valfmt = (char *)malloc(21*sizeof(char));
    if ( *Valfmt == NULL ) IOHBTerminate("Insufficient memory for Valfmt.");
    readHB_header(in_file, Title, Key, Type, M, N, nonzeros, &Nrhs,
                  Ptrfmt, Indfmt, (*Valfmt), Rhsfmt,
                  &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);
    std::fclose(in_file);
    ParseRfmt(*Valfmt,&Valperline,&Valwidth,&Valprec,&Valflag);

        *colptr = (int *)malloc((*N+1)*sizeof(int));
        if ( *colptr == NULL ) IOHBTerminate("Insufficient memory for colptr.\n");
        *rowind = (int *)malloc(*nonzeros*sizeof(int));
        if ( *rowind == NULL ) IOHBTerminate("Insufficient memory for rowind.\n");
        if ( Type[0] == 'C' ) {
/*
   std::fprintf(stderr, "Warning: Reading complex data from HB file %s.\n",filename);
   std::fprintf(stderr, "         Real and imaginary parts will be interlaced in val[].\n");
*/
           /* Malloc enough space for real AND imaginary parts of val[] */
           *val = (char *)malloc(*nonzeros*Valwidth*sizeof(char)*2);
           if ( *val == NULL ) IOHBTerminate("Insufficient memory for val.\n");
        } else {
           if ( Type[0] != 'P' ) {
             /* Malloc enough space for real array val[] */
             *val = (char *)malloc(*nonzeros*Valwidth*sizeof(char));
             if ( *val == NULL ) IOHBTerminate("Insufficient memory for val.\n");
           }
        }  /* No val[] space needed if pattern only */
        return readHB_mat_char(filename, *colptr, *rowind, *val, *Valfmt);

}

int readHB_aux_char(const char* filename, const char AuxType, char b[])
{
/****************************************************************************/
/*  This function opens and reads the specified file, placing auxilary      */
/*  vector(s) of the given type (if available) in b :                       */
/*  Return value is the number of vectors successfully read.                */
/*                                                                          */
/*                AuxType = 'F'   full right-hand-side vector(s)            */
/*                AuxType = 'G'   initial Guess vector(s)                   */
/*                AuxType = 'X'   eXact solution vector(s)                  */
/*                                                                          */
/*    ----------                                                            */
/*    **CAVEAT**                                                            */
/*    ----------                                                            */
/*  Parsing real formats from Fortran is tricky, and this file reader       */
/*  does not claim to be foolproof.   It has been tested for cases when     */
/*  the real values are printed consistently and evenly spaced on each      */
/*  line, with Fixed (F), and Exponential (E or D) formats.                 */
/*                                                                          */
/*  **  If the input file does not adhere to the H/B format, the  **        */
/*  **             results will be unpredictable.                 **        */
/*                                                                          */
/****************************************************************************/
    std::FILE *in_file;
    int i,j,n,maxcol,start,stride,col,last,linel,nvecs,rhsi;
    int Nrow, Ncol, Nnzero, Nentries,Nrhs;
    int Ptrcrd, Indcrd, Valcrd, Rhscrd;
    int Rhsperline, Rhswidth, Rhsprec;
    int Rhsflag;
    char Title[73], Key[9], Type[4] = "XXX", Rhstype[4];
    char Ptrfmt[17], Indfmt[17], Valfmt[21], Rhsfmt[21];
    char line[BUFSIZ];
    char *ThisElement;

    if ((in_file = std::fopen( filename, "r")) == NULL) {
      std::fprintf(stderr,"Error: Cannot open file: %s\n",filename);
      return 0;
     }

    readHB_header(in_file, Title, Key, Type, &Nrow, &Ncol, &Nnzero, &Nrhs,
                  Ptrfmt, Indfmt, Valfmt, Rhsfmt,
                  &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);

    if (Nrhs <= 0)
    {
      std::fprintf(stderr, "Warn: Attempt to read auxillary vector(s) when none are present.\n");
      return 0;
    }
    if (Rhstype[0] != 'F' )
    {
      std::fprintf(stderr,"Warn: Attempt to read auxillary vector(s) which are not stored in Full form.\n");
      std::fprintf(stderr,"       Rhs must be specified as full. \n");
      return 0;
    }

/* If reading complex data, allow for interleaved real and imaginary values. */
    if ( Type[0] == 'C' ) {
       Nentries = 2*Nrow;
     } else {
       Nentries = Nrow;
    }

    nvecs = 1;

    if ( Rhstype[1] == 'G' ) nvecs++;
    if ( Rhstype[2] == 'X' ) nvecs++;

    if ( AuxType == 'G' && Rhstype[1] != 'G' ) {
      std::fprintf(stderr, "Warn: Attempt to read auxillary Guess vector(s) when none are present.\n");
      return 0;
    }
    if ( AuxType == 'X' && Rhstype[2] != 'X' ) {
      std::fprintf(stderr, "Warn: Attempt to read auxillary eXact solution vector(s) when none are present.\n");
      return 0;
    }

    ParseRfmt(Rhsfmt, &Rhsperline, &Rhswidth, &Rhsprec,&Rhsflag);
    maxcol = Rhsperline*Rhswidth;

/*  Lines to skip before starting to read RHS values... */
    n = Ptrcrd + Indcrd + Valcrd;

    for (i = 0; i < n; i++) {
      if (std::fgets(line, BUFSIZ, in_file) == NULL) {
        std::fprintf(stderr,"Error: Failed to read from file.\n");
        return 0;
      }
    }

/*  start  - number of initial aux vector entries to skip   */
/*           to reach first  vector requested               */
/*  stride - number of aux vector entries to skip between   */
/*           requested vectors                              */
    if ( AuxType == 'F' ) start = 0;
    else if ( AuxType == 'G' ) start = Nentries;
    else start = (nvecs-1)*Nentries;
    stride = (nvecs-1)*Nentries;

    if (std::fgets(line, BUFSIZ, in_file) == NULL) {
      std::fprintf(stderr,"Error: Failed to read from file.\n");
      return 0;
    }
    linel= std::strchr(line,'\n')-line;
    if ( std::sscanf(line,"%*s") < 0 )
       IOHBTerminate("Trilinos_Util_iohb.cpp: Null (or blank) line in auxillary vector data region of HB file.\n");
    col = 0;
/*  Skip to initial offset */

    for (i=0;i<start;i++) {
       col += Rhswidth;
       if ( col >= ( maxcol<linel?maxcol:linel ) ) {
           if (std::fgets(line, BUFSIZ, in_file) == NULL) {
             std::fprintf(stderr,"Error: Failed to read from file.\n");
             return 0;
           }
           linel= std::strchr(line,'\n')-line;
       if ( std::sscanf(line,"%*s") < 0 )
       IOHBTerminate("Trilinos_Util_iohb.cpp: Null (or blank) line in auxillary vector data region of HB file.\n");
           col = 0;
       }
    }

    if (Rhsflag == 'D')  {
      while( std::strchr(line,'D') ) *std::strchr(line,'D') = 'E';
    }
/*  Read a vector of desired type, then skip to next */
/*  repeating to fill Nrhs vectors                   */

  for (rhsi=0;rhsi<Nrhs;rhsi++) {

    for (i=0;i<Nentries;i++) {
       if ( col >= ( maxcol<linel?maxcol:linel ) ) {
           if (std::fgets(line, BUFSIZ, in_file) == NULL) {
             std::fprintf(stderr,"Error: Failed to read from file.\n");
             return 0;
           }
           linel= std::strchr(line,'\n')-line;
       if ( std::sscanf(line,"%*s") < 0 )
       IOHBTerminate("Trilinos_Util_iohb.cpp: Null (or blank) line in auxillary vector data region of HB file.\n");
           if (Rhsflag == 'D')  {
              while( std::strchr(line,'D') ) *std::strchr(line,'D') = 'E';
           }
           col = 0;
       }
       ThisElement = &b[i*Rhswidth];
       std::strncpy(ThisElement,line+col,Rhswidth);
          if ( Rhsflag != 'F' && std::strchr(ThisElement,'E') == NULL ) {
             /* insert a char prefix for exp */
             last = std::strlen(ThisElement);
             for (j=last+1;j>=0;j--) {
                ThisElement[j] = ThisElement[j-1];
                if ( ThisElement[j] == '+' || ThisElement[j] == '-' ) {
                   ThisElement[j-1] = Rhsflag;
                   break;
                }
             }
          }
       col += Rhswidth;
    }
    b+=Nentries*Rhswidth;

/*  Skip any interleaved Guess/eXact vectors */

    for (i=0;i<stride;i++) {
       col += Rhswidth;
       if ( col >= ( maxcol<linel?maxcol:linel ) ) {
           if (std::fgets(line, BUFSIZ, in_file) == NULL) {
             std::fprintf(stderr,"Error: Failed to read from file.\n");
             return 0;
           }
           linel= std::strchr(line,'\n')-line;
       if ( std::sscanf(line,"%*s") < 0 )
       IOHBTerminate("Trilinos_Util_iohb.cpp: Null (or blank) line in auxillary vector data region of HB file.\n");
           col = 0;
       }
    }

  }


    std::fclose(in_file);
    return Nrhs;
}

int readHB_newaux_char(const char* filename, const char AuxType, char** b, char** Rhsfmt)
{
    std::FILE *in_file;
    int Ptrcrd, Indcrd, Valcrd, Rhscrd;
    int Nrow,Ncol,Nnzero,Nrhs;
    int Rhsperline, Rhswidth, Rhsprec;
    int Rhsflag;
    char Title[73], Key[9], Type[4] = "XXX", Rhstype[4];
    char Ptrfmt[17], Indfmt[17], Valfmt[21];

    if ((in_file = std::fopen( filename, "r")) == NULL) {
      std::fprintf(stderr,"Error: Cannot open file: %s\n",filename);
      return 0;
     }

    *Rhsfmt = (char *)malloc(21*sizeof(char));
    if ( *Rhsfmt == NULL ) IOHBTerminate("Insufficient memory for Rhsfmt.");
    readHB_header(in_file, Title, Key, Type, &Nrow, &Ncol, &Nnzero, &Nrhs,
                  Ptrfmt, Indfmt, Valfmt, (*Rhsfmt),
                  &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, Rhstype);
     std::fclose(in_file);
        if ( Nrhs == 0 ) {
          std::fprintf(stderr,"Warn: Requested read of aux vector(s) when none are present.\n");
          return 0;
        } else {
          ParseRfmt(*Rhsfmt,&Rhsperline,&Rhswidth,&Rhsprec,&Rhsflag);
          if ( Type[0] == 'C' ) {
            std::fprintf(stderr, "Warning: Reading complex aux vector(s) from HB file %s.",filename);
            std::fprintf(stderr, "         Real and imaginary parts will be interlaced in b[].");
            *b = (char *)malloc(Nrow*Nrhs*Rhswidth*sizeof(char)*2);
            if ( *b == NULL ) IOHBTerminate("Insufficient memory for rhs.\n");
            return readHB_aux_char(filename, AuxType, *b);
          } else {
            *b = (char *)malloc(Nrow*Nrhs*Rhswidth*sizeof(char));
            if ( *b == NULL ) IOHBTerminate("Insufficient memory for rhs.\n");
            return readHB_aux_char(filename, AuxType, *b);
          }
        }
}

int writeHB_mat_char(const char* filename, int M, int N,
                        int nz, const int colptr[], const int rowind[],
                        const char val[], int Nrhs, const char rhs[],
                        const char guess[], const char exact[],
                        const char* Title, const char* Key, const char* Type,
                        char* Ptrfmt, char* Indfmt, char* Valfmt, char* Rhsfmt,
                        const char* Rhstype)
{
/****************************************************************************/
/*  The writeHB function opens the named file and writes the specified      */
/*  matrix and optional right-hand-side(s) to that file in Harwell-Boeing   */
/*  format.                                                                 */
/*                                                                          */
/*  For a description of the Harwell Boeing standard, see:                  */
/*            Duff, et al.,  ACM TOMS Vol.15, No.1, March 1989              */
/*                                                                          */
/****************************************************************************/
    std::FILE *out_file;
    int i,j,acount,linemod,entry,offset;
    int totcrd, ptrcrd, indcrd, valcrd, rhscrd;
    int nvalentries, nrhsentries;
    int Ptrperline, Ptrwidth, Indperline, Indwidth;
    int Rhsperline, Rhswidth, Rhsprec;
    int Rhsflag;
    int Valperline, Valwidth, Valprec;
    int Valflag;           /* Indicates 'E','D', or 'F' float format */
    char pformat[16],iformat[16],vformat[19],rformat[19];

    if ( Type[0] == 'C' ) {
         nvalentries = 2*nz;
         nrhsentries = 2*M;
    } else {
         nvalentries = nz;
         nrhsentries = M;
    }

    if ( filename != NULL ) {
       if ( (out_file = std::fopen( filename, "w")) == NULL ) {
         std::fprintf(stderr,"Error: Cannot open file: %s\n",filename);
         return 0;
       }
    } else out_file = stdout;

    if ( Ptrfmt == NULL ) strcpy(Ptrfmt, "(8I10)");
    ParseIfmt(Ptrfmt,&Ptrperline,&Ptrwidth);
    std::sprintf(pformat,"%%%dd",Ptrwidth);

    if ( Indfmt == NULL ) Indfmt =  Ptrfmt;
    ParseIfmt(Indfmt,&Indperline,&Indwidth);
    std::sprintf(iformat,"%%%dd",Indwidth);

    if ( Type[0] != 'P' ) {          /* Skip if pattern only  */
      if ( Valfmt == NULL ) strcpy(Valfmt, "(4E20.13)");
      ParseRfmt(Valfmt,&Valperline,&Valwidth,&Valprec,&Valflag);
      std::sprintf(vformat,"%%%ds",Valwidth);
    }

    ptrcrd = (N+1)/Ptrperline;
    if ( (N+1)%Ptrperline != 0) ptrcrd++;

    indcrd = nz/Indperline;
    if ( nz%Indperline != 0) indcrd++;

    valcrd = nvalentries/Valperline;
    if ( nvalentries%Valperline != 0) valcrd++;

    if ( Nrhs > 0 ) {
       if ( Rhsfmt == NULL ) Rhsfmt = Valfmt;
       ParseRfmt(Rhsfmt,&Rhsperline,&Rhswidth,&Rhsprec, &Rhsflag);
       std::sprintf(rformat,"%%%ds",Rhswidth);
       rhscrd = nrhsentries/Rhsperline;
       if ( nrhsentries%Rhsperline != 0) rhscrd++;
       if ( Rhstype[1] == 'G' ) rhscrd+=rhscrd;
       if ( Rhstype[2] == 'X' ) rhscrd+=rhscrd;
       rhscrd*=Nrhs;
    } else rhscrd = 0;

    totcrd = 4+ptrcrd+indcrd+valcrd+rhscrd;


/*  Print header information:  */

    std::fprintf(out_file,"%-72s%-8s\n%14d%14d%14d%14d%14d\n",Title, Key, totcrd,
            ptrcrd, indcrd, valcrd, rhscrd);
    std::fprintf(out_file,"%3s%11s%14d%14d%14d\n",Type,"          ", M, N, nz);
    std::fprintf(out_file,"%-16s%-16s%-20s", Ptrfmt, Indfmt, Valfmt);
    if ( Nrhs != 0 ) {
/*     Print Rhsfmt on fourth line and                                    */
/*           optional fifth header line for auxillary vector information: */
       std::fprintf(out_file,"%-20s\n%-14s%d\n",Rhsfmt,Rhstype,Nrhs);
    } else std::fprintf(out_file,"\n");

    offset = 1-_SP_base;  /* if base 0 storage is declared (via macro definition), */
                          /* then storage entries are offset by 1                  */

/*  Print column pointers:   */
    for (i=0;i<N+1;i++)
    {
       entry = colptr[i]+offset;
       std::fprintf(out_file,pformat,entry);
       if ( (i+1)%Ptrperline == 0 ) std::fprintf(out_file,"\n");
    }

   if ( (N+1) % Ptrperline != 0 ) std::fprintf(out_file,"\n");

/*  Print row indices:       */
    for (i=0;i<nz;i++)
    {
       entry = rowind[i]+offset;
       std::fprintf(out_file,iformat,entry);
       if ( (i+1)%Indperline == 0 ) std::fprintf(out_file,"\n");
    }

   if ( nz % Indperline != 0 ) std::fprintf(out_file,"\n");

/*  Print values:            */

    if ( Type[0] != 'P' ) {          /* Skip if pattern only  */
    for (i=0;i<nvalentries;i++)
    {
       std::fprintf(out_file,vformat,val+i*Valwidth);
       if ( (i+1)%Valperline == 0 ) std::fprintf(out_file,"\n");
    }

    if ( nvalentries % Valperline != 0 ) std::fprintf(out_file,"\n");

/*  Print right hand sides:  */
    acount = 1;
    linemod=0;
    if ( Nrhs > 0 ) {
      for (j=0;j<Nrhs;j++) {
       for (i=0;i<nrhsentries;i++)
       {
          std::fprintf(out_file,rformat,rhs+i*Rhswidth);
          if ( acount++%Rhsperline == linemod ) std::fprintf(out_file,"\n");
       }
       if ( acount%Rhsperline != linemod ) {
          std::fprintf(out_file,"\n");
          linemod = (acount-1)%Rhsperline;
       }
       if ( Rhstype[1] == 'G' ) {
         for (i=0;i<nrhsentries;i++)
         {
           std::fprintf(out_file,rformat,guess+i*Rhswidth);
           if ( acount++%Rhsperline == linemod ) std::fprintf(out_file,"\n");
         }
         if ( acount%Rhsperline != linemod ) {
            std::fprintf(out_file,"\n");
            linemod = (acount-1)%Rhsperline;
         }
       }
       if ( Rhstype[2] == 'X' ) {
         for (i=0;i<nrhsentries;i++)
         {
           std::fprintf(out_file,rformat,exact+i*Rhswidth);
           if ( acount++%Rhsperline == linemod ) std::fprintf(out_file,"\n");
         }
         if ( acount%Rhsperline != linemod ) {
            std::fprintf(out_file,"\n");
            linemod = (acount-1)%Rhsperline;
         }
       }
      }
    }

    }

    if ( std::fclose(out_file) != 0){
      std::fprintf(stderr,"Error closing file in writeHB_mat_char().\n");
      return 0;
    } else return 1;

}

int ParseIfmt(char* fmt, int* perline, int* width)
{
/*************************************************/
/*  Parse an *integer* format field to determine */
/*  width and number of elements per line.       */
/*************************************************/
    char *tmp;
    if (fmt == NULL ) {
      *perline = 0; *width = 0; return 0;
    }
    upcase(fmt);
    tmp = std::strchr(fmt,'(');
    tmp = substr(fmt,tmp - fmt + 1, std::strchr(fmt,'I') - tmp - 1);
    *perline = std::atoi(tmp);
    if (*perline == 0 ) *perline = 1 ;
    if (tmp!=NULL) free ((void *) tmp);
    tmp = std::strchr(fmt,'I');
    tmp = substr(fmt,tmp - fmt + 1, std::strchr(fmt,')') - tmp - 1);
    *width = std::atoi(tmp);
    if (tmp!=NULL) free ((void *) tmp);
    return *width;
}

int ParseRfmt(char* fmt, int* perline, int* width, int* prec, int* flag)
{
/*************************************************/
/*  Parse a *real* format field to determine     */
/*  width and number of elements per line.       */
/*  Also sets flag indicating 'E' 'F' 'P' or 'D' */
/*  format.                                      */
/*************************************************/
    char* tmp;
    char* tmp1;
    char* tmp2;
    char* tmp3;
    int len;

    if (fmt == NULL ) {
      *perline = 0;
      *width = 0;
      flag = NULL;
      return 0;
    }

    upcase(fmt);
    if (std::strchr(fmt,'(') != NULL)  fmt = std::strchr(fmt,'(');
    if (std::strchr(fmt,')') != NULL)  {
       tmp2 = std::strchr(fmt,')');
       while ( std::strchr(tmp2+1,')') != NULL ) {
          tmp2 = std::strchr(tmp2+1,')');
       }
       *(tmp2+1) = '\0';
    }
    if (std::strchr(fmt,'P') != NULL)  /* Remove any scaling factor, which */
    {                             /* affects output only, not input */
      if (std::strchr(fmt,'(') != NULL) {
        tmp = std::strchr(fmt,'P');
        if ( *(++tmp) == ',' ) tmp++;
        tmp3 = std::strchr(fmt,'(')+1;
        len = tmp-tmp3;
        tmp2 = tmp3;
        while ( *(tmp2+len) != '\0' ) {
           *tmp2=*(tmp2+len);
           tmp2++;
        }
        *(std::strchr(fmt,')')+1) = '\0';
      }
    }
    if (std::strchr(fmt,'E') != NULL) {
       *flag = 'E';
    } else if (std::strchr(fmt,'D') != NULL) {
       *flag = 'D';
    } else if (std::strchr(fmt,'F') != NULL) {
       *flag = 'F';
    } else {
      std::fprintf(stderr,"Real format %s in H/B file not supported.\n",fmt);
      return 0;
    }
    tmp = std::strchr(fmt,'(');
    tmp = substr(fmt,tmp - fmt + 1, std::strchr(fmt,*flag) - tmp - 1);
    *perline = std::atoi(tmp);
    if (*perline == 0 ) *perline = 1 ;
    if (tmp!=NULL) free ((void *) tmp);
    tmp = std::strchr(fmt,*flag);
    if ( std::strchr(fmt,'.') ) {
      tmp1 = substr( fmt, std::strchr(fmt,'.') - fmt + 1, std::strchr(fmt,')') - std::strchr(fmt,'.')-1);
      *prec = std::atoi( tmp1 );
      if (tmp1!=NULL) free ((void *) tmp1);
      tmp1 = substr(fmt,tmp - fmt + 1, std::strchr(fmt,'.') - tmp - 1);
    } else {
      tmp1 = substr(fmt,tmp - fmt + 1, std::strchr(fmt,')') - tmp - 1);
    }
    *width = std::atoi(tmp1);
    if (tmp1!=NULL) free ((void *) tmp1);
    return *width;
}

char* substr(const char* S, const int pos, const int len)
{
    int i;
    char *SubS;
    if ( (size_t)pos+len <= std::strlen(S)) {
    SubS = (char *)malloc(len+1);
    if ( SubS == NULL ) IOHBTerminate("Insufficient memory for SubS.");
    for (i=0;i<len;i++) SubS[i] = S[pos+i];
    SubS[len] = '\0';
    } else {
      SubS = NULL;
    }
    return SubS;
}

#include<cctype>
void upcase(char* S)
{
/*  Convert S to uppercase     */
    int i,len;
    len = std::strlen(S);
    for (i=0;i< len;i++)
       S[i] = std::toupper(S[i]);
}

void IOHBTerminate(const char* message)
{
   std::fprintf(stderr,"%s",message);
   std::exit(1);
}


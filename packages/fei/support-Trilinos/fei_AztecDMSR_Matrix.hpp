#ifndef _AztecDMSR_Matrix_h_
#define _AztecDMSR_Matrix_h_

#ifdef HAVE_FEI_AZTECOO

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

//
// This class is a wrapper for the Aztec DMSR matrix data structure.
//
// Important usage notes:
//
// * The 'oneBased' argument to the constructor indicates whether the
//   matrix should use 1-based indices (row numbers and column indices) in
//   the input and output arguments to its interfaces (e.g., getRow),
//   with the exception of the update list -- keep reading.
//   'oneBased' should be 1 for 1-based indices, 0 for 0-based.
//   Here's the confusing part -- the update list should contain 0-based
//   indices, regardless of the value of 'oneBased'.  That's because the
//   update list gets used internally by Aztec functions that only work
//   in 0-based numbers.
//
// * The 'rowLengths' array, input argument to the configure function,
//   must contain the lengths of each row, *NOT* including the
//   coefficient on the diagonal.
//
#include <az_aztec.h>

namespace fei_trilinos {

class Aztec_Map;
class Aztec_Vector;

class AztecDMSR_Matrix {
    
  public:
    // Constructor.
    AztecDMSR_Matrix(Aztec_Map& map, int* update,
		     bool linearDistribution=false);

    //Copy constructor
    AztecDMSR_Matrix(const AztecDMSR_Matrix& src);

    virtual ~AztecDMSR_Matrix ();

    // Mathematical functions.
    void matvec(const Aztec_Vector& x, Aztec_Vector& y) const;

    void put(double s);
    void getDiagonal(Aztec_Vector& diagVector) const;

    const Aztec_Map& getAztec_Map() const {return(amap_);};

    int rowLength(int row) const;
    
    // ... to read matrix.
    void getRow(int row, int& length, double *coefs, int *colInd) const;
    void getRow(int row, int& length, double *coefs) const;
    void getRow(int row, int& length, int *colInd) const;

    /** write access to the diagonal entry for the specified row. */
    int setDiagEntry(int row, double value);

    /** Read-only access to the diagonal entry for the specified row. */
    double getDiagEntry(int row) const;

    // ... to write matrix.
    int putRow(int row, int len, const double *coefs, 
                       const int *colInd);

    int sumIntoRow(int numRows, const int* rows,
                 int numCols, const int* cols,
                 const double* const* coefs);

    int sumIntoRow(int row, int len, const double *coefs, 
                           const int *colInd);

    int addScaledMatrix(double scalar, const AztecDMSR_Matrix& source);

    void scale(double scalar);

    /** Special direct-access pointer function.
     */
    int getOffDiagRowPointers(int row, int*& colIndices, double*& coefs,
			      int& offDiagRowLength);

    /** Special function for enforcing an essential (dirichlet) BC.
     */
    //int enforceEssentialBCs(int numBCEqns, int* bcEqns, double* bcValues,
    //		    double* rhsVector);

    //inform about structure, so that val and bindx can be allocated.
    void allocate(int *rowLengths);

    //inform about structure, including column-indices, so that val and bindx
    //can be allocated *and* so that bindx can be populated.
    void allocate(int *rowLengths, const int* const* colIndices);

    //inform that data fill is complete, so AZ_transform can be called.
    void fillComplete();

    int getTransformedEqn(int eqn) const {
      if (isFilled_) {
	if(eqn<N_update_) return( update_[orderingUpdate_[eqn]] );
	else return( external_[eqn-N_update_] );
      }
      else {
	return( eqn );
      }
    }

    bool isFilled() const {return(isFilled_);};
    void setFilled(bool flag) {isFilled_ = flag;};
    bool isAllocated() const {return(isAllocated_);};
    void setAllocated(bool flag) {isAllocated_ = flag;};

    void copyStructure(AztecDMSR_Matrix& source);

    bool readFromFile(const char *filename);
    bool writeToFile(const char *fileName) const;
    bool rowMax() const {return true;};
    double rowMax(int row) const;
 
    int getNumNonZeros() {return(nnzeros_);};

    //Aztec-specific functions:

    AZ_MATRIX* getAZ_MATRIX_PTR() {return(Amat_);};

    int* getUpdate() {return(update_);};
    int* getUpdate_index() {return(update_index_);};

  private:
    int inUpdate(int globalIndex, int& localIndex) const;
    void messageAbort(const char* mesg);
    int insert(int item, int offset, int* list, int& len, int allocLen);
    int insert(double item, int offset, double* list, int& len, int allocLen);
    void expand_array(int*& array, int& arraylen, int newlen);
    void expand_array(double*& array, int& arraylen, int newlen);

    bool isFilled_;
    bool isAllocated_;
    bool linearDistribution_;

    int localOffset_;
    int localSize_;

    Aztec_Map& amap_;

    AZ_MATRIX* Amat_;

    bool arraysAllocated_;
    double *val;
    int *bindx;
    int *rowLengths_;
    int nnzeros_; //val and bindx are of length nnzeros_+1

    int N_update_;
    int* update_;

    int* update_index_;
    int* orderingUpdate_;
    int* external_;
    int* extern_index_;
    int* data_org_;

    int* tmp_array_;
    int tmp_array_len_;
    double* dtmp_array_;
    int dtmp_array_len_;

    bool azTransformed_;
};

}//namespace fei_trilinos

#endif //HAVE_FEI_AZTECOO

#endif

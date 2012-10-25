/*--------------------------------------------------------------------*/
/*    Copyright 2007 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_Vector_Local_hpp_
#define _fei_Vector_Local_hpp_

#include <fei_iosfwd.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_Vector.hpp>

namespace fei {

class Vector_Local : public fei::Vector {
 public:
  Vector_Local(fei::SharedPtr<fei::VectorSpace> vecSpace);

  virtual ~Vector_Local();

  const char* typeName() const { return("fei::Vector_Local"); }

    /** Update 'this' = b*'this' + a*x
     */
    int update(double a,
               const fei::Vector* x,
               double b);

    /** Use data in the underlying non-overlapping decomposition to update
        any shared data in the overlapping decomposition.

        If any data is already held for the shared positions, that data will
        be replaced by the data from the 'owning' processor.
    */
    int scatterToOverlap();

    void setCommSizes();

    /** Move any shared data from the overlapping decomposition to the
        underlying non-overlapping decomposition.
    */
    int gatherFromOverlap(bool accumulate = true);

    /** Set a specified scalar throughout the vector. */
    int putScalar(double scalar);

    /** Sum values into the vector, adding to any
        that may already exist at the specified indices.
    */
    int sumIn(int numValues, const int* indices, const double* values,
              int vectorIndex=0);

    /** Copy values into the vector, overwriting any that may already exist
        at the specified indices.
    */
    int copyIn(int numValues, const int* indices, const double* values,
               int vectorIndex=0);

    /** Obtain the VectorSpace associated with this vector.
     */
    fei::SharedPtr<fei::VectorSpace> getVectorSpace() const;

    /** Set the VectorSpace associated with this vector.
     */
    void setVectorSpace(fei::SharedPtr<fei::VectorSpace> vecSpace);

    /** Sum field data into the vector, adding to any coefficients that may
        already exist at the specified locations.
        If the specified fieldID doesn't exist at one or more of the specified
        IDs, then the corresponding positions in the data array will simply
        not be used.
    */
    int sumInFieldData(int fieldID,
                       int idType,
                       int numIDs,
                       const int* IDs,
                       const double* data,
                       int vectorIndex=0);

    /** Copy field data into the vector, overwriting any coefficients that may
        already exist at the specified locations.
        If the specified fieldID doesn't exist at one or more of the specified
        IDs, then the corresponding positions in the data array will simply
        not be used.
    */
    int copyInFieldData(int fieldID,
                        int idType,
                        int numIDs,
                        const int* IDs,
                        const double* data,
                        int vectorIndex=0);

    int copyInFieldDataLocalIDs(int fieldID,
                        int idType,
                        int numIDs,
                        const int* localIDs,
                        const double* data,
                        int vectorIndex=0);

    /** Copy field data out of the vector, into the caller-allocated data
        array.
        If the specified fieldID doesn't exist at one or more of the specified
        IDs, then the corresponding positions in the data array will simply
        not be referenced.
    */
    int copyOutFieldData(int fieldID,
                         int idType,
                         int numIDs,
                         const int* IDs,
                         double* data,
                         int vectorIndex=0);

    int copyOut(int numValues, const int* indices,
                double* values, int vectorIndex=0) const;

    int writeToFile(const char* filename,
                    bool matrixMarketFormat=true);

    int writeToStream(FEI_OSTREAM& ostrm,
                      bool matrixMarketFormat=true);

    std::vector<double>& getCoefs();

 private:
  int giveToVector(int numValues, const int* indices,
                           const double* values,
                           bool sumInto, int vectorIndex);

  int assembleFieldData(int fieldID,
                       int idType,
                       int numIDs,
                       const int* IDs,
                       const double* data,
                       bool sumInto,
                       int vectorIndex);

  int assembleFieldDataLocalIDs(int fieldID,
                       int idType,
                       int numIDs,
                       const int* localIDs,
                       const double* data,
                       bool sumInto,
                       int vectorIndex);

  fei::SharedPtr<fei::VectorSpace> vecSpace_;
  std::vector<double> coefs_;
  std::map<int,int> global_to_local_;
  std::vector<int> work_indices_;
};//class Vector_Local

}//namespace fei

#endif


#ifndef _fei_Aztec_Vector_hpp_
#define _fei_Aztec_Vector_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

//
// This class provides a vector that can be used with the AztecDMSR_matrix
// and the AztecDVBR_Matrix.
//
// An important restriction to note:
//
// * An Aztec_Vector can not be constructed until AFTER the AztecDMSR_Matrix
//   (or AztecDVBR_Matrix) that it is to be used with has been completely
//   initialized and filled (e.g., A.loadComplete() has been called (which
//   means, most importantly, that AZ_transform has been called)). This is
//   because the local data array for an aztec vector must be allocated
//   with enough extra space to hold 'boundary elements' that are exchanged
//   with other processors during the calculation of a parallel matrix-vector
//   product, and we don't know how much memory that requires until after
//   AZ_transform has been called.
//
// * Also, the calling code is responsible for keeping track of any 
//   re-ordering that AZ_transform has done. i.e., Aztec_Vector is just
//   like a raw array with respect to indexing of entries. If v is an
//   instantiation of an Aztec_Vector, then v[9] literally returns the
//   entry at position 9 (the 10th entry, since indexing is 0-based).
//

class Aztec_Map;

/**==========================================================================**/
class Aztec_Vector {
  public:
    // Constructor.
    Aztec_Vector(const Aztec_Map& map, int* data_org);

    Aztec_Vector(const Aztec_Vector& source);  // copy constructor

    virtual ~Aztec_Vector ();

    Aztec_Vector* newVector() const;

    // Mathematical functions.
    double dotProd (const Aztec_Vector& y) const;
    void scale (double s);
    void linComb (const Aztec_Vector& b, double s, const Aztec_Vector& c);
    void addVec (double s, const Aztec_Vector& c);
    double norm () const;
    double norm1 () const;
 
    // operator=
    Aztec_Vector& operator = (const Aztec_Vector& rhs);
    
    // Access functions.
    double& operator [] (int index);
    const double& operator [] (int index) const;
    
    void put (double scalar);

    const double* startPointer() const {return localCoeffs_;};

    //Special function
    bool readFromFile(const char *fileName);
    bool writeToFile(const char *fileName) const;
    
  protected:
    virtual void assign(const Aztec_Vector& rhs);
    
  private:
    void checkInput();
    int inUpdate(int globalIndex, int& localIndex) const;

    const Aztec_Map& amap_;
    double *localCoeffs_;        // local vector coefficients
    int length_;
};

#endif

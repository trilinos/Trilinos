/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_trilinos_macros.hpp>
#include <fei_iostream.hpp>

#ifdef HAVE_FEI_AZTECOO

#include <assert.h>
#include <cmath>         // needed for declaration of sqrt, abs
#include <unistd.h>
#include <stdio.h>

#include <fei_mpi.h>

#ifndef FEI_SER

#define AZTEC_MPI AZTEC_MPI
#define AZ_MPI AZ_MPI
#ifndef MPI
#define MPI MPI
#endif

#endif

#include <az_blas_wrappers.h>
#include <az_aztec.h>
#include <fei_SharedPtr.hpp>
#include <fei_Aztec_Map.hpp>
#include <fei_Aztec_LSVector.hpp>

namespace fei_trilinos {

//=============================================================================
Aztec_LSVector::Aztec_LSVector(fei::SharedPtr<Aztec_Map> map, int* data_org)
 : amap_(map)
{
//  Aztec_LSVector::Aztec_LSVector -- construct a zero filled distributed vector 
//  object.

    int tmp1 = data_org[AZ_N_internal];
    int tmp2 = data_org[AZ_N_border];
    length_ = tmp1 + tmp2 + data_org[AZ_N_external];

    localCoeffs_ = new double[length_];

    this->put(0.0);
}

/**=========================================================================**/
Aztec_LSVector::~Aztec_LSVector(){
    delete [] localCoeffs_;
    localCoeffs_ = NULL;

    length_ = 0;
}

/**=========================================================================**/
Aztec_LSVector::Aztec_LSVector(const Aztec_LSVector& source)
 : amap_(source.amap_)
{

/** Aztec_LSVector::Aztec_LSVector -- copy constructor **/

   length_ = source.length_;
   localCoeffs_ = new double[length_];

//  Since virtual dispatching will not occur in a constructor, 
//  specify call explicitly for clarity.

   Aztec_LSVector::operator=(source);
}


/**=========================================================================**/
double Aztec_LSVector::dotProd (const Aztec_LSVector& y) const {

/** Aztec_LSVector::dotProd --- returns dot product of two vectors **/

   int N_update = amap_->localSize();

   double *pv = (double*)localCoeffs_;
   double *py = (double*)y.startPointer();

   double dot = AZ_gdot(N_update, pv, py, amap_->getProcConfig());

   return dot;
}


/**=========================================================================**/
void Aztec_LSVector::put(double scalar) {

/** Aztec_LSVector::put -- fills vector with the value scalar **/

   for (int i = 0; i < length_; i++)
      localCoeffs_[i] = scalar;
}

/**=========================================================================**/
void Aztec_LSVector::scale (double s) {

/** Aztec_LSVector::scale -- scales a vector by a scalar **/

   int N_update = amap_->localSize();
   int one = 1;

   DSCAL_F77(&N_update, &s, localCoeffs_, &one);

   return;
}

/**=========================================================================**/
void Aztec_LSVector::addVec (double s, const Aztec_LSVector& c) {

/** Aztec_LSVector::addVec --- add multiple of a vector:  s*c **/

   int N_update = amap_->localSize();
   int one = 1;

   double *pv = localCoeffs_;
   double *pc = (double*)c.startPointer();

   DAXPY_F77(&N_update,&s,pc,&one,pv,&one);

   return;
}

/**=========================================================================**/
double Aztec_LSVector::norm (void) const  {

/** Aztec_LSVector::norm --- returns L2 norm of vector x **/

   int N_update = amap_->localSize();

   return(AZ_gvector_norm(N_update, 2,localCoeffs_, amap_->getProcConfig()));
}

/**=========================================================================**/
double Aztec_LSVector::norm1 (void) const  {

/** Aztec_LSVector::norm --- returns 1-norm of vector x **/

   int N_update = amap_->localSize();

   return(AZ_gvector_norm(N_update, 1,localCoeffs_, amap_->getProcConfig()));
}
 
/**=========================================================================**/
double& Aztec_LSVector::operator [] (int index)  {

//  Aztec_LSVector::operator [] --- return non-const reference 

   int offset = amap_->localOffset();

   return(localCoeffs_[index - offset]);
}

/**=========================================================================**/
const double& Aztec_LSVector::operator [] (int index) const  {

// Aztec_LSVector::operator [] --- return const reference 

   int offset = amap_->localOffset();

   return(localCoeffs_[index - offset]);
}

/**=========================================================================**/
Aztec_LSVector* Aztec_LSVector::newVector() const  {

/** Aztec_LSVector::newVector --- return a pointer to uninitialized object of same
    runtime type as *this **/
    
    Aztec_LSVector* p = new Aztec_LSVector(*this);

    return p;
}

/**=========================================================================**/
Aztec_LSVector& Aztec_LSVector::operator= (const Aztec_LSVector& rhs) {

//    check for special case of a=a
   if (this != &rhs) {
      assign(rhs);
   }
        
   return(*this);
}

/**=========================================================================**/
void Aztec_LSVector::assign(const Aztec_LSVector& rhs) {

   if ((amap_->globalSize() != rhs.amap_->globalSize()) ||
      (amap_->localSize() != rhs.amap_->localSize()) ) {
      fei::console_out() << "Aztec_LSVector::assign: ERROR, incompatible maps."
           << " Aborting assignment." << FEI_ENDL;
      return;
   }

   int N_update = amap_->localSize();
   double *pr = (double*)rhs.startPointer();

   for(int i=0; i<N_update; ++i) {
     localCoeffs_[i] = pr[i];
   }

   return;    
}

/**=========================================================================**/
bool Aztec_LSVector::readFromFile(const char *fileName) {
//
//For now, this function will just use a simple (not very scalable)
//strategy of having all processors read their own piece of the vector
//from the file.
//
//Important note: This function assumes that the equation numbers in
//the file are 1-based, and converts them to 0-based.
//
   int globalSize = amap_->globalSize();

   int i=0, nn, nnz;
   double value;
   bool binaryData;

   char line[128];

   if (fileName == NULL) {
      fei::console_out() << "Aztec_LSVector::readFromFile: ERROR, NULL fileName." << FEI_ENDL;
      return(false);
   }

   if (strstr(fileName, ".txt") != NULL) {
      binaryData = false;
   }
   else {
      fei::console_out() << "Aztec_LSVector::readFromFile: fileName doesn't contain "
           << "'.txt', assuming binary data..." << FEI_ENDL;
      binaryData = true;
   }

   FILE *file = fopen(fileName,"r");
   if (!file){
      fei::console_out() << "Aztec_LSVector: Error opening " << fileName << FEI_ENDL;
      return false;
   }

   if (binaryData) {
      fread((char *)&nn,sizeof(int),1,file);
      fread((char *)&nnz,sizeof(int),1,file);
   }
   else {
      do {
         fgets(line,128,file);
      } while(strchr(line,'%'));
      sscanf(line,"%d",&nn);
   }
   if (nn != globalSize) {
      fei::console_out() << "ERROR in Aztec_LSVector::readFromFile." << FEI_ENDL;
      fei::console_out() << "Vector in file has wrong dimension." << FEI_ENDL;
      fei::console_out() << "amap_->globalSize():" << globalSize << " file n:" << nn << FEI_ENDL;
      return(false);
   }

   int start = amap_->localOffset();
   int end = start + amap_->localSize() - 1;

   while (!feof(file)) {
      if(binaryData) {
         fread((char *)&i,sizeof(int),1,file);
         fread((char *)&value,sizeof(double),1,file);
      }
      else {
         fgets(line,128,file);
         sscanf(line,"%d %le",&i,&value);
      }
      if(feof(file))break;

      if ((start <= i) && (i <= end)) {
         localCoeffs_[i - start] = value;
      }
   } //end of 'while(!feof)

   return true;
}   

/**=========================================================================**/
bool Aztec_LSVector::writeToFile(const char *fileName) const {
   int i,p;
   int N_update = amap_->localSize();
   int start = amap_->localOffset();
   int numProcs = amap_->getProcConfig()[AZ_N_procs];
   int localRank = amap_->getProcConfig()[AZ_node];
   int masterRank = 0;
   MPI_Comm thisComm = amap_->getCommunicator();

   for(p=0; p<numProcs; p++){

      //A barrier inside the loop so each processor waits its turn.
      MPI_Barrier(thisComm);

      if (p == localRank){
         FILE *file = NULL;

         if (masterRank == localRank){
            //This is the master processor, open a new file.
            file = fopen(fileName,"w");

            //Write the vector dimension n into the file.
            fprintf(file,"%d\n",amap_->globalSize());
         }
         else {
            //This is not the master proc, open file for appending
            file = fopen(fileName,"a");
         }

         //Now loop over the local portion of the vector.
         for(i=0; i<N_update; i++) {
            fprintf(file,"%d %20.13e\n",start + i, localCoeffs_[i]);
         }

         fclose(file);
      }
   }

   return(true);
}

}//namespace fei_trilinos

#endif
//HAVE_FEI_AZTECOO

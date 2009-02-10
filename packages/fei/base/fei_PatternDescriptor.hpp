#ifndef _PatternDescriptor_hpp_
#define _PatternDescriptor_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <vector>

/** PatternDescriptor holds the information that the FEI implementation
  needs to know about an 'access pattern'. An
  access pattern is used to access sets of coefficients in the matrix or rhs
  vector through the FEI. These patterns can consist of rectangular 'regions'
  described by sets of nodes, for accessing those coefficients that correspond
  to nodal degrees-of-freedom, or by sets of fields associated with an element-
  degree-of-freedom, or a mixture of nodal and elemental degrees-of-freedom.
  
  Here is the information held in a PatternDescriptor:
  
     Pattern identifier
     number of 'row' entities
     list containing the number-of-fields-per-row
     table containing list of fieldIDs for each row-entity
     number of 'col' entities per row
     list containing the number-of-fields-per-col
     table containing list of fieldIDs for each col-entity
     interleaveStrategy
  
  Note: the tables of fieldIDs (for rows and cols) consist of lists
  of std::vector<int> objects. The length of the ith std::vector<int> in this list is
  the number-of-fields that correspond to the ith row of the pattern. The 
  number of std::vector<int>s in the list is the number of rows in the pattern.
  The list of std::vector<int>s is allocated when the function 'setNumRowIDs' is
  called (and similarly for setNumColIDsPerRow).
  
  Important Note!!! Don't change the size of the std::vector<int>s returned by the
  functions 'getNumFieldsPerRow' and 'getNumFieldsPerCol'. The lengths
  of these arrays must be the value set via 'setNumRowIDs' and
  'setNumColIDs', in order for the fieldIDs tables to be sized consistently.
*/

class PatternDescriptor {
 public:
   PatternDescriptor();
   virtual ~PatternDescriptor();

   int getPatternID() {return(patternID_);}
   void setPatternID(int ptrnID) {patternID_ = ptrnID;}

   int getNumRowIDs() {return(fieldsPerRow_.size());}
   int setNumRowIDs(int num);

   int getNumColIDsPerRow() {return(fieldsPerCol_.size());}
   int setNumColIDsPerRow(int num);

   std::vector<int>& getNumFieldsPerRow() {return(fieldsPerRow_);}
   std::vector<int>* getRowFieldIDs() {return(rowFieldIDs_);}

   std::vector<int>& getNumFieldsPerCol() {return(fieldsPerCol_);}
   std::vector<int>* getColFieldIDs() {return(colFieldIDs_);}

   int getInterleaveStrategy() const {return(interleaveStrategy_);}
   void setInterleaveStrategy(int strat) {interleaveStrategy_ = strat;}

 private:
   PatternDescriptor(const PatternDescriptor& src)
     : patternID_(0), fieldsPerRow_(), rowFieldIDs_(NULL),
     fieldsPerCol_(), colFieldIDs_(NULL), interleaveStrategy_(0) {}

   PatternDescriptor& operator=(const PatternDescriptor& src)
     {
       return(*this);
     }

   int patternID_;

   std::vector<int> fieldsPerRow_;
   std::vector<int>* rowFieldIDs_;

   std::vector<int> fieldsPerCol_;
   std::vector<int>* colFieldIDs_;

   int interleaveStrategy_;
};

#endif


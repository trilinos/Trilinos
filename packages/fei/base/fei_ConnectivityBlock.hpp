/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_ConnectivityBlock_hpp_
#define _fei_ConnectivityBlock_hpp_

#include <fei_macros.hpp>

#include <map>
#include <vector>

namespace fei {
  class Pattern;

  /**
   class to hold attributes of a connectivity-block (for example,
   an element-block). This class can handle
   connectivities for any kind of mesh-objects, not just elements,
   though elements are of course the most common.
  */
  class ConnectivityBlock {
  public:
    /** constructor */
    ConnectivityBlock(int blockID,
                      fei::Pattern* pattern,
                      int numConnectivities);
    /** constructor */
    ConnectivityBlock(int blockID,
                      fei::Pattern* rowpattern, fei::Pattern* colpattern,
                      int numConnectivities);
    /** constructor */
    ConnectivityBlock(int numRowIDs,
                      const int* rowIDs,
                      const int* rowOffsets,
                      bool offsets_are_lengths = false);

    /** constructor */
    ConnectivityBlock(int fieldID,
                      int numRowIDs,
                      const int* rowIDs,
                      const int* rowOffsets,
                      bool offsets_are_lengths = false);

    /** destructor */
    virtual ~ConnectivityBlock();

    /** get block-identifier */
    int getBlockID() const { return(blockID_); }

    /** get pattern that defines the layout of dofs in the
      row-dimension for block's contributions */
    const fei::Pattern* getRowPattern() const { return(pattern_); }

    /** get pattern that defines the layout of dofs in the
      row-dimension for block's contributions */
    fei::Pattern* getRowPattern() { return(pattern_); }

    void setRowPattern(fei::Pattern* pattern) { pattern_ = pattern; }

    /** get pattern that defines the layout of dofs in the
      column-dimension for block's contributions. probably null
     if this block is made up of symmetric contributions. */
    const fei::Pattern* getColPattern() const { return(colPattern_); }

    /** get pattern that defines the layout of dofs in the
      column-dimension for block's contributions. probably null
     if this block is made up of symmetric contributions. */
    fei::Pattern* getColPattern() { return(colPattern_); }

    void setColPattern(fei::Pattern* pattern) { colPattern_ = pattern; }

    /** get map of connectivity-ids with associated offsets
    */
    const std::map<int,int>& getConnectivityIDs() const { return( connIDsOffsetMap_ ); }

    /** get map of connectivity-ids with associated offsets
    */
    std::map<int,int>& getConnectivityIDs() { return( connIDsOffsetMap_ ); }

    /** get vector of connectivity-offsets. Only available if this
      object was constructed using constructor 3 or 4. Power users only.
    */
    std::vector<int>& getConnectivityOffsets()
      { return(connectivityOffsets_); }

    /** get array of row-connectivities */
    std::vector<int>& getRowConnectivities()
      { return(connectivities_); }

    /** get array of column-connectivities */
    std::vector<int>& getColConnectivities()
      { return(colConnectivities_); }

    /** get row-connectivity for a specified ID */
    const int* getRowConnectivity(int ID) const;
    /** get column-connectivity for a specified ID */
    const int* getColConnectivity(int ID) const;
    /** get row-connectivity for a specified ID */
    int* getRowConnectivity(int ID);
    /** get column-connectivity for a specified ID */
    int* getColConnectivity(int ID);

    /** query whether block is symmetric */
    bool isSymmetric() const { return( isSymmetric_ ); }

    /** implementation detail for power-users */
    void setIsDiagonal(bool flag) { isDiagonal_ = flag; }
    /** implementation detail for power-users */
    bool isDiagonal() const { return( isDiagonal_ ); }

    /** query whether block has a field-id */
    bool haveFieldID()
      { return( haveFieldID_ ); }

    /** return block's field-id */
    int fieldID()
      { return( fieldID_ ); }

  private:
    int blockID_;
    fei::Pattern* pattern_;
    fei::Pattern* colPattern_;
    bool isSymmetric_;
    bool isDiagonal_;

    std::map<int,int> connIDsOffsetMap_;

    std::vector<int> connectivityOffsets_;

    int numRecordsPerConnectivity_;
    std::vector<int> connectivities_;
    int numRecordsPerColConnectivity_;
    std::vector<int> colConnectivities_;

    int fieldID_;
    bool haveFieldID_;

  };//class ConnectivityBlock
} //namespace fei

#endif // _fei_ConnectivityBlock_hpp_


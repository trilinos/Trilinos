/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#ifndef _fei_ConnectivityBlock_hpp_
#define _fei_ConnectivityBlock_hpp_

#include <fei_macros.hpp>
#include <fei_MapType.hpp>
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

    /** get data structure of connectivity-ids with associated offsets
    	 */
    MapIntInt& getConnectivityIDs(){return(connIDsOffsetMap_);}

    /** get data structure of connectivity-ids with associated offsets
    	 */
    const MapIntInt& getConnectivityIDs()const {return(connIDsOffsetMap_);}

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

    MapIntInt connIDsOffsetMap_;
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


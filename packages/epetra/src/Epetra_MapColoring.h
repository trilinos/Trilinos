
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef EPETRA_MAPCOLORING_H
#define EPETRA_MAPCOLORING_H

//! Epetra_MapColoring: A class for coloring Epetra_Map and Epetra_BlockMap objects.

/*! This class allows the user to associate an integer value, i.e., a color, to each element of 
    an existing Epetra_Map or Epetra_BlockMap object.  Colors may be assigned at construction, or 
    via set methods.  Any elements that are not explicitly assigned a color are assigned the color 
    0 (integer zero).  Color information may be accessed in three basic ways:
  <ol> 
  <li> By local element ID (LID) - Returns the color of a specified LID, where the LID is associated
       with the Epetra_Map or BlockMap that was passed in to the Epetra_MapColoring constructor.
  <li> By global element ID (GID) - Returns the color of the specified GID.  There two methods
       for accessing GIDs, one assumes the request is for GIDs owned by the calling processor,
       the second allows arbitrary requested for GIDs, as long as the GID is defined on some processor
       for the Epetra_Map or Epetra_BlockMap.
  <li> By color groups - LIDs and GIDs are grouped by color so that all elements of a given color can
       be accessed.
  </ol>

*/

#include "Epetra_BlockMap.h"

class Epetra_MapColoring {
    
  public:

  //@{ \name Constructors/destructors.
  //! Epetra_MapColoring constructor.
  /*!
    \param In
            ColorList - Array of dimension Map.NumMyElements() containing the list of colors
            that should be assigned the map elements on this processor. If this argument is
	    set to 0 (zero), all elements will initially be assigned color 0 (zero).  Element
	    colors can be modified by using methods described below.
    
    \param In
            Map - An Epetra_Map or Epetra_BlockMap (Note: Epetra_BlockMap is a base class of
	    Epetra_Map, so either can be passed in to this constructor.

    \return Pointer to a Epetra_MapColoring object.

  */ 
  Epetra_MapColoring(int * ColorList, const Epetra_BlockMap& Map);

  //! Epetra_MapColoring copy constructor.
  
  Epetra_MapColoring(const Epetra_MapColoring& map);
  
  //! Epetra_MapColoring destructor.
  
  virtual ~Epetra_MapColoring(void);
  //@}
  
  //@{ \name Set Color methods.
  //! Set the color of the given GID to Color, for locally owned GIDs only.
  /*! 
    \param In
            NumIDs - Number of GIDs to color.
    \param In
            GIDList - List of GIDs to color.
    \param In
            ColorList -List of colors such that ColorList[i] should be the color assigned to 
	    GIDList[i].

   \return Returns 0 if no errors, -1 if one of the GIDs in GIDList is not part of the map.
  */ 
  int SetRemoteGIDColors(int NumIDs, const int *GIDList, const int *ColorList);

  //! LID element color assignment method.
  /*! Allows color assignment of ith LID: colormap[i] = color
    \return MapColor(LID).
  */
  int& operator [] (int LID) { return Color_[LID];};

  //! GID element color assignment method, Note:  Valid only for GIDs owned by calling processor.
  /*! Allows color assignment of specified GID \e only if the GID is owned by map on 
      the calling processor. If you are unsure about the ownership, check by using the MyGID()
      method on the map object.
    \return MapColor(GID).
  */
  int& operator () (int GID) {return Color_[Map_->LID(GID)];};
  //@}
  
  //@{ \name Local/Global color accessor methods.
  //! LID element color access method.
  /*! Returns color  of ith LID: colormap[i] = color
    \return MapColor[LID].
  */
  const int& operator [] (int LID) const { return Color_[LID];};

  //! GID element color assignment method, Note:  Valid only for GIDs owned by calling processor.
  /*! Allows color assignment of specified GID \e only if the GID is owned by map on 
      the calling processor. If you are unsure about the ownership, check by using the MyGID()
      method on the map object.
    \return MapColor(GID).
  */
    const int& operator () const (int GID) {return Color_[Map_->LID(GID)];};
  //@}
  
  //@{ \name Access by color methods.
    //! Returns number of colors, computed to be MaxColorID - MinColorID + 1.
    int NumColors() const {return(NumColors_);};

    //! Returns minimum ColorID, least integer value of all color values.
    int MinColorID() const {return(MinColorID_);};

    //! Returns maximum ColorID, largest integer value of all color values.
    int MaxColorID() const {return(MaxColorID_);};

    //! Returns number of local map elements having specified ColorID
    int NumElementsInColor(int ColorID) const {return(NumElements_[ColorID-MinColorID_]);};

    //! Returns pointer to array of Map LIDs associated with the specified color.
    /*! Returns a pointer to a list of Map LIDs associated with the specified color. 
        This is a purely local list with no information about other processors.  If there
	are no LIDs associated with the specified color, the pointer is set to zero.
    */
    const int * ColorLIDList(int ColorID) const;
  //@}

 private:

    Epetra_BlockMap * Map_;
    int NumColors_;
    int MinColorID_;
    int MaxColorID_;
    int * NumElements_;
    int * Color_;
    
};

#endif /* EPETRA_MAPCOLORING_H */

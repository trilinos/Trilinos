// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_LAPACK_DATAOUTPUT_H
#define LOCA_LAPACK_DATAOUTPUT_H

#include <fstream>

#include "LOCA_Abstract_DataOutput.H"

namespace LOCA {
  namespace LAPACK {
    class Group;
  }
}

namespace LOCA {
  namespace LAPACK {

    /*!
      \brief Implements the LOCA::Abstract::DataOutput interface for LAPACK groups.
    */
    /*!
      LAPACK group data is saved by writing the solution vector and parameters
      to a file stream in ASCII form with all solution vector components and 
      parameters on one line followed by a new line character.
    */
    class DataOutput : public LOCA::Abstract::DataOutput {

    public: 
      //! Constructor
      DataOutput(fstream& fs);

      //! Destructor
      virtual ~DataOutput();
      
      //! Assignment operator
      virtual Abstract::DataOutput& operator = (const Abstract::DataOutput& source);
      //! Assignment operator
      virtual DataOutput& operator = (const DataOutput& source);

      //! Save group data to a file
      virtual void saveGroupData(const Abstract::Group& grp);

      //! Save group data to a file
      virtual void saveGroupData(const Group& grp);

    protected:
      //! File stream to write data to
      fstream& file;
    };
  } // namespace LAPACK
} // namespace LOCA

#endif

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

#ifndef LOCA_BLAS_DATAOUTPUT_H
#define LOCA_BLAS_DATAOUTPUT_H

#include <fstream>

#include "LOCA_Abstract_DataOutput.H"

namespace LOCA {
  namespace BLAS {
    class Group;
  }
}

namespace LOCA {
  namespace BLAS {

    class DataOutput : public LOCA::Abstract::DataOutput {

    public: 
      //! Constructor
      DataOutput(fstream& fs);

      //! Destructor
      virtual ~DataOutput();
      
      virtual Abstract::DataOutput& operator = (const Abstract::DataOutput& source);
      virtual DataOutput& operator = (const DataOutput& source);
      virtual void saveGroupData(const Abstract::Group& grp);
      virtual void saveGroupData(const Group& grp);

    protected:
      fstream& file;
    };
  } // namespace BLAS
} // namespace LOCA

#endif

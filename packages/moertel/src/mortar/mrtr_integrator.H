/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Glen Hansen (Glen.Hansen@inl.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*!
 * \file mrtr_integrator.H
 *
 * \class MOERTEL::Integrator
 *
 * \brief A class to perform integration of mass matrices on the overlap of
          2 segments in 1D and 2D
 *
 * \date Last update do Doxygen: 20-March-06
 *
 */
#ifndef MOERTEL_INTEGRATOR_H
#define MOERTEL_INTEGRATOR_H

#include <ctime>
#include <iostream>
#include <iomanip>
#include <vector>

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_RefCountPtr.hpp"

#include "mrtr_overlap.H"

// ----------   User Defined Includes   ----------

/*!
\brief MOERTEL: namespace of the Moertel package

The Moertel package depends on \ref Epetra, \ref EpetraExt, \ref Teuchos,
\ref Amesos, \ref ML and \ref AztecOO:<br>
Use at least the following lines in the configure of Trilinos:<br>
\code
--enable-moertel 
--enable-epetra 
--enable-epetraext
--enable-teuchos 
--enable-ml
--enable-aztecoo --enable-aztecoo-teuchos 
--enable-amesos
\endcode

*/
namespace MOERTEL
{

// forward declarations
class Interface;
class Segment;
class Node;

/*!
\class Integrator

\brief <b> A class to perform Gaussian integration and assembly of mass matrices on the overlap of
          2 segments in 1D and 2D </b>



\author Glen Hansen (Glen.Hansen@inl.gov)

*/
class Integrator 
{
public:
  
  // @{ \name Constructors and destructors

  /*!
  \brief Constructor
  
  Constructs an instance of this class.<br>
  Note that this is \b not a collective call as overlaps are integrated in parallel by
  individual processes.
  
  \param ngp : number of gaussian points to be used
  \param oneD : flag indicating whether 1D or 2D regions shall be integrated
  \param outlevel : Level of output information written to stdout ( 0 - 10 )
  */
  explicit Integrator(int ngp, bool oneD, int outlevel);
  
  /*!
  \brief Destructor
  */
  virtual ~Integrator();
  
  //@}
  // @{ \name Public members

  /*!
  \brief Return the level of output written to stdout ( 0 - 10 )
  
  */
  int OutLevel() { return outputlevel_; }
  
  /*!
  \brief Return number of Gaussian integration points to be used
  
  */
  inline int Ngp() { return ngp_; }
  
  /*!
  \brief Return coordinate of a specific Gaussian point in 1D segment
  
  \param gp : Number of Gaussian point to get the coordinate for
  */
  inline double Coordinate(int gp) { return coords_[gp]; }

  /*!
  \brief Return coordinates of a specific Gaussian point in 2D segment
  
  \param gp : Number of Gaussian point to get the coordinate for
  */
  inline double* Coordinate(int* gp) { return &coords_[*gp*2]; }
  
  /*!
  \brief Return weight for given Gaussian point
  
  \param gp : Number of Gaussian point to get the coordinate for
  */
  inline double Weight(int gp) { return weights_[gp]; }


  /*!
  \brief Integrate mass matrix 'M' on a 1D overlap between a slave and a master segment

  Integrate over overlap the trace space shape function of the mortar side times
  the Lagrange multiplier space function of the slave side
  
  \param sseg : Slave Segment
  \param sxia : lower bound of overlap in slave segment coordinates
  \param sxib : upper bound of overlap in slave segment coordinates
  \param mseg : Mortar Segment
  \param mxia : lower bound of overlap in mortar segment coordinates
  \param mxib : upper bound of overlap in mortar segment coordinates
  */
  Epetra_SerialDenseMatrix* Integrate(MOERTEL::Segment& sseg, double sxia, double sxib,
                                      MOERTEL::Segment& mseg, double mxia, double mxib);

  /*!
  \brief Assemble integration result '-M' into global matrix 'M'
  
  \param inter : Interface
  \param sseg : Slave Segment
  \param mseg : Mortar Segment
  \param M : global sparse matrix 'M'
  \param Mdense : local dense matrix from integration of overlap between sseg and mseg
  */
  bool Assemble(MOERTEL::Interface& inter, MOERTEL::Segment& sseg, MOERTEL::Segment& mseg, 
                Epetra_CrsMatrix& M, Epetra_SerialDenseMatrix& Mdense);


  /*!
  \brief Integrate mass matrix 'D' on a 1D slave segment overlap

  Integrate over overlap the trace space shape function of the slave side times
  the Lagrange multiplier space function of the slave side
  
  \param sseg : Slave Segment
  \param sxia : lower bound of overlap in slave segment coordinates
  \param sxib : upper bound of overlap in slave segment coordinates
  */
  Epetra_SerialDenseMatrix* Integrate(MOERTEL::Segment& sseg, double sxia, double sxib);

  /*!
  \brief Assemble integration result 'D' into global matrix 'D'
  
  \param inter : Interface
  \param sseg : Slave Segment
  \param D : global sparse matrix 'D'
  \param Ddense : local dense matrix from integration of overlap between sseg and mseg
  */
  bool Assemble(MOERTEL::Interface& inter, MOERTEL::Segment& sseg,  
                Epetra_CrsMatrix& D, Epetra_SerialDenseMatrix& Ddense);


  /*!
  \brief Integrate modification to mass matrix 'M' on a 1D overlap between a slave and a master segment

  Integrate over overlap the modification of the trace space shape function of the mortar side times
  the Lagrange multiplier space function of the slave side.
  This modification due to B. wohlmuth improves behaviour for curved interfaces 
  
  \param sseg : Slave Segment
  \param sxia : lower bound of overlap in slave segment coordinates
  \param sxib : upper bound of overlap in slave segment coordinates
  \param mseg : Mortar Segment
  \param mxia : lower bound of overlap in mortar segment coordinates
  \param mxib : upper bound of overlap in mortar segment coordinates
  */
  Epetra_SerialDenseMatrix* Integrate_2D_Mmod(MOERTEL::Segment& sseg, double sxia, double sxib,
                                              MOERTEL::Segment& mseg, double mxia, double mxib);

  // Assemble the result -Mdense from the integration above into M
  /*!
  \brief Assemble modification integration result '-DELTA_M' into global matrix 'M'
  
  \param inter : Interface
  \param sseg : Slave Segment
  \param mseg : Mortar Segment
  \param M : global sparse matrix 'M'
  \param Mmod : local dense matrix from integration of modification on overlap between sseg and mseg
  */
#if 0
  bool Assemble_2D_Mod(MOERTEL::Interface& inter, MOERTEL::Segment& sseg, MOERTEL::Segment& mseg, 
                       Epetra_CrsMatrix& M, Epetra_SerialDenseMatrix& Mmod);
#endif
  bool Assemble_2D_Mod(MOERTEL::Interface& inter, MOERTEL::Segment& sseg, MOERTEL::Segment& mseg, 
                       Epetra_SerialDenseMatrix& Mmod);


  /*!
  \brief Integrate a 2D triangle overlap segment, master part M and slave part D

  Integrate a triangle which is part of the discretization of a overlap polygon between
  a slave and a mortar segment.
  
  \param actseg : triangle to be integrated over
  \param sseg : Slave Segment
  \param mseg : Mortar Segment
  \param Ddense : Dense matrix holding integration result 'D'
  \param Mdense : Dense matrix holding integration result 'M'
  \param overlap : the overlap class holds all neccessary function values
  \param eps : Threshold, skip triangles with an area < eps
  */
  bool Integrate(Teuchos::RCP<MOERTEL::Segment> actseg,
                 MOERTEL::Segment& sseg, 
                 MOERTEL::Segment& mseg,
                 Epetra_SerialDenseMatrix** Ddense, 
                 Epetra_SerialDenseMatrix** Mdense, 
                 MOERTEL::Overlap& overlap, double eps,
                 bool exactvalues);

  /*!
  \brief Assemble integration result 'D' into Node (2D interfaces only)
  
  \param inter : Interface
  \param sseg : Slave Segment
  \param Ddense : local dense matrix from integration of overlap
  */
  bool Assemble(MOERTEL::Interface& inter,MOERTEL::Segment& sseg,Epetra_SerialDenseMatrix& Ddense);

  /*!
  \brief Assemble integration result 'M' into Node (2D interfaces only)
  
  \param inter : Interface
  \param sseg : Slave Segment
  \param mseg : Mortar Segment
  \param Mdense : local dense matrix from integration of overlap
  */
  bool Assemble(MOERTEL::Interface& inter,MOERTEL::Segment& sseg,MOERTEL::Segment& mseg,
                Epetra_SerialDenseMatrix& Mdense);

  //@}

private:  
  // don't want = operator
  Integrator operator = (const Integrator& old);
  // don't want copy-ctor
  Integrator(MOERTEL::Integrator& old);

private:

  bool           oneD_;    // dimension of the integration, true if 1-dimensional
  int            ngp_;     // order of the integration
  int            outputlevel_; // output level (0-10)
  std::vector<double> coords_;  // vector holding gaussian point coordinates
  std::vector<double> weights_; // vector holding gaussian point weights


};

} // namespace MOERTEL

#endif // MOERTEL_INTEGRATOR_H

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
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
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
 * \class MoertelT::Integrator
 *
 * \brief A class to perform integration of mass matrices on the overlap of
          2 segments in 1D and 2D
 *
 * \date Last update do Doxygen: 20-March-06
 *
 */
#ifndef MOERTEL_INTEGRATORT_H
#define MOERTEL_INTEGRATORT_H

#include <ctime>
#include <iostream>
#include <iomanip>
#include <vector>

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_RefCountPtr.hpp"

#include "mrtr_overlap.hpp"

// ----------   User Defined Includes   ----------

/*!
\brief MoertelT: namespace of the Moertel package

The Moertel package depends on \ref Tpetra, \ref Teuchos,
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
namespace MoertelT
{

// forward declarations
template <class ST,
          class LO,
          class GO,
          class N >
class InterfaceT;
class Segment;
class Node;

/*!
\class Integrator

\brief <b> A class to perform Gaussian integration and assembly of mass matrices on the overlap of
          2 segments in 1D and 2D </b>



\author Glen Hansen (gahanse@sandia.gov)

*/
template <class ST,
          class LO,
          class GO,
          class N >
class IntegratorT
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
  explicit IntegratorT(int ngp, bool oneD, int outlevel);
  
  /*!
  \brief Destructor
  */
  virtual ~IntegratorT();
  
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
  Teuchos::SerialDenseMatrix<LO, ST>* Integrate(MOERTEL::Segment& sseg, double sxia, double sxib,
                                      MOERTEL::Segment& mseg, double mxia, double mxib);

  /*!
  \brief Assemble integration result '-M' into global matrix 'M'
  
  \param inter : Interface
  \param sseg : Slave Segment
  \param mseg : Mortar Segment
  \param M : global sparse matrix 'M'
  \param Mdense : local dense matrix from integration of overlap between sseg and mseg
  */
  bool Assemble(MoertelT::InterfaceT<ST, LO, GO, N>& inter, MOERTEL::Segment& sseg, 
                Tpetra::CrsMatrix<ST, LO, GO, N>& D, Teuchos::SerialDenseMatrix<LO, ST>& Ddense);

  /*!
  \brief Integrate mass matrix 'D' on a 1D slave segment overlap

  Integrate over overlap the trace space shape function of the slave side times
  the Lagrange multiplier space function of the slave side
  
  \param sseg : Slave Segment
  \param sxia : lower bound of overlap in slave segment coordinates
  \param sxib : upper bound of overlap in slave segment coordinates
  */
  Teuchos::SerialDenseMatrix<LO, ST>* Integrate(MOERTEL::Segment& sseg, double sxia, double sxib);

  /*!
  \brief Assemble integration result 'D' into global matrix 'D'
  
  \param inter : Interface
  \param sseg : Slave Segment
  \param D : global sparse matrix 'D'
  \param Ddense : local dense matrix from integration of overlap between sseg and mseg
  */
  bool Assemble(MoertelT::InterfaceT<ST, LO, GO, N>& inter, 
                MOERTEL::Segment& sseg,  MOERTEL::Segment& mseg,
                Tpetra::CrsMatrix<ST, LO, GO, N>& D, Teuchos::SerialDenseMatrix<LO, ST>& Ddense);


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
  Teuchos::SerialDenseMatrix<LO, ST>* Integrate_2D_Mmod(MOERTEL::Segment& sseg, double sxia, double sxib,
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
  bool Assemble_2D_Mod(MoertelT::Interface& inter, MOERTEL::Segment& sseg, MOERTEL::Segment& mseg, 
                       Tpetra_CrsMatrix& M, Teuchos::SerialDenseMatrix<LocalOrdinal, Scalar>& Mmod);
#endif
  bool Assemble_2D_Mod(MoertelT::InterfaceT<ST, LO, GO, N>& inter, MOERTEL::Segment& sseg, MOERTEL::Segment& mseg, 
                       Teuchos::SerialDenseMatrix<LO, ST>& Mmod);


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
                 Teuchos::SerialDenseMatrix<LO, ST>** Ddense, 
                 Teuchos::SerialDenseMatrix<LO, ST>** Mdense, 
                 MOERTEL::Overlap<MoertelT::InterfaceT<ST, LO, GO, N> >& overlap, double eps,
                 bool exactvalues);

  /*!
  \brief Assemble integration result 'D' into Node (2D interfaces only)
  
  \param inter : Interface
  \param sseg : Slave Segment
  \param Ddense : local dense matrix from integration of overlap
  */
  bool Assemble(MoertelT::InterfaceT<ST, LO, GO, N>& inter,
     MOERTEL::Segment& sseg,Teuchos::SerialDenseMatrix<LO, ST>& Ddense);

  /*!
  \brief Assemble integration result 'M' into Node (2D interfaces only)
  
  \param inter : Interface
  \param sseg : Slave Segment
  \param mseg : Mortar Segment
  \param Mdense : local dense matrix from integration of overlap
  */
  bool Assemble(MoertelT::InterfaceT<ST, LO, GO, N>& inter,MOERTEL::Segment& sseg,MOERTEL::Segment& mseg,
                Teuchos::SerialDenseMatrix<LO, ST>& Mdense);

  //@}

private:  
  // don't want = operator
  IntegratorT operator = (const IntegratorT<ST, LO, GO, N>& old);
  // don't want copy-ctor
  IntegratorT(MoertelT::IntegratorT<ST, LO, GO, N>& old);

private:

  bool           oneD_;    // dimension of the integration, true if 1-dimensional
  int            ngp_;     // order of the integration
  int            outputlevel_; // output level (0-10)
  std::vector<double> coords_;  // vector holding gaussian point coordinates
  std::vector<double> weights_; // vector holding gaussian point weights


};

} // namespace MoertelT

#ifndef HAVE_MOERTEL_EXPLICIT_INSTANTIATION
#include "Moertel_IntegratorT.hpp"
#include "Moertel_IntegratorT_Def.hpp"
#endif

#endif // MOERTEL_INTEGRATORT_H

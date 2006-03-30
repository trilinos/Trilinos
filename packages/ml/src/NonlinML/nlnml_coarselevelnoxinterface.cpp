/*
#@HEADER
# ************************************************************************
#
#               ML: A Multilevel Preconditioner Package
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
// ML-headers
#include "ml_common.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

// ----------   Includes   ----------
#include <ctime>
#include <iostream>
#include "nlnml_coarselevelnoxinterface.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
NLNML::NLNML_CoarseLevelNoxInterface::NLNML_CoarseLevelNoxInterface()
{

  return;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_CoarseLevelNoxInterface::recreate()
{

  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
NLNML::NLNML_CoarseLevelNoxInterface::~NLNML_CoarseLevelNoxInterface()
{
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate nonlinear function (public, derived)             m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_CoarseLevelNoxInterface::computeF(
                                 const Epetra_Vector& x, Epetra_Vector& F, 
			         const FillType fillFlag)
{
  return true;
}

/*----------------------------------------------------------------------*
 |  restrict from fine to this level (public)                 m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_CoarseLevelNoxInterface::restrict_fine_to_this()
{
  return;
}

/*----------------------------------------------------------------------*
 |  prolongate from this level to fine (public)               m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_CoarseLevelNoxInterface::prolong_this_to_fine()
{
  return;
}

/*----------------------------------------------------------------------*
 |  restrict from this to next coarser level (public)         m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_CoarseLevelNoxInterface::restrict_to_next_coarser_level()
{
  return;
}


/*----------------------------------------------------------------------*
 |  prolongate from next coarser level to this level (public) m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_CoarseLevelNoxInterface::prolong_to_this_level()
{
  return;
}

/*----------------------------------------------------------------------*
 |  set ptr to all prolongation operators (public)            m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_CoarseLevelNoxInterface::setP()
{
  return;
}


/*----------------------------------------------------------------------*
 |  set modified system (public)                              m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_CoarseLevelNoxInterface::setModifiedSystem()
{
  return;
}

/*----------------------------------------------------------------------*
 |  make application apply all constraints (public)           m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_CoarseLevelNoxInterface::ApplyAllConstraints()
{
  return;
}


/*----------------------------------------------------------------------*
 |  return this levels blockmap (public)                      m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_CoarseLevelNoxInterface::BlockMap()
{
  return;
}









#endif

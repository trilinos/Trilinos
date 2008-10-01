// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include "Dimension.hpp"

const char * Dim::name() const 
{ static const char n[] = "Dim" ; return n ; }
const Dim & Dim::tag() 
{ static const Dim myself ; return myself ; }

const char * QuadPoint::name() const 
{ static const char n[] = "QuadPoint" ; return n ; }
const QuadPoint & QuadPoint::tag() 
{ static const QuadPoint myself ; return myself ; }

const char * Node::name() const 
{ static const char n[] = "Node" ; return n ; }
const Node & Node::tag() 
{ static const Node myself ; return myself ; }

const char * Point::name() const 
{ static const char n[] = "Point" ; return n ; }
const Point & Point::tag() 
{ static const Point myself ; return myself ; }

const char * Cell::name() const 
{ static const char n[] = "Cell" ; return n ; }
const Cell & Cell::tag() 
{ static const Cell myself ; return myself ; }

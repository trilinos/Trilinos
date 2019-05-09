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
 * \file mrtr_segment_bilinearquad.H
 *
 * \class MOERTEL::Segment_BiLinearQuad
 *
 * \brief A class to define a 4-noded quadrilateral 2D Segment
 *
 * \date Last update do Doxygen: 20-March-06
 *
 */
#ifndef MOERTEL_SEGMENT_BILINEARQUAD_H
#define MOERTEL_SEGMENT_BILINEARQUAD_H

#include <ctime>
#include <iostream>
#include <map>

#include "mrtr_segment.H"
#include "mrtr_functions.H"
#include "mrtr_utils.H"


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

/*!
\class Segment_BiLinearQuad

\brief <b> A class to define a 4-noded quadrilateral 2D Segment</b>

This class defines a 4-noded linear 2D quadrilateral interface surface segment.

<b>Important:</b><br>
Upon construction, the user must assure that the nodes attached to this segment are
given in counter-clockwise order such that the outward normal to the domain
points out from the screen:<br>
<pre>
        Node 3     Node 2
          o----------o
          |          |
          |          |           
          |          |                   <-------| 
          |          |            domain surface | 
          o----------o                      ------  
        Node 0     Node 1

</pre>
The reason for this is that the outward normal to the interface is implicitly defined by
the orientation of the segments. There is no way for the MOERTEL package to check the
orientation of the segments as they are passed in by the user and most obscure results
will be produced when the orientation is uncorrect!

<b>Important:</b><br>
There is currently no full support for quadrilateral interface discretizations. However,
when quads are added to a \ref MOERTEL::Interface they will be split into
2 triangles internally that are then used to perform the integration. The orientation of
the 2 triangles resulting from 1 quad is consistent with the orientation of the quad.

The \ref MOERTEL::Segment_BiLinearQuad class supports the std::ostream& operator <<

\author Glen Hansen (gahanse@sandia.gov)

*/
class  Segment_BiLinearQuad : public MOERTEL::Segment
{
public:
  
  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  \param Id : A unique positive Segment id. Does not need to be continous among segments
  \param nnode : Number of nodes this segment is attached to
  \param nodeId : Pointer to vector length nnode holding unique positive 
                  node ids of nodes this segment is attached to <b> in counterclock-wise
                  orientation looking at the interface such that the outward normal to the
                  domain points out from the screen</b>
  \param outlevel : Level of output to stdout to be generated by this class (0-10)
  */
  explicit Segment_BiLinearQuad(int id, int nnode, int* nodeId, int out);

  explicit Segment_BiLinearQuad(int id, const std::vector<int>&, int out);
  
  /*!
  \brief Empty Constructor
  
  To be used for communicating Segment_BiLinearTri classes among processes together with 
  \ref Pack and \ref UnPack
  */
  explicit Segment_BiLinearQuad(int out);
  
  /*!
  \brief Copy Constructor

  Makes a deep copy of a Segment_BiLinearQuad
  
  */
  Segment_BiLinearQuad(MOERTEL::Segment_BiLinearQuad& old);
  
  /*!
  \brief Destructor

  */
  virtual ~Segment_BiLinearQuad();

  /*!
  \brief Deep copy this instance of Segment_BiLinearQuad and return pointer to the copy
  
  */
  MOERTEL::Segment* Clone();

  //@}
  // @{ \name Public members

  /*!
  \brief not implemented
  
  */
  int* Pack(int* size);
  
  /*!
  \brief not implemented
  
  */
  bool UnPack(int* pack);
  
  
  /*!
  \brief not implemented
  
  */
  double Area();

  /*!
  \brief not implemented
  
  */
  double* BuildNormal(double* xi);
  
  /*!
  \brief not implemented
  
  */
  double Metric(double* xi, double g[], double G[][3]);

  /*!
  \brief not implemented
  
  */
  bool LocalCoordinatesOfNode(int lid, double* xi);
  
  //@}

protected:  

  // don't want = operator
  Segment_BiLinearQuad operator = (const Segment& old);

};

} // namespace MOERTEL

// << operator
std::ostream& operator << (std::ostream& os, const MOERTEL::Segment_BiLinearQuad& seg);


#endif // MOERTEL_SEGMENT_BILINEARQUAD_H

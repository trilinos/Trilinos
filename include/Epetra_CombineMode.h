/*
//@HEADER
// ************************************************************************
//
//               Epetra: Linear Algebra Services Package
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef EPETRA_COMBINEMODE_H
#define EPETRA_COMBINEMODE_H

#if defined(Epetra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Epetra package is deprecated"
#endif
#endif


/*! \file Epetra_CombineMode.h
    \brief Epetra_Combine Mode enumerable type
 */

/*! \enum Epetra_CombineMode
    If set to Add, components on the receiving processor will be added
    together.    If set to Zero, off-processor components will be ignored.
    If set to Insert, off-processor components will replace existing
    components on the receiving processor.  If set to InsertAdd, off-processor components
    will replace existing components, but multiple off-processor contributions will be added.
    If set to Average, off-processor components will be averaged with
    existing components on the receiving processor. (Recursive Binary Average)
    If set to AbsMax, magnitudes of off-processor components will be maxed
    with magnitudes of existing components of the receiving processor.
    { V = Supported by Epetra_Vector and Epetra_MultiVector,
      M = Supported by Epetra_CrsMatrix and Epetra_VbrMatrix }
*/

enum Epetra_CombineMode {Add,    /*!< Components on the receiving processor
                                     will be added together. (V,M) */
                        Zero,   /*!< Off-processor components will be
                                     ignored. (V,M) */
                        Insert, /*!< Off-processor components will
                                     be inserted into locations on
                                     receiving processor replacing existing values. (V,M) */
                        InsertAdd, /*!< Off-processor components will
                                     be inserted into locations on
                                     receiving processor replacing existing values. (V,M) */
                        Average,/*!< Off-processor components will be
                                     averaged with existing components
                                     on the receiving processor. (V) */
                        Epetra_Max,    /*!< Off-processor components will be
                                     maxed with existing components
                                     on the receiving processor. (V) */
                        Epetra_Min,    /*!< Off-processor components will be
                                     min'ed with existing components
                                     on the receiving processor. (V) */
                        AbsMax, /*!< Magnitudes of Off-processor components will be
                                     maxed with magnitudes of existing components
                                     on the receiving processor. (V) */
                        AbsMin, /*!< Magnitudes of Off-processor components will be
                                     min'ed with magnitudes of existing components
                                     on the receiving processor. (V) */
                        Epetra_AddLocalAlso    /*!< Like Add but local components are also added */
                        };

#endif // EPETRA_COMBINEMODE_H

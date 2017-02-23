// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef   __mySourceTerm_hpp__
#define   __mySourceTerm_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// C++
#include <string>

// Panzer
#include "Panzer_Dimension.hpp"
#include "Panzer_FieldLibrary.hpp"

// Phalanx
#include "Phalanx_Evaluator_Derived.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"
#include "Phalanx_FieldManager.hpp"

/**
 *  \brief Description.                                                          // JMG:  Fill this out.                          
 *                                                                               //                                               
 *  Detailed description.                                                        //                                               
 *                                                                               //                                               
 *  A source for the poisson equation that results in the solution               //                                               
 *  (with homogeneous dirichlet conditions) \f$sin(2\pi x)sin(2\pi y)\f$.        //                                               
 */
template<typename EvalT, typename Traits>
class MySourceTerm
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{
  public:

    /**
     * 	\brief Description.                                                      // JMG:  Fill this out.                          
     *                                                                           //                                               
     * 	Detailed description.                                                    //                                               
     *                                                                           //                                               
     * 	\param[?] name Description.                                              //                                               
     * 	\param[?] ir   Description.                                              //                                               
     */
    MySourceTerm(
      const std::string&             name,
      const panzer::IntegrationRule& ir);

    /**
     * 	\brief Description.                                                      // JMG:  Fill this out.                          
     *                                                                           //                                               
     * 	Detailed description.                                                    //                                               
     *                                                                           //                                               
     * 	\param[?] d  Description.                                                //                                               
     * 	\param[?] fm Description.                                                //                                               
     */
    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    /**
     * 	\brief Description.                                                      // JMG:  Fill this out.                          
     *                                                                           //                                               
     * 	Detailed description.                                                    //                                               
     *                                                                           //                                               
     * 	\param[?] d Description.                                                 //                                               
     */
    void
    evaluateFields(
      typename Traits::EvalData d);

  private:
    typedef typename EvalT::ScalarT ScalarT;

    /**
     * 	\brief Description.                                                      // JMG:  Fill this out.                          
     *                                                                           //                                               
     * 	Detailed description.                                                    //                                               
     */
    PHX::MDField<ScalarT, panzer::Cell, panzer::Point> result;

    /**
     * 	\brief Description.                                                      // JMG:  Fill this out.                          
     *                                                                           //                                               
     * 	Detailed description.                                                    //                                               
     */
    int irDegree_;

    /**
     * 	\brief Description.                                                      // JMG:  Fill this out.                          
     *                                                                           //                                               
     * 	Detailed description.                                                    //                                               
     */
    int irIndex_;
}; // end of class MySourceTerm

#endif // __mySourceTerm_hpp__

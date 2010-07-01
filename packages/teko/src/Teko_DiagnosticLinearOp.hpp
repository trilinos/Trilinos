/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER
*/

/** \file Teko_DiagnosticLinearOp.hpp
  * 
  * File that implements the a linear operator for printing out
  * diagnostics.
  */

#ifndef __Teko_DiagnosticLinearOp_hpp__
#define __Teko_DiagnosticLinearOp_hpp__

#include <iostream>

#include "Teko_Utilities.hpp"
#include "Teko_ImplicitLinearOp.hpp"

#include "Teuchos_Time.hpp"

namespace Teko {

/** \brief This linear operator prints diagnostics about operator
 *         application and creation times. It is useful for debugging
 *         problems and determining bottle necks.
 */
class DiagnosticLinearOp : public ImplicitLinearOp {
public:
   /** \brief This constructor explicitly takes the linear operator
     *        that needs to be wrapped and a string for output that describes
     *        the diagnostics.
     */
   DiagnosticLinearOp(const Teuchos::RCP<std::ostream> & ostrm,const ModifiableLinearOp & A,const std::string & diagnosticString);

   /** \brief This constructor explicitly takes the linear operator
     *        that needs to be wrapped and a string for output that describes
     *        the diagnostics.
     */
   DiagnosticLinearOp(const Teuchos::RCP<std::ostream> & ostrm,const LinearOp & fwdOp,const ModifiableLinearOp & A,const std::string & diagnosticString);

   /** \brief Destructor prints out timing information about this operator.
     */
   virtual ~DiagnosticLinearOp();

   //! \name Inherited methods from Thyra::LinearOpBase
   //@{

   /** @brief Range space of this operator */
   virtual VectorSpace range() const { return wrapOpA_->range(); }

   /** @brief Domain space of this operator */
   virtual VectorSpace domain() const { return wrapOpA_->range(); }

   /** @brief Perform a matrix vector multiply with this operator. 
     *
     * The <code>apply</code> function takes one vector as input 
     * and applies the inverse \f$ LDU \f$ decomposition. The result
     * is returned in \f$y\f$. If this operator is reprsented as \f$M\f$ then
     * \f$ y = \alpha M x + \beta y \f$ (ignoring conjugation!).
     *
     * @param[in]     x 
     * @param[in,out] y 
     * @param[in]     alpha (default=1)
     * @param[in]     beta  (default=0)
     */
   virtual void implicitApply(const MultiVector & x, MultiVector & y,
              const double alpha = 1.0, const double beta = 0.0) const;
   //@}

   virtual void describe(Teuchos::FancyOStream & out_arg,
                         const Teuchos::EVerbosityLevel verbLevel) const
   { wrapOpA_->describe(out_arg,verbLevel); }

   int numApplications() const { return timer_.numCalls(); }
   double totalTime() const { return timer_.totalElapsedTime(); }

   ModifiableLinearOp getModifiableOp() const
   { return wrapOpA_; } 

   LinearOp getLinearOp() const
   { return wrapOpA_; } 
 
   void setForwardOp(const Teko::LinearOp & lo)
   { fwdOp_ = lo; }

protected:
   // fundamental operators to use
   Teuchos::RCP<std::ostream> outputStream_;
   ModifiableLinearOp wrapOpA_;      ///< inverse of \f$ S \f$
   Teko::LinearOp fwdOp_;      ///< inverse of \f$ S \f$
   std::string diagString_;
 
   mutable Teuchos::Time timer_;

private:
   // hide me!
   DiagnosticLinearOp();
   DiagnosticLinearOp(const DiagnosticLinearOp &);
};

/** \brief Constructor method for building <code>DiagnosticLinearOp</code>.
  *
  * Constructor method for building <code>DiagnosticLinearOp</code>.
  *
  * \param[in] os     Output stream to print diagnostics to
  * \param[in] A      Operator to be wrapped
  * \param[in] label  String for outputing with diagnostics
  *
  * \returns A linear operator that wrapping A that will print diagnostics
  *          on descruction.
  * 
  * \relates LU2x2InverseOp
  */
inline ModifiableLinearOp createDiagnosticLinearOp(const Teuchos::RCP<std::ostream> & os,const ModifiableLinearOp & A,const std::string & label)
{
   return Teuchos::rcp(new DiagnosticLinearOp(os,A,label));
}

/** \brief Constructor method for building <code>DiagnosticLinearOp</code>.
  *
  * Constructor method for building <code>DiagnosticLinearOp</code>.
  *
  * \param[in] os     Output stream to print diagnostics to
  * \param[in] fwdOp  Forward operator to compute residual with
  * \param[in] A      Operator to be wrapped
  * \param[in] label  String for outputing with diagnostics
  *
  * \returns A linear operator that wrapping A that will print diagnostics
  *          on descruction.
  * 
  * \relates LU2x2InverseOp
  */
inline ModifiableLinearOp createDiagnosticLinearOp(const Teuchos::RCP<std::ostream> & os,const Teko::LinearOp & fwdOp,const ModifiableLinearOp & A,const std::string & label)
{
   return Teuchos::rcp(new DiagnosticLinearOp(os,fwdOp,A,label));
}

} // end namespace Teko

#endif	

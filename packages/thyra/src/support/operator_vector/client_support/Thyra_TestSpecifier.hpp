/* @HEADER@ */
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
/* @HEADER@ */


#ifndef THYRA_TESTSPECIFIER_HPP
#define THYRA_TESTSPECIFIER_HPP

namespace Thyra
{
  /** */
  template <class Scalar>
  class TestSpecifier
  {
  public:
    /** \brief Local typedef for promoted scalar magnitude */
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

    /** 
     * 
     */
    TestSpecifier(bool doTest, 
                  const ScalarMag& warningTol,
                  const ScalarMag& errorTol)
      : doTest_(doTest),
        warningTol_(warningTol),
        errorTol_(errorTol)
    {;}

    /** */
    bool doTest() const {return doTest_;}

    /** */
    const ScalarMag& warningTol() const {return warningTol_;}

    /** */
    const ScalarMag& errorTol() const {return errorTol_;}

  private:

    bool doTest_;

    ScalarMag warningTol_;

    ScalarMag errorTol_;
  };
}
#endif

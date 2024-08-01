// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EIGENVALUESORT_STRATEGIES_H
#define LOCA_EIGENVALUESORT_STRATEGIES_H

#include "Teuchos_RCP.hpp"

#include "NOX_Abstract_Group.H"   // for ReturnType

// Forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace LOCA {
  class GlobalData;
}

namespace LOCA {

  //! Collection of eigenvalue sorting strategies
  namespace EigenvalueSort {

    //! Abstract interface for eigenvalue sorting strategies
    /*!
     * AbstractStrategy defines an abstract interface for sorting eigenvalues.
     * It is used by LOCA::Eigensolver strategies to ensure the desired
     * eigenvalues are printed/saved.
     *
     * The interface defines two pure virtual methods, sort(), to sort
     * the eigenvalues, optionally returning a permutation vector determining
     * how th eigenvalues were sorted.  There is one version for real-only
     * eigenvalues and one version for complex.  Derived classes should
     * implement these methods for a particular sorting strategy.
     * Constructors for derived classes should be of the form:
     *
     * \code
     * class Derived : public AbstractStrategy {
     * public:
     *   Derived(
     *      const Teuchos::RCP<LOCA::GlobalData>& global_data,
     *      const Teuchos::RCP<Teuchos::ParameterList>& eigenParams);
     *   ...
     * };
     * \endcode
     *
     * where \c global_data is the LOCA global data object and \c eigenParams
     * is the eigensolver parameter list
     * (see LOCA::Eigensolver::AbstractStrategy).  In addition to any
     * parameters for the chosen sorting method, this list should contain
     * the parameter "Sorting Order" giving the name of the sorting
     * strategy.
     *
     * This class and its children follow the Strategy pattern as defined
     * in Erich Gamma, et al. "Design Patterns:  Elements of Reusable
     * Object-Oriented Software." Addison Wesley, Boston, MA, 1995.
     */
    class AbstractStrategy {

    public:

      //! Constructor
      AbstractStrategy() {}

      //! Destructor
      virtual ~AbstractStrategy() {}

      //! Sort real eigenvalues optionally returning a permutation vector
      /*!
       * \param n [in] Number of eigenvalues
       * \param evals [in/out] Array of length n containing the eigenvalues
       * to be sorted.
       * \param perm [out] Vector of length n to store the permutation
       * (optional).
       * \return Returns the status of the sorting routine
       */
      virtual NOX::Abstract::Group::ReturnType
      sort(int n, double* evals, std::vector<int>* perm = NULL) const = 0;

      //! Sort complex eigenvalues optionally returning a permutation vector
      /*!
       * \param n [in] Number of eigenvalues
       * \param r_evals [in/out] Array of length n containing the real part
       * of the eigenvalues to be sorted.
       * \param i_evals [in/out] Array of length n containing the imaginary
       * part of the eigenvalues to be sorted.
       * \param perm [out] Vector of length n to store the permutation
       * (optional).
       * \return Returns the status of the sorting routine
       */
      virtual NOX::Abstract::Group::ReturnType
      sort(int n, double* r_evals, double* i_evals = NULL,
       std::vector<int>* perm = NULL) const = 0;

    private:

      //! Private to prohibit copying
      AbstractStrategy(const AbstractStrategy&);

      //! Private to prohibit copying
      AbstractStrategy& operator = (const AbstractStrategy&);

    }; // Class AbstractStrategy

    //! Largest-magnitude sorting strategy
    /*!
     * Sorts eigenvalues in decreasing order according to their magnitude.
     * This method requires no parameters in the \c eigenParams argument
     * to the constructor
     */
    class LargestMagnitude : public AbstractStrategy {

    public:

      //! Constructor
      /*!
       * \param global_data [in] Global data object
       * \param eigenParams [in] Eigensolver parameters.
       */
      LargestMagnitude(
         const Teuchos::RCP<LOCA::GlobalData>& /* global_data */,
         const Teuchos::RCP<Teuchos::ParameterList>& /* eigenParams */) {}

      //! Destructor
      ~LargestMagnitude() {}

      //! Sort real eigenvalues
      virtual NOX::Abstract::Group::ReturnType
      sort(int n, double* evals, std::vector<int>* perm = NULL) const;

      //! Sort complex eigenvalues
      virtual NOX::Abstract::Group::ReturnType
      sort(int n, double* r_evals, double* i_evals,
       std::vector<int>* perm = NULL) const;

    private:

      //! Private to prohibit copying
      LargestMagnitude(const LargestMagnitude&);

      //! Private to prohibit copying
      LargestMagnitude& operator = (const LargestMagnitude&);

    }; // Class LargestMagnitude

    //! Smallest-magnitude sorting strategy
    /*!
     * Sorts eigenvalues in increasing order according to their magnitude.
     * This method requires no parameters in the \c eigenParams argument
     * to the constructor
     */
    class SmallestMagnitude : public AbstractStrategy {

    public:

      //! Constructor
      /*!
       * \param global_data [in] Global data object
       * \param eigenParams [in] Eigensolver parameters.
       */
      SmallestMagnitude(
         const Teuchos::RCP<LOCA::GlobalData>& /* global_data */,
         const Teuchos::RCP<Teuchos::ParameterList>& /* eigenParams */) {}

      //! Destructor
      ~SmallestMagnitude() {}

      //! Sort real eigenvalues
      virtual NOX::Abstract::Group::ReturnType
      sort(int n, double* evals, std::vector<int>* perm = NULL) const;

      //! Sort complex eigenvalues
      virtual NOX::Abstract::Group::ReturnType
      sort(int n, double* r_evals, double* i_evals,
       std::vector<int>* perm = NULL) const;

    private:

      //! Private to prohibit copying
      SmallestMagnitude(const SmallestMagnitude&);

      //! Private to prohibit copying
      SmallestMagnitude& operator = (const SmallestMagnitude&);

    }; // Class SmallestMagnitude

    //! Largest-real sorting strategy
    /*!
     * Sorts eigenvalues in decreasing order according to their real part.
     * This method requires no parameters in the \c eigenParams argument
     * to the constructor
     */
    class LargestReal : public AbstractStrategy {

    public:

      //! Constructor
      /*!
       * \param global_data [in] Global data object
       * \param eigenParams [in] Eigensolver parameters.
       */
      LargestReal(
         const Teuchos::RCP<LOCA::GlobalData>& /* global_data */,
         const Teuchos::RCP<Teuchos::ParameterList>& /* eigenParams */) {}

      //! Destructor
      ~LargestReal() {}

      //! Sort real eigenvalues
      virtual NOX::Abstract::Group::ReturnType
      sort(int n, double* evals, std::vector<int>* perm = NULL) const;

      //! Sort complex eigenvalues
      virtual NOX::Abstract::Group::ReturnType
      sort(int n, double* r_evals, double* i_evals,
       std::vector<int>* perm = NULL) const;

    private:

      //! Private to prohibit copying
      LargestReal(const LargestReal&);

      //! Private to prohibit copying
      LargestReal& operator = (const LargestReal&);

    }; // Class LargestReal

    //! Smallest-real sorting strategy
    /*!
     * Sorts eigenvalues in increasing order according to their real part.
     * This method requires no parameters in the \c eigenParams argument
     * to the constructor
     */
    class SmallestReal : public AbstractStrategy {

    public:

      //! Constructor
      /*!
       * \param global_data [in] Global data object
       * \param eigenParams [in] Eigensolver parameters.
       */
      SmallestReal(
         const Teuchos::RCP<LOCA::GlobalData>& /* global_data */,
         const Teuchos::RCP<Teuchos::ParameterList>& /* eigenParams */) {}

      //! Destructor
      ~SmallestReal() {}

      //! Sort real eigenvalues
      virtual NOX::Abstract::Group::ReturnType
      sort(int n, double* evals, std::vector<int>* perm = NULL) const;

      //! Sort complex eigenvalues
      virtual NOX::Abstract::Group::ReturnType
      sort(int n, double* r_evals, double* i_evals,
       std::vector<int>* perm = NULL) const;

    private:

      //! Private to prohibit copying
      SmallestReal(const SmallestReal&);

      //! Private to prohibit copying
      SmallestReal& operator = (const SmallestReal&);

    }; // Class SmallestReal

    //! Largest-imaginary sorting strategy
    /*!
     * Sorts eigenvalues in decreasing order according to their imaginary part.
     * This method requires no parameters in the \c eigenParams argument
     * to the constructor
     */
    class LargestImaginary : public AbstractStrategy {

    public:

      //! Constructor
      /*!
       * \param global_data [in] Global data object
       * \param eigenParams [in] Eigensolver parameters.
       */
      LargestImaginary(
         const Teuchos::RCP<LOCA::GlobalData>& /* global_data */,
         const Teuchos::RCP<Teuchos::ParameterList>& /* eigenParams */) {}

      //! Destructor
      ~LargestImaginary() {}

      //! Sort real eigenvalues
      virtual NOX::Abstract::Group::ReturnType
      sort(int n, double* evals, std::vector<int>* perm = NULL) const;

      //! Sort complex eigenvalues
      virtual NOX::Abstract::Group::ReturnType
      sort(int n, double* r_evals, double* i_evals,
       std::vector<int>* perm = NULL) const;

    private:

      //! Private to prohibit copying
      LargestImaginary(const LargestImaginary&);

      //! Private to prohibit copying
      LargestImaginary& operator = (const LargestImaginary&);

    }; // Class LargestImaginary

    //! Smallest-imaginary sorting strategy
    /*!
     * Sorts eigenvalues in increasing order according to their imaginary part.
     * This method requires no parameters in the \c eigenParams argument
     * to the constructor
     */
    class SmallestImaginary : public AbstractStrategy {

    public:

      //! Constructor
      /*!
       * \param global_data [in] Global data object
       * \param eigenParams [in] Eigensolver parameters.
       */
      SmallestImaginary(
         const Teuchos::RCP<LOCA::GlobalData>& /* global_data */,
         const Teuchos::RCP<Teuchos::ParameterList>& /* eigenParams */) {}

      //! Destructor
      ~SmallestImaginary() {}

      //! Sort real eigenvalues
      virtual NOX::Abstract::Group::ReturnType
      sort(int n, double* evals, std::vector<int>* perm = NULL) const;

      //! Sort complex eigenvalues
      virtual NOX::Abstract::Group::ReturnType
      sort(int n, double* r_evals, double* i_evals,
       std::vector<int>* perm = NULL) const;

    private:

      //! Private to prohibit copying
      SmallestImaginary(const SmallestImaginary&);

      //! Private to prohibit copying
      SmallestImaginary& operator = (const SmallestImaginary&);

    }; // Class SmallestImaginary

    //! Largest-Real Cayley sorting strategy
    /*!
     * Sorts eigenvalues in decreasing order according to the real part of
     * their inverse-Cayley transformation.
     * This method references the "CayleyPole" and "CayleyZero" parameters
     * in the eigensolver parameter list.
     */
    class LargestRealInverseCayley : public AbstractStrategy {

    public:

      //! Constructor
      /*!
       * \param global_data [in] Global data object
       * \param eigenParams [in] Eigensolver parameters.
       */
      LargestRealInverseCayley(
         const Teuchos::RCP<LOCA::GlobalData>& global_data,
         const Teuchos::RCP<Teuchos::ParameterList>& eigenParams);

      //! Destructor
      ~LargestRealInverseCayley() {}

      //! Sort real eigenvalues
      virtual NOX::Abstract::Group::ReturnType
      sort(int n, double* evals, std::vector<int>* perm = NULL) const;

      //! Sort complex eigenvalues
      virtual NOX::Abstract::Group::ReturnType
      sort(int n, double* r_evals, double* i_evals,
       std::vector<int>* perm = NULL) const;

    protected:

      //! Compute the part of the inverse of the Cayley transformation
      double realLambda(double er, double ei) const;

    private:

      //! Private to prohibit copying
      LargestRealInverseCayley(const LargestRealInverseCayley&);

      //! Private to prohibit copying
      LargestRealInverseCayley& operator = (const LargestRealInverseCayley&);

    protected:

      //! Cayley pole
      double sigma;

      //! Cayley zero
      double mu;

    }; // Class SmallestImaginary

  } // Namespace EigenvalueSort

} // Namespace LOCA

#endif // LOCA_EIGENVALUESORT_STRATEGIES_H

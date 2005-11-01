/*! 
 * \file AnasaziProjectedEigenproblem.hpp
 *
 * \brief General class to construct and solve the projected, dense eigenproblem.
 */

#ifndef ANASAZI_PROJECTED_EIGENPROBLEM_HPP
#define ANASAZI_PROJECTED_EIGENPROBLEM_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziMultiVec.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Anasazi {

/*!
\class ProjectedEigenproblem

\brief General approach to construct and solve the projected eigenproblem.

ProjectedEigenproblem is a pure virtual class that defines the interface for
every concrete ProjectedEigenproblem. Currently, one concrete implementation
is given by trait ProjectedEigenproblem<int, double, MV>.

This class defines a generic approach that can be used to construct and
solve the small, dense eigenproblem that occurs in several projection
methods, like Davison and Jacobi-Davidson. Extracting an eigenpair from a subspace is an important problem in tis own,
considering the vast number of methods that rely on this ability.
Unfortunately, extracting <I>optimal</I> eigenpairs from a subspace may
beocome quite an involved task, mainly depending on the requirements the
eigenpairs have to meet. 

The assumptions are as follows.
The (distributed, sparse) eigenproblem to be solved is
\f[
A V = \Lambda B V,
\f]
where \f$A\f$ and \f$B\f$ are two operators, \f$\Lambda\f$ is a diagonal
matrix containing the eigenvalues and \f$V\f$ is a matrix containing the
eigenvectors. A method that constructs a search space \f$\mathcal{U}\f$ is adopted; this method requires the solution of the reduced problem
\f[
\mathcal{U}^H A \mathcal{U} Z = \Theta \mathcal{U}^H B \mathcal{U} Z,
\f]
where \f$Z\f$ and \f$\Theta\f$ are the eigenvalues and eigenvectors of the reduced problem, respectively. This reduced problem can be written as
\f[
A_\mathcal{U} Z = \Theta B_\mathcal{U}
\f]
and it represents the original problem in the subspace defined by the search
space. This class
constructs the matrices \f$A_\mathcal{U}\f$ and \f$B_\mathcal{U}\f$, and solves
the above equation using optimized methods, typically LAPACK routines. The user has to provide the the three spaces
\f$\mathcal{U}\f$, \f$\mathcal{A U}\f$ and \f$\mathcal{B U}\f$ are defined by the algorithm. 

The general usage is as follows. The class is templated with an OrdinalType
(for example, \c int), a ScalarType (for example, \c double), and an
Anasazi::MultiVec<OrdinalType, ScalarType>. 

First, one has to instantiate an object,
\verbatim
ProjectedEigenproblem<int, ScalarType, MV> PE(MaxSize);
\endverbatim
Then, pointers to already allocated Teuchos::SerialDenseMatrix's for \f$\Theta\f$ and \f$Z\f$ must be passed,
\verbatim
PE.SetTheta(&theta);
PE.SetZ(&Z);
\endverbatim
These Teuchos::SerialDenseMatrix objects will contain the computed eigenpairs.

The matrix type can be specified using
\verbatim
PE.SetMatrixType("Symmetric");
\endverbatim
Other available types are: \c "Hermitian" or \c "General".

Every time a vector (or a multi-vector) is added to the search spaces,
one has to update the ProjectedEigenproblem components by using
\verbatim
PE.Add(IncrementU, IncrementA, IncrementB);
\endverbatim
At this point, one can extract eigenpair approximations from the space by
simply calling
\verbatim
PE.Extract()
\endverbatim
The eigenpair is returned in the \c theta and \c Z matrix defined by the user.

\warning Still very incomplete...

\date Last updated on 26-Oct-05.

\author Oscar Chinellato (ETHZ/ICOS) and Marzio Sala (ETHZ/COLAB)
*/

template <class OrdinalType, class ScalarType, class MV>
class ProjectedEigenproblem 
{
public:
  
  //! Constructor for a given matrix type and maximum search space size.
  /*!
   * 
   * \param MatrixType - (In) can assume one of the following values:
   *                     "Symmetric", "Hermitian", "General".
   *
   * \param maxSize - (In) Maximum size of the projected eigenproblem.
   */
  ProjectedEigenproblem(const string& MatrixType, int maxSize){};

  //! Destructor.
  virtual ~ProjectedEigenproblem(){};
    
  //! Adds given vectors to the search spaces. 
  /*! Adds input vectors \c U, \c AU, and \c BU to the three search space.
   * The class internally stores shallow copies.
   */
  void Add(const MV &U, const MV &AU, const MV &BU); 

  //! Extracts the eigenpairs from the reduced system.
  void Extract(); 

  void Rotate(const Teuchos::SerialDenseMatrix<OrdinalType,ScalarType> &Q); 

  //! \fixme Move in the constuctor??
  void SetTheta(const std::vector<ScalarType> theta); 

  //! \fixme Move in the constuctor??
  void SetZ(const Teuchos::SerialDenseMatrix<OrdinalType,ScalarType> &Z); 
};

/*!
\brief Specialization for <int, double> of ProjectedEigenproblem;
see the documentation for Anasazi::ProjectedEigenproblem for more details.

\author Oscar Chinellato (ETHZ/ICOS) and Marzio Sala (ETHZ/COLAB)

\date Last updated on 01-Nov-05

*/
template <class MV>
class ProjectedEigenproblem<int, double, MV> { 
public:

  enum MatrixType {GENERAL=0, SYMMETRIC, HERMITIAN};

  //! Constructor for a given matrix type and maximum search space size.
  /*! \fixme: NOW ONLY SYMMETRIC PROBLEMS!
   */
  ProjectedEigenproblem(const string& MatrixType, int maxSize)
  {
    _A.shape(maxSize, maxSize);
    _Awork.shape(maxSize, maxSize);
    _work.resize(3*maxSize);
    _actualSize = 0;
    _maxSize = maxSize;
    if (MatrixType == "General")
      _MT = GENERAL;
    else if (MatrixType == "Symmetric")
      _MT = SYMMETRIC;
    else if (MatrixType == "Hermitian")
      _MT = HERMITIAN;
    else
      throw("Error, value of MatrixType is not correct");
  }

  //! Destructor.
  virtual ~ProjectedEigenproblem(){}

  // Adds the input vectors to the internally stored search space, using shallow copies.
  void Add(const MV &U, const MV &AU, const MV &BU)
  {
    int i, j;
    int oldSize;
    std::vector<int>    list(1);
    std::vector<double> Aij(1);
	
    // Retrieve the pointers of the new entries in U, AU, BU
    oldSize = _actualSize;
    for(i=0; i < U.GetNumberVecs(); ++i){
      list[0] = i;

      _U[_actualSize]  = U.CloneView(list);
      _AU[_actualSize] = AU.CloneView(list);
      _BU[_actualSize] = BU.CloneView(list);
      _actualSize++;
    }


    // Update _A (i.e. Atilde)
    for(i=0; i < _actualSize; ++i){
      for (j=0; j < U.GetNumberVecs(); ++j){
	_U[i]->MvDot((*_AU[oldSize+j]), &Aij);
	_A(i, oldSize + j) = Aij[0];
      }
    }
  }

  //! Extracts the eigenpairs from the current projected problem using LAPACK.
  void Extract()
  {
    (*_Z) = _A;    
    _lapack.SYEV('V', 'U', _actualSize, _Z->values(), _maxSize, &((*_theta)[0]), &(_work[0]), 3*_maxSize, &_info);
  }

  void Rotate(const Teuchos::SerialDenseMatrix<int,double> &Q)      // During a restart, this can be simplified
  {
    int i, j, ret;

    // Mirror the matrix _A                     
    for(i=1; i<_actualSize; ++i){
      for(j=0; j<i; ++j){
	_A(i,j) = _A(j,i);
      }
    }

    _A.reshape(Q.numRows(), Q.numRows());
    _Awork.shape(_A.numRows(), Q.numCols());
    ret = _Awork.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, _A, Q, 0.0); assert(ret == 0);
    _A.shape(Q.numCols(),Q.numCols());
    ret = _A.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, Q, _Awork, 0.0); assert(ret == 0);
    _A.reshape(_maxSize, _maxSize);
    _actualSize = Q.numCols();    
  }

  void SetTheta(std::vector<double> *theta)  
  {
    _theta = theta;
  }
 
  void SetZ(Teuchos::SerialDenseMatrix<int,double> *Z)
  {
    _Z = Z;
  }

  void Print()
  {
    // FIXME: more details in a more comprehensible form!
    cout << _A << endl;
  }

private:
  Teuchos::SerialDenseMatrix<int, double> _A, _Awork;
  std::vector<double> _work;
  int _info;
  int _actualSize;
  std::vector<double> *_theta;
  Teuchos::SerialDenseMatrix<int,double> *_Z;   
  int _maxSize;
  MatrixType  _MT;
  std::map<int, const MV* > _U;
  std::map<int, const MV* > _AU;
  std::map<int, const MV* > _BU;
  Teuchos::LAPACK<int, double> _lapack;
};

} // namespace Anasazi

#endif // ANASAZI_PROJECTED_EIGENPROBLEM_HPP

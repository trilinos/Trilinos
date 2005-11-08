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

\brief Pure virtual class to manage the construction and the solution of
projected eigenproblems.

ProjectedEigenproblem is a pure virtual class that defines the interface for
every concrete ProjectedEigenproblem. Currently, one concrete implementation
is given by trait ProjectedEigenproblem<int, double, MV> (with MV<double>).

ProjectedEigenproblem encapsulates all it is needed to construct and solve the small, dense eigenproblem that occers in several projection
methods, like Davison and Jacobi-Davidson. 

The assumptions are as follows.
The (distributed, sparse) eigenproblem to be solved is
\f[
A X = \Lambda B X, \vspace{2cm} (1)
\f]
where \f$A\f$ and \f$B\f$ are two operators, \f$\Lambda\f$ is a diagonal
matrix containing the eigenvalues and \f$X\f$ is a matrix containing the
eigenvectors. A method that constructs a search space \f$\mathcal{U}\f$ is
adopted. Given \f$\mathcal{U}\f$, an approximate solution of (1) can be computed by solving the reduced eigenproblem
\f[
\begin{tabular}{rcl}
$\mathcal{U}^H A \mathcal{U} \, Z$ & = & $\Theta \, \mathcal{U}^H B \mathcal{U} \, Z$ \\
$A_\mathcal{U} \, Z$ & = & $\Theta \, B_\mathcal{U}$ \\
\end{tabular}
\f]
where \f$Z\f$ and \f$\Theta\f$ are the eigenvalues and eigenvectors of the
reduced problem, respectively. This eigenproblem is a projection of the
original problem in the subspace \f$\mathcal{U}\f$, and it is small, dense,
and serial. Class ProjectedEigenproblem offers a set of capabilities to efficiently solve the projected problem, typically resorting to LAPACK routines.

The general usage is as follows. The class is templated with an OrdinalType
(for example, \c int), a ScalarType (for example, \c double), and an
Anasazi::MultiVec<ScalarType>. The search space
\f$\mathcal{U}\f$ and the auxiliary spaces
\f$\mathcal{U}_A = A \, \mathcal{U}\f$ and 
\f$\mathcal{U}_B = B \, \mathcal{U}\f$ must be provided by the user.

First, one has to instantiate an object,
\verbatim
string MatrixType = "Symmetric"; // type of matrix
OrdinalType MaxSize = 16;        // maximum number of eigenvalues to compute
ProjectedEigenproblem<OrdinalType, ScalarType, MV> PE(MatrixType, MaxSize);
\endverbatim
Other available matrix types are: \c "Hermitian" or \c "General".

Then, pointers to already allocated Teuchos::SerialDenseMatrix's for \f$\Theta\f$ and \f$Z\f$ must be passed,
\verbatim
PE.SetTheta(&theta);
PE.SetZ(&Z);
\endverbatim
These Teuchos::SerialDenseMatrix objects will contain the computed eigenpairs.

Every time a vector (or a multi-vector) is added to the search spaces,
one has to update the ProjectedEigenproblem components by using
\verbatim
PE.Add(IncrementU, IncrementAU, IncrementBU);
\endverbatim
Note that only the additional vectors are passed, <I>not</I> the entire spaces.
At this point, one can extract eigenpair approximations from the space by
simply calling
\verbatim
PE.Extract()
\endverbatim
The eigenpair is returned in the \c theta and \c Z matrix defined by the user.

\note When the search space \f$\mathcal{U}\f$ is rotated (using a dense matrix
<TT>Q</TT>), one should propagate this rotation to the projected eigenproblem
as well, using method <TT>Rotate(Q)</TT>.


\date Last updated on 01-Nov-05.

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

  //! Applies a given rotation \c Q to the projected eigenproblem.
  void Rotate(const Teuchos::SerialDenseMatrix<OrdinalType,ScalarType> &Q); 

  //! Sets the vector that will contain the computed eigenvalues.
  void SetTheta(const std::vector<ScalarType>* theta); 

  //! Sets the matrix that will contain the computed eigenvectors.
  void SetZ(const Teuchos::SerialDenseMatrix<OrdinalType,ScalarType>* Z); 
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

      _U[_actualSize]  = MVT::CloneView( U, list);
      _AU[_actualSize] = MVT::CloneView( AU, list);
      _BU[_actualSize] = MVT::CloneView( BU, list);
      _actualSize++;
    }


    // Update _A (i.e. Atilde)
    for(i=0; i < _actualSize; ++i){
      for (j=0; j < U.GetNumberVecs(); ++j){
	MVT::MvDot( *_U[i], *_AU[oldSize+j], &Aij);
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

  //! Sets the vector that will contain the computed eigenvalues.
  void SetTheta(std::vector<double> *theta)  
  {
    _theta = theta;
  }
 
  //! Sets the dense matrix that will contain the computed eigenvectors.
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
  std::map<int, Teuchos::RefCountPtr<const MV> > _U;
  std::map<int, Teuchos::RefCountPtr<const MV> > _AU;
  std::map<int, Teuchos::RefCountPtr<const MV> > _BU;
  Teuchos::LAPACK<int, double> _lapack;

  // typedef for Anasazi::MultiVecTraits class, which should be used for any
  // interaction with any multivectors.
  typedef Anasazi::MultiVecTraits<double,MV> MVT;
};

#if 0
template <class MV>
class ProjectedEigenproblem<int, complex<double>, MV>
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
  ProjectedEigenproblem(const string& MatrixType, int maxSize)
  {}

  //! Destructor.
  virtual ~ProjectedEigenproblem()
  {}
    
  //! Adds given vectors to the search spaces. 
  /*! Adds input vectors \c U, \c AU, and \c BU to the three search space.
   * The class internally stores shallow copies.
   */
  void Add(const MV &U, const MV &AU, const MV &BU)
  {}

  //! Extracts the eigenpairs from the reduced system.
  void Extract()
  {}

  //! Applies a given rotation \c Q to the projected eigenproblem.
  void Rotate(const Teuchos::SerialDenseMatrix<int,complex<double> > &Q)
  {}

  //! Sets the vector that will contain the computed eigenvalues.
  void SetTheta(const std::vector<complex<double> > *theta)
  {}

  //! Sets the matrix that will contain the computed eigenvectors.
  void SetZ(const Teuchos::SerialDenseMatrix<int,complex<double> > *Z)
  {}
  //
  // typedef for Anasazi::MultiVecTraits class, which should be used for any
  // interaction with any multivectors.
  typedef Anasazi::MultiVecTraits<complex<double>,MV> MVT;
};
#endif

/*!
\brief Specialization for <int, complex<double> > of ProjectedEigenproblem;
see the documentation for Anasazi::ProjectedEigenproblem for more details.

\author Oscar Chinellato (ETHZ/ICOS) and Marzio Sala (ETHZ/COLAB)

\date Last updated on 01-Nov-05

*/
template <class MV>
class ProjectedEigenproblem<int, complex<double>, MV> { 
public:

  enum MatrixType {GENERAL=0, SYMMETRIC, HERMITIAN};

  //! Constructor for a given matrix type and maximum search space size.
  /*! \fixme: NOW ONLY SYMMETRIC PROBLEMS!
   */
  ProjectedEigenproblem(const string& MatrixType, int maxSize)
  {
    _A.shape(maxSize, maxSize);
    _A2.shape(maxSize, maxSize);
    _Awork.shape(maxSize, maxSize);
    _work.resize(3*maxSize);
    _rwork.resize(2 * _maxSize);

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

    if (MatrixType != "Symmetric")
    {
      cerr << "Not implemented yet..." << endl;
      exit(0);
    }
  }

  //! Destructor.
  virtual ~ProjectedEigenproblem(){}

  // Adds the input vectors to the internally stored search space, using shallow copies.
  void Add(const MV &U, const MV &AU, const MV &BU)
  {
    int i, j;
    int oldSize;
    std::vector<int>    list(1);
    std::vector<complex<double> > Aij(1);
	
    // Retrieve the pointers of the new entries in U, AU, BU
    oldSize = _actualSize;
    for(i=0; i < U.GetNumberVecs(); ++i){
      list[0] = i;

      _U[_actualSize]  = MVT::CloneView( U, list);
      _AU[_actualSize] = MVT::CloneView( AU, list);
      _BU[_actualSize] = MVT::CloneView( BU, list);
      _actualSize++;
    }


    // Update _A (i.e. Atilde)
    for(i=0; i < _actualSize; ++i){
      for (j=0; j < U.GetNumberVecs(); ++j){
	MVT::MvDot( *_U[i], *_AU[oldSize+j], &Aij);
	_A(i, oldSize + j) = Aij[0];
      }
    }
  }

  //! Extracts the eigenpairs from the current projected problem using LAPACK.
  void Extract()
  {
    // Mirror the matrix _A                     
    for(int i=1; i<_actualSize; ++i){
      for(int j=0; j<i; ++j){
	_A(i,j) = _A(j,i);
      }
    }

    _A2 = _A;
    _lapack.GEEV('N', 'V', _actualSize, _A2.values(), _maxSize, 
                 &((*_theta)[0]), 0, 1, _Z->values(), _maxSize, 
                 &(_work[0]), 3*_maxSize, &(_rwork[0]), &_info);
    //for (int i = 0 ; i < _actualSize ; ++i)

  }

  void Rotate(const Teuchos::SerialDenseMatrix<int,complex<double> > &Q)      // During a restart, this can be simplified
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

  //! Sets the vector that will contain the computed eigenvalues.
  void SetTheta(std::vector<complex<double> > *theta)  
  {
    _theta = theta;
  }
 
  //! Sets the dense matrix that will contain the computed eigenvectors.
  void SetZ(Teuchos::SerialDenseMatrix<int,complex<double> > *Z)
  {
    _Z = Z;
  }

  void Print()
  {

    // FIXME: more details in a more comprehensible form!
    cout << _A << endl;
  }

private:
  Teuchos::SerialDenseMatrix<int, complex<double> > _A, _A2, _Awork;
  std::vector<complex<double> > _work;
  std::vector<double> _rwork; // this is double
  int _info;
  int _actualSize;
  std::vector<complex<double> > *_theta;
  Teuchos::SerialDenseMatrix<int,complex<double> > *_Z;   
  int _maxSize;
  MatrixType  _MT;
  std::map<int, Teuchos::RefCountPtr<const MV> > _U;
  std::map<int, Teuchos::RefCountPtr<const MV> > _AU;
  std::map<int, Teuchos::RefCountPtr<const MV> > _BU;
  Teuchos::LAPACK<int, complex<double> > _lapack;

  // typedef for Anasazi::MultiVecTraits class, which should be used for any
  // interaction with any multivectors.
  typedef Anasazi::MultiVecTraits<complex<double> ,MV> MVT;
};
} // namespace Anasazi

#endif // ANASAZI_PROJECTED_EIGENPROBLEM_HPP

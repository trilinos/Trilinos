#ifndef MUEMEX_TYPES_H
#define MUEMEX_TYPES_H

#include "Teuchos_ParameterList.hpp"
#include <string>
#include <complex>
#include <stdexcept>
#include "MueLu.hpp"
#include "MueLu_EpetraOperator.hpp"
#include "MueLu_TpetraOperator.hpp"
#include "MueLu_Hierarchy_decl.hpp"
#include "MueLu_CreateEpetraPreconditioner.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Tpetra_CrsMatrix_decl.hpp"
#include "Xpetra_EpetraCrsMatrix.hpp"

#if !defined(HAVE_MUELU_MATLAB) || !defined(HAVE_MUELU_EPETRA) || !defined(HAVE_MUELU_TPETRA)
#error "Muemex types require MATLAB, Epetra and Tpetra."
#else
#include "mex.h"

//Useful global typedefs for MueMex

enum MUEMEX_TYPE
{
	INT,
	DOUBLE,
	STRING,
	COMPLEX,
	XPETRA_ORDINAL_VECTOR,
	TPETRA_MULTIVECTOR_DOUBLE,
	TPETRA_MULTIVECTOR_COMPLEX,
	TPETRA_MATRIX_DOUBLE,
	TPETRA_MATRIX_COMPLEX,
	XPETRA_MATRIX_DOUBLE,
	XPETRA_MATRIX_COMPLEX,
	XPETRA_MULTIVECTOR_DOUBLE,
	XPETRA_MULTIVECTOR_COMPLEX,
	EPETRA_CRSMATRIX,
	EPETRA_MULTIVECTOR
};

typedef Tpetra::Vector<>::node_type mm_node_t;
typedef Tpetra::Vector<>::local_ordinal_type mm_LocalOrd;
typedef Tpetra::Vector<>::global_ordinal_type mm_GlobalOrd;
typedef Tpetra::Map<> muemex_map_type;
typedef Tpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_CrsMatrix_double;
typedef std::complex<double> complex_t;
typedef Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_CrsMatrix_complex;
typedef MueLu::TpetraOperator<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_operator_real;
typedef MueLu::TpetraOperator<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_operator_complex;
typedef Xpetra::Vector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_ordinal_vector;
typedef Xpetra::Matrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_Matrix_double;
typedef Xpetra::Matrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_Matrix_complex;
typedef Xpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_MultiVector_double;
typedef Xpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_MultiVector_complex;
typedef MueLu::Hierarchy<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Hierarchy_double;
typedef MueLu::Hierarchy<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Hierarchy_complex;

/* Static utility functions */
template<typename Scalar>
Teuchos::RCP<Tpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> loadTpetraMV(const mxArray* mxa);
//create a sparse array in Matlab
template<typename Scalar>
mxArray* createMatlabSparse(int numRows, int numCols, int nnz);
//create an ordinal (int32) vector in Matlab
mxArray* createMatlabLOVector(Teuchos::RCP<Xpetra_ordinal_vector> vec);
//copy a sparse Xpetra matrix (double or complex) to Matlab
template<typename Scalar>
mxArray* saveMatrixToMatlab(Teuchos::RCP<Xpetra::Matrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mat);
template<typename Scalar>
mxArray* createMatlabMultiVector(int numRows, int numCols);
template<typename Scalar>
mxArray* saveTpetraMV(Teuchos::RCP<Tpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mv);
template<typename Scalar>
mxArray* saveMultiVectorToMatlab(Teuchos::RCP<Xpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mv);
template<typename Scalar>
Teuchos::RCP<Xpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> loadXpetraMV(const mxArray* mxa);
template<typename Scalar>
void fillMatlabArray(Scalar* array, const mxArray* mxa, int n);
//set up Tpetra matrix from MATLAB array
template<typename Scalar>
Teuchos::RCP<Tpetra::CrsMatrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> tpetraLoadMatrix(const mxArray* mxa);
//same as above but for Xpetra
template<typename Scalar>
Teuchos::RCP<Xpetra::Matrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> xpetraLoadMatrix(const mxArray* mxa);
//Get an Epetra_MultiVector from MATLAB array
Teuchos::RCP<Epetra_MultiVector> loadEpetraMV(const mxArray* mxa);
//Save an Epetra MV to MATLAB array
mxArray* saveEpetraMV(Teuchos::RCP<Epetra_MultiVector> mv);
//Load an Epetra matrix from MATLAB array
Teuchos::RCP<Epetra_CrsMatrix> epetraLoadMatrix(const mxArray* mxa);
//Get an int from a MATLAB double or int input
int parseInt(const mxArray* mxa);
//Save an Epetra matrix to MATLAB array
mxArray* saveEpetraMatrix(Teuchos::RCP<Epetra_CrsMatrix> mat);
//Load an ordinal vector
Teuchos::RCP<Xpetra_ordinal_vector> loadLOVector(const mxArray* mxa);

class MuemexArg
{
	public:
		MuemexArg(MUEMEX_TYPE type);
		MUEMEX_TYPE type;
};

template<typename T>
class MuemexData : public MuemexArg
{
	public:
		MuemexData(T& data, MUEMEX_TYPE type);	//Construct from pre-existing data, to pass to MATLAB.
		MuemexData(const mxArray* mxa); //Construct from MATLAB array, to get from MATLAB.
		~MuemexData();
		mxArray* convertToMatlab(); //Create a MATLAB object and copy this data to it
		T& getData();				//Set and get methods
		void setData(T& data);
	private:
		T data;
};

#endif //HAVE_MUELU_MATLAB error handler
#endif //MUEMEX_TYPES_H guard

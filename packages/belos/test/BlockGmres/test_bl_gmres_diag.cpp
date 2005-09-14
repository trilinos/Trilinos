

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblemManager.hpp>
#include <BelosOutputManager.hpp>
#include <BelosStatusTestMaxIters.hpp>
#include <BelosStatusTestMaxRestarts.hpp>
#include <BelosStatusTestResNorm.hpp>
#include <BelosStatusTestCombo.hpp>
#include <BelosEpetraOperator.hpp>
#include <BelosPetraInterface.hpp>
#include <BelosBlockGmres.hpp>

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#include <Teuchos_Time.hpp>

#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#include <mpi.h>
#else
#include <Epetra_SerialComm.h>
#endif

using std::vector;
using namespace Belos;

//************************************************************************************************

class Vector_Operator
{
  public:

    Vector_Operator(unsigned m, unsigned n) : m(m), n(n) { y.resize(m, 0.0); };

    virtual ~Vector_Operator() {};

    virtual vector<double> & operator () (const vector<double> &x) = 0;

    int size (int dim) const { return (dim == 1) ? m : n; };

  protected:

    unsigned m, n;        // an (m x n) operator 
    vector<double> y;
};

//************************************************************************************************

class Diagonal_Operator : public Vector_Operator
{
  public:

    Diagonal_Operator(unsigned n, double v) : Vector_Operator(n, n), v(v) { };

    ~Diagonal_Operator() { };

    vector<double> & operator () (const vector<double> &x)
    {
        for (unsigned i=0; i < m; ++i) y[i] = v*x[i];  // NOTE: square operator!

        return y;
    };

  private:

    double v;
};

//************************************************************************************************

class Diagonal_Operator_2 : public Vector_Operator
{
  public:

    Diagonal_Operator_2(unsigned n, double v) : Vector_Operator(n, n), v(v) { };

    ~Diagonal_Operator_2() { };

    vector<double> & operator () (const vector<double> &x)
    {
        for (unsigned i=0; i < m; ++i) y[i] = (i+1)*v*x[i];  // NOTE: square operator!

        return y;
    };

  private:

    double v;
};

//************************************************************************************************

class Composed_Operator : public Vector_Operator
{
  public:

    Composed_Operator(unsigned n, Vector_Operator *pA, Vector_Operator *pB);

    virtual ~Composed_Operator() {};

    virtual vector<double> & operator () (const vector<double> &x);

  private:

    Vector_Operator *pA; 
    Vector_Operator *pB; 
};

Composed_Operator::Composed_Operator(unsigned n, 
                                     Vector_Operator *pA, 
                                     Vector_Operator *pB) 
    : Vector_Operator(n, n), pA(pA), pB(pB) 
{
}

vector<double> & Composed_Operator::operator () (const vector<double> &x)
{
    Vector_Operator &A(*pA);
    Vector_Operator &B(*pB);

    y = A(B(x));

    return y;
}

//************************************************************************************************

class Trilinos_Interface : public Epetra_Operator
{
  public:

    Trilinos_Interface(Vector_Operator   *pA,
                       const Epetra_Comm *pComm,
                       const Epetra_Map  *pMap)
        : pA(pA), pComm(pComm), pMap(pMap) 
    {
    };

    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
    {
        return(Apply(X,Y));  // No inverse
    };

    virtual ~Trilinos_Interface() {};

    const char * Label() const {return("Trilinos_Interface, an Epetra_Operator implementation");}; 
  
    bool UseTranspose() const {return(use_transpose);};      // always set to false
      
    int SetUseTranspose(bool UseTranspose_) { use_transpose = false; return(-1); };

    bool HasNormInf() const {return(false);};                // cannot return inf-norm

    double NormInf() const {return(0.0);};                   

    virtual const Epetra_Comm & Comm() const {return *pComm; }
      
    virtual const Epetra_Map & OperatorDomainMap() const {return *pMap; };
      
    virtual const Epetra_Map & OperatorRangeMap() const {return *pMap; };

  private:

    Vector_Operator   *pA;
    const Epetra_Comm *pComm;
    const Epetra_Map  *pMap;

    bool use_transpose;
};

int Trilinos_Interface::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
    vector<double> x(X.MyLength());
    X.ExtractCopy(&x[0], X.MyLength());

    Vector_Operator &A(*pA);
    x = A(x);

    for (int i = 0; i < Y.MyLength(); ++i) 
        Y.ReplaceMyValue(i, 0, x[i]);

    return(0);
}

//************************************************************************************************

class Iterative_Inverse_Operator : public Vector_Operator
{
  public:

    Iterative_Inverse_Operator(unsigned n, 
                               Vector_Operator *pA, 
                               bool print);              

    virtual ~Iterative_Inverse_Operator() 
    {
       delete pE;    delete pComm; delete pMap;
       delete pPE;   delete pPX;   delete pPB;   delete pProb;
       delete test1; delete test2; delete test3; delete test4;
       delete pTest; delete pOM;   delete pBelos;
    }

    virtual vector<double> & operator () (const vector<double> &b);

  private:

    Vector_Operator *pA;       // operator which will be inverted 
                               // supplies a matrix vector multiply
    const bool print;

    Epetra_Comm  *pComm;
    Epetra_Map      *pMap;
    Epetra_Operator *pE;

    PetraMat<double> *pPE;   
    PetraVec<double> *pPX;   
    PetraVec<double> *pPB;
    StatusTestMaxIters<double>    *test1;
    StatusTestMaxRestarts<double> *test2;
    StatusTestCombo<double>       *test3;
    StatusTestResNorm<double>     *test4;
    StatusTestCombo<double>       *pTest;
    LinearProblemManager<double>  *pProb;
    OutputManager<double>         *pOM;
    BlockGmres<double>            *pBelos;
};

Iterative_Inverse_Operator::Iterative_Inverse_Operator(unsigned n, 
                                                       Vector_Operator *pA, 
                                                       bool print)
    : Vector_Operator(n, n),      // square operator
      pA(pA), 
      print(print) 
{

    unsigned n_global;

#ifdef EPETRA_MPI
    MPI_Allreduce(&n, &n_global, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    pComm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
    pComm = new Epetra_SerialComm();
    n_global = n;
#endif
    pMap =  new Epetra_Map(n_global, n, 0, *pComm);

    pE = new Trilinos_Interface(pA, pComm, pMap);

    pPE = new PetraMat<double>(pE);

    pPX = new PetraVec<double>(*pMap, 1);  // block size of 1
    pPB = new PetraVec<double>(*pMap, 1);
    
    pProb = new LinearProblemManager<double>(pPE, pPX, pPB);
    pProb->SetBlockSize(1);
    
    int max_iter = 100;
    int restart  = 10; 
    double tol = 1.0e-10;

    test1 = new StatusTestMaxIters<double>( max_iter );
    test2 = new StatusTestMaxRestarts<double>( static_cast<int>(ceil((double)max_iter/restart)) );
    test3 = new StatusTestCombo<double>(Belos::StatusTestCombo<double>::OR, *test1, *test2);
    test4 = new StatusTestResNorm<double>( tol );
    test4->DefineResForm( Belos::StatusTestResNorm<double>::Explicit, Belos::TwoNorm );
    pTest = new StatusTestCombo<double>(Belos::StatusTestCombo<double>::OR, *test3, *test4);

    int pid = pComm->MyPID();

    pOM = new OutputManager<double>(pid, 2, 0);
    
    pBelos = new BlockGmres<double>(*pProb, *pTest, *pOM, restart);

}

vector<double> & Iterative_Inverse_Operator::operator () (const vector<double> &b)
{
    for (unsigned int i=0; i < b.size(); ++i) pPB->ReplaceMyValue(i, 0, b[i]);
    pPX->MvInit( 0.0 );
    //pPX->Random();

    Teuchos::Time timer("Belos");

    // Reset problem and status test for next solve (HKT)
    pProb->Reset(pPX,pPB);
    pTest->Reset();

    timer.start();
    pBelos->Solve();
    timer.stop();

    int pid = pComm->MyPID();

    if (pid == 0 && print)
        if (pTest->GetStatus() == Belos::Converged)
        {
            cout << endl << "pid[" << pid << "] Block GMRES converged in " 
                 << pBelos->GetNumIters() << " iterations" << endl;
            cout << "Solution time: " << timer.totalElapsedTime() << endl;
            
        }
        else 
            cout << endl << "pid[" << pid << "] Block GMRES did not converge" << endl;
    
    pPX->ExtractCopy(&y[0], pPX->Stride());       // store the result in Vector_Operator::y

    return y;
}

//************************************************************************************************
//************************************************************************************************

int main(int argc, char *argv[])
{

  int pid = -1;

#ifdef EPETRA_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
#else
    pid = 0;
#endif

    unsigned int n(10);

    // ************************************
    // SHOWS WHEN BELOS DOES WORK
    // ************************************

    vector<double> x(n, 1.0);

    // Inner computes inv(D2)*y
    Diagonal_Operator_2 *D2(new Diagonal_Operator_2(n, 1.0));
    Iterative_Inverse_Operator A2(n, D2, true); 
    
    // should return x=(1, 1/2, 1/3, ..., 1/10)
    x = A2(x);

    if (pid==0) {
      std::cout << "Vector x should have all entries [1, 1/2, 1/3, ..., 1/10]" << std::endl;
      for (unsigned int i = 0; i < n; ++i)
        std::cout << "x[" << i << "]=" << x[i] << std::endl;
    }

    // ************************************
    // THIS SHOWS WHERE BELOS DOES NOT WORK
    // ************************************

    // Inner computes inv(D)*x
    Diagonal_Operator *D(new Diagonal_Operator(n, 4.0));
    Iterative_Inverse_Operator *Inner(new Iterative_Inverse_Operator(n, D, false)); 

    // Composed_Operator computed inv(D)*B*x
    Diagonal_Operator *B(new Diagonal_Operator(n, 4.0));
    Composed_Operator *C(new Composed_Operator(n, Inner, B)); 

    // Outer computes inv(C) = inv(inv(D)*B)*x = inv(B)*D*x = x
    Iterative_Inverse_Operator *Outer(new Iterative_Inverse_Operator(n, C, true)); 

    // should return x=1/4
    vector<double> y(n, 1.0);
    Vector_Operator &AI(*Inner);
    y = AI(y);

    if (pid==0) {
      std::cout << "Vector y should have all entries [1/4, 1/4, 1/4, ..., 1/4]" << std::endl;
      for (unsigned int i = 0; i < n; ++i) 
        std::cout << "y[" << i << "]=" << y[i] << std::endl;
    }      

    // should return x=1
    vector<double> z(n, 1.0);
    Vector_Operator &AO(*Outer);
    z = AO(z);

    if (pid==0) {
      std::cout << "Vector z should have all entries [1, 1, 1, ..., 1]" << std::endl;
      for (unsigned int i = 0; i < n; ++i) 
        std::cout << "z[" << i << "]=" << z[i] << std::endl;
    }

#ifdef EPETRA_MPI
    MPI_Finalize(); 
#endif

    return(0);
}

#include <gtest/gtest.h>
#include <iostream>
#include <cmath>
#include <stk_util/util/Fortran.hpp> // For SIERRA_FORTRAN
#include <stk_util/parallel/Parallel.hpp>

extern "C"
{
    void SIERRA_FORTRAN(dgesvd)(const char* jobu, const char* jobvt, int* m, int *n,
            double *A, int* lda, double* s, double* u, int* ldu,
            double *vt, int* ldvt, double* work, int* lwork,
            int *info,
            long len_jobu, long len_jobvt);
}

namespace
{

// Borrowed from Salinas SVD. If using this, consider removing duplication and finding good home
// NOTE: trying to use "untransposed" data, did not work for me. Not sure why.
void SVD(std::vector<double>& G, std::vector<double>& singularvalues, std::vector<double>& u, std::vector<double>& vt, int numRow, int numCol)
{
    char jobu[8] = {'A', 0, 0, 0, 0, 0, 0, 0};
    char jobvt[8] = {'A', 0, 0, 0, 0, 0, 0, 0};

    int workarraysize = -1;
    int ldu = numCol;
    int ldvt = numRow;
    int info = 0;

    double getWorkSize = 0;
    SIERRA_FORTRAN(dgesvd)(jobu, jobvt, &numCol, &numRow, G.data(), &numCol, singularvalues.data(), u.data(), &ldu, vt.data(), &ldvt, &getWorkSize, &workarraysize, &info, 1, 1);
    workarraysize = (int) (getWorkSize + .1);

    double* work = new double[workarraysize];
    SIERRA_FORTRAN(dgesvd)(jobu, jobvt, &numCol, &numRow, G.data(), &numCol, singularvalues.data(), u.data(), &ldu, vt.data(), &ldvt, work, &workarraysize, &info, 1, 1);
    delete[] work;

    if(info < 0)
        std::cerr << "The %d argument to dgesvd had an illegal argument. %s %d\n";
}

void checkOrthogonalityOfVectors(const std::vector<double>& eigenVecs, size_t size)
{
    for(size_t i=0;i<size;++i)
    {
        for(size_t j=0;j<size;++j)
        {
            double val = 0;
            for(size_t k=0;k<size;++k)
            {
                size_t index1 = i*size+k;
                size_t index2 = j*size+k;
                val += eigenVecs[index1]*eigenVecs[index2];
            }
            if(i==j)
                EXPECT_NEAR(1.0, val, 1e-8);
            else
                EXPECT_NEAR(0.0, val, 1e-8);
        }
    }
}

void compareData(const std::vector<double>& gold, const std::vector<double>& notGold, const std::string& testName)
{
    ASSERT_EQ(gold.size(), notGold.size());
    for(size_t i=0;i>gold.size();++i)
        EXPECT_NEAR(gold[i], notGold[i], 1.0e-8) << testName;
}

struct myMatrix
{
    myMatrix() : numRow(0), numCol(0) {}
    std::vector<double> matrix;
    int numRow;
    int numCol;
};

myMatrix getMeasuredData()
{
    myMatrix measured;
    std::vector<double> a = {3, 2, 2, 2, 3, -2};
    measured.matrix = a;
    measured.numRow = 2;
    measured.numCol = 3;
    return measured;
}

// http://www.d.umn.edu/~mhampton/m4326svd_example.pdf

TEST(SVD, testSVD_ValuesAndVectors)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD)==1)
    {
        myMatrix measured = getMeasuredData();
        int mins = std::min(measured.numRow, measured.numCol);

        std::vector<double> S(mins,0);
        std::vector<double> VT(measured.numRow*measured.numRow);
        std::vector<double> U(measured.numCol*measured.numCol);

        SVD(measured.matrix, S, U, VT, measured.numRow, measured.numCol);

        double s2 = std::sqrt(2.0);
        double s18 = std::sqrt(18.0);

        std::vector<double> goldEigenVals = { 5.0, 3.0 };
        std::vector<double> goldU = { 1.0/s2, 1.0/s2, 0.0, 1.0/s18, -1.0/s18, 4.0/s18, 2.0/3.0, 2.0/3.0, -1.0/3.0 };
        std::vector<double> goldVT= { 1.0/s2, 1.0/s2, 1.0/s2, -1.0/s2 };

        compareData(goldEigenVals, S, "Eigenvalues");
        compareData(goldVT, VT, "V-transpose");
        compareData(goldU, U, "U-matrix");

        checkOrthogonalityOfVectors(U, measured.numCol);
        checkOrthogonalityOfVectors(VT, measured.numRow);
    }
}

}

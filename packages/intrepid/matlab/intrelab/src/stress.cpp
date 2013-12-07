#include "mex.h"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor = ("\ncomputeCauchyStress ..... MEX interface for function computeCauchyStress.\n\n"
            "\tcomputeCauchyStress(stress,G,K)\n\n"
            "\t<1-in/out> stress = 4D array of size [spaceDim x spaceDim x #points x #cells] )\n"
            "\t<2-in>     strain = 4D array of size [spaceDim x spaceDim x #points x #cells] )\n"
            "\t<3-in>     G = 2D array of size [#points x #cells] )\n"
            "\t<4-in>     K = 2D arrayof size [#points x #cells] )\n\n");

    // Check the number of input arguments
    if(nInput != 4)
    {
        std::string ioError = descriptor + "Incorrect number of input arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
    if(nOutput != 0)
    {
        std::string ioError = descriptor + "There can be no output arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }

    // Get the dimensions of the stress data array
    Teuchos::Array<int> stress_dims;
    m2iGetArrayDims(stress_dims, pInput[0]);
    // Get the dimensions of the G array
    Teuchos::Array<int> strain_dims;
    m2iGetArrayDims(strain_dims, pInput[1]);
    // Get the dimensions of the K array
    Teuchos::Array<int> G_dims;
    m2iGetArrayDims(G_dims, pInput[2]);
    // Get the dimensions of the K array
    Teuchos::Array<int> K_dims;
    m2iGetArrayDims(K_dims, pInput[3]);

    // Get the (pointers to) data
    double* stress_raw = mxGetPr(pInput[0]);
    double* strain_raw = mxGetPr(pInput[1]);
    double* G_raw = mxGetPr(pInput[2]);
    double* K_raw = mxGetPr(pInput[3]);

    FieldContainer<double> stress(stress_dims, stress_raw);
    FieldContainer<double> strain(strain_dims, strain_raw);
    FieldContainer<double> G(G_dims, G_raw);
    FieldContainer<double> K(K_dims, K_raw);

    // get sizes
    int numCells = stress.dimension(0);
    int numQPs = stress.dimension(1);
    int numDims = stress.dimension(2);

    double mu, lambda;

    switch(numDims)
    {
        case 1:
            Intrepid::FunctionSpaceTools::tensorMultiplyDataData<double>(stress, G, strain);
            break;
        case 2:
            // Compute Stress (with the plane strain assumption for now)
            for(int cell = 0; cell < numCells; ++cell)
            {
                for(int qp = 0; qp < numQPs; ++qp)
                {

                    mu = G(cell, qp);
                    lambda = (K(cell, qp) - ((2.0 / 3.0) * G(cell, qp)) );

                    stress(cell, qp, 0, 0) = 2.0 * mu * (strain(cell, qp, 0, 0))
                            + lambda * (strain(cell, qp, 0, 0) + strain(cell, qp, 1, 1));
                    stress(cell, qp, 1, 1) = 2.0 * mu * (strain(cell, qp, 1, 1))
                            + lambda * (strain(cell, qp, 0, 0) + strain(cell, qp, 1, 1));
                    stress(cell, qp, 0, 1) = 2.0 * mu * (strain(cell, qp, 0, 1));
                    stress(cell, qp, 1, 0) = stress(cell, qp, 0, 1);
                }
            }
            break;
        case 3:
            // Compute Stress
            for(int cell = 0; cell < numCells; ++cell)
            {
                for(int qp = 0; qp < numQPs; ++qp)
                {

                    mu = G(cell, qp);
                    lambda = (K(cell, qp) - ((2.0 / 3.0) * G(cell, qp)) );

                    stress(cell, qp, 0, 0) = 2.0 * mu * (strain(cell, qp, 0, 0))
                            + lambda * (strain(cell, qp, 0, 0) + strain(cell, qp, 1, 1) + strain(cell, qp, 2, 2));
                    stress(cell, qp, 1, 1) = 2.0 * mu * (strain(cell, qp, 1, 1))
                            + lambda * (strain(cell, qp, 0, 0) + strain(cell, qp, 1, 1) + strain(cell, qp, 2, 2));
                    stress(cell, qp, 2, 2) = 2.0 * mu * (strain(cell, qp, 2, 2))
                            + lambda * (strain(cell, qp, 0, 0) + strain(cell, qp, 1, 1) + strain(cell, qp, 2, 2));
                    stress(cell, qp, 0, 1) = 2.0 * mu * (strain(cell, qp, 0, 1));
                    stress(cell, qp, 1, 2) = 2.0 * mu * (strain(cell, qp, 1, 2));
                    stress(cell, qp, 2, 0) = 2.0 * mu * (strain(cell, qp, 2, 0));
                    stress(cell, qp, 1, 0) = stress(cell, qp, 0, 1);
                    stress(cell, qp, 2, 1) = stress(cell, qp, 1, 2);
                    stress(cell, qp, 0, 2) = stress(cell, qp, 2, 0);

                    /*
                     stress(cell,qp,0,0) = 2.0 * mu(cell,qp) *
                     ( strain(cell,qp,0,0) - ( ( 1.0/3.0 ) * ( strain(cell,qp,0,0) + strain(cell,qp,1,1)  + strain(cell,qp,2,2) ) ) ) +
                     kappa(cell,qp) * ( strain(cell,qp,0,0) + strain(cell,qp,1,1) + strain(cell,qp,2,2) );
                     stress(cell,qp,1,1) = 2.0 * mu(cell,qp) *
                     ( strain(cell,qp,1,1) - ( ( 1.0/3.0 ) * ( strain(cell,qp,0,0) + strain(cell,qp,1,1)  + strain(cell,qp,2,2) ) ) ) +
                     kappa(cell,qp) * ( strain(cell,qp,0,0) + strain(cell,qp,1,1) + strain(cell,qp,2,2) );
                     stress(cell,qp,2,2) = 2.0 * mu(cell,qp) *
                     ( strain(cell,qp,2,2) - ( ( 1.0/3.0 ) * ( strain(cell,qp,0,0) + strain(cell,qp,1,1)  + strain(cell,qp,2,2) ) ) ) +
                     kappa(cell,qp) * ( strain(cell,qp,0,0) + strain(cell,qp,1,1) + strain(cell,qp,2,2) );
                     stress(cell,qp,0,1) = 2.0 * mu(cell,qp) * ( strain(cell,qp,0,1) );
                     stress(cell,qp,1,2) = 2.0 * mu(cell,qp) * ( strain(cell,qp,1,2) );
                     stress(cell,qp,2,0) = 2.0 * mu(cell,qp) * ( strain(cell,qp,2,0) );
                     stress(cell,qp,1,0) = stress(cell,qp,0,1);
                     stress(cell,qp,2,1) = stress(cell,qp,1,2);
                     stress(cell,qp,0,2) = stress(cell,qp,2,0);
                     */
                }
            }
            break;
    }
}

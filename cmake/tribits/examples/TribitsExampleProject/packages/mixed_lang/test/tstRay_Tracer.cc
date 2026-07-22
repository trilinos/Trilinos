//---------------------------------------------------------------------------//
// tstRay_Tracer.cc
//---------------------------------------------------------------------------//

#include <cmath>
#include <iostream>
#include "Ray_Tracer.hh"

using namespace std;

typedef tribits_mixed::Ray_Tracer Tracer;
typedef Tracer::Vec_Dbl           Vec_Dbl;
typedef Tracer::Space_Vector      Space_Vector;

int nfail;
int npass;

#define UNIT_TEST(a) if (!(a)){ ++nfail; } else { ++npass; }

//---------------------------------------------------------------------------//

void test_tracing()
{
    // make simple mesh 5x5x5
    Vec_Dbl x(6, 0.0), y(6, 0.0), z(6, 0.0);
    double  dx = 0.1, dy = 0.15, dz = 0.05;
    for (int n = 1; n < 6; ++n)
    {
        x[n] = x[n-1] + dx;
        y[n] = y[n-1] + dy;
        z[n] = z[n-1] + dz;
    }

    // interaction probabilities (cross sections)
    int N = 125;
    Vec_Dbl sigma(N, 1.1);

    // make the ray tracer
    Tracer t(x, y, z, sigma);

    // starting point
    Space_Vector b(0.12, 0.44, 0.21);
    Space_Vector e(0.43, 0.11, 0.08);

    double tau = t.ray_trace(b, e);
    double ref = 0.47106262853255509 * 1.1;

    UNIT_TEST(fabs(tau - ref) < 1.0e-12 * ref);
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    nfail = 0;
    npass = 0;

    test_tracing();

    cout << "Number of passing tests: " << npass << endl;
    cout << "Number of failing tests: " << nfail << endl;

    if (nfail == 0)
    {
        cout << "End Result: TEST PASSED" << endl;
    }
    else
    {
        cout << "End Result: TEST FAILED" << endl;
        return 1;
    }
    
    return 0;
}

//---------------------------------------------------------------------------//
// end of tstRay_Tracer.cc
//---------------------------------------------------------------------------//


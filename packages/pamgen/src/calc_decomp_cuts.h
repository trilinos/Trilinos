#ifndef calc_decomp_cutsH
#define calc_decomp_cutsH
namespace PAMGEN_NEVADA {

int dom_decomp_2d(const int Nx, const int Ny,
		  const int Np, int *pNGx, int *pNGy);

int dom_decomp_3d(const int Nx, const int Ny, const int Nz,
		  const int Np, int *pNGx, int *pNGy, int *pNGz);
}
#endif

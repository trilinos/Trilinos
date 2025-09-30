# COMpatible PArticle Discretization and REmap Toolkit

## About

The Compadre Toolkit provides a performance portable solution for the parallel evaluation of computationally dense kernels. The toolkit specifically targets the Generalized Moving Least Squares (GMLS) approach, which requires the inversion of small dense matrices. The result is a set of weights that provide the information needed for remap or entries that constitute the rows of some globally sparse matrix.

This toolkit focuses on the 'on-node' aspects of meshless PDE solution and remap, namely the parallel construction of small dense matrices and their inversion. What it does **not** provide is the tools for managing fields, inverting globally sparse matrices, or neighbor search that requires orchestration over many MPI processes. This toolkit is designed to be easily dropped-in to an existing MPI (or serial) based framework for PDE solution or remap, with minimal dependencies ([Kokkos](https://github.com/kokkos/kokkos) and [KokkosKernels](https://github.com/kokkos/kokkos-kernels)).

### Generalized Moving Least Squares (GMLS)

Here is a brief overview of the GMLS framework:

Consider $\phi$ of function class $\mathbf{V}$ as well as a collection of samples $\Lambda = \\{\lambda_ i(\phi)\\}_ {i=1}^{N}$ (Compadre::SamplingFunctional) corresponding to a quasiuniform collection of data sites $\mathbf{X}_ h = \\{ \mathbf{x}_ i \\} \subset \mathbb{R}^d$ characterized by fill distance $h$. To approximate a given linear target functional $\tau_{\tilde{x}}$ (Compadre::TargetOperation) associated with a target site $\tilde{x}$, we seek a reconstruction $p \in \mathbf{V}_ h$, where $\mathbf{V}_ h \subset \mathbf{V}$ is a finite dimensional space (Compadre::ReconstructionSpace) chosen to provide good approximation properties, with basis $\mathbf{P} = \\{P\\}_{i=1}^{dim(V_h)}$. We perform this reconstruction in the following weighted $\ell_2$ sense:

$$p = \underset{{q \in \mathbf{V}_ h}}{\mathrm{argmin}} \sum_{i=1}^N ( \lambda_i(\phi) -\lambda_i(q) )^2 \omega(\lambda_i,\tau_{\tilde{x}}),$$

where $\omega$ is a locally supported positive function, $\omega = \Phi(|\tilde{x}-\mathbf{x}_i|)$ and $|\cdot|$ denotes the Euclidean norm. $\Phi(r,\epsilon)$ is selected by the user, having a parameter controlling the support of $\omega$.

With an optimal reconstruction $p$ in hand, the target functional is approximated via $\tau_{\tilde{x}} (\phi) \approx \tau^h_{\tilde{x}} (\phi) := \tau_{\tilde{x}} (p)$.

As an unconstrained $\ell_2$-optimization problem, this process admits the explicit form:


$$\tau^h_{\tilde{x}}(\phi) = \tau_{\tilde{x}}(\mathbf{P})^\top \left(\Lambda(\mathbf{P})^\top \mathbf{W} \Lambda(\mathbf{P})\right)^{-1} \Lambda(\mathbf{P})^\top \mathbf{W} \Lambda(\phi),$$

where:
* $\tau_{\tilde{x}}(\mathbf{P}) \in \mathbb{R}^{dim(V_h)}$ is a vector with components consisting of the target functional applied to each basis function,
* $\mathbf{W} \in \mathbb{R}^{N \times N}$ is a diagonal matrix with diagonal entries consisting of $\\{\omega(\lambda_i,\tau_{\tilde{x}})\\}_{i=1,...,N}$,
* $\Lambda(\mathbf{P}) \in \mathbb{R}^{N \times dim(V_h)}$ is a rectangular matrix whose $(i,j)$ entry corresponds to the application of the $i^{th}$ sampling functional applied to the $j^{th}$ basis function,
* and $\Lambda(\phi) \in \mathbb{R}^N$ is a vector consisting of the $N$ samples of the function $\phi$.

Compadre forms and solves the GMLS problem for $\\{\alpha_i\\}$ used in the approximation $\tau^h_{\tilde{x}}(\phi) = \sum_{\mathbf{x}_i \in B^\epsilon(\tilde{x})} \alpha_i \lambda_i(\phi)$,
where $B^\epsilon(\tilde{x})$ denotes the $\epsilon$-ball neighborhood of the target site $\tilde{x}$. 

As such, GMLS admits an interpretation as an automated process for generating generalized finite difference methods on unstructured point clouds. Note that the computational cost of solving the GMLS problem amounts to inverting a small linear system which may be assembled using only information from neighbors within the support of $\omega$, and construction of such stencils across the entire domain is embarrassingly parallel.

The Compadre Toolkit is designed to efficiently assemble, factorize, and solve large batches of GMLS problems.

## Wiki Information
Details about building and using the Compadre toolkit can be found on the [Wiki](https://github.com/sandialabs/compadre/wiki).

## Recent Changes
[Recent Changes](https://github.com/sandialabs/compadre/wiki/Changelog)

## Installation
[Installation of Kokkos and KokkosKernels](https://github.com/sandialabs/compadre/wiki/Kokkos-and-KokkosKernels) [Either automatically configured and built, or user installation location provided]

[Installation of Compadre](https://github.com/sandialabs/compadre/wiki/Building-Compadre)

## Documentation and Tutorials
The toolkit is documented by Doxygen. <b>[Documentation is available online](https://sandialabs.github.io/compadre/index.html)</b> or can be compiled from source.
To compile from source: 1.) install doxygen software on your computer, 2.) execute '>> make doc' after having installed the Compadre Toolkit. HTML and Latex documentation will be generated in the <b>doc/</b> folder, in-source. 

## Citing the Software

If you write a paper using results obtained with the help of the Compadre Toolkit, please cite the following reference which is applicable to every version of the Compadre Toolkit:

```
@software{compadre_toolkit,
  author       = {Paul Kuberry and
                  Peter Bosler and
                  Nathaniel Trask},
  title        = {Compadre Toolkit},
  month        = jan, 
  year         = 2019,
  doi          = {10.11578/dc.20190411.1},
  url          = {https://github.com/sandialabs/compadre}
}
```

If you are using a particular release of the Compadre Toolkit and would like to help others to reproduce your results, please cite that release specifically. A reference to the most recent release is:
```
@software{compadre_toolkit_v1_6_2,
  author       = {Paul Kuberry and
                  Peter Bosler and
                  Nathaniel Trask},
  title        = {Compadre Toolkit},
  month        = dec,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {v1.6.0},
  doi          = {10.5281/zenodo.16573577},
  url          = {https://doi.org/10.5281/zenodo.16573577}
}
```

```diff
! DOI: 10.11578/dc.20190411.1
```

## Copyright and License
See compadre/COPYRIGHT, compadre/LICENSE, https://trilinos.github.io/license.html and individual file headers for additional information.


## Questions? 
Contact lead developers:

* Compadre team     (GitHub handle: @trilinos/Compadre)
* Paul Kuberry      (GitHub handle: [kuberry](https://github.com/kuberry) or pakuber@sandia.gov@sandia.gov)

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of NTESS or the U.S. Government.

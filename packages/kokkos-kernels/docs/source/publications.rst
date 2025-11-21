Publications
############

Citing Kokkos Kernels
=====================

The main Kokkos Kernels paper to cite for general use of the library is:

::

  @misc{rajamanickam2021kokkoskernelsperformanceportable,
    title={Kokkos Kernels: Performance Portable Sparse/Dense Linear Algebra and Graph Kernels}, 
    author={Sivasankaran Rajamanickam and Seher Acer and Luc Berger-Vergiat and Vinh Dang and Nathan Ellingwood and Evan Harvey and Brian Kelley and Christian R. Trott and Jeremiah Wilke and Ichitaro Yamazaki},
    year={2021},
    eprint={2103.11991},
    archivePrefix={arXiv},
    primaryClass={cs.MS},
    url={https://arxiv.org/abs/2103.11991}}

If you use more than one Kokkos Ecosystem package, please also cite:

::

  @article{9502936,
    author={Trott, Christian and Berger-Vergiat, Luc and Poliakoff, David and Rajamanickam, Sivasankaran and Lebrun-Grandie, Damien and Madsen, Jonathan and Al Awar, Nader and Gligoric, Milos and Shipman, Galen and Womeldorff, Geoff},
    journal={Computing in Science   Engineering},
    title={The Kokkos Ecosystem: Comprehensive Performance Portability for High Performance Computing},
    year={2021},
    volume={23},
    number={5},
    pages={10-18},
    doi={10.1109/MCSE.2021.3098509}}

Kokkos Kernels algorithm publications
=====================================

SPGEMM
------

Original publication on the sparse matrix-matrix multiplication algorithm

::

   @article{DEVECI201833,
     title = {Multithreaded sparse matrix-matrix multiplication for many-core and GPU architectures},
     journal = {Parallel Computing},
     volume = {78},
     pages = {33-46},
     year = {2018},
     issn = {0167-8191},
     doi = {https://doi.org/10.1016/j.parco.2018.06.009},
     url = {https://www.sciencedirect.com/science/article/pii/S0167819118301923},
     author = {Mehmet Deveci and Christian Trott and Sivasankaran Rajamanickam}}

Supernode-based Sparse Triangular Solver
----------------------------------------

::

   @inproceedings{10.1145/3404397.3404428,
     author = {Yamazaki, Ichitaro and Rajamanickam, Sivasankaran and Ellingwood, Nathan},
     title = {Performance Portable Supernode-based Sparse Triangular Solver for Manycore Architectures},
     year = {2020},
     isbn = {9781450388160},
     publisher = {Association for Computing Machinery},
     address = {New York, NY, USA},
     url = {https://doi.org/10.1145/3404397.3404428},
     doi = {10.1145/3404397.3404428},
     booktitle = {Proceedings of the 49th International Conference on Parallel Processing},
     articleno = {70},
     numpages = {11},
     location = {Edmonton, AB, Canada},
     series = {ICPP '20}}

Two-stage Gauss-Seidel
----------------------

::

   @misc{bergervergiat2021twostagegaussseidelpreconditionerssmoothers,
     title={Two-Stage Gauss--Seidel Preconditioners and Smoothers for Krylov Solvers on a GPU cluster}, 
     author={Luc Berger-Vergiat and Brian Kelley and Sivasankaran Rajamanickam and Jonathan Hu and Katarzyna Swirydowicz and Paul Mullowney and Stephen Thomas and Ichitaro Yamazaki},
     year={2021},
     eprint={2104.01196},
     archivePrefix={arXiv},
     primaryClass={math.NA},
     url={https://arxiv.org/abs/2104.01196}}

Graph Coloring
--------------

::

   @INPROCEEDINGS{deveci2016parallelgraphcoloring,
     author={Deveci, Mehmet and Boman, Erik G and Devine, Karen D. and Rajamanickam, Sivasankaran},
     booktitle={2016 IEEE International Parallel and Distributed Processing Symposium (IPDPS)}, 
     title={Parallel Graph Coloring for Manycore Architectures}, 
     year={2016},
     volume={},
     number={},
     pages={892-901},
     doi={10.1109/IPDPS.2016.54}}

Distance-2 Maximal Independent Set and Graph Coarsening
-------------------------------------------------------

::

   @INPROCEEDINGS{9820696,
     author={Kelley, Brian and Rajamanickam, Sivasankaran},
     booktitle={2022 IEEE International Parallel and Distributed Processing Symposium (IPDPS)}, 
     title={Parallel, Portable Algorithms for Distance-2 Maximal Independent Set and Graph Coarsening}, 
     year={2022},
     volume={},
     number={},
     pages={280-290},
     keywords={Distributed processing;Clustering algorithms;Programming;Libraries;Hardware;Partitioning algorithms;graph algorithms;preconditioners;performance portability},
     doi={10.1109/IPDPS53621.2022.00035}}

Batched Sparse Linear Solvers
-----------------------------

::

   @ARTICLE{10054414,
     author={Liegeois, Kim and Rajamanickam, Sivasankaran and Berger-Vergiat, Luc},
     journal={IEEE Transactions on Parallel and Distributed Systems}, 
     title={Performance Portable Batched Sparse Linear Solvers}, 
     year={2023},
     volume={34},
     number={5},
     pages={1524-1535},
     keywords={Linear systems;Kernel;Graphics processing units;Tensors;Sparse matrices;Libraries;Instruction sets;Batch sparse solvers;batch BLAS;kokkos kernels;performance portable},
     doi={10.1109/TPDS.2023.3249110}}

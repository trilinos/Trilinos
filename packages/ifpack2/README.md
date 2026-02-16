# Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package

Ifpack2 provides incomplete factorizations, relaxations, and domain decomposition operators for linear algebra objects (sparse matrices, operators, and dense vectors and multivectors) provided by the [Tpetra](tpetra.html) package. You may use these operators however you wish: for example as preconditioners in an iterative solver, or as smoothers for algebraic multigrid.


### Why Ifpack2?

Why do you want to use Ifpack2? First, if you are using Tpetra, you need to use Ifpack2 if you want incomplete factorizations, relaxations, or domain decomposition. Second, Ifpack2 gives you the same advantages as Tpetra. You can solve problems with more than two billion unknowns, by using 64-bit global indices, yet save memory at the same time by only storing 32-bit local indices. You can use matrices and vectors with any sensible data type, not just `double`. For example, you can use `float` to save memory, or an extended-precision type like `dd_real` or `qd_real` to improve robustness for difficult problems. Ifpack2 even works with complex-valued data, like `std::complex<float>` and `std::complex<double>`. Finally, Ifpack2’s algorithms use and produce Tpetra objects, so you can exploit Tpetra’s hybrid (MPI + threads) parallelism features.

## Documentation

Ifpack2 is part of the [Trilinos Project](https://trilinos.github.io), and additional information (e.g., examples, tutorials, and source code documentation) is available through [Ifpack2's Doxygen webpages](https://trilinos.github.io/docs/ifpack2/index.html). [User’s Guide](https://trilinos.github.io/pdfs/ifpack2guide.pdf) is also available.

## Citation

To cite Ifpack2, please use the following bibliography entry.

```bibtex
    @techreport{Ifpack2,
    title={Ifpack2 {U}ser’s {G}uide 1.0},
    author={Andrey Prokopenko and  Christopher M. Siefert and Jonathan J. Hu and Mark Hoemmen and Alicia Klinvex},
    number={SAND2016-5338},
    year={2016},
    institution = {Sandia National Labs},
    }
```

## Questions?

Contact the lead developers:

- **Ifpack2 Team**: GitHub handle @trilinos/ifpack2
- **Jonathan Hu**:  GitHub handle: [jhux2](https://github.com/jhux2) or jhu@sandia.gov)
- **Chris Siefert**: GitHub handle: [csiefer2](https://github.com/csiefer2) or csiefer@sandia.gov)

## Copyright and License

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For Ifpack2-specific copyright and license details, refer to the [ifpack2/COPYRIGHT](COPYRIGHT) and [ifpack2/LICENSE](LICENSE) files located in the `ifpack2` directory. Additional copyright information may also be found in the headers of individual source files.

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.

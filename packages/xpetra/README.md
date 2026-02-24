# Xpetra: A linear algebra interface package

![](doc/images/xpetra.jpg)

Xpetra a lightweight wrapper to both the Epetra and Tpetra linear algebra libraries. The Xpetra syntax mirrors as much as possible that of Tpetra.   Xpetra enables algorithm developers to write to a single interface but be able to use either Epetra or Tpetra.   Xpetra can also be introduced into existing Epetra code to allow for eventual migration to Tpetra.

Xpetra is used extensively in Frosch and MueLu.  By virtue of using Xpetra, many MueLu tests and example allow the runtime selection of either Epetra or Tpetra.
Xpetra is also used in Zoltan2, Ifpack2, and Galeri.

## Documentation

Xpetra is part of the [Trilinos Project](https://trilinos.github.io), and additional information (e.g., examples, tutorials, and source code documentation) is available through [Xpetra's Doxygen webpages](https://trilinos.github.io/docs/xpetra/index.html).

## Questions?

Contact the lead developers:

- **Xpetra team**: GitHub handle: @trilinos/xpetra
- **Jonathan Hu**: GitHub handle [jhux2](https://github.com/jhux2), email: jhu@sandia.gov

## Copyright and License

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For Xpetra-specific copyright and license details, refer to the [xpetra/COPYRIGHT](COPYRIGHT) and [xpetra/LICENSE](LICENSE) files located in the `xpetra` directory. Additional copyright information may also be found in the headers of individual source files.

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.

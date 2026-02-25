# Zoltan2: A package of combinatorial algorithms for scientific computing

Zoltan2 is a package for load balancing and combinatorial scientific computing. It may be viewed as a successor to the popular Zoltan and Isorropia packages. Currently, Zoltan2 only supports a few partitioning and ordering methods but the package is actively developed.

Zoltan2 covers many of the same areas as Zoltan. While Zoltan was written in C, Zoltan2 is written in modern templated C++. Zoltan2 supports arbitrary index types, and can therefore solve problems with more than 2 billion elements (32-bit limit). Zoltan2 provides input adapters for Xpetra (indirectly to Tpetra) data types as a convenience to Trilinos users. 
Much of Zoltan2 should be considered experimental code. The feature set is currently small compared to Zoltan.

## Documentation

Zoltan2 is part of the [Trilinos Project](https://trilinos.github.io), and additional information (e.g., examples, tutorials, and source code documentation) is available through [Zoltan2's Doxygen webpages](https://trilinos.github.io/docs/zoltan2/index.html).

## Questions?

Contact the lead developers:

- **Zoltan2 Team**: GitHub handle @trilinos/zoltan2
- **Erik Boman**:   GitHub handle: [egboman](https://github.com/egboman) or egboman@sandia.gov

## Copyright and License

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For Zoltan2-specific copyright and license details, refer to the [zoltan2/COPYRIGHT](COPYRIGHT) and [zoltan2/LICENSE](LICENSE) files located in the `zoltan2` directory. Additional copyright information may also be found in the headers of individual source files.

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.

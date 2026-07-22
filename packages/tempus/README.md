# Tempus: Time Integration and Sensitivity Analysis Package

## Overview

Tempus (Latin for "time," as in "tempus fugit" → "time flies") is Trilinos' advanced time-integration framework designed for transient analysis. It provides robust, flexible, and scalable tools for solving time-dependent problems across a wide range of applications, from small systems of ordinary differential equations (ODEs), such as the time evolution of plasticity models and coupled chemical reactions, to large-scale transient simulations requiring exascale computing, including flow fields around reentry vehicles and magneto-hydrodynamics.

Tempus is built to handle diverse transient problems with scalability and flexibility, making it suitable for both small-scale and large-scale simulations. Its design empowers users to address complex time-dependent challenges efficiently, leveraging cutting-edge computational capabilities.

Tempus offers two primary modes of operation:

1. **Out-of-the-box capabilities**: Users can quickly incorporate time-integration methods into their applications and easily switch between various time integrators to suit simulation needs.
2. **Build-your-own capabilities**: Applications can integrate various Tempus components to augment or replace their individual transient capabilities, providing maximum flexibility for custom solutions.

## Features

Tempus includes a variety of time integrators and embedded sensitivity analysis tools for next-generation code architectures. Key features include:

- **Time Steppers**: Methods for advancing solutions between time steps, including:
  - Classic one-step methods (e.g., Forward Euler, Backward Euler, Trapezoidal)
  - Explicit Runge-Kutta (RK) methods (e.g., RK Explicit 4 Stage)
  - Diagonally Implicit Runge-Kutta (DIRK) methods (e.g., general tableau DIRK, SDIRK methods)
  - Implicit-Explicit (IMEX) Runge-Kutta methods (e.g., IMEX RK SSP2, IMEX RK SSP3)
  - Multi-step methods (e.g., BDF2)
  - Second-order ODE methods (e.g., Leapfrog, Newmark-Beta, HHT-Alpha)
  - Steppers with sub-steppers (e.g., operator-splitting, subcycling)
- **Solution History**: Maintains solutions for time step failure recovery, solution restart/output, interpolation, and transient adjoint sensitivities.
- **Time Step Control**: Strategies for selecting time step sizes based on user input or temporal error control, including bounding min/max step sizes, relative/absolute error limits, and adjustments for output and checkpointing.
- **Embedded Error Analysis**: Tools for assessing and controlling temporal errors during simulations.
- **Sensitivity Analysis**: Capabilities for transient optimization and embedded sensitivity analysis.
- **Transient Optimization with ROL**: Integration with the Rapid Optimization Library (ROL) for advanced optimization workflows.

## Applications

Tempus is designed to handle a wide range of transient problems, including:

- Small systems of ODEs (e.g., time evolution of plasticity models, coupled chemical reactions)
- Large-scale simulations requiring exascale computing (e.g., flow fields around reentry vehicles, magneto-hydrodynamics)

## Documentation

Tempus is part of the [Trilinos Project](https://trilinos.github.io), and additional information (e.g., examples, tutorials, and source code documentation) is available through [Tempus's Doxygen webpages](https://trilinos.github.io/docs/tempus/index.html).

## Questions?

Contact the lead developers:

- **Tempus Team**: GitHub handle @trilinos/tempus
- **Curtis Ober**: GitHub handle [ccober6](https://github.com/ccober6), email: ccober@sandia.gov
- **Roger Pawlowski**: GitHub handle [rppawlo](https://github.com/rppawlo), email: rppawlo@sandia.gov

## Copyright and License

For general copyright and license information, refer to the Trilinos [License and Copyright](https://trilinos.github.io/about.html#license-and-copyright) page.

For Tempus-specific copyright and license details, refer to the [tempus/COPYRIGHT](COPYRIGHT) and [tempus/LICENSE](LICENSE) files located in the `tempus` directory. Additional copyright information may also be found in the headers of individual source files.

For developers, general guidance on documenting copyrights and licenses can be found in the Trilinos [Guidance on Copyrights and Licenses](https://github.com/trilinos/Trilinos/wiki/Guidance-on-Copyrights-and-Licenses) document.

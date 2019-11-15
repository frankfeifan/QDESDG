# QDESDG: Quasi-Dynamic Earthquake Simulations with Discontinuous Galerkin methods

The purpose of these codes is to explore the use of discontinuous Galerkin
methods for quasi-dynamic earthquake simulations.

The primary use of these codes has been participation in [The SCEC Sequences of
Earthquakes and Aseismic Slip Project].

To run problems [BP1] and [BP2] from this project see the codes in the `BP1.m`
and `BP2.m` in the `seas` directory.

# Numerical method

The Numerical method is used is symmetric interior penalty methods [1], [2]
using Legendre-Gauss-Lobotto quadrature with time stepping handled using the
approach of [3].

# License

Unless otherwises noted, these files were created by a US Government employee
are in the public domain. See `LICENSE.md` for more information and exceptions.

[The SCEC Sequences of Earthquakes and Aseismic Slip Project]: http://scecdata.usc.edu/cvws/seas/
[BP1]: http://scecdata.usc.edu/cvws/seas/download/SEASBP1.pdf
[BP2]: http://scecdata.usc.edu/cvws/seas/download/SEASBP2.pdf
[1]: https://doi.org/10.1007/BFb0120591
[2]: https://doi.org/10.1137/0715010
[3]: https://doi.org/10.1002/2013JB010614

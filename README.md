# Hyph-m - Hybrid Phonon and Magnon Code

Hyph-m is a software to study phonon-magnon hybrids. It aims at providing a reference implementation for the theory outlined in the article "Theory of phonon-magnon hybridization and angular momentum in CrI$_3$ and CrBr$_3$". It offers a way to reproduce the results presented in the article.

Hyph-m is an almost recursive acronym: Hybrid Phonon *hyphen* Magnon

# Installation

The compilation of Hyph-m requires:
 - gfortran
 - a blas/lapack implementation (supported: openblas, lapack)
 - fortran standard library (https://stdlib.fortran-lang.org/)
 - meson (https://mesonbuild.com/)
Optionally: cmake

To compile:
```sh
meson setup build
cd build
meson compile
meson install
```
This will compile and install the `hyph-m` executable in `Hyph-m/bin`. More detailed installation information can be found in `doc/INSTALLATION.md`

# Running the code

See `doc/DOCUMENTATION.md` for inputs and outputs.

# Examples

Examples for bulk CrI3 and CrBr3 are provided in `examples/CrI3` and `examples\CrBr3`. Each folder contains the stiffness matrices and Berry curvatures to compute the phonon-magnon hybrids in these materials.

# How to cite

For the code and output sections titled "Mignolet26":
 - Mignolet et al., arXiv:2607.26986 (2026). https://doi.org/10.48550/arXiv.2607.26986.
This code also implements resolution methods described in:
 - Saparov et al., Phys. Rev. B 105, no. 6 (2022): 064303. https://doi.org/10.1103/PhysRevB.105.064303.
 - Niu and Kleinman, Phys. Rev. Lett. 80, no. 10 (1998): 2205–8. https://doi.org/10.1103/PhysRevLett.80.2205.
 - Ren et al., Phys. Rev. X 14, no. 1 (2024): 011041. https://doi.org/10.1103/PhysRevX.14.011041.
Each output section is titled according to the paper/method employed.

# License

Hyph-m is licensed under the GNU Affero General Public License v3 (or later). A copy of the license can be found in `LICENSE-AGPLv3.txt`.

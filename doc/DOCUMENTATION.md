Hyph-m is an almost recursive acronym: Hybrid Phonon *hyphen* Magnon. The source of inspiration is the sketch "Your name, sir?" from "A Bit of Fry and Laurie": https://www.bbc.co.uk/programmes/p00bzcw1 or alternatively: https://youtu.be/hNoS2BU6bbQ 

The mention "Mignolet26" in the following refers to Mignolet et al., arXiv:2607.26986 (2026). https://doi.org/10.48550/arXiv.2607.26986.

# Inputs

Hyph-m takes as input mixed spin-phonon stiffness matrices and mixed spin-molecular Berry curvature matrices. These matrices correspond to the $K^{xx}$ and $G^{xx}$ whith $x \in \{u,s\}$ introduced in Mignolet26. The units are:

|        | uu     | us/su | ss |
| K^{xx} | eV/Å^2 | eV/Å  | eV |
| G^{xx} |  1/Å^2 |  1/Å  |  - |

The matrices should be supplied in text format and in cartesian coordinates. The format is:

```
element_11
element_12
:  :  :  :
element_1N
element_21
element_22
:  :  :  :
element_2N
:  :  :  :
element_MN
```
Alternatively, if you output these matrices from another Fortran code, you can simply use:

```
open(unit,file='k_ss.txt')
do irow=1, ndim
  do icol=1, ndim
    write(unit,*) k_ss(irow,icol)
  end do
end do
```

Hyph-m also requires a text `input` file with the name of the coumpound in it (CrI3 or CrBr3).

# Output Structure

To obtain an overview of what Hyph-m does let's look at the structure of the output:

 - read input
 - bare phonons
   - phonon frequencies
   - phonon angular momentum
 - phonons with molecular Berry curvature (Saparov22)
   - deactivated for the moment as to reduce the amount of output (can be reactived by uncommenting in the code)
 - phonons with molecular Berry curvature (Mignolet26)
   - phonon frequencies (correction w.r.t. bare freq.)
   - phonon angular momentum
   - net phonon angular momentum
   - temperature dependence of net phonon angular momentum
   - atomic contribution to each mode
 - bare magnons (Niu98)
   - magnon frequencies
   - magnon eigenvectors
 - phonon-magnon hybrids (Ren24)
   - frequencies
   - phononic and magnonic character (currently broken feature)
 - phonon-magnon hybrids (Mignolet26)
   - frequencies | corr. w.r.t. bare freq | corr. w.r.t. uncoupled freq
   - hybrid freq | bare freq | uncoupled freq
   - phononic, magnonic and coupled character of each hybrid mode
   - phonon angular momentum
   - magnon angular momentum
   - decomposition of the phononic character into atomic contributions

# How to add a material

I did not take the time to implement proper input routines. The materials constants are hardcoded in the code itself (lattice parameters, atomic species, positions and masses, and magnetic atoms). If you want to study a material other than CrI3/CrBr3, add its parameters at the end of the list. You could also implement proper input routines, in which case do not hesitate to open a pull request to contribute back to Hyph-m.


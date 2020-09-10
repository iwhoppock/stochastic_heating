# Stochastic Heating

Test-particle simulation of stochastic heating of particles interacting with a spectrum of Kinetic Alfvén waves. 
This implements the Boris pusher to move the particles through phase space. 

## Overview

See these papers for backgroud information:

* [Perpendicular Ion Heating by Low-Frequency Alfvén-Wave Turbulence in the Solar Wind ](https://arxiv.org/pdf/1001.2069.pdf)
* [Stochastic Proton Heating by Kinetic-Alfvén-Wave Turbulence in Moderately High-β Plasmas](https://arxiv.org/pdf/1811.08873.pdf)

Others include:

* [Perpendicular Ion Heating by Reduced Magnetohydrodynamic Turbulence](https://arxiv.org/pdf/1309.0742.pdf)
* [Radial evolution of stochastic heating in low-β solar wind](https://arxiv.org/pdf/1905.13355.pdf)
* [The Enhancement of Proton Stochastic Heating in the near-Sun Solar Wind](https://arxiv.org/pdf/1912.02653.pdf)
* [Hybrid-Kinetic Simulations of Ion Heating in Alfvénic Turbulence](https://arxiv.org/pdf/1901.11028.pdf)

Lastly, much of the code uses analytical expressions of Kinetic Alfvén Waves, which are phenomenally treated in this 1999 work by Joe Hollweg:

* [Kinetic Alfvén wave revisited](https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/1998JA900132)

### Dependencies

* This runs with OpenMP; this is not configured with MPI, so you must run it on a shared memory system
* OpenMP can be turned off by disabling the command around the particle loop starting, i.e., `#pragma omp parallel for`. This is only recommended for local, small runs. 

## Technical Overview

### Running Locally

For large number of particles, run this on Trillian (Cray XE6m-200):
```
autoreconf -i
./configure CC=cc
make
qsub sh.qsub
```
`sh.qsub` is the file for the PBS (batch scheduler). Locally, this is not a concern, it looks like this:
```
#! /bin/bash
#PBS -l nodes=1:ppn=32
#PBS -j oe

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=32
aprun -d 32 ./sh_sim
```

`configure.ac` can be uncommented in various locations (self explanatory) to include a fortran compiler. This must be noted with the configure command (above).

# RC Kit: ISC and IC

This source distribution supports `SUB=5` (ISC). The CPU and GPU executables use the same Fortran sources. The GPU build additionally compiles `cuda_fastpath.cuf` and enables `USE_CUDA`.

Build on Linux with GNU Make:

```sh
make cpu
make gpu
```

The CPU build requires GNU Fortran, BLAS, and LAPACK. The GPU build requires NVIDIA HPC SDK (`nvfortran`), CUDA, cuBLAS, cuSOLVER, BLAS, and LAPACK. The default GPU architecture is `cc70`; select another architecture using `make gpu GPU_ARCH=cc80` as appropriate for the target GPU.

The executables are `RC_Kit_CPU.exe` and `RC_Kit_GPU.exe`. Object and module files are kept separately in `build/cpu` and `build/gpu`, so both builds can coexist.

To use NVIDIA Fortran for the CPU build:

```sh
make cpu CPU_FC=nvfortran CPU_FLAGS='-O3 -Mpreprocess' CPU_FIXED_FLAGS='-O3 -Mpreprocess -Mextend' CPU_MODFLAGS='-module '
```

To use Intel Fortran for the CPU build:

```sh
make cpu CPU_FC=ifort CPU_FLAGS='-O3 -fpp -extend-source' CPU_FIXED_FLAGS='-O3 -fpp -extend-source' CPU_MODFLAGS='-module '
```

Run from the directory containing the quantum-chemistry input data:

```sh
./RC_Kit_GPU.exe EXAMPLE_ISC_FULL.input EXAMPLE.log
./RC_Kit_CPU.exe EXAMPLE_IC.input EXAMPLE.log
```

The second argument is optional. With a specified output name, the Correlation Function file is named `EXAMPLE.log_TC.TXT`. The program reports its start time, end time, and elapsed time.

Energies and SOC values are in Hartree, temperature is in K, `TIME_STEP` is in ps, and `DAMPING_WIDTH` is in cm^-1. The integration endpoint is `INTEGRAL_STEPS * TIME_STEP`.

Cartesian NAC values are derivative couplings in atomic units, without multiplication by an electronic energy gap. ISC intermediate-state NAC files contain the lower-index-to-higher-index electronic pairs, ordered by increasing intermediate-state index with the target omitted. 
SOC and SOC derivatives must have consistent electronic phases and are used in their supplied orientation.

ISC SOC files include all states of the opposite spin, including the target. Full-SVC SOC-derivative files also include the target block: `3*N_ATOMS * NUM_STATES_SAME_SPIN_FINAL` complex values in `FN_SOC_DERIV_INIT_FINALS`, and `3*N_ATOMS * NUM_STATES_SAME_SPIN_INIT` complex values in `FN_SOC_DERIV_FINAL_INITS`.
Intermediate-state NAC files instead have `3*N_ATOMS * (NUM_STATES_SAME_SPIN_INIT-1)` or `3*N_ATOMS * (NUM_STATES_SAME_SPIN_FINAL-1)` real values.

ISC TC output preserves `SI_BIDIRECTIONAL_V1`: each active component has separate complex positive- and negative-time columns. The energy-gap phase and partition-function normalization are already included. `G_0_0` excludes the SOC modulus squared; the total integrand supplies that coefficient, rate weights, and damping. Integrate `INTEGRAND_TOTAL` over the positive time grid in atomic units without an additional factor of two. IC TC output retains its time, real damped correlation function, and optional damping columns; its reported rate includes both time directions.

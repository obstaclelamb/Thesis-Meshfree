# Thesis-Meshfree
Using meshfree scheme to solve post-surgery soft tissue prediction

# Overview
This repository contains MATLAB code for 3D analysis using the Reproducing Kernel Particle Method (RKPM). It currently includes various examples and functionalities for solving different types of problems.

# Features

## Solvers:

* Linear Patch Test
* Beam Problem
* Plate Problem
* Poisson Problem
* Tensile Test
* FEA Model

## Core RKPM Functions:

* Shape function generation
* Numerical integration
* Matrix assembly
* Boundary condition implementation
* Post-processing and visualization

## Utilities:
* Input file handling
* Convergence analysis
* Progress bar display
* Requirements
* MATLAB (tested with R2023a)

# Usage

## Clone the repository:

```{Bash}
git clone https://github.com/obstaclelamb/Thesis-Meshfree.git
```

## Navigate to the desired example directory:

```{Bash}
cd RKPM2D/01_LinearPatchTest
```

## Run the `MAIN.m` script:

```{Matlab}
MAIN
```

Examples
* 01_LinearPatchTest: Verifies the basic implementation of RKPM by performing a linear patch test.
* 02_BeamProblem: Analyzes the behavior of a beam under different loading conditions.
* 03_PlateProblem: Solves plate bending problems using RKPM.
* 04_PoissonProblem_src: Addresses the Poisson equation in 2D.
* 05_FEA_Model: Provides a framework for general finite element analysis using RKPM.
* 06_TensileTest: Simulates a tensile test on a specimen.
* PlotRKShapeFunction: Generates plots of RKPM shape functions.

# Customization

Input Files: Modify the getInput.m file in each example directory to adjust problem parameters, material properties, boundary conditions, and discretization settings.
Core Functions: Customize the core RKPM functions in the MAIN_PROGRAM directory to implement different formulations, kernel functions, or integration schemes.
Contributing

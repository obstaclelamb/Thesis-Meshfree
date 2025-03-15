# Thesis-Meshfree
Using meshfree scheme to solve post-surgery soft tissue prediction

# Overview
This repository contains MATLAB code for 3D analysis using the Reproducing Kernel Particle Method (RKPM). It currently includes various examples and functionalities for solving different types of problems.

# Usage

## Clone the repository:

```{Bash}
git clone https://github.com/obstaclelamb/Thesis-Meshfree.git
```

## Open MATLAB

Directly double click the `MAIN.m` is recommanded, since browse the file through the MATLAB terminal is kind of troublesome.

## Run the `MAIN.m` script:

Type the following in the terminal (when `MAIN.m` is in the path (current folder)):

```{Matlab}
MAIN
%  XD
```
or just press `F5` when you looking at `MAIN.m`. This file is already contain the **improtant 3 function**.

## Run the demostration of RKPM 1-dimensional explaination (when `MAIN.m` is in the path (current folder)):

```
addpath('RK1Dcode_modified')
RK1D_demo
```

## Run other

Examples (not included)
* 01_LinearPatchTest: Verifies the basic implementation of RKPM by performing a linear patch test.
* 02_BeamProblem: Analyzes the behavior of a beam under different loading conditions.
* 03_PlateProblem: Solves plate bending problems using RKPM.
* 04_PoissonProblem_src: Addresses the Poisson equation in 2D.
* 05_FEA_Model: Provides a framework for general finite element analysis using RKPM.
* 06_TensileTest: Simulates a tensile test on a specimen.
* PlotRKShapeFunction: Generates plots of RKPM shape functions.

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



# Output

## Domain setup

![](https://github.com/obstaclelamb/Thesis-Meshfree/blob/main/fig/NodalRepresentativeDomain.jpg)
![](https://github.com/obstaclelamb/Thesis-Meshfree/blob/main/fig/SupportandNeighbors.jpg)

## Solutions

![Deformed](https://github.com/obstaclelamb/Thesis-Meshfree/blob/main/fig/DeformedConfiguration.jpg)

### displacement

![](https://github.com/obstaclelamb/Thesis-Meshfree/blob/main/fig/u1.jpg)
![](https://github.com/obstaclelamb/Thesis-Meshfree/blob/main/fig/u2.jpg)


### stress

![](https://github.com/obstaclelamb/Thesis-Meshfree/blob/main/fig/sigma11.jpg)
![](https://github.com/obstaclelamb/Thesis-Meshfree/blob/main/fig/sigma12.jpg)
![](https://github.com/obstaclelamb/Thesis-Meshfree/blob/main/fig/sigma22.jpg)

### strain

![](https://github.com/obstaclelamb/Thesis-Meshfree/blob/main/fig/epsilon11.jpg)
![](https://github.com/obstaclelamb/Thesis-Meshfree/blob/main/fig/epsilon12.jpg)
![](https://github.com/obstaclelamb/Thesis-Meshfree/blob/main/fig/epsilon22.jpg)


# Customization

Input Files: Modify the getInput.m file in each example directory to adjust problem parameters, material properties, boundary conditions, and discretization settings.
Core Functions: Customize the core RKPM functions in the MAIN_PROGRAM directory to implement different formulations, kernel functions, or integration schemes.
Contributing

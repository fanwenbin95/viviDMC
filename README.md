# viviDMC

A fortran program `viviDMC` for diffusion Monte Carlo (DMC) 
using vibrational initialization. 

- Author : Wenbin FAN (wbfan21@m.fudan.edu.cn)
- Advisors : Prof. Donghui ZHANG, Prof. Xin XU @FDU

## Function

- sampling from vibrational analysis
- traditional propagation without weights

## Prerequisites

Intel OneAPI Fortran Compiler, Intel MKL

## Usage

1. expose your PES to `interface.f90`
2. provide one or more initial geometries on `EQ.xyz`
3. modify the parameter in `desk.inp`
4. make and run this program

# IO

## output

The `Eref` will printed to screen. 
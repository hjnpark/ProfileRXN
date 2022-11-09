Starting with two molecular structure (reactant and product), ProfileRXN will first optimize them and interpolate a trajectory between the optimized geometries. NEB calculation will then be used to guide the trajectory towards the minimum energy pathway.
Finally, a guessed transition state structure from the NEB calculation will be optimized and energies of reactant, product, and transition state will be written as a json file.

Author: Heejune Park

Contact Email: heepark@ucdavis.edu

## Quick Set Up
### 1. Creating a conda environment and installing ProfileRXN

Commands below will create a conda environment and install ProfileRXN.
 ```shell
git clone git@github.com:hjnpark/ProfileRXN
cd ProfileRXN
conda update conda
conda env create -f ProfileRXN_env.yaml
conda activate ProfileRXN
pip install -e .
```
Setting up the environment might take some time.

### 2. Installing QCFractal and geomeTRIC

Install `next` branch of QCFractal:
[https://github.com/MolSSI/QCFractal/tree/next](https://github.com/MolSSI/QCFractal/tree/next)

Following commands will install the QCFractal's `next` branch.

```shell
git clone -b next git@github.com:MolSSI/QCFractal.git
cd QCFractal
pip install -e qcportal/ -e qcfractalcompute/ -e qcfractal/
```

Install `neb` branch of geomeTRIC:
[https://github.com/hjnpark/geomeTRIC/tree/neb](https://github.com/hjnpark/geomeTRIC/tree/neb)

```shell
git clone -b neb git@github.com:hjnpark/geomeTRIC.git
cd geomeTRIC
pip install -e .
```

Done!

## User Guide

`profile-rxn -h` will show all the available arguments. Input xyz file should have two geometries. If it has more than two, the first and last frame will be used.
`exmaples/bash.sh` shows an example command.

## TO DO
1. Perform vibrational frequency analysis and calculate Gibbs free energies.
2. IRC step needs to be added.
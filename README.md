isoRMSD  
=======
- isoRMSD can calculate RMSD between 2 conformers with different atom names.
- The atom order of molecules also don't need to be same.
- Hydration or not is OK.
- This script is for small molecule only. It will be very slow for protein.

This script is base on RDKit cookbook -- [RMSD Calculation between N molecules](http://www.rdkit.org/docs/Cookbook.html).


Howto
-----
You need RDKit to run isoRMSD.py.
RDKit is a collection of cheminformatics and machine-learning software written in C++ and Python.
The easiest way to installing RDKit is using the excellent conda package manager in Anaconda python.  
And then get RDKit for running isoRMSD.py.  

### 1. Download and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.continuum.io/)
Install command is simple!
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# All enter 'yes' is OK.
```
### 2. Install [RDkit](http://rdkit.org/)
```
conda create -n rmsd -c conda-forge rdkit
conda activate rmsd

# Or conda install -c rdkit rdkit
```
### 3. Usage
```
cd example
python ../isoRMSD.py mol1.pdb mol2.pdb rmsd.txt
# Best_RMSD: 0.666
# Best_Not_Fit_RMSD: 1.703
```
isoRMSD.py will output two RMSD, one is fitted, another is no fit.  
Not fit  RMSD mean no change in molecules coordinates.  

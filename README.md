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

1. Download and install a new python from [Anaconda](https://www.continuum.io/). python 2 or 3 is OK
Install command is simple!
```
bash Anaconda2-4.1.1-Linux-x86_64.sh
# All enter 'yes' is OK.
```
2. Download and install new python package [RDkit](http://rdkit.org/).
```
conda install -c rdkit rdkit
```
3. Usage:
```
python isoRMSD.py mol1.pdb mol2.pdb rmsd.txt
```
isoRMSD.py will output two RMSD, one is fitted, another is no fit.  
Not fit  RMSD mean no change in molecules coordinates.  

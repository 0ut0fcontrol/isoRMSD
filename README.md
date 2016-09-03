*isoRMSD*  
========
- isoRMSD can calculates RMSD between 2 conformers with different atom names.
- The atom order of molecules also don't need to be same.
- Hydration or not is OK.
- This script is for small molecule only. It will be very slow for protein.

This script base on RDKit cookbook -- [RMSD Calculation between N molecules](http://www.rdkit.org/docs/Cookbook.html)



Howto
-----
You need rdkit python module to run isoRMSD.py.  
The easiest way to get rdkit is to install Anaconda python.  
And then get rdkit for run isoRMSD.py.  

1. Download and install a new python form [Anaconda](https://www.continuum.io/). python 2 or 3 is OK
install command is simple!
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

# You need rdkit python module to run isoRMSD.py
# The easier way to get rdkit is to install Anaconda python 
# And then get rdkit for run isoRMSD.py



# 1.Download and install a new python 2.7 [Anaconda Version](https://www.continuum.io/).

bash Anaconda2-4.1.1-Linux-x86_64.sh
# All enter 'yes' is OK.

# 2. Download and install new python package RDkit from Anaconda User rdkit.
# [RDKit: Open-Source Cheminformatics Software](http://rdkit.org/)

conda install -c rdkit rdkit

# 3. Usage:

python isoRMSD.py mol1.pdb mol2.pdb rmsd.txt

# isoRMSD.py will output two RMSD, one is fitted, another is no fit.
# not fit  RMSD mean no change in molecules coordinates.  

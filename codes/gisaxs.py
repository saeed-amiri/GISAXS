"""
Modeling GISAXS experiments for a 2 dimentional model.
First:
    Read the LAMMPS datafile
    calculate the poits in the reciprocal space
    ...
"""
import sys
import read_lmp_data as relamp
import reciprocal as recip

fname: str = sys.argv[1]  # LAMMPS datafile to work with

data = relamp.ReadData(fname)
qspace = recip.MakeReciprocal(data.Atoms_df)

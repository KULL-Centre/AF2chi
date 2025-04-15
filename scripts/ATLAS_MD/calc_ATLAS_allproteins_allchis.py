import os
import numpy as np
import mdtraj as md
import json
import glob

def save_json(filename, object):
    with open(filename, 'w') as f:
        json.dump(object, f)

no_chi1_restype = ['GLY', 'ALA', 'PRO']
output_dir = 'allproteins_ATLAS_allchis'
ATLAS_pdb_file = 'pdb_list.dat'
n_chi_aa = {'CYS': 1,
 'ASP': 2,
 'GLU': 3,
 'PHE': 2,
 'HIS': 2,
 'ILE': 2,
 'LYS': 4,
 'LEU': 2,
 'MET': 3,
 'ASN': 2,
 'PRO': 1,
 'GLN': 3,
 'ARG': 4,
 'SER': 1,
 'THR': 1,
 'VAL': 1,
 'TRP': 2,
 'TYR': 2}

os.system(f'mkdir {output_dir}')
pdbs = np.genfromtxt(ATLAS_pdb_file, dtype=str)
for pdb in pdbs:

    print(f'Starting {pdb}')

    #download and unzip traj from ATLAS
    os.system(f'wget https://www.dsimb.inserm.fr/ATLAS/database/ATLAS/{pdb}/{pdb}_protein.zip --quiet')
    os.system(f'unzip -d {pdb}_tmp {pdb}_protein.zip ')
    
    #calculate chi1 angles for all three replica simulations and make dict with results
    chi1_residues = {}
    chi2_residues = {}
    chi3_residues = {}
    chi4_residues = {}

    for replica in range(1,4):
        ensemble = md.load(glob.glob(f'{pdb}_tmp/{pdb}_prod_R{replica}_fit.xtc')[0], top=glob.glob(f'{pdb}_tmp/{pdb}.pdb')[0])
        for residue in ensemble.top.residues:
            restype = str(residue)[:3]
            resi = int(str(residue)[3:])-1
            if restype not in no_chi1_restype:
                resatoms = ensemble.top.select(f'resi {resi} and resname {restype}')
                res_slice = ensemble.atom_slice(resatoms)
                chi1 = md.compute_chi1(res_slice)[1]
                if replica == 1:
                    chi1_residues[str(residue)] = np.reshape(chi1, len(chi1)).tolist()
                else:
                    chi1_residues[str(residue)] = [*chi1_residues[str(residue)], *np.reshape(chi1, len(chi1)).tolist()]

                if n_chi_aa[restype] >= 2:
                    chi2 = md.compute_chi2(res_slice)[1]
                    if replica == 1:
                        chi2_residues[str(residue)] = np.reshape(chi2, len(chi2)).tolist()
                    else:
                        chi2_residues[str(residue)] = [*chi2_residues[str(residue)], *np.reshape(chi2, len(chi2)).tolist()]

                if n_chi_aa[restype] >=3:
                    chi3 = md.compute_chi3(res_slice)[1]
                    if replica == 1:
                        chi3_residues[str(residue)] = np.reshape(chi3, len(chi3)).tolist()
                    else:
                        chi3_residues[str(residue)] = [*chi3_residues[str(residue)], *np.reshape(chi3, len(chi3)).tolist()]

                if n_chi_aa[restype] >=4:
                    chi4 = md.compute_chi4(res_slice)[1]
                    if replica == 1:
                        chi4_residues[str(residue)] = np.reshape(chi4, len(chi4)).tolist()
                    else:
                        chi4_residues[str(residue)] = [*chi4_residues[str(residue)], *np.reshape(chi4, len(chi4)).tolist()]

    all_chis = {}

    all_chis['chi1'] = chi1_residues
    all_chis['chi2'] = chi2_residues
    all_chis['chi3'] = chi3_residues
    all_chis['chi4'] = chi4_residues

    save_json(f'{output_dir}/{pdb}.json', all_chis)

    #remove temporary files
    os.system(f'rm -r {pdb}_tmp {pdb}_protein.zip')


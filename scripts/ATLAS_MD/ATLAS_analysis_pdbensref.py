import os
import glob
import numpy as np
import pandas as pd
from scipy.stats import circmean
from scipy.spatial.distance import jensenshannon
import json
import mdtraj as md
import sys

pdbens_files = glob.glob('ATLAS_pdbensembles/*.json')
uniprot_ids = []
for file in pdbens_files:
    uniprot_ids.append(file[19:25])

ATLAS_info_file = '2022_06_13_ATLAS_info.tsv'

pdbs = []
pdb2uniprot= {}
ATLAS_info = np.genfromtxt(ATLAS_info_file, dtype=str, skip_header=1, usecols=(0,5))
for uniprot_id in uniprot_ids:
    for info in ATLAS_info:
        if uniprot_id == info[1]:
            pdb = info[0]
            pdbs.append(pdb)
            pdb2uniprot[pdb] = uniprot_id

x=[i for i in range(5,360,10)]
x_range=[i for i in range(0,365,10)]

chi_type = str(sys.argv[1])  

def adjust_period_angle(angle):
    while angle < 0:
        angle += 360
    while angle >= 360:
        angle -= 360
    return angle

def save_json(filename, object):
    with open(filename, 'w') as f:
        json.dump(object, f)

def load_json(file):
    with open(file, 'r') as f:
        loaded_json = json.load(f)
    return loaded_json

def rot_pop_rebinned_to_three(input_distr):
    assert len(input_distr) == 36, "not correct length"
    w = np.array([np.sum(input_distr[0:12]),np.sum(input_distr[12:24]),np.sum(input_distr[24:36])])
    w /= np.sum(w)

    #Set very small weights to 0 to avoid numerical errors
    w[w<1e-8] = 0.0
    w /= np.sum(w)
    
    return w

top8000_prior=pd.read_csv('../extra_data/Top8000_rebinned_all_chi_distributions.csv')

subfolders = [ f.path for f in os.scandir('af2chi_runs/') if f.is_dir()]

js_af2_proteins = {}
js_top8000_proteins = {}
js_prior_proteins = {}
js_re_proteins = {}
js_md_proteins = {}

js_af2_fixed_proteins = {}
js_top8000_fixed_proteins = {}
js_prior_fixed_proteins = {}
js_re_fixed_proteins = {}
js_md_fixed_proteins = {}

js_af2_multi_proteins = {}
js_top8000_multi_proteins = {}
js_prior_multi_proteins = {}
js_re_multi_proteins = {}
js_md_multi_proteins = {}

for dir in subfolders:
    pdb = dir[-6:]

    if pdb not in pdbs:
        continue

    print(f'Starting {dir}')

    reweighted_distributions_allchis = load_json(glob.glob(f'{dir}/output/*_rank_001_sc_pops_fitted.json')[0])
    prior_distributions_allchis = load_json(glob.glob(f'{dir}/output/*_rank_001_sc_pops_priori.json')[0])
    md_ensemble_allchis = load_json(f'allproteins_ATLAS_allchis_jsons/{pdb}.json')
    pdb_ensemble_allchis = load_json(f'ATLAS_pdbensembles/{pdb2uniprot[pdb]}_pdb_chis_pops.json')

    #Get af2 structure angles
    af2_structure = md.load(glob.glob(f'{dir}/output/*unrelaxed_rank_001_alphafold2*.pdb')[0])
    chi1_residues = {}
    chi2_residues = {}
    chi3_residues = {}
    chi4_residues = {}

    for residue in af2_structure.top.residues:

        restype = str(residue)[:3]
        resi = int(str(residue)[3:])-1
        resatoms = af2_structure.top.select(f'resi {resi} and resname {restype}')
        res_slice = af2_structure.atom_slice(resatoms)
        chi1 = md.compute_chi1(res_slice)[1][0]
        chi2 = md.compute_chi2(res_slice)[1][0]
        chi3 = md.compute_chi3(res_slice)[1][0]
        chi4 = md.compute_chi4(res_slice)[1][0]

        if len(chi1) != 0:
            chi1_residues[str(residue)] = chi1[0]
        if len(chi2) != 0:
            chi2_residues[str(residue)] = chi2[0]
        if len(chi3) != 0:
            chi3_residues[str(residue)] = chi3[0]
        if len(chi4) != 0:
            chi4_residues[str(residue)] = chi4[0]

    af2_structures_allchis = {}
    af2_structures_allchis['chi1'] = chi1_residues
    af2_structures_allchis['chi2'] = chi2_residues
    af2_structures_allchis['chi3'] = chi3_residues
    af2_structures_allchis['chi4'] = chi4_residues

    reweighted_distributions = reweighted_distributions_allchis[chi_type]
    prior_distributions = prior_distributions_allchis[chi_type]
    af2_structures_angles = af2_structures_allchis[chi_type]
    md_ensemble_angles = md_ensemble_allchis[chi_type]
    pdb_ensemble_angles = pdb_ensemble_allchis[chi_type]

    md_ensemble_angles_distr={}
    for key in md_ensemble_angles.keys():

        #print('md angles', md_ensemble_angles[key])

        md_degr_conv = [np.rad2deg(circmean(ang,np.deg2rad(-5),np.deg2rad(355))) for ang in md_ensemble_angles[key]]
        #print('md angles conv', md_degr_conv)

        md_ensemble_angles_distr[key]=np.histogram(md_degr_conv,bins=x_range)[0]/len(md_degr_conv)

    pdb_ensemble_distributions={}
    for key in pdb_ensemble_angles.keys():
       
        #Filter out residues with no structures
        if len(pdb_ensemble_angles[key]) >= 10:

            #print('pdbens angles', pdb_ensemble_angles[key])

            pdbens_degr_conv = [np.rad2deg(circmean(ang,np.deg2rad(-5),np.deg2rad(355))) for ang in pdb_ensemble_angles[key]]
            #print('pdbens angles conv', pdbens_degr_conv)

            pdb_ensemble_distributions[key]=np.histogram(pdbens_degr_conv,bins=x_range)[0]/len(pdbens_degr_conv)
            
    af2_distributions={}
    for residue in reweighted_distributions.keys():
        restype = residue[:3]
         
        #print('af struc angles', af2_structures_angles[residue])

        angle_structure = np.rad2deg(af2_structures_angles[residue])#[0] ##chi1 HAVE TO been in degrees not radians
        angle_structure = adjust_period_angle(angle_structure)
        
        #print('af struc angles conv', angle_structure)

        af2_distributions[residue] = np.histogram(angle_structure,bins=x_range)[0]

    js_re=[]
    js_prior=[]
    js_top8000=[]
    js_af2=[]
    js_md=[]
    i=0
    for idx,res in enumerate(reweighted_distributions.keys()):
        if res in pdb_ensemble_distributions.keys() and res in md_ensemble_angles_distr.keys():

            js_re.append(jensenshannon(reweighted_distributions[res],pdb_ensemble_distributions[res]))
            js_prior.append(jensenshannon(prior_distributions[res],pdb_ensemble_distributions[res]))
            js_top8000.append(jensenshannon(top8000_prior[f'{res[:3]}_{chi_type}'],pdb_ensemble_distributions[res]))
            js_af2.append(jensenshannon(af2_distributions[res],pdb_ensemble_distributions[res]))
            js_md.append(jensenshannon(md_ensemble_angles_distr[res],pdb_ensemble_distributions[res]))

    is_fixed={}
    
    for idx,res in enumerate(reweighted_distributions.keys()):
        if res in pdb_ensemble_distributions.keys() and res in md_ensemble_angles_distr.keys():
            red_pop=rot_pop_rebinned_to_three(pdb_ensemble_distributions[res])
            if np.count_nonzero(red_pop)>1:
                is_fixed[res]=False
            else:
                is_fixed[res]=True
                                            
    js_re_fixed=[]
    js_prior_fixed=[]
    js_top8000_fixed=[]
    js_af2_fixed=[]
    js_md_fixed=[]
    
    js_re_multi=[]
    js_prior_multi=[]
    js_top8000_multi=[]
    js_af2_multi=[]
    js_md_multi=[]
    
    for idx,res in enumerate(reweighted_distributions.keys()):
        if res in pdb_ensemble_distributions.keys() and res in md_ensemble_angles_distr.keys():
            if is_fixed[res]:
                js_re_fixed.append(jensenshannon(reweighted_distributions[res],pdb_ensemble_distributions[res]))
                js_prior_fixed.append(jensenshannon(prior_distributions[res],pdb_ensemble_distributions[res]))
                js_top8000_fixed.append(jensenshannon(top8000_prior[f'{res[:3]}_{chi_type}'],pdb_ensemble_distributions[res]))
                js_af2_fixed.append(jensenshannon(af2_distributions[res],pdb_ensemble_distributions[res]))
                js_md_fixed.append(jensenshannon(md_ensemble_angles_distr[res],pdb_ensemble_distributions[res]))
            else:
                js_re_multi.append(jensenshannon(reweighted_distributions[res],pdb_ensemble_distributions[res]))
                js_prior_multi.append(jensenshannon(prior_distributions[res],pdb_ensemble_distributions[res]))
                js_top8000_multi.append(jensenshannon(top8000_prior[f'{res[:3]}_{chi_type}'],pdb_ensemble_distributions[res]))
                js_af2_multi.append(jensenshannon(af2_distributions[res],pdb_ensemble_distributions[res]))
                js_md_multi.append(jensenshannon(md_ensemble_angles_distr[res],pdb_ensemble_distributions[res]))

    js_af2_proteins[pdb] = js_af2
    js_top8000_proteins[pdb] = js_top8000
    js_prior_proteins[pdb] = js_prior
    js_re_proteins[pdb] = js_re
    js_md_proteins[pdb] = js_md

    js_af2_fixed_proteins[pdb] = js_af2_fixed
    js_top8000_fixed_proteins[pdb] = js_top8000_fixed
    js_prior_fixed_proteins[pdb] = js_prior_fixed
    js_re_fixed_proteins[pdb] = js_re_fixed
    js_md_fixed_proteins[pdb] = js_md_fixed

    js_af2_multi_proteins[pdb] = js_af2_multi
    js_top8000_multi_proteins[pdb] = js_top8000_multi
    js_prior_multi_proteins[pdb] = js_prior_multi
    js_re_multi_proteins[pdb] = js_re_multi
    js_md_multi_proteins[pdb] = js_md_multi

js_results_dict = {
    'af2':js_af2_proteins,
    'top8000':js_top8000_proteins,
    'prior':js_prior_proteins,
    'reweighted':js_re_proteins,
    'md':js_md_proteins 
}
save_json(f'output_js_results/ATLAS_proteins_js_pdbensemble_results_{chi_type}.json', js_results_dict)

js_fixed_results_dict = {
    'af2':js_af2_fixed_proteins,
    'top8000':js_top8000_fixed_proteins,
    'prior':js_prior_fixed_proteins,
    'reweighted':js_re_fixed_proteins,
    'md':js_md_fixed_proteins  
}
save_json(f'output_js_results/ATLAS_proteins_js_pdbensemble_fixedresis_results_{chi_type}.json', js_fixed_results_dict)

js_multi_results_dict = {
    'af2':js_af2_multi_proteins,
    'top8000':js_top8000_multi_proteins,
    'prior':js_prior_multi_proteins,
    'reweighted':js_re_multi_proteins,
    'md':js_md_multi_proteins
}
save_json(f'output_js_results/ATLAS_proteins_js_pdbensemble_multiresis_results_{chi_type}.json', js_multi_results_dict)


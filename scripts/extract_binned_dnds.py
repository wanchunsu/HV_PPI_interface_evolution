# Extract binned dN/dS values from the codeml outputs

import argparse
import json
import os.path as osp
import os

def extract_dnds(res_fi_with_path):
	with open(res_fi_with_path) as rf:
		for line in rf:
			if line.startswith('omega (dN/dS) =  '):
				dn_ds = float(line.split()[-1])
				#print(res_fi_with_path, dn_ds)
	return dn_ds

def extract_tstv(res_fi_with_path):
	with open(res_fi_with_path) as rf:
		for line in rf:
			if line.startswith('kappa (ts/tv) =  '):
				ts_tv = float(line.split()[-1])
				#print(res_fi_with_path, dn_ds)
	return ts_tv

def extract_SE_for_binned_dnds(res_fi_with_path): #get the standard error
	with open(res_fi_with_path) as rf:
		for line in rf:
			if line.startswith("SEs for parameters:"):
				line_of_SEs = next(rf)
				SE_for_dnds = float(line_of_SEs.strip("\n").split()[-1])
	return SE_for_dnds
            

def main():
	script_dir = osp.dirname(__file__)
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', '--file_type') #codeml_unrooted or codeml_unrooted_ds_neq_0
	args = parser.parse_args()

	fi_type = args.file_type

	binned_dnds_folder = osp.join(script_dir, '..', 'data', 'binned_dnds')

	residue_files = [fi for fi in os.listdir(binned_dnds_folder) if fi.endswith(fi_type + '.out')]

	dnds_extracted_fi = osp.join(script_dir, '..', 'data', 'binned_dnds', fi_type + '.binned_dnds.json')

	dnds_SE_extracted_fi = osp.join(script_dir, '..', 'data', 'binned_dnds', fi_type + '.binned_dnds_SEs.json')

	tstv_extracted_fi = osp.join(script_dir, '..', 'data', 'binned_dnds', fi_type + '.binned_tstv.json')

	# rate of nonsynonymous/synonymous substitutions
	dict_of_dnds = {'exo': 0.0, 'endo': 0.0, 'mimicry': 0.0, 'surface': 0.0, 'buried': 0.0, 'all_exo': 0.0, 'all_endo': 0.0} #initialize
	for res_fi in residue_files:
		residue_type = res_fi.split(".")[0]
		res_fi_with_path = osp.join(binned_dnds_folder, res_fi)
		dnds_val = extract_dnds(res_fi_with_path)
		dict_of_dnds[residue_type] = dnds_val

	# SE for dN/dS values
	dict_of_dnds_SEs = {'exo': 0.0, 'endo': 0.0, 'mimicry': 0.0, 'surface': 0.0, 'buried': 0.0, 'all_exo': 0.0, 'all_endo': 0.0} #initialize
	for res_fi in residue_files:
		residue_type = res_fi.split(".")[0]
		res_fi_with_path = osp.join(binned_dnds_folder, res_fi)
		SE_val = extract_SE_for_binned_dnds(res_fi_with_path)
		dict_of_dnds_SEs[residue_type] = SE_val

	#rate of transitions/transversions
	dict_of_tstv = {'exo': 0.0, 'endo': 0.0, 'mimicry': 0.0, 'surface': 0.0, 'buried': 0.0, 'all_exo': 0.0, 'all_endo': 0.0} #initialize
	for res_fi in residue_files:
		residue_type = res_fi.split(".")[0]
		res_fi_with_path = osp.join(binned_dnds_folder, res_fi)
		tstv_val = extract_tstv(res_fi_with_path)
		dict_of_tstv[residue_type] = tstv_val

	with open(dnds_extracted_fi, 'w') as dnds_o:
		json.dump(dict_of_dnds, dnds_o)

	with open(dnds_SE_extracted_fi, 'w') as dnds_se_o:
		json.dump(dict_of_dnds_SEs, dnds_se_o)

	with open(tstv_extracted_fi, 'w') as tstv_o:
		json.dump(dict_of_tstv, tstv_o)


if __name__ == '__main__':
	main()
import sys
from datetime import datetime
from multiprocessing import Pool
import multiprocessing as multi
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
rdBase.DisableLog('rdApp.warning')
	
def calc_desclist( smi ):
	# MW, tPSA, AROM, HBA, HBD, ROTB, LogP, RING
	desclist = ['MolWt', 'TPSA', 'NumAromaticRings', 'NumHAcceptors',\
				 'NumHDonors', 'NumRotatableBonds', 'MolLogP', 'RingCount']
	calculator = MoleculeDescriptors.MolecularDescriptorCalculator( desclist )

	#Calculation
	res = []
	try:
		mol = Chem.MolFromSmiles(smi)
		descs = calculator.CalcDescriptors(mol)
		for num in range(0, len(desclist)):
			if descs[num]:
				res.append(str(descs[num]))
			else:
				res.append("")
		resrow = "\t".join(res)
	except Exception:
		ex, ms, tb = sys.exc_info()
		sys.stderr.write("ERROR : " + str(ms) + "\n")

	return resrow

def run_calc_multi(inpfile, num_cpus=multi.cpu_count()):
	with open(inpfile, 'r') as inF:
		print("Calculation started")
		counter = 0
		smilist = []
		for l in inF:
			if counter != 0:
				smi = l.rstrip()
				smilist.append(smi)
			counter += 1

		print("Compounds processing :", len(smilist))

		# Multiprocessing
		p = Pool(num_cpus)
		result = p.map(calc_desclist, smilist)
		p.close()

		# Output preparation
		pref = inpfile.split(".")[0]
		with open(pref + '_calc.tsv', 'w') as out:
			out.write("SMILES	MW	TPSA	AROM	HBA	" +\
				"HBD	ROTB	LOGP	RING\n")
			
			for smi, line in zip(smilist, result):
				out.write(smi + "\t" + line + "\n")


if __name__ == "__main__":
	inF = sys.argv[1]
	print(multi.cpu_count())
	if inF:
		for numcpu in [4, 2, 1]:
			stime = datetime.now()
			run_calc_multi(inF, numcpu)
			etime = datetime.now()
			delta = etime - stime	
			print("Num of CPUs: {}, calculation time: {}".format(numcpu, delta))
	else:
		print("Please specify the input file.")

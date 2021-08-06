import sys, os
from Gff import GffGenes
from collections import OrderedDict
from Bio import SeqIO

def classify_fold(inGff3, inSeq, outPre=None):
	if outPre is None:
		outPre = os.path.basename(inSeq)
	d_seqs = {rc.id:rc.seq for rc in SeqIO.parse(inSeq, 'fasta')}
	d_position = OrderedDict()
	d_count ={}
	d_handles = {}

	for rc in GffGenes(inGff3):
		rc.classify_fold(d_seqs, d_position)
	all_type = 'ALL'
	d_count = {all_type: 0}
	d_handles = {}
		
	for pos, types in d_position.iteritems():
		for type in types:
			if type not in d_handles:
				outfile = '{}.{}.pos'.format(outPre, type)
				d_handles[type] = open(outfile, 'w')
			handle = d_handles[type]
			pos.write(handle)
			try: d_count[type] += 1
			except KeyError: d_count[type] = 1
		d_count[all_type] += 1
	for type,count in sorted(d_count.items()):
		print '{}\t{}'.format(type,count)
	for handle in d_handles.values():
		handle.close()

def main():
	classify_fold(inGff3=sys.argv[1], inSeq=sys.argv[2])

if __name__ == '__main__':
	main()

#from collections import namedtuple

#Region = namedtuple('Region', ['chrom', 'start', 'end'])	# 1-based
#Position = namedtuple('Position', ['chrom', 'pos'])	# 1-based
class Region():	# 1-based
	def __init__(self, chrom=None, start=None, end=None):	# 1-based
		self.chrom = chrom
		if start >end:
			start, end = end, start
		self.start = start if start is None else int(start)
		self.end   = end if end is None else int(end)
	def __str__(self):
		return '{chrom}:{start}-{end}'.format(chrom=self.chrom, start=self.start, end=self.end)
	def __repr__(self):
		return self.__format__("Region(chrom='{chrom}', start={start}, end={end})")
	def __format__(self, spec=None):
		if spec is None:
			spec = '{chrom}\t{start}\t{end}'
		return spec.format(chrom=self.chrom, pos=self.pos)

	def __eq__(self, other):
		if self.key == other.key:
			return True
		return False
	def __hash__(self):
		return hash(self.key)
	def __len__(self):
		return self.end - self.start + 1
	@property
	def key(self):
		return (self.chrom, self.start, self.end)
	@property
	def positions(self):
		return [Position(chrom=self.chrom, pos=pos) for pos in range(self.start, self.end+1)]
	def overlaps(self, other):
		return max(0, min(self.end, other.end)-max(self.start, other.start))
class Regions():	# parser
	def __init__(self):
		pass
	@classmethod
	def subtract(cls, regions1, regions2):	# regions1 >= regions2
		subtracts = set(regions1) - set(regions2)
		subtracts2 = set(regions2) - set(regions1)
		if not subtracts2:
			return cls.sort_regions(subtracts)
		subsubtracts = set([])
		for region1 in subtracts:
			for region2 in subtracts2:
				if cls.overlap(region1, region2):
					if region2.start > region1.start:
						subsubtracts.add(Region(start=region1.start, end=region2.start-1))
					if region1.end > region2.end:
						subsubtracts.add(Region(start=region2.end+1, end=region1.end))
					break
			else:
				subsubtracts.add(region1)
		return cls.sort_regions(subsubtracts)

	@classmethod
	def overlap(cls, region1, region2):
		return max(0, min(region1.end, region2.end)-max(region1.start, region2.start))
	@classmethod
	def sort_regions(cls, regions):
		return sorted(regions, key=lambda x: x.start)
	@classmethod
	def get_introns(cls, exon_regions):
		intron_regions = []
		exon_regions = sorted(exon_regions, key=lambda x:x.start)
		last_end = exon_regions[0].end
		for exon_region in exon_regions[1:]:
			intron_regions += [Region(chrom=exon_region.chrom, start=last_end+1, end=exon_region.start-1) ]
			last_end = exon_region.end
		return intron_regions
	@classmethod
	def get_splicing(cls, exon_regions, span=2):
		splicing_regions = []
		exon_regions = sorted(exon_regions, key=lambda x:x.start)
		for i, exon_region in enumerate(exon_regions):
			if  i < len(exon_regions) - 1:	# end
				splicing_regions += [Region(chrom=exon_region.chrom, 
						start=exon_region.end-span, end=exon_region.end+span)]
			if i > 0:	# start
				splicing_regions += [Region(chrom=exon_region.chrom,
                        start=exon_region.start-span, end=exon_region.start+span)]
		return splicing_regions	# regions

class Position():
	def __init__(self, chrom=None, pos=None):
		self.chrom = chrom
		self.pos = pos if pos is None else int(pos)
	def __add__(self, number):
		return Position(chrom=self.chrom, pos=self.pos+number)
	def __sub__(self, other):
		if self.chrom != other.chrom:
			return None
		return self.pos - other.pos
	def __lt__(self, other):
		if self.chrom == other.chrom and self.pos < other.pos:
			return True
		return False
	def __le__(self, other):
		if self.__lt__(other) or self.__eq__(other):
			return True
		return False
	def __gt__(self, other):
		if self.chrom == other.chrom and self.pos > other.pos:
			return True
		return False
	def __ge__(self, other):
		if self.__gt__(other) or self.__eq__(other):
			return True
		return False
	def __eq__(self, other):
		if self.key == other.key:
			return True
		return False
	def __hash__(self):
		return hash(self.key)
	def __str__(self):
		return '{chrom}:{pos}'.format(chrom=self.chrom, pos=self.pos)
	def __repr__(self):
		return self.__format__("Position(chrom='{chrom}', pos={pos})")
	def __format__(self, spec=None):
		if spec is None:
			spec = '{chrom}\t{pos}'
		elif spec == "":  # format_spec defaults to ""
			return str(self)
		return spec.format(chrom=self.chrom, pos=self.pos)
	@property
	def key(self):
		return (self.chrom, self.pos)
	def write(self, fout):
		print >> fout, self.__format__()

class Positions():	# parser
	def __init__(self, posfile=None, chrom_col=0, pos_col=1, based=1):
		self.posfile = posfile
		self.chrom_col = chrom_col
		self.pos_col = pos_col
		self.based = based
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in open(self.posfile):
			temp = line.strip().split('\t')
			chrom, pos = temp[self.chrom_col], temp[self.pos_col]
			yield Position(chrom=chrom, pos=pos)
	def to_regions(self, positons=None):
		if positons is None:
			positons  = list(self)
		positons = sorted(positons, key=lambda x:(x.chrom, x.pos))
		last_pos = positons[0]
		adj_positons = [[last_pos]]
		for pos in positons[1:]:
			if pos - last_pos == 1:
				adj_positons[-1] += [pos]
			else:
				adj_positons += [[pos]]
		for cluster in adj_positons:
			chrom = cluster[0].chrom
			start = cluster[0].pos
			end = cluster[-1].pos
			yield Region(chrom, start, end)

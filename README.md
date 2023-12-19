Get 0-fold and 4-fold degenerate sites from Genome sequences (`*.fasta`) and Gene annotations (standard `*.gff3`).
### Installation ###
Dependencies:
+ [python 2.7](https://www.python.org/)
    + [biopython](https://biopython.org/)
    + [networkx](http://networkx.github.io/)
    + [lazy_property](https://github.com/jackmaney/lazy-property)
	+ quickly install by `pip2 install biopython<=1.76 networkx<2.0 lazy-property`

### Quick Start ###
```
git clone https://github.com/zhangrengang/degeneracy
cd degeneracy
python2 get_degeneracy.py test/test.gff3 test/test.fa
python2 get_degeneracy.py test/test.gff3 test/test.fa out_prefix
```

### Output ###
```
out_prefix.codon1-fold0.pos
out_prefix.codon2-fold0.pos
out_prefix.codon3-fold4.pos
```
The postions are in 1-based coordinates.

### Note ###
1. It concatenates the multiple CDS regions of a transcript and then identify 0d/4d positions, 
so the results include these positions in splicing sites, 
but it is slower than the methods based on the one-by-one single CDS region which have to omit the splicing sites. 
However, there are likely more evolutionary constraints in the splicing sites 
and it may be better to exclude the splicing sites.

2. It sorts the multiple CDS regions, so it is not suitable to trans-splicing genes.

Get 0-fold and 4-fold degeneracy sites from Genome sequences (`*.fasta`) and Gene annotations (standard `*.gff3`).
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

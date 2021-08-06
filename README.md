Get 0-fold and 4-fold degeneracy sites from Genome sequences (`*.fasta`) and Gene annotations (`*.gff3`).
### Installation ###
Dependencies:
+ [python 2.7](https://www.python.org/)
    + [biopython](https://biopython.org/): quickly install by `pip2 install biopython<=1.76`
    + [networkx](http://networkx.github.io/): quickly install by `pip2 install networkx<2.0`
    + [lazy_property](https://github.com/jackmaney/lazy-property): quickly install by `pip2 install lazy-property`

### Quick Start ###
```
git clone https://github.com/zhangrengang/degeneracy
cd degeneracy
python get_degeneracy.py test/test.gff3 test/test.fa
```

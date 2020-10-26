# BITS: miscellaneous BioInformatics ToolS

Basically bioinformatics-related, reusable Python codes as a package. Developed in Python 3.7, but most of the codes should work with any Python >= 3.

## How to install

```bash
$ git clone https://github.com/yoshihikosuzuki/bits
$ cd bits
$ python3 setup.py install
```

## How to use

Import and use any function/class you like:

```python
> from bits.seq.util import revcomp_seq
> print(revcomp_seq('aactg'))
```

```text
cagtt
```



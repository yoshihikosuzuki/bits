# BITS: miscellaneous BioInformatics ToolS

Basically bioinformatics-related, reusable Python codes as a package. Developed in Python 3.7, but most of the codes should work with any Python >= 3.

## How to install

```bash
$ git clone https://github.com/yoshihikosuzuki/bits
$ cd bits
$ python3 setup.py install
```

## How to use

BITS has two submodules: `bits.seq` and `bits.util`. The former is for DNA sequence manipulations and the latter is for other, general purposes in Python programming.

Import and use any function/class you like:

```python
> from bits.seq import revcomp_seq
> print(revcomp_seq('aactg'))
```

```text
cagtt
```

## List of functions



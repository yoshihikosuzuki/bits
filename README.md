# BITS: miscellaneous BioInformatics ToolS

Basically bioinformatics-related, reusable Python3 codes as a package.

## Requirements

- Python packages listed in `setup.cfg` (automatically installed)
- [Seqkit](https://bioinf.shenwei.me/seqkit/) (if you use `bits.seq.load_fast[a|q]` with `range=` option)
- [Edlib](https://github.com/Martinsos/edlib) (if you use `bits.seq.EdlibRunner`)
- [Gepard](https://github.com/univieCUBE/gepard) (if you use `bits.seq.DotPlot`)
- [DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB) (if you use `bits.seq.load_db` etc)

## How to install

```bash
$ git clone https://github.com/yoshihikosuzuki/bits
$ cd bits
$ pip install .
```

## How to use

BITS has two submodules:

- `bits.seq` for DNA sequence manipulations, and
- `bits.util` for other, general purposes in Python programming.

Import and use any function/class as follows:

```python
import bits.seq as bs
print(bs.revcomp_seq('aactg'))
```

```text
> cagtt
```

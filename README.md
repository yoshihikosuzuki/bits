# BITS: miscellaneous BioInformatics ToolS

Basically bioinformatics-related, reusable Python3 codes as a package.

## How to install

```bash
git clone https://github.com/yoshihikosuzuki/bits
cd bits
pip install .
```

To run unit tests,

```bash
python -m unittest
```

## Optional requirements

If you use `bits.seq.load_fast[a|q]` with `range=` option,

- [Seqkit](https://bioinf.shenwei.me/seqkit/)

If you use `bits.seq.DotPlot`,

- [Gepard](https://github.com/univieCUBE/gepard), or
- [modified Gepard](https://github.com/yoshihikosuzuki/gepard) (recommended)

If you use `bits.seq.load_db` etc.,

- [DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB)

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

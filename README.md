Here I put some reusable (bioinformatics-related) codes as a python3 package.



## How to install

```bash
$ git clone https://github.com/yoshihikosuzuki/BITS
$ cd BITS
$ python3 setup.py install
```



## How to use

In detail, **see [this Jupyter Notebook](https://nbviewer.jupyter.org/github/yoshihikosuzuki/BITS/blob/master/docs/usage.ipynb)** providing description of the modules and some examples.



### Basic usage

Import what you want:

```python
from BITS.seq.utils import revcomp
```

And use it:

```python
print(revcomp('aactg'))
```

Then you will get the result:

```ini
cagtt
```



### Short description on modules

* `BITS.clustering`: Classes for clustering of numerical and sequence data.
* `BITS.plot`: Utility functions for Plotly, enabling much simpler codes.
* `BITS.seq`: Alignment using edlib, dot plot using Gepard, and functions for CIGAR string and other sequence manipulations.
* `BITS.util`: Others including utilities for pickle, subprocess, pyinterval, shell command execution, and job scheduler.
Here I put some reusable (bioinformatics-related) codes as a python3 package.

## 0. Installation

```bash
$ git clone --recursive https://github.com/yoshihikosuzuki/BITS
$ cd BITS
$ python3 setup.py install
```


## 1. Command-line functions

You can use these functions from a shell as follows:

```bash
$ python -m BITS.<module_name> [arguments]
```


### Convert script to job management system (SGE or SLURM)-style script

```bash
# For SGE
$ python -m BITS.sge_nize <your_script_file> [options]

# For SLURM
$ python -m BITS.slurm_nize <your_script_file> [options]
```

In `[options]`, you can specify the arguments in the table below:

| Option name |           Default value            |     Note     |
| :---------: | :--------------------------------: | :----------: |
|  job_name   |            "run_script"            |              |
|   out_log   | "sge.log" (SGE)<br>"sbatch_stdout" |              |
|   err_log   |          "sbatch_stderr"           | (only slurm) |
|   n_core    |                 1                  |              |
| time_limit  |             "24:00:00"             | (only slurm) |
| mem_per_cpu |                1024                | (only slurm) |
|  partition  |              "batch"               | (only slurm) |
|    wait     |                True                |              |

Example usage:

```bash
$ python -m BITS.sge_nize my_script.sh job_name="run_my_script" n_core=12 wait=False
```



## 2. Python-module functions and classes

You can call these functions in a python script as follows:

```python
from BITS.<module_name> import <function_name>
```

To see the details of the functions, use the help command of IPython (`<function_name>???`).


### `BITS.utils`

|       Name        |   Type   |                         description                          |
| :---------------: | :------: | :----------------------------------------------------------: |
|    run_command    | function |              General-purpose command executer.               |
|     print_log     | function |  Simple decolator for watching start and end of a function.  |
|     sge_nize      | function |                 Add headers for qsub of SGE.                 |
|    slurm_nize     | function |               Add headers for sbatch of SLURM.               |
|   NoDaemonPool    |  class   | Inherited class of Pool so that it runs as a non-daemon process and thus can have child processes. |
|     make_line     | function |      For Plotly. Create a line-shape object for Plotly.      |
|   interval_len    | function | For pyinterval (https://pyinterval.readthedocs.io/en/latest/).<br> Return the sum of the interval lengths in <intvls>. |
| subtract_interval | function | For pyinterval (https://pyinterval.readthedocs.io/en/latest/).<br> Calculate A - B, where A and B are interval objects.
    Afterwards, remove remaining intervals shorter than <length_threshold>. |


### `BITS.run`

|        Name         |   Type   |           description           |
| :-----------------: | :------: | :-----------------------------: |
|        Cigar        |  class   | For cigar string manipulations. |
|    FlattenCigar     |  class   |                                 |
|        Align        |  class   |                                 |
|      run_edlib      | function |                                 |
|     run_consed      | function |                                 |
| consed_to_consensus | function |                                 |


### `BITS.seq`

|           Name            |   Type   |                         description                          |
| :-----------------------: | :------: | :----------------------------------------------------------: |
|        load_fasta         | function |                                                              |
|      single_to_multi      | function | Cut a single sequence at every <width> bp and return as a list. |
|        save_fasta         | function | NOTE: <reads> must be a dictionary. If <sort> is True, the headers in <reads> will be sorted.<br/>    <out_type> defines the existence of newlines within the sequences (by every <width> bp). |
|          revcomp          | function |     Return the reverse complement of the given sequence.     |
|          DotPlot          |  class   |                                                              |
| extract_adapters_from_bax | function |   Extract all adapter sequences detected in PacBio reads.    |
|    consensus_adaptors     | function | Extract all adapter sequences from *.bax.h5, and then take consensus of them. |

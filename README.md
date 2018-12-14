Here I put some reusable bioinformatics-related codes as a python3 package.

## Function list

### Add headers for job management system (currently SGE and SLURM) to a script

```bash
$ python -m BITS.sge_nize <your_script> [options]
$ python -m BITS.slurm_nize <your_script> [options]
```

### `BITS.io`

### `BITS.utils`

* run_command(command)
* qsub_script(in_script_name, out_script_name, args=(), job_name="run_script", out_log="sge.log", n_core=1, run_directly=False)

### `BITS.run`

* extract_adapters_from_bax(in_bax, out_adapters)
* run_consed(in_seqs, out_consensus, variant_vector=False, variant_graph=False, variant_fraction=0.3, display_width=80)
* consensus_adaptors(in_bax, out_adapters='adapters', out_consensus='adapter.consensus.fasta')

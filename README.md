Here I put some reusable bioinformatics-related codes as a python3 package.

## List of the functions

### `BITS.io`

### `BITS.utils`

* run_command(command)
* qsub_script(in_script_name, out_script_name, args=(), job_name="run_script", out_log="sge.log", n_core=1, run_directly=False)

### `BITS.run`

* extract_adapters_from_bax(in_bax, out_adapters)
* run_consed(in_seqs, out_consensus, variant_vector=False, variant_graph=False, variant_fraction=0.3, display_width=80)
* consensus_adaptors(in_bax, out_adapters='adapters', out_consensus='adapter.consensus.fasta')

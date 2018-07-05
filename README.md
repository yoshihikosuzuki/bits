Dipsered codes make it hard to re-use them. Here I gather small bioinformatics-related codes as a python3 library as well as some shell scripts used in it.

## Function list

### `BITS.core`

* extract_adapters_from_bax(in_bax, out_adapters)
* run_consed(in_seqs, out_consensus, variant_vector=False, variant_graph=False, variant_fraction=0.3, display_width=80)
* consensus_adaptors(in_bax, out_adapters='adapters', out_consensus='adapter.consensus.fasta')

### `BITS.utils`

* run_command(command)
* qsub_script(in_script_name, out_script_name, args=(), job_name="run_script", out_log="sge.log", n_core=1, run_directly=False)
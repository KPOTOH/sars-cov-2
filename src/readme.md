## Scripts description and execution order

<!-- - "./pipeline/0.sequences-rand-N.pl" -->
- `./filter_fasta_from_dublicates_and_NNN.py`
- `./pipeline/1.bwa.sh`
- `./pipeline/2.longest-alignments.pl`
- `./pipeline/3.sam2fasta.py`
- `./count_seq_distance_to_reference.py`
- `./mulal_qc.py`
- `./pipeline/4.fasttree.sh`
- `./tree_pruner.py`
- `./tree_simplifier.sh`
- `./resolve_polytomies_in_ete3.py`
- `./pipeline/6.ancestors-reconstruction.sh`
- `./mutations_extractor_with_context.py`
- `./calculate_distances_to_closest.py`


`add_edge_level_to_table` in `add_features_to_rerions.py`
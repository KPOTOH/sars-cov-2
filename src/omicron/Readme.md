## Notes

1. While filtration of sequences omicron reference dropped due to low quality. Consequently, need to append it to first place separately.
```
> awk 'BEGIN {RS=">"} /EPI_ISL_6752027/ {print ">"$o}' data/omicron/gisaid_omicron_13-01-22.fasta > data/omicron/omicron_ref.fasta
> cat data/omicron/omicron_ref.fasta data/omicron/sequences.filtered.fasta > data/omicron/sequences.filtered.fasta
```
### Mapping

* Reads are mapped to the *P. davidsonii* reference genome with bwa mem
* Duplicates marked with samtools markdup
* Overlapping paired end reads clipped with bamutil clipOverlap
* Some basic summary statistics generated with samtools coverage and samtools stats

All of the main commands are piped to avoid the creation of many intermediate files. See [`array_mapping_pipe.sh`]()

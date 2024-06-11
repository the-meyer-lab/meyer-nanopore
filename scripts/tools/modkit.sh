# Create bash script template
#!/bin/bash
# Create array of output_files
OUTPUT_FILES=
IN_FILES=("/Data1/seq_data/AI_N2_dimelo_antiDPY27_mpx_8_19_23/pod5_pass/barcode05/basecalls/barcode05.mod_mappings.sorted.bam" )
LOG_FILE="${OUT_DIR}AI_barcode05.pileup.log"
#"/Data1/seq_data/AI_N2_dimelo_antiDPY27_mpx_8_19_23/pod5_pass/barcode06/basecalls/barcode06.mod_mappings.sorted.bam")

#modkit pileup path/to/reads.bam output/path/pileup.bed --log-filepath pileup.log

/Data1/software/modkit/modkit pileup $IN_FILES $OUT_FILES --log-filepath $LOG_FILE --filter-threshold 0.5

'''
pileup
Tabulates base modification calls across genomic positions. This command produces a bedMethyl
formatted file. Schema and description of fields can be found in the README.

Usage: modkit pileup [OPTIONS] <IN_BAM> <OUT_BED>

Arguments:
  <IN_BAM>
          Input BAM, should be sorted and have associated index available.

  <OUT_BED>
          Output file (or directory with --bedgraph option) to write results into. Specify "-" or
          "stdout" to direct output to stdout.

Options:
      --log-filepath <LOG_FILEPATH>
          Specify a file for debug logs to be written to, otherwise ignore them. Setting a file is
          recommended. (alias: log)

      --region <REGION>
          Process only the specified region of the BAM when performing pileup. Format should be
          <chrom_name>:<start>-<end> or <chrom_name>

      --max-depth <MAX_DEPTH>
          Maximum number of records to use when calculating puleup. This argument is passed to the
          pileup engine. If you have high depth data, consider increasing this value substantially.
          Must be less than 2147483647 or an error will be raised.

          [default: 8000]

  -t, --threads <THREADS>
          Number of threads to use while processing chunks concurrently.

          [default: 4]

  -i, --interval-size <INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently. Smaller interval chunk sizes
          will use less memory but incur more overhead.

          [default: 100000]

      --chunk-size <CHUNK_SIZE>
          Break contigs into chunks containing this many intervals (see `interval_size`). This
          option can be used to help prevent excessive memory usage, usually with no performance
          penalty. By default, modkit will set this value to 1.5x the number of threads specified,
          so if 4 threads are specified the chunk_size will be 6. A warning will be shown if this
          option is less than the number of threads specified.

      --suppress-progress
          Hide the progress bar.

  -n, --num-reads <NUM_READS>
          Sample this many reads when estimating the filtering threshold. Reads will be sampled
          evenly across aligned genome. If a region is specified, either with the --region option or
          the --sample-region option, then reads will be sampled evenly across the region given.
          This option is useful for large BAM files. In practice, 10-50 thousand reads is sufficient
          to estimate the model output distribution and determine the filtering threshold.

          [default: 10042]

  -f, --sampling-frac <SAMPLING_FRAC>
          Sample this fraction of the reads when estimating the filter-percentile. In practice,
          50-100 thousand reads is sufficient to estimate the model output distribution and
          determine the filtering threshold. See filtering.md for details on filtering.

      --seed <SEED>
          Set a random seed for deterministic running, the default is non-deterministic.

      --no-filtering
          Do not perform any filtering, include all mod base calls in output. See filtering.md for
          details on filtering.

  -p, --filter-percentile <FILTER_PERCENTILE>
          Filter out modified base calls where the probability of the predicted variant is below
          this confidence percentile. For example, 0.1 will filter out the 10% lowest confidence
          modification calls.

          [default: 0.1]

      --filter-threshold <FILTER_THRESHOLD>
          Specify the filter threshold globally or per-base. Global filter threshold can be
          specified with by a decimal number (e.g. 0.75). Per-base thresholds can be specified by
          colon-separated values, for example C:0.75 specifies a threshold value of 0.75 for
          cytosine modification calls. Additional per-base thresholds can be specified by repeating
          the option: for example --filter-threshold C:0.75 --filter-threshold A:0.70 or specify a
          single base option and a default for all other bases with: --filter-threshold A:0.70
          --filter-threshold 0.9 will specify a threshold value of 0.70 for adenosine and 0.9 for
          all other base modification calls.

      --mod-thresholds <MOD_THRESHOLDS>
          Specify a passing threshold to use for a base modification, independent of the threshold
          for the primary sequence base or the default. For example, to set the pass threshold for
          5hmC to 0.8 use `--mod-threshold h:0.8`. The pass threshold will still be estimated as
          usual and used for canonical cytosine and 5mC unless the `--filter-threshold` option is
          also passed. See the online documentation for more details.

      --sample-region <SAMPLE_REGION>
          Specify a region for sampling reads from when estimating the threshold probability. If
          this option is not provided, but --region is provided, the genomic interval passed to
          --region will be used. Format should be <chrom_name>:<start>-<end> or <chrom_name>.

      --sampling-interval-size <SAMPLING_INTERVAL_SIZE>
          Interval chunk size in base pairs to process concurrently when estimating the threshold
          probability, can be larger than the pileup processing interval.

          [default: 1000000]

      --include-bed <INCLUDE_BED>
          BED file that will restrict threshold estimation and pileup results to positions
          overlapping intervals in the file. (alias: include-positions)

      --include-unmapped
          Include unmapped base modifications when estimating the pass threshold.

      --ignore <IGNORE>
          Ignore a modified base class  _in_situ_ by redistributing base modification probability
          equally across other options. For example, if collapsing 'h', with 'm' and canonical
          options, half of the probability of 'h' will be added to both 'm' and 'C'. A full
          description of the methods can be found in collapse.md.

      --force-allow-implicit
          Force allow implicit-canonical mode. By default modkit does not allow pileup with the
          implicit mode ('.', or silent). The `update-tags` subcommand is provided to update tags to
          the new mode. This option allows the interpretation of implicit mode tags: residues
          without modified base probability will be interpreted as being the non-modified base.

      --motif <MOTIF> <MOTIF>
          Output pileup counts for only sequence motifs provided. The first argument should be the
          sequence motif and the second argument is the 0-based offset to the base to pileup base
          modification counts for. For example: --motif CGCG 0 indicates to pileup counts for the
          first C on the top strand and the last C (complement to G) on the bottom strand. The --cpg
          argument is short hand for --motif CG 0.

          This argument can be passed multiple times. When more than one motif is used, the
          resulting output BED file will indicate the motif in the "name" field as
          <mod_code>,<motif>,<offset>. For example, given `--motif CGCG 2 --motif CG 0` there will
          be output lines with name fields such as "m,CG,0" and "m,CGCG,2". To use this option with
          `--combine-strands`, all motifs must be reverse-complement palindromic or an error will be
          raised.

      --cpg
          Only output counts at CpG motifs. Requires a reference sequence to be provided.

  -r, --ref <REFERENCE_FASTA>
          Reference sequence in FASTA format. Required for CpG motif filtering.

  -k, --mask
          Respect soft masking in the reference FASTA.

      --preset <PRESET>
          Optional preset options for specific applications. traditional: Prepares bedMethyl
          analogous to that generated from other technologies for the analysis of 5mC modified
          bases. Shorthand for --cpg --combine-strands --ignore h.

          [possible values: traditional]

      --combine-mods
          Combine base modification calls, all counts of modified bases are summed together. See
          collapse.md for details.

      --combine-strands
          When performing CpG analysis, sum the counts from the positive and negative strands into
          the counts for the positive strand.

      --edge-filter <EDGE_FILTER>
          Discard base modification calls that are this many bases from the start or the end of the
          read. For example, a value of 10 will require that the base modification is at least the
          11th base or 11 bases from the end.

      --only-tabs
          For bedMethyl output, separate columns with only tabs. The default is to use tabs for the
          first 10 fields and spaces thereafter. The default behavior is more likely to be
          compatible with genome viewers. Enabling this option may make it easier to parse the
          output with tabular data handlers that expect a single kind of separator.

      --bedgraph
          Output bedGraph format, see https://genome.ucsc.edu/goldenPath/help/bedgraph.html. For
          this setting, specify a directory for output files to be make in. Two files for each
          modification will be produced, one for the positive strand and one for the negative
          strand. So for 5mC (m) and 5hmC (h) there will be 4 files produced.

      --prefix <PREFIX>
          Prefix to prepend on bedgraph output file names. Without this option the files will be
          <mod_code>_<strand>.bedgraph.

      --partition-tag <PARTITION_TAG>
          Partition output into multiple bedMethyl files based on tag-value pairs. The output will
          be multiple bedMethyl files with the format
          `<prefix>_<tag_value_1>_<tag_value_2>_<tag_value_n>.bed` prefix is optional and set with
          the `--prefix` flag.

  -h, --help
          Print help information (use `-h` for a summary).
'''
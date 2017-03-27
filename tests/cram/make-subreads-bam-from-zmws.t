Test pbsv align

Set up
  $ . $TESTDIR/setup.sh
  $ fofn=$DATDIR/na12878_bam.fofn 
  $ txt=$DATDIR/na12878_subset_read_names.txt

  $ make-subreads-bam-from-zmws $txt $fofn subset

Test make-subreads-bam-from-zmws

Set up
  $ . $TESTDIR/setup.sh
  $ fofn=$DATDIR/na12878_bam.fofn 
  $ txt=$DATDIR/na12878_subset_read_names.txt
  $ O=$OUTDIR/cram-make-subreads-bam-from-zmws

  $ make-subreads-bam-from-zmws $txt $fofn $O --dry_run 1>/dev/null 2>/dev/null && echo $?
  0

  $ rm -rf $O.subreads.bam $O.subreadset.xml
  $ make-subreads-bam-from-zmws $txt $fofn $O 1>/dev/null 2>/dev/null && echo $?
  0
  $ ls $O.subreads.bam 1>/dev/null && echo $?
  0

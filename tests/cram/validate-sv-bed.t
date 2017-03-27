Test pbsv align

Set up
  $ . $TESTDIR/setup.sh
  $ bed=$DATDIR/test.bed

  $ validate-sv-bed $txt $bed $bed
  all_calls=11
  good_calls=11
  false_positives=0

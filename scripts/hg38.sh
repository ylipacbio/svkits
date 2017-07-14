hg=/pbi/analysis/smrtportal/beta/userdata/references/hg38/sequence/hg38.fasta
lens=50,100,200,500,1000,2000,5000
numbers=50,125,200,250,200,125,50 
time -v add-indels-to-fasta $hg hg.ins.fasta hg.ins.bed --insert --lens $lens --numbers $numbers --name_pattern '^chr(.|..)$' 2>hg.ins.log
time -v add-indels-to-fasta $hg hg.del.fasta hg.del.bed --delete --lens $lens --numbers $numbers --name_pattern '^chr(.|..)$' 2>hg.del.log

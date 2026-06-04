#Peak files are here:
#https://regmedsrv1.wustl.edu/Public_SPACE/bmiao/Public_html/TaRGET/Result/DAR/motif/

for file in `ls *.txt`
  do

    dir=$(echo $file|sed 's/\.txt//')
    echo $file
    echo $dir
    /L_Space/miao/TaRGET/homer/bin/findMotifsGenome.pl $file mm10 $dir -size given -preparsedDir preparsed

done

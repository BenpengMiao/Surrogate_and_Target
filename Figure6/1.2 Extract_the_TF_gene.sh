ls *filter.txt|sed 's/.txt//' > folder.txt

cat folder.txt|while read fold
do
 cd ${fold}
 grep "target=" homerResults.html |sed 's/<\/TD><TD>//'|sed 's/<\/TD><TD>/\t/g'|sed 's/\//\t/g'|cut -f 1-4,6 > a1
 grep "target=" homerResults.html |sed 's/<\/TD><TD>//'|sed 's/<\/TD><TD>/\t/g'|cut -f 6|sed 's/<BR/\t/'|cut -f 1|perl -F/ -ane 'chomp(@F); my $a=pop @F; print $a."\n";'|cut -d "(" -f 2|sed 's/)//' > a2
 grep "target=" homerResults.html |sed 's/<\/TD><TD>//'|sed 's/<\/TD><TD>/\t/g'|cut -f 6|sed 's/<BR/\t/'|cut -f 2|cut -d "/" -f 3|cut -d "\"" -f 1 > a3
 paste a1 a2 a3 > a
 awk '$1 < 1e-11' a|awk '$6 >= 0.8' > b
 rm a
 rm a1
 rm a2
 rm a3

 cat b|while read l
 do
   echo $fold >> k
   file=$(echo $l|cut -d " " -f 7) 
   grep "<H4>" homerResults/$file|cut -d "/" -f 1|cut -d "(" -f 1|sed 's/\s*//g'|sed 's/<H4>//' > g1
   grep "Score" homerResults/$file|sed 's/<\/TD><TD>/\t/g'|cut -f 2|cut -d "<" -f 1 > g2
   paste g1 g2 > g
   awk '$2 >= 0.8' g|cut -f 1|perl -e 'my %hash;foreach(<>){chomp;my $a=lc($_);$a=ucfirst($a);$hash{$a}=1;} my @d;foreach(keys %hash){push @d,$_;}my $b=join(",",@d);print $b."\n";' >> gene
 done

 paste b gene k > motif_gene.txt
 rm b
 rm gene
 rm g1
 rm g2
 rm g
 rm k
 
 cd ../ 
 
 cat */motif_gene.txt > motif_gene.txt

done





#!/usr/bin/perl

#use Cwd;

$dir = "I:/F/tumor_immune_interaction/computation/A_webserver/symbolChangeData";
chdir $dir;

#$dir2 = getcwd;
#print $dir2;

open(IN, "gene2accession") or die "$!";

open(IN, "gene2accession") or die "$!";
open(OUT,">>human_symbol2ID2.txt") or die "$!";
open(OUT2,">>mouse_symbol2ID2.txt") or die "$!";

while(defined($line=<IN>)){
   chomp $line;
   @gene = split(/\t/,$line);

   if($gene[0] eq 9606){
       $hash{$gene[1]} = $gene[15]."\t".$gene[0];
   }
   if($gene[0] eq 10090){
       $hash2{$gene[1]} = $gene[15]."\t".$gene[0];
   }
}

foreach $v (keys %hash){
    print OUT $v."\t".$hash{$v}."\n";
}

foreach $v (keys %hash2){
    print OUT2 $v."\t".$hash2{$v}."\n";
}

close(IN);
close(OUT);
close(OUT2);



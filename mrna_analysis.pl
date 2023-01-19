#!/usr/bin/perl

$opt = shift(@ARGV); 

$refseq = 'GGGACTAGTATTTCAGTTGCAGTATAATATAGTTAAGTAGTTAGAAGTATTGTAATATTCCCATGGGTGACTACAAAGACCATGACGGTGATTATAAAGATCATGACATCGATTACAAGGATGACGATGACAAGCTCgaagacgccaaaaacataaagaaaggcccggcgccattctatccgctggaagatggaaccgctggagagcaactgcataaggctatgaagagatacgccctggttcctggaacaattgcttttacagatgcacatatcgaggtggacatcacttacgctgagtacttcgaaatgtccgttcggttggcagaagctatgaaacgatatgggctgaatacaaatcacagaatcgtcgtatgcagtgaaaactctcttcaattctttatgccggtgttgggcgcgttatttatcggagttgcagttgcgcccgcgaacgacatttataatgaacgtgaattgctcaacagtatgggcatttcgcagcctaccgtggtgttcgtttccaaaaaggggttgcaaaaaattttgaacgtgcaaaaaaagctcccaatcatccaaaaaattattatcatggattctaaaacggattaccagggatttcagtcgatgtacacgttcgtcacatctcatctacctcccggttttaatgaatacgattttgtgccagagtccttcgatagggacaagacaattgcactgatcatgaactcctctggatctactggtctgcctaaaggtgtcgctctgcctcatagaactgcctgcgtgagattctcgcatgccagagatcctatttttcggcaatcaaatcattccggatactgcgattttaagtgttgttccattccatcacggttttggaatgtttactacactcggatatttgatatgtggatttcgagtcgtcttaatgtatagatttgaagaagagctgtttctgaggagccttcaggattacaagattcaaagtgcgctgctggtgccaaccctattctccttcttcgccaaaagcactctgattgacaaatacgatttatctaatttacacgaaattgcttctggtggcgctcccctctctaaggaagtcggggaagcggttgccaagaggttccatctgccaggtatcaggcaaggatatgggctcactgagactacatcagctattctgattacacccgagggggatgataaaccgggcgcggtcggtaaagttgttccattttttgaagcgaaggttgtggatctggataccgggaaaacgctgggcgttaatcaaagaggcgaactgtgtgtgagaggtcctatgattatgtccggttatgtaaacaatccggaagcgaccaacgccttgattgacaaggatggatggctacattctggagacatagcttactgggacgaagacgaacacttcttcatcgttgaccgcctgaagtctctgattaagtacaaaggctatcaggtggctcccgctgaattggaatccatcttgctccaacaccccaacatcttcgacgcaggtgtcgcaggtcttcccgacgatgacgccggtgaacttcccgccgccgttgttgttttggagcacggaaagacgatgacggaaaaagagatcgtggattacgtcgccagtcaagtaacaaccgcgaaaaagttgcgcggaggagttgtgtttgtggacgaagtaccgaaaggtcttaccggaaaactcgacgcaagaaaaatcagagagatcctcataaaggccaagaagggcggaaagatcgccgtgTAATTAACATAAGCTAGCTACCCATACGATGTTCCAGATTACGCTCTCGAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA';
$refseq = uc($refseq); 

if ($opt eq "help" || $opt eq ''){
    print STDERR "\nUsage: \"mod_mrna_analysis.pl <running option> [additional arguments]\"\n\n";
    print STDERR "Running options:\n";
    print STDERR "   help - prints out this help message\n";
    print STDERR "   split - splits the SAM alignment file to a collection full-matching and a collection of non-full matching aligned reads SAM file. Additional argument is the input SAM file.\n";
    print STDERR "   simpleQChist - calculate average QC score distribution of a sam file. Additional argument is the input SAM file.\n";
    print STDERR "   qcfilter - filter reads in the sam file based on average QC score. Additional arguments: 1. input SAM file, 2. QC score threshold\n";
    print STDERR "   mapmutations - maps insertions and deletions on the sequence and will create datafiles for plotting these in R. Additional argument is the input SAM file.\n\n";
}

if ($opt eq "split"){
    $infile = shift(@ARGV); 
    $outfile1 = $outfile2 = $infile;
    $outfile1 =~ s/\.sam/_fullmatch\.sam/;
    $outfile2 =~ s/\.sam/_nofullmatch\.sam/; 
    open(IN, $infile); 
    open(O1, ">".$outfile1);
    open(O2, ">".$outfile2);
    $cnt1 = $cnt2 = 0;
    while(<IN>){
        if ($_ =~ /^@/){ 
            print O1 $_; 
            print O2 $_;
        }
        else { 
            @fl = split(/\t/, $_);
            if ($fl[5] =~ /^\d+M$/ || $fl[5] =~ /^\d+S\d+M$/ || $fl[5] =~ /^\d+M\d+S$/ || $fl[5] =~ /^\d+S\d+M\d+S$/){ 
                $cnt1++;
                print O1 $_;
            }
            else { 
                $cnt2++;
                print O2 $_;
            }
        }
    }
    close(IN);
    close(O1);
    close(O2);
    print STDERR "Number of matching reads: ".$cnt1.", number of non-matching reads: ".$cnt2."\n";
}

if ($opt eq "simpleQChist"){
    open(IN, shift(@ARGV));
    for $i (0..40){ 
        $qscores[$i] = 0; 
    }
    while(<IN>){
        if (!/^@/){ 
            @fl = split(/\t/, $_); 
            $avg_qscore = qcScore($fl[10]); 
            $avg_qscore = int($avg_qscore+0.5);
            $qscores[$avg_qscore]++; 
        }
    }
    close(IN);
    for $i (0..40){ 
        print $i."\t".$qscores[$i]."\n"; 
    }
}

if ($opt eq "qcfilter"){
    open(IN, shift(@ARGV));
    $tr = shift(@ARGV);
    while(<IN>){
        if (/^@/){ 
            print $_;
        }
        else { 
            @fl = split(/\t/, $_);
            $avg_qscore = qcScore($fl[10]); 
            if ($avg_qscore >= $tr){ 
                print $_;
            }
        }
    }
    close(IN);
}

if ($opt eq "mapmutations"){
    @seqarr = split(/|/, $refseq); 
    for $pos (1..scalar(@seqarr)){
        $data{$pos}{'seq'} = $seqarr[$pos-1]; 
        $data{$pos}{'N'} = 0;
        $data{$pos}{'I'} = 0;
        $data{$pos}{'D'} = 0;
    }
    open(IN, shift(@ARGV));
    while(<IN>){
        if (!/^@/){
            @fl = split(/\t/, $_);
            if ($fl[1] == 16){ 
                @regs = $fl[5] =~ /(\d+[SMIDN])/g; # 
                $current_pos = $fl[3]; 
                for $r (@regs){ 
                    if ($r =~ /[SM]/){ 
                        $r =~ s/[SM]//;
                        $current_pos = $current_pos + $r;
                    }
                    else{
                        ($type) = $r =~ /([IDN])/;
                        $r =~ s/[IDN]//;
                        for $i (1..$r){
                            $data{($current_pos+$i-1)}{$type}++;
                        }
                        $current_pos = $current_pos + $r;
                    }
                }
            }
        }
    }
    print "POS\tREF\tNUM_DEL\tNUM_INS\n";
    for $i (1..scalar(@seqarr)){
        print $i."\t".$data{$i}{'seq'}."\t".$data{$i}{'D'}."\t".$data{$i}{'I'}."\n";
    }
    close(IN);
}


sub qcScore {
    if (length($_[0]) == 1){ 
        return (ord($_[0])-33); 
    }
    if (length($_[0]) > 1){
        my @chr = split(/|/, $_[0]); 
        my $sum = 0;
        for my $c (@chr){
            $sum += ord($c)-33; 
        }
        return ($sum / length($_[0]));
    }
}

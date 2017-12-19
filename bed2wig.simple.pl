#!/usr/bin/perl -w

my $temp_dir = '/scratch/tmp';
my $randn=substr(rand(),2,6);

my $usage = 'bed2wig.c.pl bed_file > bedgraph_no_header_no_orientation_file
';

my $nn = 'n';
my $or = 'n';
my $clean = 'y';
my $col4 = 'y';
my $prefix='';

open(FF,$ARGV[scalar(@ARGV)-1]) or die "cant open $ARGV[scalar(@ARGV)-2]\n";

# make temp directory, sort BED file by orient, chr, start
my $tdir = $temp_dir."/"."temp".$randn;
`mkdir $tdir`;

# split up multiple-exon beds into multiple records
open(SPLIT,">$tdir/split.bed") or die "nah uh-uh man, $tdir/split.bed no good for open\n";
while(my $r = <FF>) {  #split it

    chomp($r);
    my @MULT = split("\t",$r);
	
	my @LEN;
	my @STARTS;
    unless($MULT[10]=~/,/){$MULT[10]=$MULT[10].",";}
	unless($MULT[11]=~/,/){$MULT[11]=$MULT[11].",";}
	# remove dependency on 12 column requirement:
	# if 12 col entry
	if(scalar(@MULT)==12){	
		if(check12($r)==0){next;} # formatting issue in this line
        @LEN = split(",",$MULT[10]);
        @STARTS = split(",",$MULT[11]);
    
	}
	else { # no splicing, use cols 2 and 3 to determine start and length
	    if(check3($r)==0){next;} # formatting issue in this line
		my $ll = $MULT[2]-$MULT[1];
		unshift(@LEN,$ll);
		unshift(@STARTS,0);		
	}
    my $genome_s = $MULT[1];
    
    for (my $i=0; $i<scalar(@STARTS); ++$i) {
    	
		if( ($col4 eq 'y')&(scalar(@MULT)>3)) {
		    # for number of times sequence exists in bed file
			# if have fewer than 6 columns, fill in with + orient
    		for (my $k=0; $k<$MULT[3]; ++$k) { 
    		
				if(scalar(@MULT)>=6){				
        				print SPLIT $MULT[0], "\t", $genome_s+$STARTS[$i], "\t", $genome_s+$LEN[$i]+$STARTS[$i], "\t", $MULT[3], "\t", $MULT[4], "\t", $MULT[5], "\n";
           		}
				else { print SPLIT $MULT[0], "\t", $genome_s+$STARTS[$i], "\t", $genome_s+$LEN[$i]+$STARTS[$i], "\t", "1", "\t", "0", "\t", "+", "\n"; }
    		} # for repeat sequence mappings (column 4 in BED file)	
		}
		else { # each entry is single read
		
				if(scalar(@MULT)>=6){				
        				print SPLIT $MULT[0], "\t", $genome_s+$STARTS[$i], "\t", $genome_s+$LEN[$i]+$STARTS[$i], "\t", $MULT[3], "\t", $MULT[4], "\t", $MULT[5], "\n";
           		}
				else { print SPLIT $MULT[0], "\t", $genome_s+$STARTS[$i], "\t", $genome_s+$LEN[$i]+$STARTS[$i], "\t", "1", "\t", "0", "\t", "+", "\n"; }
    		}
	} # for all starts

}
close SPLIT;
# sort the single exon/bed file
if($or eq 'y') {
`sort -k 6,6 -k 1,1 -k 2,2n -k 3,3n -S 4G -T $temp_dir $tdir/split.bed > $tdir/sorted.bed`;
}
else {
`sort -k 1,1 -k 2,2n -k 3,3n -S 4G -T $temp_dir $tdir/split.bed > $tdir/sorted.bed`;
}

print STDERR "parsing\n";
# go through sorted file
my $oldid = 0;
my @recs = ();

open(F,"$tdir/sorted.bed") or die "nope no $tdir/sorted.bed\n";

# initialize on first record
my $r;
my $chr;
my @a = ();
$r = <F>;
$r =~ /(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\n/;
my $old_chr = $1;
my $first_s = $2;
my $e = $3;
my $max_e = $e;
my $old_orient = $6;
my $orient;
seek(F,0,0);
my $s;
my $sss;
my $eee;
my $c = 0;
my $track;

while ($r = <F>) {

    $r =~ /(.*?)\t(.*?)\t(.*?)\t/;
    
    $chr = $1;
    $s = $2;
    $e = $3;

    # if new chromo, dump last cluster
    unless ($chr eq $old_chr) {
    
    	   print STDERR "working on $chr\n"; 
    	   
    	   if ( $a[1][0] ) {
    	   
        	   ($track,$sss,$eee) = get_track(@a);   
        	   print_track($track,$old_chr,$a[0][0]); 
        	   @a=(); 
        	   $a[0][0]= $s; 
        	   $a[0][1]=$e; 
        	   $max_e = $e; 
        	   $c=1;
    	   
    	   }
    	   
    	    else {
       
           		print $old_chr, "\t", $a[0][0], "\t", $a[0][1], "\t", "1", "\n";  
           		@a=(); $a[0][0]= $s; $a[0][1]=$e; $max_e = $e; $c=1;
    		
       		}
    	   
    	   $old_chr = $chr; 
    	   next;

    }#chr switch
    $old_chr = $chr;
	
    # if new hyper-read, dump to track
    if ($max_e < $s) {
    
       if ( $a[1][0] ) {
           ($track,$sss,$eee) = get_track(@a);
           print_track($track,$chr,$a[0][0]); 
           @a=(); 
    	   $a[0][0]= $s; 
           $a[0][1]=$e;  
           $max_e = $e; $c=1;
       }
       
       else {
       
       		print $chr, "\t", ($a[0][0]), "\t", $a[0][1], "\t", "1", "\n";  
       		@a=(); $a[0][0]= $s; $a[0][1]=$e; $max_e = $e; $c=1;
    		
       }
       
    }
    
    else { $a[$c][0]=$s; $a[$c][1]=$e; if ($e>$max_e) {$max_e = $e;} ++$c;}

} # while sorted.bed file

# dump last cluster
($track,$sss,$eee) = get_track(@a);   
print_track($track,$old_chr,$a[0][0]); 
 # or = n

close F;

if($clean eq 'y') {
#clean up a bit
`rm -rf $tdir`;
}

########################  SUBS  #############

sub get_track {

    my(@a) = @_;
    my $ret = '';
	
    
    # find max e
    
    my $ee = 0;
    for (my $i=0; $i<scalar(@a); ++$i) {
    	
    	if($a[$i][1] > $ee) {$ee = $a[$i][1];}
    	
    }
    
    my $ss = $a[0][0];
    
    # make track
    for (my $j=$ss; $j<$ee; ++$j) {
    	
    	my $c = 0;
    	for (my $ii = 0; $ii < scalar(@a); ++$ii) {
    	
    		if ( ($j >= $a[$ii][0]) & ($j < $a[$ii][1]) ) { ++$c; }
    		
    	}
		
    	unless($c == 0) {$ret = $ret."\n".$c;}
    	if ($c == 0) {print "no man\n";}
    }
    
    $ret =~ s/\n//;
    
    return($ret,$ss,$ee);

}

sub print_track {

    my($tt,$chr,$ss)=@_;
        
    my @t = split("\n",$tt);
    
    my $old_nt = $t[0];
    my $old_ss = $ss;
    
	my $i=0;
    for ($i=0; $i<scalar(@t); ++$i) {
    	
    	unless ( $t[$i] eq $old_nt ) { 
    	
    		   print $chr, "\t", ($old_ss), "\t", ($ss+$i), "\t", $old_nt, "\n";
    		   $old_ss = ($ss+$i);
    		   
    	} # unless
		
		$old_nt = $t[$i];
    
    } #for 
	print $chr, "\t", ($old_ss), "\t", ($ss+$i), "\t", $old_nt, "\n";
}

sub check12 {
my($rr)=@_;

my @a = split("\t",$rr);
if( 
#($a[0] =~ /^chr[MmXxYyUn0-9]+(_random)?$/) &
($a[1] =~ /^[0-9]+$/) &
($a[2] =~ /^[0-9]+$/) &
($a[5] =~ /^[\+\-]$/) &
($a[6] =~ /^[0-9]+$/) &
($a[7] =~ /^[0-9]+$/) &
($a[8] =~ /^[0-9\,]+$/) &
($a[9] =~ /^[0-9]+$/) &
($a[10] =~ /^[\,0-9]+$/) &
($a[11] =~ /^[\,0-9]+$/) ) 
{return 1;}
else {return 0;}

}

sub check3 {
my($rr)=@_;

my @a = split("\t",$rr);
if( 
#($a[0] =~ /^chr[MmXxYyUn0-9]+(_random)?$/) &
($a[1] =~ /^[0-9]+$/) &
($a[2] =~ /^[0-9]+$/) )
{return 1;}
else {return 0;}

}

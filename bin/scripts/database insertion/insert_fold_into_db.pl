#!/usr/bin/perl -w 
use strict ; 
use 5.20.1 ; 
use DBI;
my $dbfile = "structure_database.db"  ;
my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","");
my $scopt = shift || die "need scope_total.txt as arg 1 " ; 
open(SCOP,"<./$scopt") ;
my $in = shift || die "need input hhr file as arg 2 " ; 
my %scop_info  ; 
while(<SCOP>){
    chomp ; 
    $_ =~ s/\r//gi ; 
    my @dat = split(/\|/, $_ ); 
    $scop_info{$dat[0]}{$dat[1]} = $dat[2] ; 
}
close SCOP ; 
open(HHR, "<$in") ; 
my @hhr = <HHR> ; 
my %scops ;
my $uid ;
my %seen ; 
foreach (@hhr)
{
	chomp  $_;
    if ($_ =~ m/^Query/)
    {
        $_ =~ s/Query\s+// ; 
        $uid = $_; 
    }
	next unless $_ =~ m/^\s+[0-9]+/gi ; 
	chomp $_ ;
	my @data = split ('\s+', $_) ;
	$data[3] =~ m/([a-z]\.[0-9]+\.[0-9]+\.[0-9]+)/gi; 
	my $fold = $1 ; 
    if (exists $seen{$fold})
    {
        next ;
    }
    else 
    {
        $seen{$fold} = 1 ;
    }
    my $prob ; 
	my $pvalue ; 
	my $evalue ; 
    my $col ; 
    my $hhmlen ; 
    my $start ; 
    my $stop ; 
	for ( my $i =  0 ;  $i < scalar @data ; $i++ ) 
	{
		my $info = $data[$i] ; 
		if ($info =~ m/^[0-9]+\.[0-9]$/gi ) 
		{
			$prob = $info ;
            $evalue = lc($data[$i+1]) ; 
            $pvalue = lc($data[$i+2]) ; 
            $col = $data[$i+5] ; 
            my $range = $data[$i+6]; 
            my $rang2 = $data[$i+7] ;
            $range =~ m/([0-9]+)-([0-9]+)/ ; 
            $start = $1 ; 
            $stop = $2 ; 
            $hhmlen = $data[$i+8] ;  
            if( ! defined $hhmlen )
            {
                $rang2 =~ m/[0-9]+-[0-9]+\(([0-9]+)\)/ ; 
                $hhmlen = $1  ; 
            }
            $hhmlen =~ s/\(//gi ;
            $hhmlen =~ s/\)//gi ;        
			last; 
		} 
	}
    my $cov = $col / $hhmlen ; 
    my $pos_in_hhr = 0; 
    foreach my $line (@hhr)
    {
        chomp $line ; 
        $pos_in_hhr++ ;
        if ($line =~ m/^No\s+$data[1]/)
        {
            last ; 
        } 
    }
    $pos_in_hhr = $pos_in_hhr + 1 ; 
    my $content  = $hhr[$pos_in_hhr] ; 
    my @content = split(' ',  $content  ) ;
    my @fold_components = split ( '\.' , $fold ) ;   
    my $class = $fold_components[0] ; 
    my $foldi = $fold_components[0] . "." . $fold_components[1] ; 
    my $sfam = $fold_components[0] . "." . $fold_components[1] . "." . $fold_components[2]  ;
    my $fam = $fold_components[0] . "." . $fold_components[1] . "." . $fold_components[2]  . "." . $fold_components[3]  ; 
    my $classd = $scop_info{'class'}{$class}; 
    my $foldid = $scop_info{'folds'}{$foldi};
    my $sfamd = $scop_info{'superfamilies'}{$sfam} ;
    my $famd = $scop_info{'families'}{$fam} ;
    $content[0] =~ s/^[A-Z][a-z]+=// ; #hhpred probability 
    $content[2] =~ s/^[A-Z][a-z]+=// ; #hhpred score 
    $content[3] =~ s/^[A-Z][a-z]+_[a-z]+=// ; #hhpred cols
    $content[4] =~ s/^[A-Z][a-z]+=// ; #hhpred %identity 
    $content[4] =~ s/\%// ; #hhpred %identity
    $content[5] =~ s/^[A-Z][a-z]+=// ;  #hhpred similarity  
    # say join (" ", 
    #     $uid, 
    #     $data[1],
    #     $class, 
    #     $classd, 
    #     $foldi,
    #     $foldid, 
    #     $sfam , 
    #     $sfamd , 
    #     $fam,  
    #     $famd, 
    #     $evalue, 
    #     $pvalue, 
    #     $cov ,  
    #     $content[0], 
    #     $content[2] , 
    #     $content[3] , 
    #     $content[4] , 
    #     $content[5] ) ; 
    my $final = [
        $uid, 
        $famd,
        $class,
        $classd,
        $foldi,
        $foldid, 
        $sfam , 
        $sfamd , 
        $fam,  
        $famd, 
        $content[0],
        $evalue, 
        $pvalue, 
        $cov ,  
        $content[2] , 
        $content[3] , 
        $content[4] , 
        $content[5] ,
        $start , 
        $stop 
    ] ;
    my $sth = $dbh->prepare("INSERT INTO fold (      
        uid ,                                          
        scop_code ,                                                                
        class ,                                        
        class_desc ,                              
        fold ,                         
        fold_desc ,                           
        superfam ,                      
        superfam_desc  , 
        fam , 
        fam_desc ,                                     
        prob ,                                         
        evalue ,                                      
        pvalue ,                                        
        coverage ,                                     
        score , 
        columns, 
        percent_identity , 
        similarity , 
        start, 
        stop                                        
        )                                              
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");   
    my @final2 = @$final ;   
    #say join ",", @final2  ; 
    $sth->execute(@final2);
}


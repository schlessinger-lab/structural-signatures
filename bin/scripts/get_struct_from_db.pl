#!/usr/bin/perl -w
use strict ; 
use 5.10.1 ; 
use DBI ;
use  Parallel::ForkManager ;
###helps 
my $optional = "####optional minimal thresholds for idenitifying folds####\n" .
    " length: arg 6 | default:30\n probablity: arg 7 | default:50\n overlapp coverage: arg 8 | default:0.8\n pvalue: arg 9 | default:1e-5\n" .
    " evalue: arg 10 | default:1e-5\n percent idenitity: arg 11 | default:0.3\n coverage with template: arg 12 | default:0.3\n" ;
my $help = "    need structure_database.db as arg 1
    need type: [domain] or [fold] or [both] as arg 2
    need protein name type: uniprot id: [uid] or  tgene name: [gn] as arg 3
    need ParentChildTreeFile.txt as arg 4
    need input genelist as arg 5
    need output file name as arg 6
    need available uniprot ids as arg 7" ; 
##inputs 
my $dbfile = shift || die "$help" ; 
my $type = shift || die "need type: [domain] or [fold] or [both] as arg 2 " ; 
my $name = shift || die "need protein name type:\n\tuniprot id: [uid] or\n\tgene name: [gn]\nas arg 3 " ; 
my $interpro_tree = shift || die "need ParentChildTreeFile.txt as arg 4" ; 
my $input = shift || die "need genelist as arg 5" ; 
my $fn = shift || die "need output file name as arg 6" ;
my $avail_uid = shift || die "need available uniprot ids as arg 7" ; 
my $remove_dup = shift || die "remove duplicate proteins from gene names? yes or no as arg 8" ; 
open(AVAI, "<$avail_uid") ;
my %avail ;
while(<AVAI>)
{
    chomp ; 
    $_ =~ s/\r|\n//g;
    $avail{$_} = 1; 
}
close(AVAI) ; 
##database 
my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","");
###optional parameters 
my $len = shift @ARGV || 30 ; #template length threshold arg 7
my $prob = shift @ARGV || 50 ; #templateprobablity threshold arg 8
my $overlappcov = shift @ARGV || 0.8 ;  #overlapp coverage threshold arg 9
my $pvalue = shift @ARGV || 1e-5 ; #template pvalue threshold arg 10
my $evalue = shift @ARGV || 1e-5 ; #template evalue threshold arg 11
my $pid = shift @ARGV || 0.3  ; #template percent identity threshold arg 12  
my $cov = shift @ARGV || .3 ; #template percent identity threshold arg 13 
my $maxproc = 25 ;
##mpi 
my $pm = new Parallel::ForkManager(4);
$pm->set_max_procs($maxproc);
open(IN, "<$input") ;
my @names ; 
while(<IN>)
{
    chomp ;  
    push @names , $_ ; 
} 
my @uid ; 
#####Convert names into uid 
if ($name=~ "gn")
{
    my $query = "select uid from uniprot_info where gname like ?" ; 
    @uid = get_uid($query)  ;  
    say STDERR "Converted Gene Names to UniprotIDs" ; 
}
elsif ($name =~ "uid")
{
    @uid = @names ; 
}
if ( $remove_dup eq "yes") 
{
    my %uid = map {$_ => 1} @uid ;
    @uid = keys %uid ; 
}
my @uidf ; 
##checks if uids are available ; 
foreach (@uid)
{
    chomp ;
    $_ =~ s/\r|\n//g; 
    if (exists $avail{$_})
    {
        push @uidf , $_ ;
    } 
}

say "Found: " , scalar(@uidf) , " out of " , scalar(@names)  , " in database" ;
@uid = @uidf ; 
open(STRUCT_INFO, ">>$fn.scop.structure.info.csv");
open(DOMAIN_INFO, ">>$fn.ipr.info.csv"); 
open(FOUND, ">$fn.found.genes") ; 
say FOUND scalar(@uid);
close (FOUND) ; 
#####Get structure counts 
if ($type =~ "domain" )
{
    say STDERR "Counting Domain Freqencies" ; 
    open( DOMAIN , ">>$fn.domain.cnt") ;
    domain_analysis(); 
}
elsif ( $type =~ "fold" )
{
    say STDERR "Counting Scop Freqencies" ;
    open( CL , ">>$fn.scop.class.cnt") ;
    open( F , ">>$fn.scop.fold.cnt") ;
    open( S , ">>$fn.scop.superfam.cnt") ;
    open( FA , ">>$fn.scop.family.cnt") ;
    fold_analysis(); 
}
elsif ( $type =~ "both" ) 
{
    say STDERR "Counting Domain Freqencies" ; 
    open( DOMAIN , ">>$fn.domain.cnt") ;
    domain_analysis(); 
    say STDERR "Counting Scop Freqencies" ;
    open( CL , ">>$fn.scop.class.cnt") ;
    open( F , ">>$fn.scop.fold.cnt") ;
    open( S , ">>$fn.scop.superfam.cnt") ;
    open( FA , ">>$fn.scop.family.cnt") ;
    fold_analysis(); 
}
else 
{
    die "need type: [domain] or [fold] or [both] as arg 2" ;
}
#####subroutines and functions#####

sub get_uid 
{
    my $query = shift ; 
    open (TMP , ">>$fn.uid.converted"); 
    say "Converting genenames to uniprot ids" ;
    foreach (@names )
    {  
        $pm->start and next;
        my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","","");
        my $sth = $dbh->prepare($query) ; 
        my $q = "%" . $_ . "%" ; 
        $sth->execute($q); 
        my $uid = $sth->fetchrow_array() ; 
        if (defined $uid )
        {   
            say TMP $uid ; 
        }
        else 
        {
            print STDERR "$_ was not found in the database\n" ; 
        }
        $pm->finish;
    }
    $pm->wait_all_children;
    open(UID, "<./$fn.uid.converted") ; 
    my @uid =  <UID> ; 
    #system("rm ./$fn.uid.converted") ;
    return @uid ; 
}

sub fold_analysis 
{
    my $sth1 = $dbh->prepare("select count(uid) from fold where uid like ? "); 
    my $sth2 = $dbh->prepare("select * from fold where uid like ? ") ; 
    print "Performing Fold analysis with the following minimum threshold values for assigning folds:\n" . 
    "\tlength: $len, probablity: $prob, overlap coverage: $overlappcov, pvalue: $pvalue\n" .
    "\tevalue: $evalue, percent idenitity: $pid, coverage with template: $cov\n" ; 
    my %famcnts ;
    my $cccc = 1 ;  
    my @identified_folds ; 
    foreach (@uid)
    {
        
        print  "extracting $cccc of " , scalar @uid ,  " proteins\n" ;
        $cccc++; 
        $_ =~ s/\r|\n//gi ; 
        $sth1->execute($_) ; 
        my $cnt_uid = $sth1->fetchrow_array() ; #number of folds
        next unless $cnt_uid > 0 ; #each fold row in db 
        $sth2->execute($_) ;
        my %folds ;  
        while ( my @row = $sth2->fetchrow_array())
        {
            if ( not defined $row[7])
            {
                say STDERR "Some fold were not idenitified for uid: $_ " ; 
                next ; 
            }
            else 
            {
		my $line_info =join( ',', @row[2,3,4,5,6,7,8,9,10] ) ;
		my $line = $_ .  "," .  $line_info ; 
		say STRUCT_INFO $line ; 
                $folds{$row[8]}{'class'} = $row[2] ; 
                $folds{$row[8]}{'class_d'} = $row[3] ; 
                $folds{$row[8]}{'fold'} = $row[4] ; 
                $folds{$row[8]}{'fold_d'} = $row[5] ; 
                $folds{$row[8]}{'sfam'} = $row[6] ; 
                $folds{$row[8]}{'sfam_d'} = $row[7] ; 
                $folds{$row[8]}{'fam'} = $row[8] ; 
                $folds{$row[8]}{'fam_d'} = $row[9] ; 
                $folds{$row[8]}{'prob'} = $row[10] ;
                $folds{$row[8]}{'ev'} = $row[11] ;
                $folds{$row[8]}{'pv'} = $row[12] ;
                $folds{$row[8]}{'cov'} = $row[13] ;
                $folds{$row[8]}{'score'} = $row[14] ;
                $folds{$row[8]}{'col'} = $row[15] ; #length/columns in common with template
                $folds{$row[8]}{'pid'} = $row[16] ; #percident
                $folds{$row[8]}{'start'} = $row[17] ;
                $folds{$row[8]}{'stop'} = $row[18] ;
                $folds{$row[8]}{'sim'} = $row[19] ; #similarity 
            }
        }
        #say $_ ;
        my %final = eval_fold(\%folds) ; 
        foreach my $scop ( keys %final)
        {
            if ($scop =~ "No-Fold" )
            {
                #say "No-Fold" ;
                $famcnts{'class'}{'No-Fold'}{'cnt'}++ ;  
                $famcnts{'class'}{'No-Fold'}{'desc'} = "No-Fold" ;
                $famcnts{'fold'}{'No-Fold'}{'cnt'}++ ;  
                $famcnts{'fold'}{'No-Fold'}{'desc'} = "No-Fold" ;
                $famcnts{'sfam'}{'No-Fold'}{'cnt'}++ ;  
                $famcnts{'sfam'}{'No-Fold'}{'desc'} = "No-Fold" ;
                $famcnts{'fam'}{'No-Fold'}{'cnt'}++ ;  
                $famcnts{'fam'}{'No-Fold'}{'desc'} = "No-Fold" ;  
                next ;
            } 
            #say "$scop\n" . "\t" . $folds{$scop}{'class'} . " " . $folds{$scop}{'class_d'} . "\n" . "\t\t". $folds{$scop}{'fold'} . " " . $folds{$scop}{'fold_d'} . "\n" . 
            #"\t\t\t" . $folds{$scop}{'sfam'} . " " .  $folds{$scop}{'sfam_d'} . "\n" ."\t\t\t\t"  . $folds{$scop}{'fam'} . " " . $folds{$scop}{'fam_d'} ;
            my $cl = $folds{$scop}{'class'} ; 
            my $cld = $folds{$scop}{'class_d'} ; 
            my $fo = $folds{$scop}{'fold'} ; 
            my $fod = $folds{$scop}{'fold_d'} ; 
            my $sf = $folds{$scop}{'sfam'} ; 
            my $sfd = $folds{$scop}{'sfam_d'} ; 
            my $fa = $folds{$scop}{'fam'}  ;
            my $fad = $folds{$scop}{'fam_d'} ; 
            $famcnts{'class'}{$cl}{'desc'} = $cld ;
            $famcnts{'class'}{$cl}{'cnt'}++ ;
            $famcnts{'fold'}{$fo}{'desc'} = $fod ;
            $famcnts{'fold'}{$fo}{'cnt'}++ ;
            $famcnts{'sfam'}{$sf}{'desc'} = $sfd ;
            $famcnts{'sfam'}{$sf}{'cnt'}++ ;
            $famcnts{'fam'}{$fa}{'desc'} = $fad ;
            $famcnts{'fam'}{$fa}{'cnt'}++ ;
        }
    }
    foreach my $type (keys %famcnts) 
    {
        if ($type =~ m/^class/)
        {
            foreach my $feat ( keys %{$famcnts{$type}})
            {
                say CL  $feat . "," . $famcnts{$type}{$feat}{'cnt'} , "," ,$famcnts{$type}{$feat}{'desc'} ; 
            }
        }
        if ($type =~ m/^fold/)
        {
            foreach my $feat ( keys %{$famcnts{$type}})
            {
                say F $feat . "," . $famcnts{$type}{$feat}{'cnt'} , "," ,$famcnts{$type}{$feat}{'desc'} ; 
            }
        }
        if ($type =~ m/^sfam/)
        {
            foreach my $feat ( keys %{$famcnts{$type}})
            {
                say S $feat . "," . $famcnts{$type}{$feat}{'cnt'} , "," ,$famcnts{$type}{$feat}{'desc'} ; 
            }
        }
        if ($type =~ m/^fam/)
        {
            foreach my $feat ( keys %{$famcnts{$type}})
            {
                say FA $feat . "," . $famcnts{$type}{$feat}{'cnt'} , ",", $famcnts{$type}{$feat}{'desc'} ; 
            }
        }        
    }
}
sub eval_fold 
{
    my $fold = shift @_ ; 
    my %hash = %$fold ;
    my %final ; 
    my @fold ; 
    my %filteredfolds ; 
    my $start ; 
    my $stop ; 
    # first order hash by longest fold 
    # my $overlappcov = shift @ARGV || 0.5 ;  #overlapp coverage threshold arg 8
    # $folds{$row[8]}{'start'} = $row[17] ;
    # $folds{$row[8]}{'stop'} = $row[18] ;
    foreach my $fold ( sort { $hash{$b}{'col'} <=>  $hash{$a}{'col'}}  keys %hash)
    {  
        ##check minimum thresholds: 
	if ( not defined $hash{$fold}{'stop'} or not defined  $hash{$fold}{'start'}  )
	{
	    next ; 
	}
        my $aln = $hash{$fold}{'stop'} - $hash{$fold}{'start'}  ; 
        next unless ($aln >= $len && $hash{$fold}{'prob'} >= $prob && $hash{$fold}{'ev'} <= $evalue && $hash{$fold}{'pv'} <= $pvalue && $hash{$fold}{'pid'} >= $pid && $hash{$fold}{'cov'}  >= $cov ); 
        push @fold , $fold ; 
        #say "\t$fold " . "len: " .$aln . " start " . $hash{$fold}{'start'}  . " stop " . $hash{$fold}{'stop'}   ; 
    }
    if ( scalar @fold == 0 )
    { 
        #say "\t No-Fold" ;
        $final{'No-Fold'} = 1; 
        return ( %final) ;  
        last ; 
    }
    for ( my $i = 0;  $i < scalar @fold ; $i++  ) 
    {   
        $start = $hash{$fold[$i]}{'start'}; 
        $stop = $hash{$fold[$i]}{'stop'} ; 
        next if ( exists $filteredfolds{$fold[$i]} ) ; 
        for (my $ii = 0 ;   $ii < scalar @fold ; $ii++ ) 
        {
            next if ( $ii == $i ) ; 
            my $qstart = $hash{$fold[$ii]}{'start'}; 
            my $qstop = $hash{$fold[$ii]}{'stop'} ; 
            #check for complete overlaps ; 
            if ($start <= $qstart && $stop >= $qstop )
            {
                $filteredfolds{$fold[$ii]} =1 ;
                #say "\t" . $fold[$ii] . " is within " . $fold[$i] ; 
                #say "\t#$fold[$i] was chosen as rep "; 
                next ;
            }
            elsif  ( $start <= $qstart && $stop >= $qstart  )
            {
                #say "\t\t" . $fold[$ii] . " is overlapped with " . $fold[$i]  ;
                my $overcovA = ( $stop - $qstart ) / ( $stop - $start ) ; 
                my $overcovB = ( $stop - $qstart ) / ( $qstop - $qstart ) ;
                if ( $overcovA > $overlappcov || $overcovB > $overlappcov )
                {
                    my $tlen = $stop - $start ; 
                    my $qlen = $qstop - $qstart ; 
                    if ( $tlen > $qlen )
                    {
                        $filteredfolds{$fold[$ii]} =1 ;
                        #say "\t\t\t#$fold[$i] was chosen as  rep $tlen vs $qlen ";
                        next ; 
                    }
                    elsif ($tlen == $qlen)
                    {
                        if ($hash{$fold[$i]}{'prob'} >= $hash{$fold[$ii]}{'prob'} )
                        {
                            #say "\t\t\t#$fold[$i] was chosen as rep  prob " . $hash{$fold[$i]}{'prob'} .  " vs " . $hash{$fold[$ii]}{'prob'} ;
                            $filteredfolds{$fold[$ii]} =1 ;
                            next ; 
                        }
                        else 
                        {
                            $filteredfolds{$fold[$i]} =1 ;
                            #say "\t\t\t#$fold[$ii] was chosen as rep  prob " . $hash{$fold[$ii]}{'prob'} .  " vs " . $hash{$fold[$i]}{'prob'} ;
                            next ;
                        }
                    }
                    else 
                    {
                        $filteredfolds{$fold[$i]} =1 ;
                        #say "\t\t\t#$fold[$ii] was chosen as rep  $tlen vs $qlen ";
                        next ; 
                    }
                }
                else 
                {
                    #say "\t\t#Both $fold[$i] and $fold[$ii] were chosen as rep" ;
                }
            }
            elsif ( $start >= $qstart && $start <= $qstop )
            {
                #say "\t\t" . $fold[$ii] . " is overlapped with " . $fold[$i]  ;
                my $overcovA = ( $qstop - $start ) / ( $qstop - $qstart ) ; 
                my $overcovB = ( $qstop - $start ) / ( $stop - $start ) ;
                if ( $overcovA > $overlappcov || $overcovB > $overlappcov )
                {
                    my $tlen = $stop - $start ; 
                    my $qlen = $qstop - $qstart ; 
                    if ( $tlen > $qlen )
                    {
                        $filteredfolds{$fold[$ii]} =1 ;
                        #say "\t\t\t#$fold[$i] was chosen as rep $tlen vs $qlen ";
                        next ; 
                    }
                    else 
                    {
                        $filteredfolds{$fold[$i]} =1 ;
                        #say "\t\t\t#$fold[$ii] was chosen as rep $tlen vs $qlen ";
                        next ; 
                    }
                }
                else 
                {
                    #say "\t\t#Both $fold[$i] and $fold[$ii] were chosen as rep" ;
                    next ; 
                }
            }
        } 
        next if ( exists $filteredfolds{$fold[$i]} ) ;
        #say "$fold[$i] was chosen as rep" ;
        $final{$fold[$i]} = 1;  ;
    }
    return (%final) ; 
}

sub domain_analysis 
{
    open(INTER, "<$interpro_tree") ; 
    my %interpro_roots ; 
    my %interpro_convert_roots ; ##take sub families/domains and convert to root 
    my $current_root ; 
    while(<INTER>)
    {
        chomp ; 
        if ($_ =~ m/^\-/ )
        {
            my @dat = split ("::" , $_ ); 
            $dat[0] =~ s/\-//gi ; 
            $interpro_convert_roots{$dat[0]} = $current_root ; 
        }
        else 
        {
            my @dat = split ("::" , $_ ); 
            $current_root = $dat[0] ;
            $interpro_roots{$current_root} = 1; 
        }
    }
    my %domain_count ; 
    my $sth = $dbh->prepare("select 
        uid,
        interpro_subfamily , 
        interpro_family, 
        interpro_subfamily_desc ,
        interpro_family_desc  from domain where uid like ? ") ; 
        my $cccc =1 ; 
    foreach (@uid)
    {
        #say "extracting $cccc of " , scalar @uid ,  " proteins" ;
        $cccc++; 
        $_ =~ s/\r|\n//gi ; 
        $sth->execute($_); 
        my ($uid, $idoma , $ifam , $idomad , $ifamd ) = $sth->fetchrow_array() ; 
        if ( not defined $uid )
        {
            $domain_count{'No-Domain'}++ ; 
            next ;
        }
	say DOMAIN_INFO $uid ."," .$idoma . ",".  $ifam . "," . $idomad . ",". $ifamd ; 
        if ( $idoma =~ m/NULL/ )
        {
            if ( $ifam =~ m/NULL/ )
            {   
                $domain_count{'No-Domain'}++ ; 
            }
            else 
            {
                my @root ; 
                my @sp = split(';', $ifam) ; 
                foreach my $ipr (@sp) 
                {     
                    if (exists $interpro_roots{$ipr} ) 
                    {
                        push (@root , $ipr ) ; 
                    } 
                    else 
                    {
                        if (exists $interpro_convert_roots{$ipr})
                        {
                            push (@root , $interpro_convert_roots{$ipr}) ; 
                        }
                        else
                        {
                            push (@root , $ipr ) ;
                        }
                    }   
                }
                my %roots = map {$_ => 1} @root ; # one unique domain per protein 
                #say join("\t", $uid, keys %roots) ;
                foreach my $root (keys %roots)
                {
                    $domain_count{$root}++ ;
                }
            }
        }
        else 
        {
            my @root ; 
            my @sp = split(';', $idoma) ; 
            foreach my $ipr (@sp) 
            {
                if (exists $interpro_roots{$ipr} ) 
                {
                    push (@root , $ipr ) ; 
                } 
                else 
                {
                    if (exists $interpro_convert_roots{$ipr})
                    {
                        push (@root , $interpro_convert_roots{$ipr}) ; 
                    }
                    else
                    {
                        push (@root , $ipr ) ;
                    }
                }
            }
            my %roots = map {$_ => 1} @root ; # one unique domain per protein 
            #say join("\t", $uid, keys %roots) ;
            foreach my $root (keys %roots)
            {
                $domain_count{$root}++ ;
            }
        }
    }
    foreach (keys %domain_count)
    {
        say DOMAIN "$_," . $domain_count{$_}; 
    }
}

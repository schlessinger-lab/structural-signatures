#!/usr/bin/perl -w 
use strict ; 
use 5.20.1 ; 
use DBI;
my $dbfile = "structure_database.db"  ; 
my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile","",""); 
open(F, "<./human_proteome_uniprot.12.9.2017.tab") ; 
open(IPR, "<./entry.list");
my %ipr ; 
my $fam ; 
my $subfam ; 
while(<IPR>)
{
    next if $. == 1 ; 
    my @dat = split ("\t", $_) ; 
    next unless $dat[1] =~ "Domain" || $dat[1] =~ "Family"; 
    $ipr{$dat[0]}{$dat[1]} = $dat[2]; 
}
#populate_uniprot_info(); ##run if uniprot info is empty 
#populate_domain();

sub populate_domain 
{
    while(<F>)
    {
        chomp ; 
        $_ =~ s/\r//gi ; 
        $_ =~ s/\"//gi ; 
        next if $. == 1 ;     
        print "working on line $. \n" ;
        my @dat = split ("\t", $_ ) ; 
        my @final ; 
        push @final , $dat[0] ; 
        push @final , @dat[51..53] ;
        my @fams ; 
        my @domains; 
        my @fam_desc ; 
        my @domain_desc ; 
        if( not defined ($dat[73]) )  
        {   
            next; 
        }
        else 
        {
            $dat[73] =~ s/\"//gi ; 
            $dat[73] =~ s/\r\n//gi ;
            my @ipr = split ("," ,$dat[73] ); 
            foreach my $i (@ipr )
            {
                if ( exists $ipr{$i} )
                {
                    my @type = keys $ipr{$i} ; 
                    if ($type[0] =~ "Family"){
                        push @fams , $i ;
                        my $d = $ipr{$i}{$type[0]} ; 
                        $d =~ s/,/-/gi ; 
                        $d =~ s/\n//gi ; 
                        push @fam_desc , $d  ;
                    } 
                    else {
                        push @domains , $i ; 
                        my $d = $ipr{$i}{$type[0]} ; 
                        $d =~ s/,/-/gi ; 
                        $d =~ s/\n//gi ; 
                        push @domain_desc , $d ;
                    }
                }
            }
        }
        if ( scalar @fams > 0 ) 
        {
            my $fams = join ";" , @fams ; 
            my $famd = join ";" , @fam_desc ; 
            push @final , $fams ; 
            push @final , $famd ; 
            if ( scalar @domains > 0 )
            {
                my $do = join ";" , @domains ; 
                my $dod = join ";", @domain_desc ; 
                push @final, $do ; 
                push @final, $dod ; 
            }
            else 
            {
                push @final , ( "NULL" ,"NULL" ) ; 
            }
        }
        else 
        {
            push @final , ( "NULL" ,"NULL")  ; 
            if ( scalar @domains > 0 )
            {
                my $do = join ";" , @domains ; 
                my $dod = join ";", @domain_desc ; 
                push @final, $do ; 
                push @final, $dod ; 
            }
            else 
            {
                push @final , ( "NULL" ,"NULL" ) ; 
            }
        }
        push @final , @dat[74..79] ; 
        my @final2 ; 
        for ( my $i =0 ; $i < scalar @final ; $i++ ) 
        {
            if( length ($final[$i]) > 0 )  
            {   
                $final[$i] =~ s/\"//gi ; 
                $final[$i] =~ s/\'//gi ; 
                push @final2 , $final[$i] ; 
                next; 
            }
            else 
            {
                push @final2 , "NULL" ; 
            }
        }
        my $sth = $dbh->prepare("INSERT INTO domain (
            uid , 
            domain_desc , 
            domain_prosite , 
            motif , 
            interpro_family , 
            interpro_family_desc , 
            interpro_subfamily , 
            interpro_subfamily_desc , 
            panther  , 
            pfam , 
            prosite ,
            smart , 
            superfam , 
            tigrfam 
            ) 
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
        $sth->execute(@final2); 
    }
}



sub populate_uniprot_info 
{
    while (<F>)
    {
        chomp ; 
        $_ =~ s/\r//gi ; 
        next if $. == 1 ; 
        print "working on line $. \n" ;
        my @dat = split ("\t", $_ ) ; 
        my @final ; 
        push @final, @dat[0..2] ; 
        $dat[3] =~ s/\s+/,/gi ;
        push @final , $dat[3] ; 
        push @final , $dat[5] ;
        push @final , @dat[10..11] ; 
        push @final , @dat[13..15] ;  
        push @final , @dat[17..22] ; 
        push @final , @dat[24..34] ; 
        push @final , @dat[36..49] ; 
        push @final , @dat[54..79] ;
        my @final2 ; 
        for ( my $i =0 ; $i < scalar @final ; $i++ ) 
        {
            if( length ($final[$i]) > 0 )  
            {   
                $final[$i] =~ s/\"//gi ; 
                $final[$i] =~ s/\'//gi ; 
                push @final2 , $final[$i] ; 
                next; 
            }
            else 
            {
                push @final2 , "NULL" ; 
            }
        }
        my $sth = $dbh->prepare("INSERT INTO uniprot_info (
            uid , 
            entry_name , 
            protein_name , 
            gname  , 
            length  , 
            sequence , 
            active_site , 
            catalytic_activity , 
            cofactor , 
            dna_binding , 
            enzyme_regulation , 
            prot_funct , 
            metal_binding ,
            nucleotide_binding , 
            pathway , 
            pH_dependence , 
            temperature , 
            binding_sites , 
            kinetics , 
            tissue  , 
            induction , 
            development , 
            go_bp , 
            go_cc , 
            go_loc , 
            go_molfunct , 
            disease , 
            mutagenesis , 
            pharma , 
            intramembrane , 
            cellular_loc , 
            topological_domain , 
            transmembrane_helix , 
            three_d_structures , 
            beta_strands , 
            alpha_helix ,
            loops_turn , 
            date_of_creation , 
            date_of_modification , 
            version , 
            coiled_coil , 
            family ,
            region , 
            repeats , 
            sequence_similarity , 
            zinc_fing , 
            embl_id , 
            ref_seq , 
            unigene , 
            pir_id , 
            ccds_id , 
            pdb_ids , 
            disprot_id , 
            chembl_id , 
            dbsnp_ids , 
            ensembl_ids , 
            geneid_id , 
            genedb_id , 
            kegg , 
            gene_three_d_id ,
            interproscan , 
            panther , 
            pfam , 
            prosite , 
            smart , 
            superfam , 
            tigrfam 
        ) 
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)");
        $sth->execute(@final2);
        #exit ; 
    }
}
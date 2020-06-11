# structural-signatures-2.0
### Leaner, faster, better,

## Structural Signatures Generator! (now with more emojis 🎉)

To start, untar the tar.gz file in the database folder and run structural-signatures-2.0.sh.

Next export the `homed` variable to be the installation directory of structural signatures, add it to your `.bashrc`. 

If you are in the installation directory and want it automatically done, you can use the following command on linux systems (Ubuntu, WSL), not OSX. 

`$ pwd | sed -E "s/(.+)/export homed='\1'/gi"  >> ~/.bashrc`

Otherwise just export the installation directory, using the following command, changing the `/path/to/structural-signatures-2.0/` to your installation path: 

`export homed='/path/to/structural-signatures-2.0/'  >> ~/.bashrc`

## Dependancies
R >=3.0

Rscript

Perl >=5.22.8

Perl DBI

Perl Parrallel Fork Manager

sqlite >3.0

## Running structural-signatures-2.0.sh 

Running `structural-signatures-2.0.sh` without any arguments will give the following: 

![](https://i.imgur.com/EVFfFpO.png  "log" )

To run, specify the `-i`, `-t`, `n` and `-o` arguments. A quick example can be run using the following code (provided that all the dependancies are met): 

`./structural-signatures-2.0.sh -i apoptosis.gene.list.head -t both -n gn -o test` 

This will generate the following output if successful: 
![](https://i.imgur.com/seaJzyN.png )

And will generate the following output: (***note: outputs will change given how you specify the `-t` option*** )

![](https://i.imgur.com/sX6xwYt.png) 

Where: 
1) `*-enrichment.csv` files are the structural signatures output 
2) `*.info.csv` contain information about what structures were found for which proteins 
  (ipr.info.csv has domain information and scop.structure.info has SCOPe structure information) 
3) `*.uid.converted` contains the uniprot identifiers for the input gene list 

## Adding custom backgrounds to compare signatures against. 

Structural signatures generates pvalues and fold changes by comparing the counted structures against the structure counts from the human proteome defined by SwissProt (~20,000 proteins)

You may want to compute structural enrichment against a background that is a subset of human proteome from SwissProt. 

To do this first you can run `structural-signatures-2.0.sh -t both` on the genelist you want to generate the background for. 

Next for each of the structral signature files (`*-enrichments.csv`) run the following command to get the structure counts: 

`cut -f1,2 -d","` **`filename`**`-`**`structure`**`-enrichment.csv > `**`name`**`.background.`**`[ipr/scop]`**`.`**`[domain/family/superfam/fold]`**`.csv`

Where:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**filename** 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;is the output file name from structural signatures. 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**structure** 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;is the output structure type: *domain*, *fold*, *family*, *superfam* or *fold*

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**name** 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;is the name of the background you give 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**[ipr/scop]** 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;if structure type is domain use ***ipr*** 


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;if fold, family, superfamily or fold use ***scop***

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**[domain/family/superfam/fold]** 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;assign for the specific input structure type you are working with 

Next make a directory here: `bin/files/backgrounds/**name**`  where **name** is the name of the backgrounds you assigned earlier (___everything must have the same name___). Move all the files to that directory. 

To use the custom background just specify the `-r` option to `structural-signatures-2.0.sh` and input the **name** of the database. Structural-signatures will verify if the counts are in the correct format. 

For examples of how you should format your custom background counts look at the `bin/files/backgrounds/human_proteome/` directory files. 

## To Do: 
- [ ] Make an installation script 
- [ ] Backport `structural-signatures-1.0 features`(disordered regions, 2D protein feature enrichment, etc)
- [ ] Multiprocessing of protein structure extraction (currently runs one at a time) 
- [ ] Dockerize 
- [ ] Add various error tests 

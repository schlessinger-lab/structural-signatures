#!/usr/bin/Rscript  
args <- commandArgs(TRUE)
counted.frequencies = read.table(args[1], sep = ",", header = F, ) 
background.frequencies = read.table(args[2], sep = ",", header =F )
number.proteins = as.numeric(args[3] )
output = as.character(args[4])
type = as.character(args[5])
ip = c()
c = c()
bg.c = c()
np = c()
cpr = c()
pv = c()
pco = c()
folc = c()
count.proteome = 20230 ##number of proteins in reviewed proteome 
for ( i in 1:nrow(counted.frequencies)){
    ipr = as.character(counted.frequencies[i,1])
    cnt = as.numeric(counted.frequencies[i,2])
    wo = as.numeric(number.proteins - cnt )
    bg.cnt = as.numeric(background.frequencies[which(as.character(background.frequencies$V1) == ipr ) , 2] )
    if ( length(bg.cnt) == 0 ) 
    {
	next ;
    } 
    wo.bg.cnt = count.proteome - bg.cnt
    mat = as.data.frame(cbind(c(cnt,wo),c(bg.cnt,wo.bg.cnt)) )
    mat.ft = fisher.test(mat , alternative = "greater" )
    pvalue = mat.ft$p.value
    pvalue.corrected = .05 / nrow(counted.frequencies)
    freq = cnt/number.proteins 
    bg.freq = bg.cnt/count.proteome 
    fc = log10(freq/bg.freq)
    ip = c(ip,ipr)
    c = c(c, cnt)
    bg.c = c(bg.c, bg.cnt )
    np = c(np, number.proteins )
    cpr = c(cpr, count.proteome )
    pv = c(pv,pvalue)
    pco = c(pco, pvalue.corrected)
    folc = c(folc, fc )
#    print(paste(ipr, cnt , bg.cnt , number.proteins , count.proteome ,pvalue, pvalue.corrected, fc  , sep = " ")) 
}
fdr = as.numeric(as.character(p.adjust(as.numeric(as.character(pv)), method = "fdr")))
se = seq(1,length(ip), by = 1 )
ou = rep(output,length(ip))
df = as.data.frame(cbind(se,as.character(ip),c,bg.c,np, cpr, pv, fdr,  pco,  folc, ou )) 
df = df[order(fdr),]
#row.names(df) = df$se 
df$se = NULL
#df$fdr = format(df$fdr, scientific = FALSE) 
names(df) = c("structure", "counts_observed", "background_counts", "number_of_genes_in_set",  "total_number_proteins_proteome" , "pvalue" ,  "fdr" , "bonforroni_cutoff" , "log_fold_change", "name")
write.table(file = paste0("./", output, "-" , type , "-enrichments.csv" ), x= df,  sep = "," , row.names = F, quote = F ) 
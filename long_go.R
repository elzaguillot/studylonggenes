## short script to retrieve gos with long and short transcripts in average

## first retrieve from biomart ensembl, short and long genes with their go
library(biomaRt)
mart = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl", mart = mart)
allgenes <- getBM(attributes=c('start_position','end_position',"transcript_length",'ensembl_gene_id',"go_id"), mart = ensembl)

## then filter go with minimum size 10, measure average transcript lengths and order them
library(dplyr)
gotostudy <- allgenes %>% group_by(go_id) %>% summarize(n=n(),mean_size=mean(transcript_length,na.rm=T))  %>% filter(n>10) %>% arrange(desc(mean_size)) 
head(gotostudy)[,1]
tail(gotostudy)[,1]

## retrieve go description from GO.db
library(GO.db)
getdescription<-function(x)
{
  return(Term(GOTERM[x]))
}

## extract longest and shortest average length go
ngo=dim(gotostudy)[1]
longgo=gotostudy[(ngo-10):ngo,]%>% mutate(description=getdescription(go_id))
shortgo=gotostudy[1:10,]%>% mutate(description=getdescription(go_id))

## collate and output in a table
write.table(rbind(longgo,shortgo),file="short_and_long_go.txt")

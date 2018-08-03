#!/bin/bash

WORKDIR=~/CADcausalitytest
GENENAME="$(cat $1)"
GENE=$(pwd)/$1
CHR=$2
DISTANCE=$3
TRESHOLD=$4
NELSON=UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt
CC4D=~/CADcausalitytest/CARDIOGRAMC4D
EXPR=$WORKDIR/HCASMC_expr
GENOTYPES=$WORKDIR/HCASMC_genotypes
VCF=$WORKDIR/HCASMC_genotypes/vcf
REV=$EXPR/reverse

if [ ! -d $WORKDIR ]
then
mkdir $WORKDIR
fi

if [ ! -d $CC4D ]
then
mkdir $CC4D
fi

cd $CC4D

if [ ! -f $NELSON ]
then
wget https://www.dropbox.com/s/hjmkk5rjj0tfswy/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt
echo "Downloaded CARDIOGRAM plus C4D data, Nelson et al."
fi


awk '{if ($3=="'"$CHR"'") print $0}' UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt > UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.chr$CHR.txt

#write R script to get ENSEMBL id, needs biomaRt in R

echo "#!/usr/bin/Rscript
library(biomaRt)
listMarts(host=\"grch37.ensembl.org\")
ensembl = useMart(\"ENSEMBL_MART_ENSEMBL\",dataset=\"hsapiens_gene_ensembl\", host=\"grch37.ensembl.org\")
id_merge = getBM(attributes=c(\"ensembl_gene_id\",\"external_gene_name\",\"transcript_start\",\"transcript_end\"),mart=ensembl)
write.table(id_merge, file=\"id_merge.txt\", sep = \"\t\", quote =F, col.names=F, row.names=F)
" > script.r

#run R script

chmod 775 script.r
./script.r

grep "	$GENENAME	" id_merge.txt > gene_id_merge

START="$(awk 'BEGIN{min=100000000000000000}{if($3<min){min=$3;out=$3}}END{print out}' gene_id_merge)"
END="$(awk 'BEGIN{max=0}{if($4>max){max=$4;out=$4}}END{print out}' gene_id_merge)"

echo "Composite gene start:" $START
echo "Composite gene end:" $END
echo "Extending with the distance:" $DISTANCE "bp"

awk 'NR==FNR {if ($4>="'"$START"'"-"'"$DISTANCE"'" && $4<="'"$END"'"+"'"$DISTANCE"'") print $0}' UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.chr$CHR.txt  > chr$CHR.$GENENAME.region.txt
 
cut -f2,5,10 chr$CHR.$GENENAME.region.txt > SNP_effect.alele_pval.txt

echo "Selecting p-value threshold" $TRESHOLD
awk -v var="$TRESHOLD" '{if ($3<=var) print $0}' SNP_effect.alele_pval.txt > SNP_effect.alele_pval.threshold.txt

#create genotype files

if [ ! -d $GENOTYPES ]
then
mkdir $GENOTYPES
fi


cd $GENOTYPES


if [ ! -d $VCF ]
then 
mkdir $VCF
cd $VCF
wget https://www.dropbox.com/s/nnytxlbx1v0gh8y/phased_and_imputed.tar
echo "Unpacking genome vcf files..."

tar -xvf phased_and_imputed.tar
gunzip phased_and_imputed*
fi

cd $VCF

grep CHROM phased_and_imputed.chr$CHR.vcf > HEADER.txt
awk '{$1=$2=$3=$4=$5=$6=$7=$8=$9=""; print $0}' HEADER.txt > HEADER.txt.cut
#cat HEADER.txt.cut


while read -r a b c; do
        grep "$a	" phased_and_imputed.chr$CHR.vcf > SNP.txt
	#cat SNP.txt
        #clean genotype files
        sed -i -E 's/:[0-9.,]*//g' SNP.txt 
	#cat SNP.txt
       
#grab reference and alternative aleles
        REF="$(awk '{printf $4}' SNP.txt)"
        ALT="$(awk '{printf $5}' SNP.txt)"
	#echo $REF
	#echo $ALT

        EFFECT=$b
        #echo $EFFECT
	
	sed -i "s/0|0/$REF$REF/g" SNP.txt
	sed -i -E "s/0\|(1|2)/$REF$ALT/g" SNP.txt
        sed -i -E "s/(1|2)\|0/$REF$ALT/g" SNP.txt
        sed -i "s/1|1/$ALT$ALT/g" SNP.txt
	#cat SNP.txt
        
        awk '{$1=$2=$3=$4=$5=$6=$7=$8=$9=""; print $0}' SNP.txt > SNP.txt.cut 
        #cat SNP.txt.cut

	grep -v -P "[ATCG,]{3,}" SNP.txt.cut > SNP.txt.cut2
	sed -i "s/"$EFFECT"/1/g" SNP.txt.cut2
        sed -i -E "s/[ATGC]/0/g" SNP.txt.cut2
        
#cat SNP.txt.cut

        cat HEADER.txt.cut SNP.txt.cut2 > GENOTYPES.$a.txt
	#cat GENOTYPES.$a.txt

done < $CC4D/SNP_effect.alele_pval.threshold.txt

cat GENOTYPES.rs* > GENOTYPES.combined

rm GENOTYPES.rs*

grep -v "1020301 102901" GENOTYPES.combined > GENOTYPES.combined.even
#awk 'NR%2==0' GENOTYPES.combined > GENOTYPES.combined.even

sed -i 's/11/2/g' GENOTYPES.combined.even
sed -i -E 's/(10|01)/1/g' GENOTYPES.combined.even
sed -i 's/00/0/g' GENOTYPES.combined.even

cat HEADER.txt.cut GENOTYPES.combined.even > GENOTYPES.combined.even.HEADER

####awk 'BEGIN{print "count", "lineNum"}{print gsub(/t/,"") "\t" NR}' file
####awk -F'|' -v fld=2 'BEGIN{print "count", "lineNum"}{print gsub(/t/,"",$fld) "\t" NR}' file

grep -v "1020301 102901" GENOTYPES.combined > GENOTYPES.combined.even.b
cat HEADER.txt.cut GENOTYPES.combined.even.b > GENOTYPES.combined.even.b.header

if [ ! -d $EXPR ]
then
mkdir $EXPR
fi

cd $EXPR

#check if reverse folder with expression levels counted on reverse strand exists, if not, download from Dropbox link

if [ ! -d $REV ]
then
wget https://www.dropbox.com/s/edm0ykexjmue5yf/reverse.zip
echo "Unpacking expression files..."
unzip reverse.zip
fi

#write R script to get ENSEMBL id, needs biomaRt in R

echo "#!/usr/bin/Rscript
library(biomaRt)
listMarts(host=\"grch37.ensembl.org\")
ensembl = useMart(\"ENSEMBL_MART_ENSEMBL\",dataset=\"hsapiens_gene_ensembl\", host=\"grch37.ensembl.org\")
id_merge = getBM(attributes=c(\"ensembl_gene_id\",\"external_gene_name\"),mart=ensembl)
write.table(id_merge, file=\"id_merge.txt\", sep = \"\t\", quote =F, col.names=F, row.names=F)
" > script.r

#run R script

chmod 775 script.r
./script.r

#Use awk to append gene names

awk 'NR==FNR {h[$2] = $1; h2[$2] = $2; next} {print h[$1]}' id_merge.txt $GENE > genename

#remove temporary files

rm id_merge.txt
rm script.r

#get gene counts for gene of interest

while read line; do
                set $line
                find . -name *gene.count | xargs grep $line > COUNTS.txt
done < genename

sed -E "s/^.\/reverse\///g" COUNTS.txt | sed -E "s/\/.*\.[0-9]*//g" > COUNTS.txt.cut

#get total gene counts per sample

find . -name *gene.count | xargs -I % awk 'BEGIN {FS = " "} ; {sum+=$2} END {print sum}' % > TOTAL.txt
find . -name *gene.count | xargs -I % echo % > SAMPLES.txt 

sed -E "s/^.\/reverse\///g" SAMPLES.txt | sed -E "s/\/.*\.[0-9a-zA-Z]*//g" > SAMPLES.txt.cut

paste SAMPLES.txt.cut TOTAL.txt > TOTALCOUNTS.txt

awk 'NR==FNR {h[$1] = $0; next} {if(h[$1]) print h[$1]"\t"$0}' COUNTS.txt.cut TOTALCOUNTS.txt > TABLE.txt

awk '{print $1 "\t" $2/$4*1000000}' TABLE.txt > TABLE.RPM.txt

cp TABLE.RPM.txt $VCF

cd $VCF

#create ggplot2 graph 

echo "#!/usr/bin/Rscript
library(\"ggplot2\")
data<-read.table (file=\"GENOTYPES.combined.even.HEADER\", head=T, check.names=F)
data<-rbind(data, Total = colSums(data))

data.tr<-t(data)
#data.sum<-rowSums(data.tr)
expr<-read.table(\"TABLE.RPM.txt\", row.names=1, check.names=F)

data.merge<-merge(x = data.tr, y = expr, by = \"row.names\", all = F)

write.table(file=\"data.merge.table.txt\", data.merge)


x<-read.table(\"GENOTYPES.combined.even.b.header\", header=T, check.names=F, colClasses='character')
#m<-matrix(, ncol = 58, nrow = length(x[,1]))
#for (i in 1:58){m[,i]<-ifelse(x[,1]==x[,i],1,0)}
#n<-matrix(,ncol= 58, nrow=1)
#for (i in 1:58){n[,i]<-ifelse(sum(m[,i])==nrow(m),1,0)}

        m<-matrix(, ncol = ncol(x), nrow = length(x[,1]))
        n<-matrix(,ncol= ncol(x), nrow=ncol(x))
        p<-matrix(,ncol= ncol(x), nrow=ncol(x))

for (j in 1:ncol(x)){
        for (i in 1:ncol(x)){m[,i]<-ifelse(x[,j]==x[,i],1,0)}
        for (i in 1:ncol(x)){
                n[j,i]<-ifelse(sum(m[,i])==nrow(m),1,0) 
        }
}

for (k in 1:ncol(x)){
	for (i in 1:ncol(x)){
		p[k,i]<-ifelse(n[k,i]==1,colnames(x)[i],NA)}
		#print(colnames(x)[i])}
}

p.dat<-as.data.frame(p)

pdf(\"correlation_matrix.pdf\")
require(lattice)
levelplot(n)
dev.off()

listing<-list()

for (l in 1:ncol(x)){
	listing[[paste0(\"k.\", l)]]<-p.dat[l,][!is.na(p.dat[l,])]
}

listing2<-unique(listing)

isEmpty <- function(x) {
    return(length(x)==0)
}

k=0

plot<-matrix(,ncol= 2, nrow=length(listing2))
for (l in 1:length(listing2)){
#listing2[[l]][1]
	for (o in 1:length(listing2[[l]])) {if (!isEmpty(data.merge[data.merge\$Row.names==listing2[[l]][o],]\$Total)) {k=1}}

		for (o in 1:length(listing2[[l]])){ 
			if (k==1 && !isEmpty(data.merge[data.merge\$Row.names==listing2[[l]][o],]\$Total)) {
				plot[l,1]<-data.merge[data.merge\$Row.names==listing2[[l]][o],]\$Total}

}
	if (k==1) {plot[l,2]<-mean(data.merge[data.merge\$Row.names %in% listing2[[l]],]\$V2)}
	k=0
}

plot2<-na.omit(plot)
plot3<-as.data.frame(plot2)

colnames(plot3)<-c(\"Total\", \"V2\")


pdf(\"output.pdf\")
ggplot(plot3, aes(Total, V2, color = Total)) +
  geom_point(shape = 16, size = 5, show.legend = FALSE, alpha = .4) +
  theme_minimal() +
  scale_color_gradient(low = \"#0091ff\", high = \"#f0650e\") +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=26))+ labs(title = \"$GENENAME Causality Test\", x=\"CAD Risk SNPs Total\", y=\"$GENENAME expression\") +
  theme(plot.title = element_text(size = rel(2))) +
  geom_smooth(method=lm) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(plot3\$V2)+max(plot3\$V2)/5))+
  annotate(x=min(plot3\$Total)+(max(plot3\$Total)-min(plot3\$Total))/5, y=max(plot3\$V2)+(max(plot3\$V2)-min(plot3\$V2))/6,label=paste(\"R = \", round(cor(plot3\$V2,plot3\$Total),2)),geom=\"text\", size=8, col=\"darkred\")


  

dev.off()
"> script.R

chmod 775 script.R
./script.R


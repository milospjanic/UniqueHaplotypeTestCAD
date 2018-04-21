# UniqueHaplotypeTestCAD

UniqueHaplotypeTestCAD is a version of GeneCausalityTest for coronary artery disease that resolves local haplotype structure and averages on samples with unique haplotype structure. The script improves correlations made by GeneCausalityTest for defining causality of a gene for CAD, especially in the local regions of strong linkage disequilibrium. The script outputs the directionality of expression change with the increasing number of risk SNPs and uses CAD GWAS data (Nelson et al.) and HCASMC eQTL data for regression analysis.

User provides as arguments:

1. gene of interest
2. chromosome where the gene is located
3. defines p-value threshold to select SNPs from Nelson et al. meta analysis for CAD 
4. defines plus/minus region to select SNPs from, starting from composite gene coordinates obtained after collapsing all gene isoforms into a composite gene.

# Usage

To run the script download the .sh file
<pre>
wget https://raw.githubusercontent.com/milospjanic/UniqueHaplotypeTestCAD/master/UniqueHaplotypeTestCAD.sh
chmod 755 UniqueHaplotypeTestCAD.sh
</pre>

Place the script in your home or any other folder. The script will create ~/CADcausalitytest as its working folder, and three subfolders: CARDIOGRAMC4D, HCASMC_expr and HCASMC_genotypes. HCASMC_expr will contain per-gene RNAseq read counts for each HCASMC sample, while HCASMC_genotypes/vcf/ contains whole genome sequencing vcf files of HCASMC samples. CARDIOGRAMC4D folder will contain summary data from Nelson et al. Script will check if all three folders and data sets are present and if not download.

In the folder CARDIOGRAMC4D, file **chrN.GENE.region.txt** will contain all the SNPs from Nelson et al. in the selected region from the tested gene. File **SNP_effect.alele_pval.threshold.txt** will contain SNPs selected with a p-value threshold:

<pre>
rs7483886	G	2.13e-06
rs11532052	T	4.91e-06
rs11226029	G	3.52e-10
rs11226031	C	2.02e-06
...
</pre>

In the HCASMC_expr folder, file **TABLE.RPM.txt** contains per-gene RNAseq read counts for each HCASMC sample. In the HCASMC_genotypes folder, file GENOTYPES.combined.even.b.header contains risk SNP matrix with counts representing phased SNPs using 0 or 1 to mark risk alleles:

<pre>
1020301 102901 1042702 1051601 1060602 10705 112201 1278 1346 1347 1369 1386 1448 1483 1497 1522 1559 1576 1587 1596 177089 1795 1923 200212 2030801 2040401 20805 2102 2109 2115 2135 2139 2161 2228 2282 2305 2356 24156 2435 2463 2477 2510 289727 2989 3003 3100203 3101801 313605 59386145 59885590 7103002 8072501 8100901 9052004 9070202 9071501 9090701 CA1401
         11 11 11 11 11 11 10 11 11 11 11 11 11 10 11 11 11 10 10 10 11 11 11 11 10 11 11 11 10 11 11 10 11 11 11 11 10 11 11 11 11 11 11 11 00 10 10 10 11 11 10 11 11 11 11 10 10 10
         11 11 11 11 11 11 10 11 11 11 11 11 11 10 10 11 11 10 10 10 11 11 11 11 10 11 11 11 10 11 11 10 11 11 11 11 10 10 11 11 11 11 11 11 00 10 10 10 11 11 10 11 11 11 11 10 10 10
         11 11 11 11 11 11 10 10 10 11 11 11 11 10 11 11 11 10 10 10 11 10 11 11 10 11 11 11 10 11 11 10 11 11 11 11 10 10 11 11 11 10 11 11 00 10 10 11 11 11 10 11 11 11 11 10 10 10
         10 11 10 11 11 11 10 11 11 10 11 11 11 10 00 00 11 10 10 00 10 10 11 00 10 11 11 11 10 10 10 10 11 11 10 11 00 10 10 11 11 11 11 11 00 10 10 10 11 11 10 10 11 11 11 10 10 10
         10 11 10 11 11 11 10 11 11 10 11 11 11 10 00 10 11 10 10 00 10 10 11 00 10 11 11 11 10 10 10 10 11 11 10 11 00 10 10 11 11 11 11 11 00 10 10 10 11 11 10 10 11 11 11 10 10 10
         10 11 10 11 11 11 10 11 11 10 11 11 11 10 10 10 11 10 10 00 10 10 11 00 10 11 11 11 10 10 10 10 11 11 10 11 00 10 10 11 11 11 11 11 00 10 10 10 11 11 10 10 11 11 11 10 10 10
         10 11 10 11 11 11 10 11 11 10 11 11 11 10 00 10 11 10 10 00 10 10 11 00 10 11 11 11 10 10 10 10 11 11 10 11 00 10 10 11 11 11 11 11 00 10 10 10 11 11 10 10 11 11 11 10 10 10
         11 11 11 11 11 11 10 11 11 11 11 11 11 10 10 10 11 10 10 10 11 11 11 11 10 11 11 11 10 11 11 10 11 11 11 11 10 10 11 11 11 11 11 11 00 10 10 10 11 11 10 11 11 11 11 10 10 10
         11 11 11 11 11 11 10 11 11 11 11 11 11 10 10 11 11 10 10 10 11 11 11 11 10 11 11 11 10 11 11 10 11 11 11 11 10 10 11 11 11 11 11 11 00 10 10 10 11 11 10 11 11 11 11 10 10 10
         10 11 10 11 11 11 10 11 11 10 11 11 11 11 00 10 11 10 10 00 10 10 11 00 10 11 11 11 10 10 10 10 11 11 10 11 00 10 10 ...

</pre>

File GENOTYPES.combined.even.HEADER contains risk SNP matrix with counts representing 0, 1 or 2 risk SNPs:

<pre>

1020301 102901 1042702 1051601 1060602 10705 112201 1278 1346 1347 1369 1386 1448 1483 1497 1522 1559 1576 1587 1596 177089 1795 1923 200212 2030801 2040401 20805 2102 2109 2115 2135 2139 2161 2228 2282 2305 2356 24156 2435 2463 2477 2510 289727 2989 3003 3100203 3101801 313605 59386145 59885590 7103002 8072501 8100901 9052004 9070202 9071501 9090701 CA1401
         0 1 1 0 2 1 0 0 0 2 1 1 0 0 0 0 1 0 1 0 2 0 1 1 1 1 0 1 0 1 1 0 2 2 1 0 1 1 2 1 2 0 1 1 0 1 2 1 2 0 1 0 0 1 0 1 1 0
         0 1 0 0 0 1 0 0 0 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 2 0 2 0 1 1 0 1 1 1 2 0 1 0 0 0 0 1 0 0
         0 1 0 0 0 1 0 0 0 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 2 0 2 0 1 1 0 1 1 1 2 0 1 0 0 0 0 1 0 0
         0 1 0 0 0 1 0 0 0 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 2 0 2 0 1 1 0 1 1 1 2 0 1 0 0 0 0 1 0 0
         0 1 0 0 0 1 0 0 0 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 2 0 2 0 1 1 0 1 1 1 2 0 1 0 0 0 0 1 0 0
         0 1 1 0 2 1 0 0 0 2 1 1 0 0 0 0 1 0 1 0 2 0 1 1 1 1 0 1 0 1 1 0 2 2 1 0 1 1 2 1 2 0 1 1 0 1 2 1 2 0 1 0 0 1 0 1 1 0
         0 1 0 0 0 1 0 0 0 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 2 0 2 0 1 1 0 1 1 1 2 0 1 0 0 0 0 1 0 0
         0 1 0 0 0 1 0 0 0 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 2 0 2 0 1 1 0 1 1 1 2 0 1 0 0 0 0 1 0 0
         0 1 0 0 0 1 0 0 0 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 2 0 2 0 1 1 0 1 1 1 2 0 1 0 0 0 0 1 0 0
</pre>

The sample correlation matrix representing the haplotype structure for the tested region is calculated as a 0/1 matrix with 1 representing full correlation between haplotypes and outputed as **correlation_matrix.pdf** in the folder **HCASMC_genotypes**.

![alt text](https://github.com/milospjanic/UniqueHaplotypeTestCAD/blob/master/correlation_matrix.png)

Running the script. Place gene name in a file gene.txt, and provide gene.txt, chromosome, distance and p-value threshold as arguments:

<pre>
./UniqueHaplotypeTestCAD.sh gene.txt 11 100000 0.00001
</pre>

# Dependences

UniqueHaplotypeTestCAD requires biomaRt installed in R. In case the biomaRt repository is down run this command to substitute the repository URL with the archive version of biomaRt:

<pre>
sed -i 's/grch37.ensembl.org/jul2016.archive.ensembl.org/g' UniqueHaplotypeTestCAD.sh
</pre>

# Examples

By resolving the local haplotype structure, UniqueHaplotypeTestCAD might improve correlation of gene expression with the increased number of risk SNPs made by GeneCausalityTest. This is particularly evident in regions of strong linkage disequilibrium where SNPs are co-inherited in blocks and passed through generations in groups. 

Here is an example of a gene that is positively assotiated with the CAD trait (i.e. with the CAD risk allele count), and the correlation that GeneCausalityTest reports is 0.25, probably lower due to the strong linkage disequilibrium in the region that is evident from the graph as samples tend to pile up in groups, due to co-inheritence of risk SNPs. UniqueHaplotypeTestCAD will resolve the local haplotype structure and improve correlation from 0.25 to 0.45:

![alt text](https://github.com/milospjanic/UniqueHaplotypeTestCAD/blob/master/gene1_cad_causality_test.png)

After resolving the local haplotype structure by UniqueHaplotypeTestCAD:

![alt text](https://github.com/milospjanic/UniqueHaplotypeTestCAD/blob/master/gene1_unique_haplotype_test.png)

Another example is a gene that is negatively assotiated with the CAD trait (i.e. with the CAD risk allele count). The correlation that GeneCausalityTest reports is 0.11, while UniqueHaplotypeTestCAD improves the correlation to 0.25 by resolving the local haplotype structure of this locus with strong linkage disequilibrium.

![alt text](https://github.com/milospjanic/UniqueHaplotypeTestCAD/blob/master/gene2_cad_causality_test.png)

After UniqueHaplotypeTestCAD correction:

![alt text](https://github.com/milospjanic/UniqueHaplotypeTestCAD/blob/master/gene2_unique_haplotype-test.png)


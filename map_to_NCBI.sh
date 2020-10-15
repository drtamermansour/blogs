#### Maps of Entrez gene ids

## Generation of a map between NCBI IDs (Entrez gene ids) to the official NCBI symbols:
if [ ! -f Homo_sapiens.gene_info ];then
  wget ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz  ## March 24, 2020
  gunzip Homo_sapiens.gene_info.gz
else echo "Homo_sapiens.gene_info file exists in the geneRepo";fi

#cat Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{print $2,$3,$5}' > Homo_sapiens.gene_info.ID_to_symbol ## 61631
head -n1 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t"}{print $2,$3,"HGNC","Ensembl","MIM",$5}' > Homo_sapiens.gene_info.ID_to_symbol #Homo_sapiens.gene_info.map
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t"}{if($6!="-")print $2,$3,$6,$5}' | awk 'BEGIN{FS=OFS="\t"}{ delete vars; split($3,a,"|");for(i in a) { n = index(a[i], ":"); if(n) { x = substr(a[i], n + 1); key = substr(a[i], 1, n - 1); val = substr(a[i], n + 1, length(x)); if(vars[key]=="")vars[key] = val;else vars[key] = vars[key]","val; } } MIM = vars["MIM"]; HGNC = vars["HGNC"]; Ensembl = vars["Ensembl"]; print $1,$2,HGNC,Ensembl,MIM,$4; }' >> Homo_sapiens.gene_info.ID_to_symbol
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t"}{if($6=="-")print $2,$3,"","","",$5}' >> Homo_sapiens.gene_info.ID_to_symbol

###############################################
#### Maps of Entrez gene ids to several databases through HGNC

## Generation of a map of NCBI info (Entrez gene ids, official NCBI symbols, Synonyms) to HGNC info (IDS, symbols, prev_symbol)
if [ ! -f hgnc_complete_set.txt ];then
  wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt # HGNC dataset from the "Statistics & download files" page
else echo "hgnc_complete_set.txt file exists in the geneRepo";fi
cat hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{if($19!="")print $19,$1,$6,$2,$9,$11,$20,$21,$22,$23,$24,$25,$26}' | sed 's/entrez_id/GeneID/' > hgnc_complete_set.simplified ## 42189
awk 'BEGIN{FS=OFS="\t";}FNR==NR{key=$1;sub($1 FS,"");a[key]=$0;next}{if(a[$1])print $0,a[$1];else{$18="";print}}' hgnc_complete_set.simplified Homo_sapiens.gene_info.ID_to_symbol > Homo_sapiens.gene_info.ID_to_symbol.ext

### Assessemnt
## How many HGNC ids lack NCBI IDs in the HGNC database
cat hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{if($19=="")print}' | wc -l ## 1812 (approved + withdrawn)
cat hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{if($19=="" && $6=="Approved")print}' | wc -l ## 49
## How many HGNC ids have discontinued NCBI IDs and thus would fail to show in "Homo_sapiens.gene_info.ID_to_symbol.ext"
awk 'BEGIN{FS=OFS="\t";}FNR==NR{a[$7]=1;next}{if(a[$2]!=1)print $0}' Homo_sapiens.gene_info.ID_to_symbol.ext hgnc_complete_set.simplified > hgnc_complete_set.simplified.unused ## 4
## Do we have HGNC non-approved records in "Homo_sapiens.gene_info.ID_to_symbol.ext"
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ext | awk -F"\t" '{if($7!="" && $8!="Approved")print}' | wc -l  ## 0
## diff between NCBI and HGNC in NCBI to HGNC mapping
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ext | awk -F"\t" '{if($3!=$7)print}' | wc -l ## 8

###############################################
#### Add ENS info to "Homo_sapiens.gene_info.ID_to_symbol.ext"
if [ ! -f gencode.v35.annotation.gtf ];then
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz
  gunzip gencode.v35.annotation.gtf.gz
else echo "gencode.v35.annotation.gtf file exists in the geneRepo";fi
echo "ens_gene_id.gene_version ens_transcript_id symbol hgnc_id havana_gene" | tr ' ' '\t' > gencode.v35.trans.annotation
cat gencode.v35.annotation.gtf | awk -F"\t" '!/#/{if($3=="transcript")print $9}' | sed 's/; /;/g' | sed 's/\"//g' | awk -F";" 'BEGIN{FS=";";OFS="\t"}{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, " "); if(n) { x = substr($i, n + 1); vars[substr($i, 1, n - 1)] = substr($i, n + 1, length(x)) } } id = vars["gene_id"]; trans = vars["transcript_id"]; name = vars["gene_name"]; hgnc = vars["hgnc_id"]; hav = vars["havana_gene"]; print id,trans,name,hgnc,hav; }' >> gencode.v35.trans.annotation

if [ ! -f gencode.v35.metadata.EntrezGene ];then
  wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.metadata.EntrezGene.gz
  gunzip gencode.v35.metadata.EntrezGene.gz
else echo "gencode.v35.metadata.EntrezGene file exists in the geneRepo";fi
head -n1 gencode.v35.trans.annotation | awk 'BEGIN{FS=OFS="\t"}{print "GeneID",$0}' > gencode.v35.trans.annotation.EntrezGene
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{ print a[$2],$0}' gencode.v35.metadata.EntrezGene <(tail -n+2 gencode.v35.trans.annotation) >> gencode.v35.trans.annotation.EntrezGene
head -n1 gencode.v35.trans.annotation.EntrezGene | awk 'BEGIN{FS=OFS="\t"}{sub(/\./,FS,$2);print}' > gencode.v35.gene.annotation.EntrezGene
tail -n+2 gencode.v35.trans.annotation.EntrezGene | awk 'BEGIN{FS=OFS="\t"}{sub(/\./,FS,$2);print}' | sort -k2,2 -t$'\t' | sort -u -k2,2 --merge -t$'\t' >> gencode.v35.gene.annotation.EntrezGene

## join by the HGNC id
awk 'BEGIN{FS=OFS="\t";}FNR==NR{if($6){key=$6;sub($6 FS,"");a[key]=$0;}next}{if(a[$7])print $0,a[$7];else{$24="";print}}' gencode.v35.gene.annotation.EntrezGene Homo_sapiens.gene_info.ID_to_symbol.ext > Homo_sapiens.gene_info.ID_to_symbol.ext2 ## 61677 ## key columns: 1.GeneID 2.Symbol 3.HGNC 4.Ensembl 7.hgnc_id 9.symbol 12.ensembl_gene_id 19.GeneID 20.ens_gene_id 23.symbol


## Find the best consensus of ENS IDs across the 3 databases and assign score to the final ID
## Score 0 indicates no ENS to NCBI match in the 3 DBs. Score 1:3 indicates one, two, or three matching relations respectively 
## Score -1 indicates agreement of 2 DBs but with disagreement with the 3rd. Finally, score -2 indicate either 2 or 3 contradictory matches  

## Approach 1
head -n1 Homo_sapiens.gene_info.ID_to_symbol.ext2 | awk -F"\t" 'BEGIN{FS=OFS="\t"}{print $0,"ens_cons","ens_score" }' > Homo_sapiens.gene_info.ID_to_symbol.ens
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ext2 |  awk -F"\t" 'BEGIN{FS=OFS="\t"}{delete a; if($4!="")split($4,a,",");else a[1]=""; \
for(i in a){$4=a[i];delete c;if(a[i])++c[a[i]];if($12)++c[$12];if($20)++c[$20];l=length(c); \
if(l==0){ens="";s=0;}else if(l==1){for(g in c){ens=g;s=c[g]}}else{ens="";s=-2;for(g in c){if(c[g]>1){ens=g;s=-1;}}} \
print $0,ens,s }}' | sort -k1,1 -k26,26nr -t$'\t' | sort -u -k1,1 --merge -t$'\t' >> Homo_sapiens.gene_info.ID_to_symbol.ens ## 61677
# summary of scores:    31 -2, 61 -1, 19974 0, 2887 1, 6474 2, 32249 3
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{print $26}' | sort | uniq -c | sort -k2,2n

head -n1 Homo_sapiens.gene_info.ID_to_symbol.ens | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$25,$26}' > Homo_sapiens.gene_info.ID_to_symbol.ens.sim
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk 'BEGIN{FS=OFS="\t"}{if($25!="")print $1,$2,$25,$26}' | sort -k1,1 -k3,3 -k4,4nr -t$'\t' | sort -u -k1,1 -k3,3 --merge -t$'\t' >> Homo_sapiens.gene_info.ID_to_symbol.ens.sim ## 41672 (41672 before dedup)

# select the best match to ENS IDs with duplicate enteries
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens.sim | awk 'BEGIN{FS=OFS="\t"}{print $3}' | sort | uniq -c | awk '{if($1>1)print}' | wc -l ## 103
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens.sim | awk 'BEGIN{FS=OFS="\t"}{print $3}' | sort | uniq -c | awk '{if($1>1)print $2}' | grep -Fwf - Homo_sapiens.gene_info.ID_to_symbol.ens2.sim | sort -k3,3 -k4,4nr -t$'\t' > dup_ens_ids ## 208
cat dup_ens_ids | sort -u -k3,3 --merge -t$'\t' > dup_ens_ids.keep ## 103
comm -23 <(sort dup_ens_ids) <(sort dup_ens_ids.keep) > dup_ens_ids.discard ## 105
head -n1 Homo_sapiens.gene_info.ID_to_symbol.ens.sim > Homo_sapiens.gene_info.ID_to_symbol.ens.sim2 ## 41567
comm -23 <(tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens.sim | sort) <(sort dup_ens_ids.discard) >> Homo_sapiens.gene_info.ID_to_symbol.ens.sim2 ## 41567
# check for NCBI with multiple ENS IDs
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens.sim2 | awk 'BEGIN{FS=OFS="\t"}{print $1}' | sort | uniq -c | awk '{if($1>1)print}' | wc -l ## 0


### Assessement of Approach 1
#NCBI ID might have multiple ens IDs but HGNC ID has max one ens ID. On the other hand, ens ID might have multiple NCBI and HGNC IDs
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($4!="" && $12!="" && $4==$12)print}' | wc -l  ## similar ens IDs after mapping in NCBI and HGNC ## 32516
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($4=="" && $12!="")print}' | wc -l ## missing ens IDs in NCBI but found in HGNC ## 6653
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($4!="" && $12=="")print}' | wc -l ## missing ens IDs in HGNC but found in NCBI ## 2458
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($4!="" && $12!="" && $4!=$12)print}' > ens_diff ## diff in ens IDs after mapping in NCBI and HGNC ## 75
cat ens_diff | awk -F"\t" '{print $4}' | awk -F"," '{if($2!="")print}' | wc -l ## records in ens_diff with multiple ens IDs in NCBI ## 0
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($12!="")print $12}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)print}' | wc -l ## ens IDs in HGNC that match > 1 record ## 3


tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($20!="" && $12!="" && $20==$12)print}' | wc -l  ## similar HGNC to ens mapping in the 2 databases ## 38504
  tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($20!="" && $12!="" && $20==$12 && $4==$12)print}' | wc -l  ## how many of those are matching ens ids in NCBI ## 32249
  tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($20!="" && $12!="" && $20==$12 && $4=="")print}' | wc -l  ## how many of those do not have ens ids in NCBI ## 6217
  tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($20!="" && $12!="" && $20==$12 && $4!="" && $4!=$12)print}' | wc -l ## how many of those are NOT matching ens ids in NCBI ## 38

tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($20=="" && $12!="")print}' | wc -l ## HGNC has a match to ens reported in HGNC but not in ens gtf ## 710
 tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($20=="" && $12!="" && $4==$12)print}' | wc -l ## how many of those are matching ens ids in NCBI ## 257
 tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($20=="" && $12!="" && $4=="")print}' | wc -l ## how many of those do not have ens ids in NCBI ## 429
 tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($20=="" && $12!="" && $4!="" && $4!=$12)print}' | wc -l ## how many of those are NOT matching ens ids in NCBI ## 24

tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($20!="" && $12=="")print}' | wc -l ## HGNC has a match to ens reported in ens gtf but not in HGNC ## 0
 tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($20!="" && $12=="" && $4==$20)print}' | wc -l ## how many of those are matching ens ids in NCBI ## 0
 tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($20!="" && $12=="" && $4=="")print}' | wc -l ## how many of those do not have matching ens ids in NCBI ## 0
 tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($20!="" && $12=="" && $4!="" && $4!=$20)print}' | wc -l ## how many of those are NOT matching ens ids in NCBI ## 0

tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens | awk -F"\t" '{if($20!="" && $12!="" && $20!=$12)print}' | wc -l  ## diff in HGNC to ens mapping between the 2 databases ## 30



## Approach 2
head -n1 Homo_sapiens.gene_info.ID_to_symbol.ext2 | awk -F"\t" 'BEGIN{FS=OFS="\t"}{print $0,"ens_cons","ens_score" }' > Homo_sapiens.gene_info.ID_to_symbol.ens2
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ext2 |  awk -F"\t" 'BEGIN{FS=OFS="\t"}{delete a; if($4!="")split($4,a,",");else a[1]=""; \
for(i in a){$4=a[i];delete c;if(a[i])++c[a[i]];if($12)++c[$12];if($20)++c[$20];l=length(c); \
if(l==0)print $0,"",0; else if(l==1){for(g in c)print $0,g,c[g];}else{for(g in c){if(c[g]>1)print $0,g,-1; else print $0,g,-2;}} \
}}' >> Homo_sapiens.gene_info.ID_to_symbol.ens2 ## 61892
# summary of scores:   184 -2,122 -1, 19974 0, 2888 1, 6474 2, 32249 3
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens2 | awk -F"\t" '{print $26}' | sort | uniq -c | sort -k2,2n

head -n1 Homo_sapiens.gene_info.ID_to_symbol.ens2 | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$25,$26}' > Homo_sapiens.gene_info.ID_to_symbol.ens2.sim
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens2 | awk 'BEGIN{FS=OFS="\t"}{if($25!="")print $1,$2,$25,$26}' | sort -k1,1 -k3,3 -k4,4nr -t$'\t' | sort -u -k1,1 -k3,3 --merge -t$'\t' >> Homo_sapiens.gene_info.ID_to_symbol.ens2.sim ## 41855 (41918 before dedup)

# select the best match to ENS IDs with duplicate enteries
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens2.sim | awk 'BEGIN{FS=OFS="\t"}{print $3}' | sort | uniq -c | awk '{if($1>1)print}' | wc -l ## 135
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens2.sim | awk 'BEGIN{FS=OFS="\t"}{print $3}' | sort | uniq -c | awk '{if($1>1)print $2}' | grep -Fwf - Homo_sapiens.gene_info.ID_to_symbol.ens2.sim | sort -k3,3 -k4,4nr -t$'\t' > dup_ens2_ids ## 272
cat dup_ens2_ids | sort -u -k3,3 --merge -t$'\t' > dup_ens2_ids.keep ## 135
comm -23 <(sort dup_ens2_ids) <(sort dup_ens2_ids.keep) > dup_ens2_ids.discard ## 137
head -n1 Homo_sapiens.gene_info.ID_to_symbol.ens2.sim > Homo_sapiens.gene_info.ID_to_symbol.ens2.sim2
comm -23 <(tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens2.sim | sort) <(sort dup_ens2_ids.discard) >> Homo_sapiens.gene_info.ID_to_symbol.ens2.sim2 ## 41718
# For NCBI with multiple ENS IDs, keep only one match with best score and isolate the ENS IDs with the lower score 
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens2.sim2 | awk 'BEGIN{FS=OFS="\t"}{print $1}' | sort | uniq -c | awk '{if($1>1)print}' | wc -l ## 120
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens2.sim2 | awk 'BEGIN{FS=OFS="\t"}{print $1}' | sort | uniq -c | awk '{if($1>1)print $2}' | grep -Fwf - Homo_sapiens.gene_info.ID_to_symbol.ens2.sim | sort -k1,1 -k4,4nr -t$'\t' > dup_ncbi2_ids ## 241
cat dup_ncbi2_ids | sort -u -k1,1 --merge -t$'\t' > dup_ncbi2_ids.keep ## 120
comm -23 <(sort dup_ncbi2_ids) <(sort dup_ncbi2_ids.keep) > dup_ncbi2_ids.discard ## 121
head -n1 Homo_sapiens.gene_info.ID_to_symbol.ens2.sim2 > Homo_sapiens.gene_info.ID_to_symbol.ens2.sim3
comm -23 <(tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens2.sim2 | sort) <(sort dup_ncbi2_ids.discard) >> Homo_sapiens.gene_info.ID_to_symbol.ens2.sim3 ## 41597
# summary of scores:   118 -2, 2 -1, 1 1
cat dup_ncbi2_ids.discard | awk -F"\t" '{print $4}' | sort | uniq -c | sort -k2,2n
# summary of scores:   30 -2, 61 -1, 2783 1, 6473 2, 32249 3
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol.ens2.sim3 | awk -F"\t" '{print $4}' | sort | uniq -c | sort -k2,2n







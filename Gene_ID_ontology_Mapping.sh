#### 1. Ambiguity of human gene symbols within the same database:

#### HGNC on 8/30/2020
## Download HGNC dataset (link in the "Statistics & download files" page)
wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt
tail -n+2 hgnc_complete_set.txt | wc -l ## 43942

## Explore
# gene IDs
tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS="\t";}{print $1}' | sort | uniq | wc -l ## 43942
# Symbols
tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS="\t";}{print $2}' | sort | uniq | wc -l ## 43942
# Approved vs Withdrawn
tail -n+2 hgnc_complete_set.txt | awk -F"\t" '{print $6}' | sort | uniq -c ## 42181 Approved, 1761 Withdrawn

## Generate a map of gene IDs to symbols and list of alias and previous symbols
cat hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$6,$9,$11}' | sed 's/"//g' > hgnc.ID_to_symbol

## Generate a map of gene IDs to symbols and each one of the alias
head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$9}' > hgnc.ID_to_symbol_to_EachAlias
tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{if($9!="")print $1,$2,$9}' | sed 's/"//g' | awk 'BEGIN{FS="\t";OFS="\n";}{split($3,a,"|");for(i in a)print $1"\t"$2"\t"a[i];}' >> hgnc.ID_to_symbol_to_EachAlias  ## 42451

head -n1 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$11}' > hgnc.ID_to_symbol_to_EachPrev
tail -n+2 hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{if($11!="")print $1,$2,$11}' | sed 's/"//g' | awk 'BEGIN{FS="\t";OFS="\n";}{split($3,a,"|");for(i in a)print $1"\t"$2"\t"a[i];}' >> hgnc.ID_to_symbol_to_EachPrev  ## 15166

## Generate a map of all ambiguous gene symbols: Genes with ambiguous alias  (the alias is ambiguous if it matches another alias, previous or current gene symbol)
# create list of unique gene symbols
tail -n+2 hgnc_complete_set.txt | awk -F"\t" '{print $2}' | sort | uniq > hgnc.Symbols  ## 43942
# create list of all alias symbols
tail -n+2 hgnc.ID_to_symbol_to_EachAlias | awk -F "\t" '{print $3}' > hgnc.Alias ## 42451
# create list of all previous symbols
tail -n+2 hgnc.ID_to_symbol_to_EachPrev | awk -F "\t" '{print $3}' > hgnc.Prev ## 15166
# Identify ambiguous alias or previous symbols
cat hgnc.Symbols hgnc.Alias hgnc.Prev | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > hgnc.Alias_amb  ## 4668
cat hgnc.Alias_amb | sed -r "s/^(.*)$/\t\1$/" | grep -f - hgnc.ID_to_symbol_to_EachAlias > hgnc.ID_to_symbol_to_EachAlias_amb_symbols ## 3108
cat hgnc.Alias_amb | sed -r "s/^(.*)$/\t\1$/" | grep -f - hgnc.ID_to_symbol_to_EachPrev > hgnc.ID_to_symbol_to_EachPrev_amb_symbols ## 811
#cat hgnc.Symbols hgnc.Alias hgnc.Prev | sort | uniq -c | awk '{if($1>1){print $0}}' | sort -k1,1nr > temp
#echo "HAP1" | sed -r "s/^(.*)$/\t\1\t/" | grep -f - hgnc.ID_to_symbol
#echo "HAP1" | sed -r "s/^(.*)$/\t\1$/" | grep -f - hgnc.ID_to_symbol_to_EachAlias
#echo "HAP1" | sed -r "s/^(.*)$/\t\1$/" | grep -f - hgnc.ID_to_symbol_to_EachPrev


#### NCBI genes
## Download
wget ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
gunzip Homo_sapiens.gene_info.gz

## Explore
# gene IDs
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{print $2}' | sort | uniq | wc -l ## 61622
# Symbols
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{print $3}' | sort | uniq | wc -l ## 61563

## Generate a map of gene IDs to symbols and list of alias
cat Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{print $2,$3,$5}' > Homo_sapiens.gene_info.ID_to_symbol

## Generate a map of all ambiguous gene symbols: A) Same gene symbol with multiple Entrez gene ids
tail -n+2 Homo_sapiens.gene_info | awk -F"\t" '{print $3}' | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Homo_sapiens.gene_info.Symbols_amb ## 41
cat Homo_sapiens.gene_info.Symbols_amb | sed -r "s/^(.*)$/\t\1\t/" | grep -f - Homo_sapiens.gene_info.ID_to_symbol > Homo_sapiens.gene_info.ID_to_symbol_amb_symbols  ## 100

## Generate a map of gene IDs to symbols and each one of the alias
head -n1  Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{print $2,$3,$5}' > Homo_sapiens.gene_info.ID_to_symbol_to_EachAlias
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t";}{if($5!="-")print $2,$3,$5}' | awk 'BEGIN{FS="\t";OFS="\n";}{split($3,a,"|");for(i in a)print $1"\t"$2"\t"a[i];}' >> Homo_sapiens.gene_info.ID_to_symbol_to_EachAlias  ## 70268

## Generate a map of all ambiguous gene symbols: B) Genes with ambiguous alias  (the alias is ambiguous if it matches another alias or another gene symbol)
# create list of unique gene symbols
tail -n+2 Homo_sapiens.gene_info | awk -F"\t" '{print $3}' | sort | uniq > Homo_sapiens.gene_info.Symbols  ## 61563
# create list of all alias symbols
tail -n+2 Homo_sapiens.gene_info.ID_to_symbol_to_EachAlias | awk -F "\t" '{print $3}' > Homo_sapiens.gene_info.Alias ## 70268
# Identify ambiguous alias
cat Homo_sapiens.gene_info.Symbols Homo_sapiens.gene_info.Alias | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > Homo_sapiens.gene_info.Alias_amb  ## 4668
cat Homo_sapiens.gene_info.Alias_amb | sed -r "s/^(.*)$/\t\1$/" | grep -f - Homo_sapiens.gene_info.ID_to_symbol_to_EachAlias > Homo_sapiens.gene_info.ID_to_symbol_to_EachAlias_amb_symbols ## 9396
#cat Homo_sapiens.gene_info.Symbols Homo_sapiens.gene_info.Alias | sort | uniq -c | awk '{if($1>1){print $0}}' | sort -k1,1nr > temp
#echo "HOX4" | sed -r "s/^(.*)$/\t\1\t/" | grep -f - Homo_sapiens.gene_info.ID_to_symbol
#echo "HOX4" | sed -r "s/^(.*)$/\t\1$/" | grep -f - Homo_sapiens.gene_info.ID_to_symbol_to_EachAlias


#### Gencode_human
## Download
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz
gunzip gencode.v35.annotation.gtf.gz

## Explore
cat gencode.v35.annotation.gtf | awk -F"\t" '!/#/{if($3=="gene")a+=1;}END{print "gene count=",a}'  ## gene count= 60656

## Generate a map of gene IDs to symbols and list of alias
echo "gene_id gene_type gene_name hgnc_id havana_gene" | tr ' ' '\t' > gencode.v35.gene.annotation
cat gencode.v35.annotation.gtf | awk -F"\t" '!/#/{if($3=="gene")print $9}' | sed 's/; /;/g' | sed 's/\"//g' | awk -F";" 'BEGIN{FS=";";OFS="\t"}{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, " "); if(n) { x = substr($i, n + 1); key = substr($i, 1, n - 1); val = substr($i, n + 1, length(x));if(vars[key]=="")vars[key] = val;else vars[key] = vars[key]","val;} } id = vars["gene_id"]; ann = vars["gene_type"]; name = vars["gene_name"]; hgnc = vars["hgnc_id"]; hav = vars["havana_gene"]; print id,ann,name,hgnc,hav; }' >> gencode.v35.gene.annotation
tail -n+2 gencode.v35.gene.annotation | awk -F"\t" '{print $3}' | wc -l ## 60656
tail -n+2 gencode.v35.gene.annotation | awk -F"\t" '{print $3}' | sort | uniq | wc -l ## 59609

## Generate a map of all ambiguous gene symbols: A) Same gene symbol with multiple Entrez gene ids
tail -n+2 gencode.v35.gene.annotation | awk -F"\t" '{print $3}' | sort | uniq -c | awk '{if($1>1){$1="";print $0}}' | sed 's/ //' > gencode.Symbols_amb ## 119
cat gencode.Symbols_amb | sed -r "s/^(.*)$/\t\1\t/" | grep -f - gencode.v35.gene.annotation > gencode.ID_to_symbol_amb_symbols  ## 1166

## compare to version 7
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_7/gencode.v7.annotation.gtf.gz (Dec 2010)
gunzip gencode.v7.annotation.gtf.gz
echo "gene_id gene_type gene_name hgnc_id havana_gene" | tr ' ' '\t' > gencode.v7.gene.annotation
cat gencode.v7.annotation.gtf | awk -F"\t" '!/#/{if($3=="gene")print $9}' | sed 's/; /;/g' | sed 's/\"//g' | awk -F";" 'BEGIN{FS=";";OFS="\t"}{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, " "); if(n) { x = substr($i, n + 1); key = substr($i, 1, n - 1); val = substr($i, n + 1, length(x));if(vars[key]=="")vars[key] = val;else vars[key] = vars[key]","val;} } id = vars["gene_id"]; ann = vars["gene_type"]; name = vars["gene_name"]; hgnc = vars["hgnc_id"]; hav = vars["havana_gene"]; print id,ann,name,hgnc,hav; }' >> gencode.v7.gene.annotation

head -n1 gencode.v35.gene.annotation | awk 'BEGIN{FS=OFS="\t"}{$1=$1 FS "version";print}' > gencode.v35
tail -n+2 gencode.v35.gene.annotation | sed 's/\./|/' | tr "|" "\t" >> gencode.v35
head -n1 gencode.v7.gene.annotation | awk 'BEGIN{FS=OFS="\t"}{$1=$1 FS "version";print}' > gencode.v7
tail -n+2 gencode.v7.gene.annotation | sed 's/\./|/' | tr "|" "\t" >> gencode.v7

tail -n+2 gencode.v35 | wc -l ## 60656
tail -n+2 gencode.v7 | wc -l  ## 51082
comm -12 <(tail -n+2 gencode.v35 | awk -F"\t" '{print $1}' | sort) <(tail -n+2 gencode.v7 | awk -F"\t" '{print $1}' | sort) | wc -l ## 43552 shared IDs
comm -23 <(tail -n+2 gencode.v35 | awk -F"\t" '{print $1}' | sort) <(tail -n+2 gencode.v7 | awk -F"\t" '{print $1}' | sort) | wc -l ## 17104 novel IDs
comm -13 <(tail -n+2 gencode.v35 | awk -F"\t" '{print $1}' | sort) <(tail -n+2 gencode.v7 | awk -F"\t" '{print $1}' | sort) | wc -l ## 7530 discontinued IDs

echo "gene_id version.v7 symbol.v7 version.v35 symbol.v35" | tr " " "\t" > gencode.v7.vs.v35
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2 FS $4;next;}{if(a[$1]!="")print $1,$2,$4,a[$1]}' <(tail -n+2 gencode.v35) <(tail -n+2 gencode.v7) >> gencode.v7.vs.v35
tail -n+2 gencode.v7.vs.v35 | awk -F"\t" '{if($3!=$5)a+=1}END{print a}' ## 22402

#####################
#### 2) Different gene annotations (NCBI vs Ensemble/Gencode vs HGNC) & completeness of cross-reference

#### HGNC on 8/30/2020
## Explore
tail -n+2 hgnc_complete_set.txt | awk -F"\t" '{if($6=="Approved" && $19!="")print $19}' | sort | uniq | wc -l ## 42110 with Entrez ID

## Generate a map of all HGNC gene IDs to dbXrefs
cat hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$19,$20,$9,$11}' > hgnc.map
tail -n+2 hgnc.map | awk -F"\t" '{a+=1;if($3!="")b+=1;if($4!="")c+=1;if($5!="")d+=1;if($6!="")e+=1;}END{print "hgnc_id=",a,"entrez_id=",b,"ensembl_id=",c,"alias_symbol=",d,"prev_symbol=",e}'  ## hgnc_id= 43942 entrez_id= 42110 ensembl_id= 39199 alias_symbol= 22034 prev_symbol= 12193
tail -n+2 hgnc.map | awk -F"\t" '{if($3!="" && $4!="")a+=1;}END{print "hgnc_id with entrez and ensembl ids=",a}'  ## hgnc_id with entrez and ensembl ids= 39167

## Generate a map of approved HGNC gene IDs to dbXrefs
cat hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{print $1,$2,$19,$20,$9,$11}' > hgnc_approved.map ## hgnc_id symbol  entrez_id       ensembl_gene_id alias_symbol    prev_symbol
cat hgnc_complete_set.txt | awk 'BEGIN{FS=OFS="\t";}{if($6=="Approved")print $1,$2,$19,$20,$9,$11}' >> hgnc_approved.map

## Stats
tail -n+2 hgnc_approved.map | awk -F"\t" '{a+=1;if($3!="")b+=1;if($4!="")c+=1;if($5!="")d+=1;if($6!="")e+=1;}END{print "hgnc_id=",a,"entrez_id=",b,"ensembl_id=",c,"alias_symbol=",d,"prev_symbol=",e}'  ## hgnc_id= 42181 entrez_id= 42110 ensembl_id= 39199 alias_symbol= 21740 prev_symbol= 12049
tail -n+2 hgnc_approved.map | awk -F"\t" '{if($3!="" && $4!="")a+=1;}END{print "hgnc_id with entrez and ensembl ids=",a}'  ## hgnc_id with entrez and ensembl ids= 39167
#tail -n+2 hgnc_approved.map | awk -F"\t" '{a[$2]=1;if($3!="")b[$3]=1;if($4!="")c[$4]=1;}END{print "uniq_symb=",length(a),"uniq_entrez=",length(b),"uniq_ens=",length(c)}' ## uniq_symb= 42181 uniq_entrez= 42110 uniq_ens= 39196


#### NCBI (Entrez gene ids from Homo_sapiens.gene_info which does not have the discontinued genes)
## Explore
## col 10 = type_of_gene, 11 = Symbol_from_nomenclature_authority, 13 = Nomenclature_status
# Explore types of genes
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{print $10}' | sort | uniq -c
# genes with non authorized names
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{if($3!=$11)print;}' > non_authority_naming #19553
# Classification of genes based on naming source
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{print $13}' | sort | uniq -c ## 19516 - && 42106 O
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{if($13=="O")print;}' > Nomenclature_status_O #42106
grep "HGNC:" Nomenclature_status_O | wc -l     ## 42106
grep "Ensembl:" Nomenclature_status_O | wc -l  ## 33165
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{if($13=="O")print $10}' | sort | uniq -c
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{if($13=="O" && $3!=$11)print}' | wc -l # 37
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS="\t";}{if($13!="O")print;}' > Nomenclature_status_NA #42106
grep "HGNC:" Nomenclature_status_NA | wc -l     ## 0
grep "Ensembl:" Nomenclature_status_NA | wc -l  ## 1898

## Generate a map of gene IDs to dbXrefs
head -n1 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t"}{print $2,$3,"HGNC","Ensembl","MIM",$5}' > Homo_sapiens.gene_info.map
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t"}{if($6!="-")print $2,$3,$6,$5}' | awk 'BEGIN{FS=OFS="\t"}{ delete vars; split($3,a,"|");for(i in a) { n = index(a[i], ":"); if(n) { x = substr(a[i], n + 1); key = substr(a[i], 1, n - 1); val = substr(a[i], n + 1, length(x)); if(vars[key]=="")vars[key] = val;else vars[key] = vars[key]","val; } } MIM = vars["MIM"]; HGNC = vars["HGNC"]; Ensembl = vars["Ensembl"]; print $1,$2,HGNC,Ensembl,MIM,$4; }' >> Homo_sapiens.gene_info.map
tail -n+2 Homo_sapiens.gene_info | awk 'BEGIN{FS=OFS="\t"}{if($6=="-")print $2,$3,"","","",$5}' >> Homo_sapiens.gene_info.map

## Stats
tail -n+2 Homo_sapiens.gene_info.map | awk -F"\t" '{a+=1;if($3!="")b+=1;if($4!="")c+=1;if($5!="")d+=1;if($6!="-")e+=1;}END{print "entrez_id=",a,"hgnc_id=",b,"ensembl_id=",c,"mim_id=",d,"alias_symbol=",e}'  ## entrez_id= 61622 hgnc_id= 42106 ensembl_id= 35063 mim_id= 17648 alias_symbol= 26915
tail -n+2 Homo_sapiens.gene_info.map | awk -F"\t" '{if($3!="" && $4!="")a+=1;}END{print "entrez_id with hgnc and ensembl ids=",a}'  ## entrez_id with hgnc and ensembl ids= 33165
#tail -n+2 Homo_sapiens.gene_info.map | awk -F"\t" '{a[$2]=1;if($3!="")b[$3]=1;if($4!="")c[$4]=1;}END{print "uniq_symb=",length(a),"uniq_hgnc=",length(b),"uniq_ens=",length(c)}' ## uniq_symb= 61563 uniq_hgnc= 42106 uniq_ens= 34966


#### Gencode_human
## Explore
#echo "gene_id gene_type gene_name hgnc_id havana_gene" | tr ' ' '\t' > gencode.v35.gene.annotation
#cat gencode.v35.annotation.gtf | awk -F"\t" '!/#/{if($3=="gene")print $9}' | sed 's/; /;/g' | sed 's/\"//g' | awk -F";" 'BEGIN{FS=";";OFS="\t"}{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, " "); if(n) { x = substr($i, n + 1); vars[substr($i, 1, n - 1)] = substr($i, n + 1, length(x)) } } id = vars["gene_id"]; ann = vars["gene_type"]; name = vars["gene_name"]; hgnc = vars["hgnc_id"]; hav = vars["havana_gene"]; print id,ann,name,hgnc,hav; }' >> gencode.v35.gene.annotation
#tail -n+2 gencode.v35.gene.annotation | awk -F"\t" '{if($4!="")print $4}' | wc -l ## 38596
#tail -n+2 gencode.v35.gene.annotation | awk -F"\t" '{if($4!="")print $4}' | sort | uniq | wc -l ## 38543

## Generate a map of approved Gencide gene IDs to dbXrefs
echo "gene_id transcript_id gene_name hgnc_id havana_gene" | tr ' ' '\t' > gencode.v35.trans.annotation
cat gencode.v35.annotation.gtf | awk -F"\t" '!/#/{if($3=="transcript")print $9}' | sed 's/; /;/g' | sed 's/\"//g' | awk -F";" 'BEGIN{FS=";";OFS="\t"}{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, " "); if(n) { x = substr($i, n + 1); vars[substr($i, 1, n - 1)] = substr($i, n + 1, length(x)) } } id = vars["gene_id"]; trans = vars["transcript_id"]; name = vars["gene_name"]; hgnc = vars["hgnc_id"]; hav = vars["havana_gene"]; print id,trans,name,hgnc,hav; }' >> gencode.v35.trans.annotation

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.metadata.EntrezGene.gz
gunzip gencode.v35.metadata.EntrezGene.gz
head -n1 gencode.v35.trans.annotation | awk 'BEGIN{FS=OFS="\t"}{print $0,"EntrezGene"}' > gencode.v35.trans2.annotation
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{ print $0, a[$2]}' gencode.v35.metadata.EntrezGene <(tail -n+2 gencode.v35.trans.annotation) >> gencode.v35.trans2.annotation

## Stats
cat gencode.v35.trans2.annotation | awk 'BEGIN{FS=OFS="\t"}{print $1,$3,$4,$6}' | uniq > gencode.v35.map
tail -n+2 gencode.v35.map | awk -F"\t" '{a+=1;if($3!="")b+=1;if($4!="")c+=1;}END{print "ensembl_id=",a,"hgnc_id=",b,"entrez_id=",c}'  ## ensembl_id= 60656 hgnc_id= 38596 entrez_id= 25600
tail -n+2 gencode.v35.map | awk -F"\t" '{if($3!="" && $4!="")a+=1;}END{print "gencode_id with hgnc and entrez  ids=",a}'  ## gencode_id with hgnc and entrez  ids= 24528
#tail -n+2 gencode.v35.map | awk -F"\t" '{a[$2]=1;if($3!="")b[$3]=1;if($4!="")c[$4]=1;}END{print "uniq_symb=",length(a),"uniq_hgnc=",length(b),"uniq_entrez=",length(c)}' ## uniq_symb= 59609 uniq_hgnc= 38543 uniq_entrez= 25532


#################################
#### 3. discrepancies between databases

## NCBI and HGNC
## merge HGNC ids and symbols into the NCBI DB to compare
head -n1 Homo_sapiens.gene_info.map | awk 'BEGIN{FS=OFS="\t"}{print $0,"hgnc-hgnc_id","hgnc-symbol"}' > NCBI.map.hgncExt
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$3]=$1 FS $2;next}{print $0,a[$1]}' <(tail -n+2 hgnc_approved.map) <(tail -n+2 Homo_sapiens.gene_info.map) >> NCBI.map.hgncExt

tail -n+2 NCBI.map.hgncExt | awk -F"\t" '{if($3!=$7)print}' ## Entrez ids have HGNC ids and symbols in NCBI different from those in HGNC database (7 records)
tail -n+2 NCBI.map.hgncExt | awk -F"\t" '{if($3!="" && $3==$7 && $2!=$8)print}' ## Entrez ids have the same HGNC ids in NCBI and HGNC database but with different symbols (37 mitochondrial genes)
comm -23 <(tail -n+2 hgnc_approved.map | awk -F"\t" '{if($3!="")print $1}' | sort)  <(tail -n+2 NCBI.map.hgncExt | awk -F"\t" '{if($7!="")print $7}' | sort) ## Entrez ids without HGNC ids in NCBI database but have HGNC ids in HGNC database (4 records)

## Gencode and HGNC
## A) In HGNC
tail -n+2 hgnc_approved.map | awk -F"\t" '{if($4!="")print $4}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)print $2}' | grep -Fwf - hgnc_approved.map ## there are 3 pairs of HGNC IDs where each pair maps to one gencode IDs
## B) In Gencode
tail -n+2 gencode.v35.map | grep -v "_PAR_Y" | awk -F"\t" '{if($3!="")print $3}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)print $2}' | grep -Fwf - gencode.v35.map | sort -k3,3 ## there are 17 pairs of gencode IDs where each pair maps to one  HGNC IDs. Each pair has the gene symbol as well. Similarly, each pair has the same Entrez ID but - unexpectedly - one member in 12 (out of the 17) pairs is missing the Entrez ID.
## merge HGNC ids and symbols into Gencode to compare
cat gencode.v35.map | awk 'BEGIN{FS=OFS="\t"}{split($1,a,".");print a[1],$2,$3,$4}' > gencode.v35.noGeneVer.map
head -n1 gencode.v35.noGeneVer.map | awk 'BEGIN{FS=OFS="\t"}{print $0,"hgnc-hgnc_id","hgnc-symbol","hgnc-hgnc_id2","hgnc-symbol2"}' > gencode.v35.map.hgncExt
awk 'BEGIN{FS=OFS="\t"}FNR==NR{if(a[$4]=="")a[$4]=$1 FS $2;else a[$4]=a[$4] FS $1 FS $2;next}{print $0,a[$1]}' <(tail -n+2 hgnc_approved.map) <(tail -n+2 gencode.v35.noGeneVer.map) >> gencode.v35.map.hgncExt
tail -n+2 gencode.v35.map.hgncExt | awk -F"\t" '{if(3!="" && $5!="" && $3!=$5 && $3!=$7)print}' ## Gencode ids have HGNC ids and symbols in Gencode different from those in HGNC database (6 records)
tail -n+2 gencode.v35.map.hgncExt | awk -F"\t" '{if($3!="" && (($3==$5 && $2!=$6) || ($3==$7 && $2!=$8)))print}' ## Gencode ids have the same HGNC ids in Gencode and HGNC database but with different symbols (48 genes)
comm -23 <(tail -n+2 hgnc_approved.map | awk -F"\t" '{if($4!="")print $1}' | sort)  <(tail -n+2 gencode.v35.map.hgncExt | awk -F"\t" 'BEGIN{OFS="\n";}{if($5!="")print $5;if($7!="")print $7}' | sort | uniq) | grep -Fwf - hgnc_approved.map ## Gencode ids without HGNC ids in Gencode database but have HGNC ids in HGNC database (84 records)

## Gencode and NCBI
## A) In NCBI: entrez_id with corresponding several Gencode IDs
grep "ENSG.*ENSG" Homo_sapiens.gene_info.map | wc -l ## 66
## B) In Gencode
tail -n+2 gencode.v35.map | grep -v "_PAR_Y" | awk -F"\t" '{if($4!="")print $4}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)print $2}' > dup.Enterz_ids
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=1;next}{if(a[$4]==1)print $1,$2,$4,$3}' dup.Enterz_ids gencode.v35.map | sort -k3,3nr  ## There are 45 sets (44 pairs and an additional set of 3) of Gencode IDs where each set maps to one Entrez IDs. In most sets, the genes have different symbols and HGNC IDs

######################################################################################
#### 4. Ambiguity of gene symbols across multiple species

# Download HGNC dataset from the "Custom downloads" page
# Select these columns: HGNC ID Approved symbol Status  Mouse genome database ID        Mouse genome database ID(supplied by MGI)       Rat genome database ID(supplied by RGD)
# Submit and save output webpage as custom_ortho.txt

## a) analysis for the column of Mouse genome database ID(supplied by MGI)
tail -n+2 custom_ortho.txt | awk 'BEGIN{FS=OFS="\t"}{if($3=="Approved" && ($4!="" || $5!=""))print $1,$4,$5}' > HGNC-to-Mouse.hgnc-custom
cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($3!="")print $3}' | wc -l ## 18206
cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($3!="")print $3}' | sort | uniq | wc -l ## 17976
cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($3!="")print $3}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)a+=$1}END{print a}' ## 413 
cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($3!="")print $3}' | sort | uniq -c | sort -k1,1nr | awk '{if($1>1)print $0}' | wc -l ## 183
cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($3!="")print $3}' | grep "," | wc -l ## 416
cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($3!="")print $1,$3}' | awk -F"," '{print NF}' | sort -nr | head  ## 15
cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($3!="")print $1,$3}' | awk -F"," '{if(NF==15) print $0}' 

## b) analysis for the column of curated Mouse genome database ID
cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($2!="")print $2}' | wc -l ## 17692
cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($2!="")print $2}' | sort | uniq | wc -l ## 17561
cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($2!="" && $2==$3)print $0}' | wc -l ## 17515 
cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($2!="" && $3=="")print $0}' | wc -l ## 3
cat HGNC-to-Mouse.hgnc-custom | awk 'BEGIN{FS=OFS="\t"}{if($2!="" && $3!="" && $2!=$3)print $0}' | wc -l ## 174 


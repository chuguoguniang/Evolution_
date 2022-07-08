# Dynamic Actin Gene Family Evolution in Primates


## 一、鉴定肌动蛋白基因 
1. Downloaded protein sequences which were limited to genes with actin domain (Pfam:PF00022) from Biomart
```
#计数其中人氨基酸序列
cat human.fa | grep '>' | wc -l
#89


#对氨基酸序列计数
faops size human.fa
#序列长度差异较大，不利于MEGA比对

#去除小于160aa的氨基酸序列后计数（仍然是31个基因）
faops size human.fa | tsv-filter --ge 2:160  |wc -l
#74去除了15个蛋白序列

#大于等于160aa的氨基酸序列输进一个文件
faops size human.fa | tsv-filter --ge 2:160 >human.160.fa
```

2. the amino acid sequences of all known actin genes were adopted as queries in local BLASTP 
```
#安装blast
brew install blast

#对人总氨基酸序列建立索引
makeblastdb -in ./sequence/Homo_sapiens.GRCh38.pep.all.fa -dbtype prot -parse_seqids -out ./aindex

#BLASTP
blastp -query ./human_ge160.fa -db ./index -evalue 1e-30 -qcov_hsp_perc 80 -outfmt 6 -num_t
hreads 4 -out result1.tsv

#查看一轮BLASTP后匹配上的氨基酸序列数及基因数：
cat result1.tsv | cut -f 2 | sort | uniq | wc -l  #80
cat result1.tsv | cut -f 2 | sort | uniq > protein_ID1.lst
cat ./Homo_sapiens.GRCh38.pep.all.fa | grep ">" | grep -f protein_ID1.lst | cut -d " " -f 4 | sort | uniq | wc -l  #31

#第二轮blastp
faops some ./Homo_sapiens.GRCh38.pep.all.fa ./protein_ID1.lst blastp2query.fa
blastp -query ./blastp2query.fa -db ./index -evalue 1e-30 -qcov_hsp_perc 80 -outfmt 6 -num_threads 4 -out result2.tsv
cat result2.tsv | cut -f 2 | sort | uniq | wc -l #80不再增加
cat result2.tsv | cut -f 2 | sort | uniq > protein_ID2.lst  
#protein_ID2.lst 与protein_ID1.lst同，所以CDD输入序列可以用blastp2query.fa

```
3. Based on the BLASTP results, all genes were verified with the conserved actin domain by researching in corresponding Conserved Domain Database (CDD) online
```
在Batch CD-Search中上传blastp2query.fa


```

以上是人类中肌动蛋白基因鉴定的流程，其余6个物种方法同上。









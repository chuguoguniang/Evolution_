# Dynamic Actin Gene Family Evolution in Primates


## 一、鉴定人肌动蛋白基因 
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
#74去除了15个氨基酸序列

#大于等于160aa的氨基酸序列输进一个文件
faops size human.fa | tsv-filter --ge 2:160 | cut -f 1 >human.160.fa
faops some human.fa human.160.fa human_ge160.fa
```

2. the amino acid sequences of all known actin genes were adopted as queries in local BLASTP 
```
#安装blast
brew install blast

#对人总氨基酸序列建立索引
makeblastdb -in ./sequence/Homo_sapiens.GRCh38.pep.all.fa -dbtype prot -parse_seqids -out ./index

#BLASTP
blastp -query ./human_ge160.fa -db ./index -evalue 1e-30 -qcov_hsp_perc 80 -outfmt 6 -num_threads 4 -out result1.tsv

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
在Batch CD-Search中上传blastp2query.fa 数据库选择NCBI_Curated --17937 PSSMs 其它参数默认
下载consize结果（比对效果佳的留下），对结果逐行查看，与actin均有关系
提取所需信息：
cat human.txt | cut -f 1 >CDD.result.tsv
cat CDD.result.tsv | cut -d ">" -f 2  >CDD.result.1.tsv #以>号分割，取>号后的也算是第二列
查找对应基因数：
cat ./Homo_sapiens.GRCh38.pep.all.fa | grep ">" | grep -f CDD.result.1.tsv | cut -d " " -f 4 | sort | uniq | wc -l #30 与文章数量一样



```

以上是人类中肌动蛋白基因鉴定的流程，其余6个物种方法同上。

## 二、鉴定黑猩猩肌动蛋白基因 
```
cat chimpanzee.fa | grep '>' | sort | uniq | wc -l
#71条氨基酸序列

faops size chimpanzee.fa | tsv-filter --ge 2:160  |wc -l
#60 去除了11个氨基酸序列

faops size chimpanzee.fa | tsv-filter --ge 2:160 | cut -f 1 >chimpanzee.lst
faops some chimpanzee.fa chimpanzee.lst chimpanzee_ge160.fa

makeblastdb -in ./sequence/Pan_troglodytes.Pan_tro_3.0.pep.all.fa -dbtype prot -parse_seqids -out ./index
blastp -query ./chimpanzee_ge160.fa -db ./index -evalue 1e-30 -qcov_hsp_perc 80 -outfmt 6 -num_threads 4 -out result1.tsv

cat result1.tsv | cut -f 2 | sort | uniq | wc -l 
#60 第一轮blastp后就不再增加
cat result1.tsv | cut -f 2 | sort | uniq > protein_ID1.lst 
#protein_ID1.lst与chimpanzee_ge160.fa一样，CDD输入序列可以用chimpanzee_ge160.fa，似乎>后含有基因名CDD识别不了，所以去掉
faops some Pan_troglodytes.Pan_tro_3.0.pep.all.fa protein_ID1.lst chimpanzee_cdd.tsv
在CDD中上传chimpanzee_cdd.tsv
下载输出结果文件：chimpanzee.cdd.result.txt
cat chimpanzee.cdd.result.txt | cut -f 1 | cut -d ">" -f 2 > chimpanzee.cdd.result.1.tsv
cat ./Pan_troglodytes.Pan_tro_3.0.pep.all.fa | grep ">" | grep -f chimpanzee.cdd.result.1.tsv | cut -d " " -f 4 | sort | uniq | wc -l #33个基因
```

## 三、鉴定大猩猩肌动蛋白基因
```
cat gorilla.fa | grep '>' | sort | uniq | wc -l
#64条氨基酸序列  Biomart下载时count 32个基因
faops size gorilla.fa | tsv-filter --ge 2:160  | wc -l 
#52条氨基酸序列，去除12条

faops size gorilla.fa | tsv-filter --ge 2:160 | cut -f 1 >gorilla.lst
faops some gorilla.fa gorilla.lst gorilla_ge160.fa

makeblastdb -in ./Gorilla_gorilla.gorGor4.pep.all.fa -dbtype prot -parse_seqids -out ./index
blastp -query ./gorilla_ge160.fa -db ./index -evalue 1e-30 -qcov_hsp_perc 80 -outfmt 6 -num_threads 4 -out result1.tsv
cat result1.tsv | cut -f 2 | sort | uniq | wc -l 
#53 第一轮blastp增加1个，所以第二轮blastp
cat result1.tsv | cut -f 2 | sort | uniq > protein_ID1.lst
cat ./Gorilla_gorilla.gorGor4.pep.all.fa | grep ">" | grep -f protein_ID1.lst | cut -d " " -f 4 | sort | uniq | wc -l  #31个基因
faops some ./Gorilla_gorilla.gorGor4.pep.all.fa ./protein_ID1.lst blastp2query.fa
blastp -query ./blastp2query.fa -db ./index -evalue 1e-30 -qcov_hsp_perc 80 -outfmt 6 -num_threads 4 -out result2.tsv
cat result2.tsv | cut -f 2 | sort | uniq | wc -l #53第二轮blastp不再增加
cat result2.tsv | cut -f 2 | sort | uniq > protein_ID2.lst  
#protein_ID2.lst 与protein_ID1.lst同，所以CDD输入序列可以用blastp2query.fa

CDD:
输入序列文件：blastp2query.fa
下载输出结果文件：gorilla.cdd.result.txt
cat gorilla.cdd.result.txt | cut -f 1 | cut -d ">" -f 2 > gorilla.cdd.result.1.tsv
cat ./Gorilla_gorilla.gorGor4.pep.all.fa | grep ">" | grep -f gorilla.cdd.result.1.tsv | cut -d " " -f 4 | sort | uniq | wc -l  #30个基因
```

## 四、鉴定红毛猩猩肌动蛋白基因
```
cat sumatran_orangutan.fa | grep ">" | wc -l 
#42条氨基酸序列 Biomart下载时count 29个基因

faops size sumatran_orangutan.fa | tsv-filter --ge 2:160  | wc -l
#40

faops size sumatran_orangutan.fa | tsv-filter --ge 2:160 | cut -f 1 >sumatran_orangutan.lst
faops some sumatran_orangutan.fa sumatran_orangutan.lst sumatran_orangutan_ge160.fa

makeblastdb -in ./Pongo_abelii.Susie_PABv2.pep.all.fa -dbtype prot -parse_seqids -out ./index
blastp -query ./sumatran_orangutan_ge160.fa -db ./index -evalue 1e-30 -qcov_hsp_perc 80 -outfmt 6 -num_threads 4 -out result1.tsv
cat result1.tsv | cut -f 2 | sort | uniq | wc -l 
#40 第一轮blastp后不增加
cat result1.tsv | cut -f 2 | sort | uniq > protein_ID1.lst 
faops some Pongo_abelii.Susie_PABv2.pep.all.fa protein_ID1.lst sumatran_orangutan_cdd.tsv

在CDD中上传sumatran_orangutan_cdd.tsv
下载输出文件：sumatran_orangutan.cdd.result.txt 逐行查看，和actin无关的看是否有和其相同的序列含有actin ，比如Q#6 - >ENSPPYP00000033681.1可能在两行出现，一行里含actin相关的，另一行不含
cat sumatran_orangutan.cdd.result.txt | cut -f 1 | cut -d ">" -f 2 > sumatran_orangutan.cdd.result.1.tsv
cat ./Pongo_abelii.Susie_PABv2.pep.all.fa | grep ">" | grep -f sumatran_orangutan.cdd.result.1.tsv | cut -d " " -f 4 | sort | uniq | wc -l  #27个基因
```

## 五、鉴定长臂猿肌动蛋白基因
```
cat gibbon.fa | grep ">" | wc -l 
#48条氨基酸序列 26个基因
faops size gibbon.fa | tsv-filter --ge 2:160  | wc -l 
#43
faops size gibbon.fa | tsv-filter --ge 2:160 | cut -f 1 >gibbon.lst
faops some gibbon.fa gibbon.lst gibbon_ge160.fa

makeblastdb -in ./Nomascus_leucogenys.Nleu_3.0.pep.all.fa -dbtype prot -parse_seqids -out ./index
blastp -query ./gibbon_ge160.fa -db ./index -evalue 1e-30 -qcov_hsp_perc 80 -outfmt 6 -num_threads 4 -out result1.tsv
cat result1.tsv | cut -f 2 | sort | uniq | wc -l 
#43 第一轮blastp后不再增加
cat result1.tsv | cut -f 2 | sort | uniq > protein_ID1.lst 
faops some Nomascus_leucogenys.Nleu_3.0.pep.all.fa protein_ID1.lst gibbon_cdd.tsv

CDD:
输入序列文件：gibbon_cdd.tsv
下载输出文件：gibbon.cdd.result.txt
cat gibbon.cdd.result.txt | cut -f 1 | cut -d ">" -f 2 > gibbon.cdd.result.1.tsv
cat ./Nomascus_leucogenys.Nleu_3.0.pep.all.fa | grep ">" | grep -f gibbon.cdd.result.1.tsv | cut -d " " -f 4 | sort | uniq | wc -l
#26个基因
```

## 六、鉴定猕猴肌动蛋白基因
```
cat macaque.fa | grep ">" | wc -l 
#55条氨基酸序列  30个基因 （看Gene stable ID version得30 而不是Gene name ，因为有些没有Gene name）
faops size macaque.fa | tsv-filter --ge 2:160  | wc -l 
#54
faops size macaque.fa | tsv-filter --ge 2:160 | cut -f 1 >macaque.lst
faops some macaque.fa macaque.lst macaque_ge160.fa

makeblastdb -in ./Macaca_mulatta.Mmul_10.pep.all.fa -dbtype prot -parse_seqids -out ./index
blastp -query ./macaque_ge160.fa -db ./index -evalue 1e-30 -qcov_hsp_perc 80 -outfmt 6 -num_threads 4 -out result1.tsv
cat result1.tsv | cut -f 2 | sort | uniq | wc -l 
#54第一轮blastp后不再增加
cat result1.tsv | cut -f 2 | sort | uniq > protein_ID1.lst 
faops some Macaca_mulatta.Mmul_10.pep.all.fa protein_ID1.lst macaque_cdd.tsv

CDD:
输入序列文件：macaque_cdd.tsv
下载输出文件：macaque.cdd.result.txt 逐行核对，删除表头上面的几行，7个物种都是如此做的
cat macaque.cdd.result.txt | cut -f 1 | cut -d ">" -f 2 > macaque.cdd.result.1.tsv
cat ./Macaca_mulatta.Mmul_10.pep.all.fa | grep ">" | grep -f macaque.cdd.result.1.tsv | cut -d " " -f 4 | sort | uniq | wc -l
#30个基因
```

## 七、鉴定白簇耳狨猴肌动蛋白基因
```
cat  marmoset.fa | grep ">" | wc -l 
#62条氨基酸序列 38个基因
faops size marmoset.fa | tsv-filter --ge 2:160  | wc -l 
#59
faops size marmoset.fa | tsv-filter --ge 2:160 | cut -f 1 >marmoset.lst
faops some marmoset.fa marmoset.lst marmoset_ge160.fa

makeblastdb -in ./Callithrix_jacchus.mCalJac1.pat.X.pep.all.fa -dbtype prot -parse_seqids -out ./index
blastp -query ./marmoset_ge160.fa -db ./index -evalue 1e-30 -qcov_hsp_perc 80 -outfmt 6 -num_threads 4 -out result1.tsv
cat result1.tsv | cut -f 2 | sort | uniq | wc -l 
#59第一轮blastp后不再增加
cat result1.tsv | cut -f 2 | sort | uniq > protein_ID1.lst
faops some Callithrix_jacchus.mCalJac1.pat.X.pep.all.fa protein_ID1.lst marmoset_cdd.tsv

CDD:
输入序列文件：marmoset_cdd.tsv
下载输出文件：marmoset.cdd.result.txt
cat marmoset.cdd.result.txt | cut -f 1 | cut -d ">" -f 2 > marmoset.cdd.result.1.tsv
cat ./Callithrix_jacchus.mCalJac1.pat.X.pep.all.fa | grep ">" | grep -f marmoset.cdd.result.1.tsv | cut -d " " -f 4 | sort | uniq | wc -l
#34个基因
```
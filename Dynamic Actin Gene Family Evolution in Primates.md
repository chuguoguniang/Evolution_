# Dynamic Actin Gene Family Evolution in Primates


## 一、鉴定肌动蛋白基因
### 1、鉴定人肌动蛋白基因 
+ Downloaded protein sequences which were limited to genes with actin domain (Pfam:PF00022) from Biomart
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

+ the amino acid sequences of all known actin genes were adopted as queries in local BLASTP 
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
+ Based on the BLASTP results, all genes were verified with the conserved actin domain by researching in corresponding Conserved Domain Database (CDD) online
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

### 2、鉴定黑猩猩肌动蛋白基因 
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

### 3、鉴定大猩猩肌动蛋白基因
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

### 4、鉴定红毛猩猩肌动蛋白基因
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

### 5、鉴定长臂猿肌动蛋白基因
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

### 6、鉴定猕猴肌动蛋白基因
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

### 7、鉴定白簇耳狨猴肌动蛋白基因
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
### 小结
+ 共鉴定7个物种210个actin gene

+ 整理7个物种actin基因的氨基酸序列:
```
1. human:
mv CDD.result.1.tsv CDD.result.1.lst
faops some Homo_sapiens.GRCh38.pep.all.fa CDD.result.1.lst human.actin.aa.fa
2. chimpanzee:
mv chimpanzee.cdd.result.1.tsv chimpanzee.cdd.result.1.lst
faops some Pan_troglodytes.Pan_tro_3.0.pep.all.fa chimpanzee.cdd.result.1.lst chimpanzee.actin.aa.fa
3. gibbon:
mv gibbon.cdd.result.1.tsv gibbon.cdd.result.1.lst
faops some Nomascus_leucogenys.Nleu_3.0.pep.all.fa gibbon.cdd.result.1.lst gibbon.actin.aa.fa
4. gorilla:
mv gorilla.cdd.result.1.tsv gorilla.cdd.result.1.lst
faops some Gorilla_gorilla.gorGor4.pep.all.fa gorilla.cdd.result.1.lst gorilla.actin.aa.fa
5. macaque:
mv macaque.cdd.result.1.tsv macaque.cdd.result.1.lst
faops some Macaca_mulatta.Mmul_10.pep.all.fa macaque.cdd.result.1.lst macaque.actin.aa.fa
6. marmoset:
mv marmoset.cdd.result.1.tsv marmoset.cdd.result.1.lst
faops some Callithrix_jacchus.mCalJac1.pat.X.pep.all.fa marmoset.cdd.result.1.lst marmoset.actin.aa.fa
7. sumatran_orangutan
mv sumatran_orangutan.cdd.result.1.tsv sumatran_orangutan.cdd.result.1.lst
faops some Pongo_abelii.Susie_PABv2.pep.all.fa sumatran_orangutan.cdd.result.1.lst sumatran_orangutan.actin.aa.fa

```
+ 合并7个物种actin基因的氨基酸序列:(进入指定文件夹/mnt/d/0~GitHub/Evolution_/MEGA/aa)
```
JOB=$(ls)
for i in $JOB; do cat $i >> aa.all.fa; done

cat aa.all.fa | grep ">" | wc -l 
#372
```

## 二、序列比对和系统发育分析
### 1、氨基酸序列比对和建树
+ 氨基酸序列比对，选ClustalW（大约8分钟）后，保存。#可以比对了，或许正是前面将小于160aa的序列去除的缘故。
+ 建树，但是找不到文章中的Jukes-Cantor模型，所以建树采用JTT 模型，NJ bootstrap 1000replicates 。其余参数默认。约六个小时出现错误，显示bootstrap只进行998replicates. 所以我调整，将bootstrap由1000replicates设到997，还是出现类似错误，似乎有3个replicates没进行，然后也没见结果。bootstrap调成500 replicates出现和白天一样的提示，并没有到设定的replicates。但是有tree文件。不过文章是用核苷酸序列建树的，所以氨基酸序列还需要建吗？


### 2、CDS序列比对和建树
+ 得到CDS序列尝试1——失败
```
生成Biomart上传所需的文件(含Gene stable ID(s) with version)
cat ./Pan_troglodytes.Pan_tro_3.0.pep.all.fa | grep ">" | grep -f chimpanzee.cdd.result.1.lst | cut -d " " -f 4 | sort | uniq | cut -d ":" -f 2 > chimpanzee.cds.lst
cat ./Nomascus_leucogenys.Nleu_3.0.pep.all.fa | grep ">" | grep -f gibbon.cdd.result.1.lst | cut -d " " -f 4 | sort | uniq | cut -d ":" -f 2 > gibbon.cds.lst
cat ./Gorilla_gorilla.gorGor4.pep.all.fa | grep ">" | grep -f gorilla.cdd.result.1.lst | cut -d " " -f 4 | sort | uniq | cut -d ":" -f 2 > gorilla.cds.lst
cat ./Homo_sapiens.GRCh38.pep.all.fa | grep ">" | grep -f CDD.result.1.lst | cut -d " " -f 4 | sort | uniq | cut -d ":" -f 2 > human.cds.lst
cat ./Macaca_mulatta.Mmul_10.pep.all.fa | grep ">" | grep -f macaque.cdd.result.1.lst | cut -d " " -f 4 | sort | uniq | cut -d ":" -f 2 > macaque.cds.lst
cat ./Callithrix_jacchus.mCalJac1.pat.X.pep.all.fa | grep ">" | grep -f marmoset.cdd.result.1.lst | cut -d " " -f 4 | sort | uniq | cut -d ":" -f 2 > marmoset.cds.lst
cat ./Pongo_abelii.Susie_PABv2.pep.all.fa | grep ">" | grep -f sumatran_orangutan.cdd.result.1.lst | cut -d " " -f 4 | sort | uniq | cut -d ":" -f 2 > sumatran_orangutan.cds.lst

Biomart下载好后，将7个物种的cds合并为一个文件用于MEGA比对
JOB=$(ls)
for i in $JOB; do cat $i >> cds.all.fa; done
cat cds.all.fa | grep ">" | wc -l 
#514 由基因下CDS，对应的蛋白是有小于160aa的，数量和aa数量372对不上，所以换成在Biomart上传含有Protein stable ID(s) with version的文件去下载CDS
```
+ 得到CDS序列尝试2——成功
```
cat human.actin.aa.fa | grep ">" | cut -d ">" -f 2 > human.cds.lst
cat chimpanzee.actin.aa.fa | grep ">" | cut -d ">" -f 2 > chimpanzee.cds.lst
cat gorilla.actin.aa.fa | grep ">" | cut -d ">" -f 2 > gorilla.cds.lst
cat marmoset.actin.aa.fa | grep ">" | cut -d ">" -f 2 > marmoset.cds.lst
cat gibbon.actin.aa.fa | grep ">" | cut -d ">" -f 2 > gibbon.cds.lst
cat macaque.actin.aa.fa | grep ">" | cut -d ">" -f 2 > macaque.cds.lst
cat sumatran_orangutan.actin.aa.fa | grep ">" | cut -d ">" -f 2 > sumatran_orangutan.cds.lst

Biomart下载好后，将7个物种的cds合并为一个文件用于MEGA比对
JOB=$(ls) #文件夹下只有这七个文件
for i in $JOB; do cat $i >> cds.all.fa; done
cat cds.all.fa | grep ">" | wc -l 
#372 此时和aa序列数是一致的，这才方便aa序列指导cds序列比对吧，或者说前面删除的小于160的aa才有意义吧
```
+ CDS序列比对
```
选ClustalW。10分钟还未比对完。后面发现出现序列差异大的提示，所以把小于1000核苷酸的序列去除后再比
faops size cds.all.fa |tsv-filter --ge 2:1000 | cut -f 1 > cds.ge1000.lst
faops some cds.all.fa cds.ge1000.lst cds.all.ge1000.fa
再次比对CDS（7个物种在一起的CDS序列），比对选ClustalW. 选align DNA 而不是align codons，因为选align codons显示终止密码子出现在翻译区
把小于1000核苷酸的cds序列去除后还是出现序列差异大。
将ClustalW换成muscle比对核苷酸序列，没有报错，但是当把比对后的核苷酸序列用于建树时，发现有名称重复，因为CDS下载时选的Gene stable ID version和基因，7个物种CDS序列合并后计数是372，但不重复的才210

所以重新下载CDS序列，名称改为基因和Protein stable ID，按文章的意思，没有基因的以Protein stable ID表示序列
Biomart下载好后，将7个物种的cds合并为一个文件用于MEGA比对
JOB=$(ls) #文件夹下只有这七个文件
for i in $JOB; do cat $i >> cds.all.fa; done
cat cds.all.fa | grep ">" | wc -l  #372
cat cds.all.fa | grep ">" | sort | uniq | wc -l  #372不重复了，这次应该可以正常比对和建树，按依据文章，改名字：由基因名的用基因名，没有的用Protein stable ID，后面加每个物种的后缀。
```
cat human.cds.fa | grep ">" >human.replace.csv
在csv中手动操作
cat human.replace.csv | tr "," "\t" > human.replace.tsv
faops replace -s human.cds.fa human.replace.tsv human.cds.1.fa (在此human.cds.fa 和 human.replace.tsv是对应的，加不加-s都一样)

cat chimpanzee.cds.fa | grep ">" >chimpanzee.replace.csv
在csv中手动操作（以|分两列，第一列是空的即没基因名的和第二列合并）
cat chimpanzee.replace.csv | tr "," "\t" > chimpanzee.replace.tsv
faops replace -s chimpanzee.cds.fa chimpanzee.replace.tsv chimpanzee.cds.1.fa

cat gorilla.cds.fa | grep ">" >gorilla.replace.csv
在csv中手动操作（以|分两列，第一列是空的即没基因名的和第二列合并）
cat gorilla.replace.csv | tr "," "\t" > gorilla.replace.tsv
faops replace -s gorilla.cds.fa gorilla.replace.tsv gorilla.cds.1.fa

cat sumatran_orangutan.cds.fa | grep ">" >sumatran_orangutan.replace.csv
在csv中手动操作（以|分两列，第一列是空的即没基因名的和第二列合并）
cat sumatran_orangutan.replace.csv | tr "," "\t" > sumatran_orangutan.replace.tsv
faops replace -s sumatran_orangutan.cds.fa sumatran_orangutan.replace.tsv sumatran_orangutan.cds.1.fa

cat gibbon.cds.fa | grep ">" >gibbon.replace.csv
在csv中手动操作（以|分两列，第一列是空的即没基因名的和第二列合并）
cat gibbon.replace.csv | tr "," "\t" > gibbon.replace.tsv
faops replace -s gibbon.cds.fa gibbon.replace.tsv gibbon.cds.1.fa

cat macaque.cds.fa | grep ">" >macaque.replace.csv
在csv中手动操作（以|分两列，第一列是空的即没基因名的和第二列合并）
cat macaque.replace.csv | tr "," "\t" > macaque.replace.tsv
faops replace -s macaque.cds.fa macaque.replace.tsv macaque.cds.1.fa

cat marmoset.cds.fa | grep ">" >marmoset.replace.csv
在csv中手动操作（以|分两列，第一列是空的即没基因名的和第二列合并）
cat marmoset.replace.csv | tr "," "\t" > marmoset.replace.tsv
faops replace -s marmoset.cds.fa marmoset.replace.tsv marmoset.cds.1.fa
cat marmoset.cds.1.fa | head -n 20 #看一下替换序列名后的效果

每个物种更改序列名后合并
JOB=$(ls) #进入一个新的文件夹，文件夹下只有这七个文件
for i in $JOB; do cat $i >> cds.all.changename.fa; done
cat cds.all.changename.fa | grep ">" | wc -l  #372
cat cds.all.changename.fa | grep ">" | sort | uniq | wc -l  #240
可见以基因名加物种有重复名字出现，所以我一个基因只保留一个最长的氨基酸序列,但由于有些序列没有基因名，所以可能出现15条不重复的序列名字其实对应10个基因。所以用曾师兄帮我写的程序后，达到不再出现重复序列名，但存在有些基因几个序列，这不妨碍比对。
打开powershell ，输入：python 程序 待处理文件（最后一行下面加>) 生成文件   #注意，程序 待处理文件 生成文件 三者都是带路径的，并且运行前为了判断长度需在待处理文件最后一行下面加>，运行结束将待处理文件和生成文件中最下面一行的>去掉。
具体地：
python C:\Users\15560\Desktop\jbs_get_the_longest_sq.py D:\0~GitHub\Evolution_\MEGA\cds\outputv2\cds.all.fa.change\human.cds.1.fa D:\0~GitHub\Evolution_\MEGA\cds\outputv2\cds.all.fa.change\human.cds.1.longest.fa
python C:\Users\15560\Desktop\jbs_get_the_longest_sq.py D:\0~GitHub\Evolution_\MEGA\cds\outputv2\cds.all.fa.change\chimpanzee.cds.1.fa D:\0~GitHub\Evolution_\MEGA\cds\outputv2\cds.all.fa.change\chimpanzee.cds.1.longest.fa
python C:\Users\15560\Desktop\jbs_get_the_longest_sq.py D:\0~GitHub\Evolution_\MEGA\cds\outputv2\cds.all.fa.change\gorilla.cds.1.fa D:\0~GitHub\Evolution_\MEGA\cds\outputv2\cds.all.fa.change\gorilla.cds.1.longest.fa
python C:\Users\15560\Desktop\jbs_get_the_longest_sq.py D:\0~GitHub\Evolution_\MEGA\cds\outputv2\cds.all.fa.change\marmoset.cds.1.fa D:\0~GitHub\Evolution_\MEGA\cds\outputv2\cds.all.fa.change\marmoset.cds.1.longest.fa
python C:\Users\15560\Desktop\jbs_get_the_longest_sq.py D:\0~GitHub\Evolution_\MEGA\cds\outputv2\cds.all.fa.change\gibbon.cds.1.fa D:\0~GitHub\Evolution_\MEGA\cds\outputv2\cds.all.fa.change\gibbon.cds.1.longest.fa
python C:\Users\15560\Desktop\jbs_get_the_longest_sq.py D:\0~GitHub\Evolution_\MEGA\cds\outputv2\cds.all.fa.change\macaque.cds.1.fa D:\0~GitHub\Evolution_\MEGA\cds\outputv2\cds.all.fa.change\macaque.cds.1.longest.fa
python C:\Users\15560\Desktop\jbs_get_the_longest_sq.py D:\0~GitHub\Evolution_\MEGA\cds\outputv2\cds.all.fa.change\sumatran_orangutan.cds.1.fa D:\0~GitHub\Evolution_\MEGA\cds\outputv2\cds.all.fa.change\sumatran_orangutan.cds.1.longest.fa

合并为一个
JOB=$(ls) #文件夹下只有这七个文件
for i in $JOB; do cat $i >> cds.all.cn.uniq.fa; done
cat cds.all.cn.uniq.fa | grep ">" | wc -l #240
cat cds.all.cn.uniq.fa | grep ">" | sort | uniq | wc -l  #240

开始对CDS多序列比对和建树
使用muscle进行多序列比对，然后neighbor-joining method with a Jukes-Cantor model   bootstrap analysis with 1000 replicates  大约5分钟建完，生成Original Tree 和 Bootstrap consensus Tree

Ka Ks计算方法：MEGA打开比对好的核苷酸序列，点击DISTANCE, 选择第一个compute pairwise distance , 替换类型同义-非同义，替换包括仅仅同义（同义和非同义分开来），其它参数默认，但记下来。导出到EXCEL，选择老的2007+的版本。另存到一个文件夹，改名ds。dn同理。然后复制dn excel。全选复制ds内容，选择性粘贴到dn，选除。去错，另存为一个文件。

建树后如何分类分组计算PI？
```
## 三、识别假基因
假基因（Pseudogenes，Pseudo-意为“假”）是一类染色体上的基因片段。假基因的序列通常与对应的基因相似，但至少是丧失了一部分功能，如基因不能表达或编码的蛋白质没有功能。 #维基百科
区分假基因和活基因有助于我们了解肌动蛋白基因家族的进化史
编码一个基因的DNA序列在基因组内完整出现一次，称为该基因的一个拷贝。
```
all of the nucleotide sequences of actin domains from seven species were employed to search in all the genomes used in this work (BLASTN).
blastn：是将给定的核酸序列与核酸数据库中的序列进行比对

1. 下载各物种基因组序列
方法：进入Ensemble 点击FTP Download ，选择物种，DNA(FASTA),选toplevel下载完整基因组序列。

2. 建库
makeblastdb -in ./Homo_sapiens.GRCh38.dna.toplevel.fa -dbtype nucl -parse_seqids -out ./index
blastn -query ./human.cds.fa -db ./index -perc_identity 35 -outfmt 6 -num_threads 4 -out out_file



chimpanzee:
gzip -d Pan_troglodytes.Pan_tro_3.0.dna.toplevel.fa.gz
makeblastdb -in ./Pan_troglodytes.Pan_tro_3.0.dna.toplevel.fa -dbtype nucl -parse_seqids -out ./index
blastn -query ./chimpanzee.cds.fa -db ./index -evalue 1e-30 -qcov_hsp_perc 80 -perc_identity 35 -outfmt 6 -num_threads 4 -out out_file
```
## 四、密码子使用
编码相同氨基酸的不同密码子被称为同义密码子，几乎在所有物种中不均等使用，呈现出基因的进化模式。
RSCU相对同义密码子使用度。以某一个同义密码子的使用次数为分子，以该密码子预期出现的次数为分母。预测出现的次数为该密码子所编码的氨基酸的所有密码子平均使用的次数。 如果密码子使用没有偏好，则该密码子的RSCU值等于1，大于1表明其使用频率相对较高。

## 五、 肌动蛋白基因表达模式
计算肌动蛋白基因表达值的变异系数（CV；SD/mean）以估计表达模式。


讨论：
一般来说，基因获得新的功能是由于自我复杂度的增加或拷贝数的变化。增加基因的长度或与其他领域的融合可以增加其复杂性，而复制则提供了一个获得新功能而不失去原有功能的机会。
对于多基因家族的分化与进化，通常采用“协同进化”和“生灭”模型来解释 。
由于肌动蛋白家族在细胞活动的各个方面都起着至关重要的作用，其相关功能不容易改变或移除。然而，肌动蛋白基因拷贝数的变化方式可能提供了另一种进化途径，以满足相互冲突的需求，即肌动蛋白被保守以维持重要功能，并在体内进化出新的功能，以帮助适应环境压力。在这种情况下，生物体可能不仅保持身体的正常工作，而且使物种从简单到复杂，从粗略到精细。我们推测，出生和死亡进化模式可能是其他高度保守的多基因家族的共同进化机制。
Actin genes interspecifically clustered that belong to the orthologous groups were highly conserved because of fundamental importance. On the contrary, complex groups contained
actin gene members that displayed copy number variation with significantly higher levels of average nucleotide divergence and Ka/Ks ratios compared to the orthologous groups.
对密码子偏好性和基因表达水平的分析表明，灵长类动物的肌动蛋白基因在物种内存在极大的差异，但在不同物种的群体内高度保守。这些结果可能可以用肌动蛋白基因家族的生灭进化过程来解释，这可能是其他高度保守的多基因家族的一般进化机制。


总结：本篇文章拿到手大概两周之久，已了解文章中的技术和用途，具体细节有待钻研。2022.07.12开始阅读一些IV型胶原蛋白演化的文章，尤其关注Gxy重复，制定方案，预计本周末开始做IV型胶原蛋白演化分析。



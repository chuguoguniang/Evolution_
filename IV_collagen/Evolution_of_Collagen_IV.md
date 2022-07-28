>工作计划时间线：  
2022.07.22-2022.08.10将如下有关胶原蛋白演化的工作完成;  
2022.08.10-2022.08.20从上述工作结果中探讨其它方面。

# 一、灵长目中IV型胶原基因的演化（参考肌动蛋白演化的标准演化方法）
## 1、确定研究对象及基因
+ 研究对象：
```
7个灵长目物种：
human (Homo sapiens)  
chimpanzee (Pan troglodytes)  #黑猩猩
gorilla (Gorilla gorilla)  #大猩猩
Sumatran orangutan (Pongo abelii)  #苏门答腊猩猩
gibbon (Nomascus leucogenys) #长臂猿 
crab-eating macaque (Macaca fascicularis) #长尾猕猴
White-tufted-ear marmoset (Callithrix jacchus)  #普通狨猴
```
+ 确定IV型胶原基因
```
在ensemble官网进入相应物种界面，搜索COL4,点击基因，去掉ncRNA、假基因、结合蛋白之类后，查看得：
7个物种均有6个IV型胶原基因,:COL4A1 COL4A2 COL4A3 COL4A4 COL4A5 COL4A6
```
## 2、多序列比对和建树  
（此步我输入MEGA的文件是CRLF，多序列比对完保存的是.meg格式进行的后续分析）
+ 下载、整理7个物种IV型胶原蛋白基因的氨基酸序列和CDS序列 
```
下载：进入ensemble官网下选相应物种，每个基因暂取最长的转录本对应的蛋白序列和CDS序列，手动下载。（其实我想用命令行通过FTP下载再保留最长转录本再提取，但一时没想好怎么用命令行完成，暂且手动下载） 
注：下载过程中，搜索Gorilla基因时Only searching Western Lowland Gorilla; 搜索Gibbon基因时Only searching Northern white-cheeked gibbon  
添加后缀：为每个物种中的序列名添加后缀，规则是学名的三个字母，如人用Hsa作为后缀;cds序列在后缀前加_cds用以与氨基酸序列区分
合并序列：以合并7个物种氨基酸序列为例：  
JOB=$(ls)
 for i in $JOB ; do cat $i >> protein.all.fa ; done  
 cat protein.all.fa | grep ">" | wc -l #42  
```
+ 氨基酸与cds序列的多序列比对及建树
```
使用muscle多序列比对，保存为mega格式 参数默认  
NJ建树 参数默认 建出的树COL4A5确实与COL4A1挺近 树的PDF文件见文件夹MEGA内的protein和cds文件夹
氨基酸序列得到的Original Tree见下图1：有个小疑问：氨基酸序列得到的Original Tree直接从MEGA导出时是两大分支，在iTOL是三个大分支，为什么？
```
![Original Tree](./iRtOOle2jbw5HH18GvcYfA.svg)  
图1——氨基酸序列得到的Original Tree
+ 分组
```
很明显，可以分为6个直系同源组（每个物种的基因出现一次），没有complex组（一个物种的一个基因多拷贝或少于一个拷贝）.
```
+ 计算Ka Ks和Π
```
1.  
计算ka ks（6个直系同源组，所以分6次计算，每次里分别算ka ks）  
整理6个直系同源组的序列
 cat cds.all.fa | grep "COL4A1" | cut -d ">" -f 2 >COL4A1.lst  
 faops some cds.all.fa COL4A1.lst ./COL4A1/COL4A1.fa
   
cat cds.all.fa | grep "COL4A2" | cut -d ">" -f 2 >COL4A2.lst  
 faops some cds.all.fa COL4A2.lst ./COL4A2/COL4A2.fa

cat cds.all.fa | grep "COL4A3" | cut -d ">" -f 2 >COL4A3.lst  
 faops some cds.all.fa COL4A3.lst ./COL4A3/COL4A3.fa

cat cds.all.fa | grep "COL4A4" | cut -d ">" -f 2 >COL4A4.lst  
 faops some cds.all.fa COL4A4.lst ./COL4A4/COL4A4.fa

cat cds.all.fa | grep "COL4A5" | cut -d ">" -f 2 >COL4A5.lst  
 faops some cds.all.fa COL4A5.lst ./COL4A5/COL4A5.fa

cat cds.all.fa | grep "COL4A6" | cut -d ">" -f 2 >COL4A6.lst  
 faops some cds.all.fa COL4A6.lst ./COL4A6/COL4A6.fa  

多序列比对，使用MUSCLE 参数默认  
使用MEGA计算Ka Ks，步骤见Dynamic Actin Gene Family Evolution in Primates.md  
结果见Table1和supplement.table  
疑问：文章Dynamic Actin Gene Family Evolution in Primates（以下简称肌动蛋白文章）的Table1的ka/ks是怎么算出的，似乎不是附表里ka/ks的平均值。另外肌动蛋白文章中说Whether in orthologous groups or
complex groups, the average 𝐾𝑎/𝐾𝑠 ratios of most groups(82.4%) are much smaller than 1 (only six groups of average 𝐾𝑎/𝐾𝑠 ratios are greater than 0.5, all of them belong to
complex groups), indicating that the actin genes code highly conserved proteins because of important functions and were under strong negative selection.从𝐾𝑎/𝐾𝑠 ratios可以看出选择压力我理解，但怎么indicating肌动蛋白基因编码高度保守的蛋白呢？  

2.
计算Π
见教程https://blog.csdn.net/qq_50637636/article/details/122612457  
使用的参数见图2  
有一个小疑问：MEGA给出Π的结果是只有对角线没有值，但导出为excel后，少了上面一部分值，其实我个人认为两个序列之间的核苷酸差异和这两个序列前后没有关系，序列1与序列2的差异和序列2与序列1的差异应该一样才是，所以没上面的结果也正常，看了两个肌动蛋白文章附表2的两个序列的PI，序列顺序不同，结果同。 
我得出的Π的结果见Evolution_/IV_collagen下的Table1和supplement.table 
```


![](./MEGA/Π/计算Π时的参数.png)  
图2

![](./MEGA/Π/MEGA上看到的Π结果.png)  
图3——MEGA上看到的Π结果
![](./MEGA/Π/导出为excel后看到的结果.png)  
图4——Π导出为excel后看到的结果  

## 3、鉴定IV型胶原假基因
可能找不到假基因。假若有假基因，由假基因推核苷酸内秉突变率。
+ BLASTN
```
#将IV型胶原基因的核苷酸序列（暂以上步下载的CDS序列做）和对应物种基因组序列比较
#chimpanzee
makeblastdb -in ./Pan_troglodytes.Pan_tro_3.0.dna.toplevel.fa -dbtype nucl -parse_seqids -out ./index
blastn -query ./chimpanzee.cds.fa -db ./index -evalue 1e-5 -outfmt 6 -num_threads 4 -out out_file1
```
#结果解读：输出格式选择 6 （--outfmt 6） ，默认输出为：qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore 。
![BLASTN_format6](./BLASTN/BLASTN_format6.png)
图5 BLASTN_format6  
  
(echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" && cat out_file1) > tem&&
    mv tem out_file1 # 添加表头   

不太清楚怎么从结果中去除已知的胶原蛋白基因。  

换另一种输出格式：blastn -query ./chimpanzee.cds.fa -db ./index -evalue 1e-5 -outfmt 0 -num_threads 4 -out out_file2



## 4、密码子使用情况
使用在线网站使用在线网站计算RSCU，网址：http://cloud.genepioneer.com:9929/#/tool/alltool/detail/214  
在线网站返回的结果log.txt里显示 start codon is wrong(ATG)，但由于后续会把起始密码子删去，所以暂不认为此错误有影响。  
将结果整理为tsv文件（全选在线网站解压出的excel的内容，粘贴到新建的的LF的tsv文件。之所以要用LF是因为tsv-filter对CRLF的文件会报错。注意在VSCode打开一个空tsv或只有一行的tsv，就算改成LF保存再打开仍是CRLF，但有两行文字时由CRLF改成LF保存再打开即为LF.）  
对人六个IV型胶原蛋白基因的cds进行RSCU分析后，用在线软件生成热图，见图6.  
![](./rscu_heatmap.png)  
图6


+ 对数据进行预处理（删除终止密码子和ATG、TGG）
```
cd RSCU
JOB=$(find ./ -maxdepth 2 -type f -name "*.tsv")
for J in $JOB;do
  echo -e "====> $J"
  tsv-filter -H --str-ne Codon:UAA $J | 
  tsv-filter -H --str-ne Codon:UAG | 
  tsv-filter -H --str-ne Codon:UGA |
  tsv-filter -H --str-ne Codon:AUG |
  tsv-filter -H --str-ne Codon:UGG > tem&&
  mv tem $J
done
```
+ 合并
```
cd within_groups/
JOB=$(ls)
tsv-select -H --fields Codon o1.tsv > merge.tsv
for J in $JOB;do
  echo -e "===> $J"
  tsv-join --filter-file $J -H --key-fields Codon --append-fields RSCU merge.tsv > tem
    mv tem merge.tsv
done
```
```
cd within_species/
JOB=$(ls)
tsv-select -H --fields Codon 0.tsv > merge.tsv
for J in $JOB;do
  echo -e "===> $J"
  tsv-join --filter-file $J -H --key-fields Codon --append-fields RSCU merge.tsv > tem
    mv tem merge.tsv
done
```
 
+ 计算方差
但我不确定方差和肌动蛋白文章中说的variation不是一个意思。正在琢磨。暂时以方差进行计算和画图。
```
#利用datamash中的svar计算(brew install datamash)
datamash --help

# 查看矩阵
head merge.tsv

# 计算
cat merge.tsv | datamash transpose | datamash --header-in --header-out svar 2-60 | datamash transpose > svar.tsv
# transpose 转置
# --header-in           first input line is column headers
# --header-out          print column headers as first line
```
+ R画图   
within_groups和within_species分开作图以及合并为簇状条形图都在下方写了过程。展示的图6是簇状条形图。
```
#within_groups
setwd("D:/0~GitHub/Evolution_/IV_collagen/RSCU/within_groups")
group_svar <- read.table("svar_within_groups.tsv",header=FALSE,sep='\t')
library(ggplot2)
p <- ggplot(group_svar,aes(V1,V2))+geom_bar(stat = 'identity')+ylab("variance of RSCU within groups") +xlab("Codons")+ theme(panel.grid = element_blank())
p

#within_species
setwd("D:/0~GitHub/Evolution_/IV_collagen/RSCU/within_species")
species_svar <- read.table("svar_within_species.tsv",header=FALSE,sep='\t')
library(ggplot2)
p <- ggplot(species_svar,aes(V1,V2))+geom_bar(stat = 'identity')+ylab("variance of RSCU within species") +xlab("Codons")+ theme(panel.grid = element_blank())
p
```
```
画簇状条形图
#整理数据
cd R #进入新文件夹
(echo -e "Codons\tvariance_of_RSCU_within_groups" && cat svar_within_groups.tsv) > tem&&
    mv tem svar_within_groups.tsv # 添加表头 
#手动创建cultivar.tsv
paste cultivar.tsv svar_within_groups.tsv >svar_within_groups_paste.tsv #按列合并，也相当于增加新列

(echo -e "Codons\tvariance_of_RSCU_within_species" && cat svar_within_species.tsv) > tem&&
    mv tem svar_within_species.tsv # 添加表头 
#手动创建cultivar.tsv
paste cultivar.tsv svar_within_species.tsv >svar_within_species_paste.tsv

#合并，按行追加 到RSCU目录
cat svar_within_species_paste.tsv >> svar_within_groups_paste.tsv
手动删中间表头
mv svar_within_groups_paste.tsv svar_within_groups_species.tsv
setwd("D:/0~GitHub/Evolution_/IV_collagen/RSCU")
svar <- read.table("svar_within_groups_species.tsv",header=TRUE,sep='\t')
library(ggplot2)
ggplot(svar, aes(x = Codons, y = variance_of_RSCU_within_groups, fill = Cultivar)) +geom_bar(position = "dodge", stat = "identity")+ theme(panel.grid = element_blank())
#见图7,可见物种间有相似的密码子使用模式，相比起来直系同源组间使用模式更多样化。
```
![图7](./RSCU/Figure1_variations_of_RSCU.jpeg)
图7 variations of RSCU

## 5、IV型胶原基因表达模式
进入http://www.ebi.ac.uk/arrayexpress，搜索E-AFMX-11，下载文件 (后续可以找更新一点的表达量数据)  
 A-AFFY-44.adf.txt为提供了探针的各种注释信息  
尝试提取表达值：  
安装affy包 BiocManager::install("affy")  
setwd("D:/0~GitHub/Evolution_/IV_collagen/EXPRESSION/")  
BiocManager::install("GEOquery")  
library(GEOquery)  
write.table(data.frame(fData(study),exprs(study)),file="expression.txt",row.names=FALSE,sep="\t") #提取表达值，但没有基因对应  

WYF师兄帮我整理了表达水平的文件，见expression.txt,该文件每行开头为基因名，每列最上是样本。  
出现的IV型胶原蛋白基因如下：
COL4A1
COL4A1
COL4A2
COL4A2
COL4A3
COL4A3
COL4A3
COL4A3
COL4A3
COL4A3
COL4A3
COL4A4
COL4A4
COL4A4
COL4A5
COL4A5
COL4A6
COL4A6
COL4A6  (即人和黑猩猩的6个IV型胶原基因都有检测到，并且不止一次，可能是不同探针。)

55个样本，其实是两个物种分别5个组织：  
C代表黑猩猩，H代表人，B,H,K,L,T分别代表脑，心，肾，肝，睾丸。CB,CH,CK,CL,CT各5个；HB,HH,HK,HL,HT各6个。所以共计55个样本。

因此求表达水平在不同组织中的平均值：
见excel表expression.xlsx  

整理好的数据见Table2.Average expression level of collagen IV genes from human and chimpanzee genomes in all the five tissues.  
整理用于作图的Table3 The coefficient of variation of the expression level values for collagen IV genes within species and within different groups in all the five tissues. 因为只研究人和黑猩猩的基因表达水平，所以Table3中Orthologous1只包括人和黑猩猩的COL4A1，其它组同理。最终用于作图的数据见Table3的Sheet2，将其整理到tsv文件，见Evolution_/IV_collagen/EXPRESSION/Figure_data.tsv

R中作簇状条形图：
```
setwd("D:/0~GitHub/Evolution_/IV_collagen/EXPRESSION/")
expression <- read.table("Figure_data.tsv",header=TRUE,sep='\t')
library(ggplot2)
ggplot(expression, aes(x = Tissues, y = Average_CV_of_expression_levels, fill = Cultivar)) +geom_bar(position = "dodge", stat = "identity")+ theme(panel.grid = element_blank())
```
结果见图8
![图8](./EXPRESSION/Figure.jpeg)  
图8  
结果表明，within_species在5个组织中平均CV均高于within_groups.表明物种内，比如人，六个基因执行不同功能，有不同的表达水平。而within_groups，比如Orthologous1里的来自人和黑猩猩的COL4A1，可能有相同的功能，平均CV较小。  


注：标准差可缩写为SD，在数学文本和方程中最常用小写希腊字母σ (sigma)表示总体标准差，或拉丁字母s表示样本标准差。——维基百科。注意计算CV(SD/mean)，在EXCEL中STDEV和AVERAGE使用时两个number之间以逗号连接而非加号。







# 二、IV型胶原中G-X-Y的演化探究（在标准演化之外的特征演化方面）
## 旁系同源(人）中GXY
+ 查看进化关系  
```
将人的6个IV型胶原基因的氨基酸序列和核苷酸序列进行多序列比对后建树，查看进化关系。  
其实第一部分标准演化的多序列中已经包含了这部分，在此只是单独做人的。  
将IV_collagen/sequence/下人6个基因的cds和protein复制到IV_collagen/PART2_GXY/Paralogs_GXY  
MEGA进行多序列比对，选muscel，参数默认  
MEGA对多序列比对完的序列建树，NJ法   
结果与第一部分标准演化没有大的区别。输出的pdf文件所在目录为Evolution_/IV_collagen/PART2_GXY/Paralogs_GXY 
```

+ 查看每个基因GXY repeat 中的的中断
```
根据InterProScan查看Collagen repeat  
COL4A5：57始-1460末 22个中断 4个G1G 7个G4G
COL4A1: 46始-1440末 21个中断 7个G1G 5个G4G
COL4A2: 61始-1336末 
COL4A3: 46始-1437末 23个中断 2个G1G 7个G4G
COL4A4: 62始-1370末 24个中断 3个G1G 5个G4G  根据InterProScan查看Collagen repeat,与某篇文章报道的27个中断不同,我又往后延了一段GXY，与文献一样了。其实repeat结束位置后面很长才有零星GXY出现。再核对COL4A4最后一个中断的判断。
COL4A6: 46始-1478末 25个中断 4个G1G 2个G4G
始和末之间并不是一直连续。

手动标注中断中

在手动标注中断的时候发现，有些GXY中的XY也有G，这些或许可以和1号位G对比研究？  
InterPro上说重复的第一个位置是甘氨酸，第二个和第三个位置可以是任何残基但通常是脯氨酸和羟脯氨酸。
中断看哪些较为古老，哪些较为新出。将CLO4A5基因出现的中断分类：一些COL4A5独有的、在6个基因中部分有、在6个基因中都有的（此类较为保守）


```
## 直系同源中GXY
用统计分析部分所选的6个等级30个物种（灵长目到北方兽类）的COL4A5基因的氨基酸序列多序列比对，比较中断。将人科作为参照，看哪些中断是古老是哪些是新的（人有而其它物种没有的可能是新的，人有而其它物种有的可能是老的），人不中断其它物种中断的位置可能不重要，因为允许在其它物种中中断。
中断的长度和位置也需要考虑，因为The location of interruptions is often similar for the three chains within a heterotrimer, while there may be different lengths of the interruptions at a given site in each of the three chains.[Thiagarajan G, Li Y, Mohs A, et al. Common interruptions in the repeating tripeptide sequence of non-fibrillar collagens: sequence analysis and structural studies on triple-helix peptide models. J Mol Biol. 2008;376(3):736-748. doi:10.1016/j.jmb.2007.11.075]

# 可选——考虑其它基因家族的演化（比如将肌动蛋白基因作为对比）





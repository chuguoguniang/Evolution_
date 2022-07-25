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
![](./MEGA/Π/MEGA上看到的Π结果.png)  
图3——MEGA上看到的Π结果
![](./MEGA/Π/导出为excel后看到的结果.png)  
图4——Π导出为excel后看到的结果  

## 3、鉴定IV型胶原假基因
## 4、密码子使用情况
使用在线网站使用在线网站计算RSCU，网址：http://cloud.genepioneer.com:9929/#/tool/alltool/detail/214
在线网站返回的结果log.txt里显示 start codon is wrong(ATG)，但由于后续会把起始密码子删去，所以暂不认为此错误有影响。
将结果整理为tsv文件（全选在线网站解压出的excel的内容，粘贴到新建的的LF的tsv文件）

大概率不会有complex组，可能找不到假基因。假若有假基因，由假基因推核苷酸内秉突变率。

# 二、IV型胶原中G-X-Y的演化探究（在标准演化之外的特征演化方面）
## 旁系同源中GXY
将人的6个IV型胶原基因的氨基酸序列和核苷酸序列进行多序列比对后建树，查看进化关系。
查看每个基因中的中断，看哪些较为古老，哪些较为新出。将CLO4A5基因出现的中断分类：一些COL4A5独有的、在6个基因中部分有、在6个基因中都有的（此类较为保守）
## 直系同源中GXY
用统计分析部分所选的6个等级30个物种（灵长目到北方兽类）的COL4A5基因的氨基酸序列多序列比对，比较中断。将人科作为参照，看哪些中断是古老是哪些是新的（人有而其它物种没有的可能是新的，人有而其它物种有的可能是老的），人不中断其它物种中断的位置可能不重要，因为允许在其它物种中中断。
中断的长度和位置也需要考虑，因为The location of interruptions is often similar for the three chains within a heterotrimer, while there may be different lengths of the interruptions at a given site in each of the three chains.[Thiagarajan G, Li Y, Mohs A, et al. Common interruptions in the repeating tripeptide sequence of non-fibrillar collagens: sequence analysis and structural studies on triple-helix peptide models. J Mol Biol. 2008;376(3):736-748. doi:10.1016/j.jmb.2007.11.075]

# 可选——考虑其它基因家族的演化（比如将肌动蛋白基因作为对比）

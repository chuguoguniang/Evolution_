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
+ 下载7个物种IV型胶原蛋白基因的氨基酸序列和CDS序列 
```
进入ensemble官网下选相应物种，每个基因暂取最长的转录本对应的蛋白序列和CDS序列，手动下载。（其实我想用命令行通过FTP下载再保留最长转录本再提取，但一时没想好怎么用命令行完成，暂且手动下载） 
注：搜索Gorilla基因时Only searching Western Lowland Gorilla; 搜索Gibbon基因时Only searching Northern white-cheeked gibbon  
```
+ 氨基酸序列多序列比对及建树
```



```


```
大概率不会有complex组，可能找不到假基因。假若有假基因，由假基因推核苷酸内秉突变率。

# 二、IV型胶原中G-X-Y的演化探究（在标准演化之外的特征演化方面）
## 旁系同源中GXY
将人的6个IV型胶原基因的氨基酸序列和核苷酸序列进行多序列比对后建树，查看进化关系。
查看每个基因中的中断，看哪些较为古老，哪些较为新出。将CLO4A5基因出现的中断分类：一些COL4A5独有的、在6个基因中部分有、在6个基因中都有的（此类较为保守）
## 直系同源中GXY
用统计分析部分所选的6个等级30个物种（灵长目到北方兽类）的COL4A5基因的氨基酸序列多序列比对，比较中断。将人科作为参照，看哪些中断是古老是哪些是新的（人有而其它物种没有的可能是新的，人有而其它物种有的可能是老的），人不中断其它物种中断的位置可能不重要，因为允许在其它物种中中断。
中断的长度和位置也需要考虑，因为The location of interruptions is often similar for the three chains within a heterotrimer, while there may be different lengths of the interruptions at a given site in each of the three chains.[Thiagarajan G, Li Y, Mohs A, et al. Common interruptions in the repeating tripeptide sequence of non-fibrillar collagens: sequence analysis and structural studies on triple-helix peptide models. J Mol Biol. 2008;376(3):736-748. doi:10.1016/j.jmb.2007.11.075]

# 可选——考虑其它基因家族的演化（比如将肌动蛋白基因作为对比）

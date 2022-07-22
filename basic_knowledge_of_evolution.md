积累的演化知识：
1. 演化的基本假设：相似序列有共同祖先。  
注: 因为一个基因要通过氨基酸（20种）突变产生相似序列概率非常低，几乎为零。（但是这个相似最小是多少程度？）
2. 内秉增长率和内秉突变率  
内秉增长率：按Andrewartha和Brich的定义（1954），内禀增长率是指具有稳定年龄结构的种群，在食物与空间不受限制，同种其它个体的密度维持在最适水平，环境中没有天敌，并在某一特定的温度、湿度、光照和食物性质的环境条件组配下，种群的最大瞬时增长率。内禀增长率反映了种群在理想状态下，生物种群的扩繁能力。--摘自百度百科  
内秉突变率：

3. blastn使用和各个参数  
blastn -help 而非blastn --help获得命令行参数的详细描述  
也可参考网上介绍：https://max.book118.com/html/2020/1231/7061146031003036.shtm   
blastn -task blastn -evalue 1e-3 -num_threads 4 -num_descriptions 10 -num_alignments 10 -outfmt 0 \-dust yes -soft_masking true \ -db S288c.fa -query ../PARS10/sce_genes.fasta -out sce_genes.blast  
-task <String, Permissible values: 'blastn' 'blastn-short' 'dc-megablast'
                'megablast' 'rmblastn' >
   Task to execute
   Default = `megablast'  
  -db <String>
   BLAST database name * Incompatible with:  subject, subject_loc  
   
   -out <File_Out, file name length < 256>
   Output file name
   Default = `-'  
   
   -soft_masking <Boolean>
   Apply filtering locations as soft masks
   Default = `true'  
   
   -evalue <Real>
   Expectation value (E) threshold for saving hits 
   Default = `10'  
     
     -num_descriptions <Integer, >=0>
   Number of database sequences to show one-line descriptions for
   Not applicable for outfmt > 4
   Default = `500'* Incompatible with:  max_target_seqs    
   
   -num_alignments <Integer, >=0>
   Number of database sequences to show alignments for
   Default = `250' * Incompatible with:  max_target_seqs  
    -outfmt <String>
   alignment view options:
     0 = Pairwise,    
      
  
  

   -dust <String>
   Filter query sequence with DUST (Format: 'yes', 'level window linker', or
   'no' to disable)
   Default = `20 64 1'
  
   -soft_masking <Boolean>
   Apply filtering locations as soft masks
   Default = `true'  
   参照王老师GitHub

4. 基因复制  
基因复制后有功能冗余，会发生新功能化，亚功能化或无功能化，这些都是加速进化的结果   
亚功能化即以前的功能分到不同基因，分别执行功能。

5. 同源：  
Orthology：直系同源，描述在不同物种中来自共同祖先的基因。orthologous直系同源的基因可能有相同的功能，也可能没有。
Paralogy：旁系同源，描述在同一物种内由于基因复制而分离的同源基因。

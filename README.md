![logo](transgram.png)
Description
================

TransGram is a genome-duided tool for recovering transcripts from long RNA-seq reads.


Prerequisites
================

  g++ with support for C++11 (e.g. 4.7.2)
  
  [xgboost][xgboost] https://github.com/dmlc/xgboost
 

# 1. Installing and test
===========================================================================
    
    Please ensure that xgboost are properly installed.
   
    (A) Download the lastest version of TransGram and make.
    
          $ git clone https://github.com/yutingsdu/TransGram 
          $ cd TransGram
          $ make release
          
    (B) Test TransGram on the demo data set.
        
        Change to TransGram/sample_test/, and type the following command:
        
          $ ./run_me.sh
          
        If you get the ** congrats, you succesfully installed TransGram.
      
        
===========================================================================

# 2. Usage 
===========================================================================
    
    TransGram v.1.0 usage:

    ** Required **
    
    
---------------------------------------------------------------------------

    ** Options **
    
    --help/-h			  : Output TransGram Help Information.

    --version/-v			  : Print current version of TransGram.

---------------------------------------------------------------------------

    ** Typical Command **
    
    A typical TransGram command might be:

    TransGram -b file.bam -o transgram_outdir

---------------------------------------------------------------------------

===========================================================================


Authors: Ting Yu and Zitong Ren designed and wrote TransGram.
 
Contact:
 
Any questions, problems, bugs are welcome and should be dumped to Ting Yu <yutingsdu@163.com>
 
Created on July 15, 2024.
[xgboost]: https://github.com/dmlc/xgboost
 
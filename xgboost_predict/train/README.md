Description
================

TransGram is a genome-guided transcriptome reconstruction tool for long RNA-seq reads.

This is a function of TransGram-v.0.0 ：You can train a new filtering model with your sample.
In situations where the default model you though is not suitable for the sample, TransGram provides a function to train a new model. To ensure that the trained model captures data features effectively, the data used for training should be well annotated.The new model can be used on the same or similar species.

Prerequisites
================

  g++ with support for C++11 (e.g. 4.7.2)
  
  Please ensure that [xgboost][xgboost], StringTie2 (v2.2.1) and cufflinks (v2.2.1) is properly installed !

# 1. Installing and test
===========================================================================

    (A) Installing xgboost, StringTie2(v2.2.1) and cufflinks(v2.2.1)
    
          $ pip install xgboost

          $ wget https://github.com/gpertea/stringtie/releases/download/v2.2.1/stringtie-2.2.1.Linux_x86_64.tar.gz
          $ tar xvfz ~/Downloads/stringtie-2.2.1.Linux_x86_64.tar.gz
          
          $ wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
          $ tar zxvf ~/Downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz


    (B) Download the lastest version of TransGram and make.
    
          $ git clone https://github.com/yutingsdu/TransGram 
          $ cd TransGram
          $ make release
          
    (C) An example training a model
        
          # train a new model using a sample1.bam and a reference genome annotation
          export PATH=~/Downloads/stringtie:$PATH
          export PATH=~/Downloads/cufflinks:$PATH
          train_model -b sample1.bam -a annotations.gtf

          # use the trained model on another sample
          ..........(sample2.bam)

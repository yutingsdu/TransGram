Description
================

TransGram is a genome-guided transcriptome reconstruction tool for long RNA-seq reads.
This is a function of TransGram: You can train a new filtering model with your sample.
In situations where the default model you though is not suitable for the sample, TransGram provides a function to train a new model. 
To ensure that the trained model captures data features effectively, the data used for training should be well annotated.The new model 
can be used on the same or similar species.

The users should provide a BAM file generated by minimap2, an annotation file to train a new model.


Prerequisites
================

Please ensure that [xgboost][xgboost], [TransGram][TransGram], [StringTie2][StringTie2], and [cufflinks][cufflinks] are properly installed 
and add the software binary directory to the PATH environment variable in the shell configuration file
!

An illustration of installing the used tools.
================

Installing xgboost
    
    $ pip install xgboost

Installing StringTie2

    $ wget https://github.com/gpertea/stringtie/releases/download/v2.2.1/stringtie-2.2.1.Linux_x86_64.tar.gz
    $ tar xvfz ~/Downloads/stringtie-2.2.1.Linux_x86_64.tar.gz

Installing cufflinks

    $ wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
    $ tar zxvf ~/Downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz

Installing TransGram.

    $ git clone https://github.com/yutingsdu/TransGram
    $ cd TransGram
    $ make release


Please don't forget to set the  Set the enviroment variable for StringTie2 and cufflinks.

OR

Type the command

export PATH=your-path/stringtie:$PATH

export PATH=your-path/cufflinks:$PATH

befor running, and /your-path/ is the directory where StringTie2 and cufflinks installed




An example for training a new model
================
    
    train_model -b sample1.bam -a annotations.gtf
    
    If everything goes well, you will see a directory "transgram_new_model"

Runing TransGram with the new model
================
    
    TransGram -b alignment.bam --ont -o transgram-outdir --CusModel transgram_new_model

[xgboost]: https://github.com/dmlc/xgboost
[TransGram]: https://github.com/yutingsdu/TransGram
[StringTie2]: https://github.com/skovaka/stringtie2 
[cufflinks]: http://cole-trapnell-lab.github.io/cufflinks/
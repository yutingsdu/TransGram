import os
import sys
import random
from tqdm import tqdm

source_dir_path = '/home/yuting/yuting/ThirdG-transcriptome-Guided'
def find_dataset(source_dir_path = '/home/yuting/yuting/ThirdG-transcriptome-Guided', name = 'R13'):
    all_list = os.listdir(source_dir_path)
    for _name in all_list:
        prefix = _name.split('_')[0]
        if prefix == name:
            file_list = os.listdir(f'{source_dir_path}/{_name}')
            for _file in file_list:
                if '.fastq' in _file:
                    return f'{source_dir_path}/{_name}/{_file}'
    return ''
    '''
    ###### test command ######
    _list = [i for i in range(1, 19)]
    for item in _list:
        x = find_dataset(source_dir_path, f'R{item}')
        print(x)
    ###### test command ######
    '''
import gzip
def read_fq(fq_path):
    reads = []
    x = open(fq_path, "r")
    if 'gz' in fq_path:
        x = gzip.open(fq_path, "rt")
    with x as f:
        while True:
            header = f.readline()
            if not header:
                break 
            seq = f.readline().strip()   
            _ = f.readline()             
            _ = f.readline()         
            reads.append(seq)       
    return reads
    '''
    ###### test command ######
    x = read_fq(f'{source_dir_path}/R13_Data_PacBio_3/SRR22522147_1.fastq')
    print(len(x), x[0])
    ###### test command ######
    '''
    
def read_fa(fa_path):
    reads = []
    x = open(fa_path, "r")
    if 'gz' in fa_path:
        x = gzip.open(fa_path, "rt")
    for line in x:
        if line[0] == '>':
            continue
        else:
            reads.append( line.strip('\n') )
    return reads

MeanLenR13 = 1974
def sample_reads_by_base_limit(read_list, base_limit, out_fa, mode = 'Uniform', len_cutoff = MeanLenR13):
    random.shuffle(read_list)  #
    total_bases = 0
    selected_reads = []

    if mode == 'Uniform':
        for read in tqdm(read_list):
            read_len = len(read)
            selected_reads.append(read)
            total_bases += read_len
            if total_bases > base_limit:
                break
            
    elif mode == 'Short':   # sampling short reads (length < mean length of all reads in raw fastq file)
        for read in tqdm(read_list):
            read_len = len(read)
            if read_len >= len_cutoff:
                continue
            selected_reads.append(read)
            total_bases += read_len
            if total_bases > base_limit:
                break
    
    elif mode == 'Long':    # sampling long reads (length >= mean length of all reads in raw fastq file)
        for read in tqdm(read_list):
            read_len = len(read)
            if read_len < len_cutoff:
                continue
            selected_reads.append(read)
            total_bases += read_len
            if total_bases > base_limit:
                break
            
    if total_bases < base_limit - 10000:
        print(f'The number of bases ({total_bases}) sampled cannot be met!')
        return
        
    _file = open(out_fa, 'w')
    for i in range( len(selected_reads)):
        _file.write( f'>{i}\n' )
        _file.write(f'{selected_reads[i]}\n')
    _file.close()
    '''
    ###### test command ######
    sample_reads_by_base_limit(R13_reads, 3_000_000_000, out_fa = 'u_sampling_R13_3GB.fa', mode = 'Uniform')
    ###### test command ######
    '''

import statistics

def get_fastq_infomation(read_lst):
    len_lst = [len(item) for item in read_lst]
    total_bases = sum( len_lst )
    ava_len = statistics.mean(len_lst)
    stdev_p = statistics.pstdev(len_lst)
    print(f'{total_bases}\t{ava_len}\t{stdev_p}')
    '''
    _list = [i for i in range(1, 19)]
    for item in _list:
        x = find_dataset(source_dir_path, f'R{item}')
        if x != '':
            read_lst = read_fq(x)
            get_fastq_infomation(read_lst)
        else:
            print(f'None\tNone\tNone')
            
    read_lst = read_fq('SRR14286070_1.fastq.gz')
    get_fastq_infomation(read_lst)
    read_lst = read_fq('SRR14286054_1.fastq.gz')
    get_fastq_infomation(read_lst)
    read_lst = read_fq('SRR14286061_1.fastq.gz')
    get_fastq_infomation(read_lst)
    
    uniform_gb_size_lst = [0.2, 0.4, 0.6, 0.8, 1, 3, 5, 7, 9, 11]
    for _size in uniform_gb_size_lst:
        fa_file = f'u_sampling_R13_{_size}GB.fa'
        read_lst = read_fa(fa_file)
        get_fastq_infomation(read_lst)
        
    for t in range(3):
        read_lst = read_fa(f's_sampling_R13_3GB_lessthan_1k_{t}.fa')
        get_fastq_infomation(read_lst)
    '''

def uniform_sampling_from_R13():
    uniform_gb_size_lst = [10,20,30]
    R13_reads = read_fq(fq_path = f'{source_dir_path}/R13_Data_PacBio_3/SRR22522147_1.fastq')
    for _size in uniform_gb_size_lst:
        real_size = _size * 1_000_000_000
        sample_reads_by_base_limit(R13_reads, real_size, out_fa = f'u_sampling_R13_{_size}GB.fa', mode = 'Uniform')
        
    '''
    uniform_sampling_from_R13()
    '''
        
def predict_all_uniform_data_use_TransGram():
    sequencinf_type = 'pacbio'
    uniform_gb_size_lst = [10,20,30]
    for _size in uniform_gb_size_lst:
        fa_file = f'u_sampling_R13_{_size}GB.fa'
        cmd = f'bash get_datatype.sh {fa_file} {sequencinf_type}'
        os.system(cmd)
        
def short_long_sampling_from_R13():
    random_count = 3
    R13_reads = read_fq(fq_path = f'{source_dir_path}/R13_Data_PacBio_3/SRR22522147_1.fastq')
    s_gb_size_lst = [3]
    for _size in s_gb_size_lst:
        real_size = _size * 1_000_000_000
        for t in range(random_count):
            sample_reads_by_base_limit(R13_reads, real_size, out_fa = f's_sampling_R13_{_size}GB_{t}.fa', mode = 'Short')
            
    for _size in s_gb_size_lst:
        real_size = _size * 1_000_000_000
        for t in range(random_count):
            sample_reads_by_base_limit(R13_reads, real_size, out_fa = f'l_sampling_R13_{_size}GB_{t}.fa', mode = 'Long')
    
    l_gb_size_lst = [1]
    for _size in l_gb_size_lst:
        real_size = _size * 1_000_000_000
        for t in range(random_count):
            sample_reads_by_base_limit(R13_reads, real_size, out_fa = f'l_sampling_R13_{_size}GB_{t}.fa', mode = 'Long')
            
    '''
    R13_reads = read_fq(fq_path = f'{source_dir_path}/R13_Data_PacBio_3/SRR22522147_1.fastq')
    for t in range(0, 3):
        sample_reads_by_base_limit(R13_reads, 3_000_000_000, out_fa = f's_sampling_R13_3GB_lessthan_1k_{t}.fa', mode = 'Short', len_cutoff=1000)
    '''
            
def get_datatype_of_TransGram(our_dir_TransGram):
    _file = f'{our_dir_TransGram}/data.info'
    with open( _file ) as f:
        lines = f.readlines()
        res = lines[-1].strip('\n')
        return res

if __name__ == "__main__":
    
    
    ######################## Explore whether coverage affects data types (alpha or beta) ########################
    uniform_sampling_from_R13()   # uniform samping
    predict_all_uniform_data_use_TransGram() # predict data type by TransGram
    ######################## Explore whether coverage affects data types (alpha or beta) ########################
    
    # ######################## Explore whether coverage affects data types (alpha or beta) ########################
    # uniform_sampling_from_R13()   # uniform samping
    # predict_all_uniform_data_use_TransGram() # predict data type by TransGram
    # ######################## Explore whether coverage affects data types (alpha or beta) ########################
    
    '''
    ######################## Explore whether read length affects data types (alpha or beta) ########################
    ######################## Explore whether read length affects data types (alpha or beta) ########################
    '''
    
    # short_long_sampling_from_R13()
    
    # for t in range(3):
    #     x = get_datatype_of_TransGram(f's_sampling_R13_3GB_lessthan_1k_{t}')
    #     print(x)
    
    # x = find_dataset(source_dir_path, f'R2')
    # print(x)
    
    
    
    
    
    
       
############################################################################
# Copyright (c) 2022-2024 Shandong University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


import re
import random

def select_transcripts(input_gtf, output_gtf,proportion):
    
    # 第一步：收集所有唯一的transcript_id
    transcript_ids = set()
    with open(output_gtf, 'w') as outfile:
        with open(input_gtf, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue  # 跳过注释行
                parts = line.strip().split('\t')
                attributes = parts[8]
                # 使用正则表达式提取transcript_id
                match = re.search(r'transcript_id\s+"([^"]+)"', attributes)
                if match:
                    transcript_ids.add(match.group(1))
                    if proportion == 75:
                        if len(transcript_ids) % 4 == 1:
                            continue
                        else:
                            outfile.write(line)
                    elif proportion == 50:
                        if len(transcript_ids) % 2 == 1:
                            continue
                        else:
                            outfile.write(line)
                    elif proportion == 25:
                        if len(transcript_ids) % 4 != 0:
                            continue
                        else:
                            outfile.write(line)

def select_transcripts(input_gtf, output_gtf, input_gtf2):
    lines_to_exclude = set()
    with open(input_gtf2, 'r') as f2:
        for line in f2:
            lines_to_exclude.add(line.strip())

    with open(output_gtf, 'w') as outfile:
        with open(input_gtf, 'r') as f:
            for line in f:
                if line.strip() not in lines_to_exclude:
                    outfile.write(line)

###################################################################################
### Example usage:

input_path = "hg38.ncbiRefSeq.gtf "
output_path = "selected_transcripts_75.gtf"
output_path_groundtruth = "selected_transcripts_25_groundtruth.gtf"

# Randomly select 75% of the transcripts from the annotation file as input.
select_transcripts(input_path, output_path, 75)
# The remaining 25% of the transcripts are used as groundtruth for evaluation.
select_transcripts(input_path, output_path_groundtruth,output_path)

###################################################################################
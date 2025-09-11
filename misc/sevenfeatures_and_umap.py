import numpy as np
from scipy.special import kl_div
import numpy as np
# from importlib_metadata import version
import umap
import umap.umap_ as umap
import matplotlib.pyplot as plt
# from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
from sklearn.neighbors import KNeighborsClassifier

def read_gene_data(gene_file):
    info_lst = []
    with open(gene_file) as f:
        lines = f.readlines()
        for line in lines:
            info = line.strip('\n').split(' ')
            info = [int(item) for item in info]
            info_lst.append( info )
    return info_lst

def cal_percentile(info_lst, num_idx = 3):
    gene_num_lst = [item[ num_idx ] for item in info_lst]
    p50 = np.percentile(gene_num_lst, 50)
    p75 = np.percentile(gene_num_lst, 75)
    return p50, p75

import pandas as pd
def cal_n50(_l):
    df = pd.DataFrame({"Read": _l})
    sorted_lengths = np.sort(df["Read"])[::-1]
    cumulative_sum = np.cumsum(sorted_lengths)
    n50_position = np.where(cumulative_sum >= 0.5 * cumulative_sum[-1])[0][0]
    _n50 = sorted_lengths[n50_position]
    return _n50

def cal_n75(_l):
    df = pd.DataFrame({"Read":_l})
    sorted_lengths = np.sort(df["Read"])[::-1]
    cumulative_sum = np.cumsum(sorted_lengths)
    n75_position = np.where(cumulative_sum >= 0.75 * cumulative_sum[-1])[0][0]
    _n75 = sorted_lengths[n75_position]
    return _n75

def calculate_features(flie,length_file,kk):
    p = []
    q = []
    d=[]
    
    with open(flie, mode='r', encoding='utf-8') as file:
        lines = file.readlines()
        for line in lines[1:-2]:
            columns = line.split()
            if len(columns) >= 3:
                p.append(int(columns[1]))
                q.append(int(columns[2]))
                d.append(float(columns[-1]))
                
    mean_depth=np.mean(d)
    depth_50=np.median(d)
    depth_75=np.percentile(d, 75)

    #calculate the average value and variance of proportion of reliable
    result = []
    m=0
    for p_val, q_val in zip(p, q):
        result.append(q_val / p_val)
        if q_val / p_val > kk:
            m=m+1
    result=np.array(result)
    mean_proportion=np.mean(result)
    variance_proportion=np.var(result)

    count = np.sum(result>kk)
    proportion_of_realible_gene=count/len(p)
    file_out.write(str(datatype_list[k]) + "," + str(m/(len(p)-m)) + "\n")
    
    print(flie,datatype_list[k],m/(len(p)-m))
    
    
    l=[]
    with open(length_file, mode='r', encoding='utf-8') as file:
        lines = file.readlines()
        for line in lines:
            columns = line.split()
            l.append(int(columns[1]))
    base_number= np.sum(l)
    mean_length=np.mean(l)
    variance_length=np.var(l)

    return [depth_50,depth_75,base_number,mean_length,variance_length]


if __name__ == "__main__":
    # data_list = ["data1_path:data1_type", "data2_path:data2_type", ...]

    file_name_datatype = {"/home/yuting/yuting/ThirdG-transcriptome-Guided/R1_Data4-stringtie/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/R2_Data7/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/R3_Data_ONT_BB_1/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/R4_Data_ONT_BB_7/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/R5_Data_ONT_BB_8/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/R6_Data_NC_R1/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/R7_Data_NM_R4_HCC827_ONT/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/R8_NM_R2_HCC827_PacBio/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"beta",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/R9_Data_NM_R3_H1975_PacBio/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"beta",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/R10_Data_NM_R7_PacBio_HCC827/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"beta",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/R11_Data_PacBio_1/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha",
                # "/home/yuting/yuting/ThirdG-transcriptome-Guided/R12_Data_PacBio_6/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/R13_Data_PacBio_3/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"beta",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/R14_Data_PacBio_4/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/R15_Data_PacBio_5/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/T1_Data_NM_R1_H1975_ONT/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/T2_Data_NM_R6_H1975_PacBio/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"beta",
                # "/home/yuting/yuting/ThirdG-transcriptome-Guided/Revision2/R1/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/Revision2/R2/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"beta",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/Revision2/R4/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"beta",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/Revision2/R15/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/Revision2/R16/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha",
                # "/home/yuting/yuting/ThirdG-transcriptome-Guided/Revision2/R4/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/Revision2/R17/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"beta",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/Revision2/R18/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"beta",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/Revision2/R19/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"beta",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/Revision2/R20/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha",
                # "/home/yuting/yuting/ThirdG-transcriptome-Guided/Revision/R57/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha",
                "/home/yuting/yuting/ThirdG-transcriptome-Guided/Revision/R7/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt":"alpha"
                # "/home/yuting/yuting/ThirdG-transcriptome-Guided/YASIM-sim2/Revision/transgram_outdir_v1_rzt-RL/ReadNum_forGene.txt"
    }

    data_list = list(file_name_datatype.keys())
    datatype_list = list(file_name_datatype.values())
    res_l_50=[]
    res_l_75=[]
    for file_n in data_list:
        #
        file_n = file_n + "ReadNum_forGene.txt"
        info_lst = read_gene_data( file_n)
        gene_num_lst = [item[3] for item in info_lst]
        res = cal_percentile(info_lst, num_idx = 3) # num_idx = 3, indicate calculate the percentile based on the last column in the file
        res_l_50.append(res[0])
        res_l_75.append(res[1])
    print(res_l_50,res_l_75)

    with open("output_0.2.txt", "w") as file_out:
        # for kk in [1/3,1/4,1/5,1/6,1/7]:
        for kk in [1/5]:
            k=0
            all_feature_list=[]
            for i in data_list:
                feature_list=calculate_features(i+"/TG-data.info",i+"/transgram.reads.length",kk)

                feature_list.append(res_l_50[k])
                feature_list.append(res_l_75[k])
                k=k+1
                all_feature_list.append(feature_list)

        #################################Apply UMAP for dimensionality reduction.#################################
        all_feature_list_normalization = np.log1p(all_feature_list) 
        umap_model = umap.UMAP(n_components=2,n_neighbors=5,random_state=42) 
        X_umap = umap_model.fit_transform(all_feature_list_normalization)

        alpha_indices = [index for index, value in enumerate(datatype_list) if value == "alpha"]
        beta_indices = [index for index, value in enumerate(datatype_list) if value == "beta"]

        alpha_points  = X_umap[alpha_indices,:] 

        beta_points =  X_umap[beta_indices,:] 

        X = np.vstack((alpha_points, beta_points))  
        y = np.array(['Alpha datasets'] * len(alpha_points) + ['Beta datasets'] * len(beta_points))

        plt.figure(figsize=(6, 6))  
        plt.scatter(alpha_points[:, 0], alpha_points[:, 1], color="#239fcc", label='Alpha datasets')
        plt.scatter(beta_points[:, 0], beta_points[:, 1], color="#D23B70", label='Beta datasets')

        plt.legend(fontsize=12)
        plt.savefig("umap.svg")


































    # path="/home/yuting/yuting/ThirdG-transcriptome-Guided/Revision2"
    # file_name = ["/R5/TransGram_all_features/ReadNum_forGene.txt",
    #              "/R6/TransGram_all_features/ReadNum_forGene.txt",
    #              "/R8/TransGram_all_features/ReadNum_forGene.txt",
    #              "/R9/TransGram_all_features/ReadNum_forGene.txt",
    #              "/R10/TransGram_all_features/ReadNum_forGene.txt",
    #              "/R11/TransGram_all_features/ReadNum_forGene.txt",
    #              "/R12/TransGram_all_features/ReadNum_forGene.txt",
    #              "/R13/TransGram_all_features/ReadNum_forGene.txt",
    #              "/R14/TransGram_all_features/ReadNum_forGene.txt"]
    

    # file_name = [path + item for item in file_name]             
    # # txt file like 'R1_Data4-stringtie/unmatch_model/TransGram_ReadNum_forGene/ReadNum_forGene.txt'
    # res_l_50=[]
    # res_l_75=[]
    # for file_n in file_name:
    #     info_lst = read_gene_data( file_n)
    #     gene_num_lst = [item[3] for item in info_lst]
    #     # res = cal_percentile(info_lst, num_idx = 3) # num_idx = 3, indicate calculate the percentile based on the last column in the file
    #     res_l_50.append(cal_n50(gene_num_lst))
    #     res_l_75.append(cal_n75(gene_num_lst))
    # print(res_l_50,res_l_75)

    # file_name = ["/home/yuting/yuting/Software/lrgasp-simulation/MySimulation_new_5M_2/simulated/Revision/transgram_outdir_newversion/ReadNum_forGene.txt",
    #              "/home/yuting/yuting/Software/lrgasp-simulation/Revision2/MySimulation_Revision_5M/simulated/transgram_outdir_v1/ReadNum_forGene.txt",
    #              "/home/yuting/yuting/ThirdG-transcriptome-Guided/SpikeIn-IsoQuant/unmatch_model/TransGram_all_features/ReadNum_forGene.txt",
    #              ""]
    # res_l_50=[]
    # res_l_75=[]
    # for file_n in file_name:
    #     info_lst = read_gene_data( file_n)
    #     gene_num_lst = [item[3] for item in info_lst]
    #     # res = cal_percentile(info_lst, num_idx = 3) # num_idx = 3, indicate calculate the percentile based on the last column in the file
    #     res_l_50.append(cal_n50(gene_num_lst))
    #     res_l_75.append(cal_n75(gene_num_lst))
    # print(res_l_50,res_l_75) 
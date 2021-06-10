import numpy
import numpy as np
import pandas as pd
import scipy
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.model_selection import KFold
import random

mator_d_p_01='../data/MATADOR/chem_protein_01_1.xlsx'
d_p_1_01='../data/data1_01.xlsx' #4547 association
d_p_1_01_scores='../data/d_p_1_01_scores.xlsx'
test_file_path='../data/test.xlsx'


''' Resource allocation algorithm -----Calculating the similarity matrix'''
def ra_similarity(x):
    ss=sum(x)
    ss[ss==0]=1
    result_1 = x /ss  # !!!!!sum by column
    y = np.array(result_1)  # Resource allocation algorithm----first propagation
    # print("first propagation：\n",y)

    '''second propagation： part 2.1'''
    result_2 = []
    y = np.squeeze(y)  # remove a dimension
    for i in range(y.shape[0]):
        if x[i, :].sum() == 0:  # If the column of the test set is empty, the division error is avoided
            continue
        result_2.append(y[i, :] / np.count_nonzero(y[i, :]))
    z = np.array(result_2)
    # print("second propagation：\n", z)

    '''second propagation： part 2.2'''
    z_result = np.zeros((z.shape[1], z.shape[1]))
    for i in range(z.shape[1]):  # !!!!
        ii = np.where(z[:, i] > 0)
        # list(ii)[0] Index with non-zero per column     sum(axis=0) Array add by row
        z_result[i] = z[list(ii)[0], :].sum(axis=0)
    return z_result


''' CF model ----- Score matrix calculation'''
def cf_model(inteartion_matrix): # inteartion_matrix (m,n)；cos_maxrix   Cosine similarity
    simi_maxrix=ra_similarity(inteartion_matrix) #(n,n)

    row_sum_matrix=np.sum(simi_maxrix,axis=1)# --shape(1,n)
    #diag()  Diagonal matrix； pinv() Pseudo inverse matrix；
    sum_diagonal_matrix=np.linalg.pinv(np.diag(row_sum_matrix))#(n,n)

    row_normalized_protein_similarity_matrix=np.dot(sum_diagonal_matrix,simi_maxrix)#(n,n)*(n,n)= Characteristic matrix
    score_matrix=np.dot(inteartion_matrix,row_normalized_protein_similarity_matrix)#(m,n)*(n,n)= Score matrix
    print('1_score_matrix:\n',score_matrix)
    score_matrix[np.where(interaction_matrix == 1)]=0
    pd.DataFrame(score_matrix).to_excel(d_p_1_01_scores,header=None,index=None)
    print('2_score_matrix:\n',score_matrix)
    return score_matrix

if __name__ == '__main__':

    interaction_matrix = pd.read_excel(d_p_1_01, header=None)
    interaction_matrix = numpy.array(interaction_matrix)
    cf_model(interaction_matrix)





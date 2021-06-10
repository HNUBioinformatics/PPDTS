import  pandas as pd
import  numpy as np



scores_path='../data/drug_protein/d_p_1_01_scores.xlsx'
max_column_values='../data/test4.xlsx'


column_index_rowsMax_Resultcolumn='../data/zhang_lnc_protein/column_index_rowsMax_Resultcolumn.xlsx'

#Replace 1 with 0 in the original score matrix
def association_1_to_0(path):
    df=pd.read_excel(path,header=None)
    df.replace(1,0,inplace=True) #Must pay attention:!!! There is no save here, just save the replacement object in the variable df
    df.to_excel(path,index=False,header=False)#Do not create a new row column index
    print("1 to 0 ok")

#Returns the index of the maximum value of each line element
# def row_max_index(path,to_path):
def row_max_index(path):
    arr=pd.read_excel(path)

    arr_max_index=arr.idxmax(axis=1)
    arr_max_vla=arr.max(axis=1)
    print("Returns the index of the maximum value of each line element：\n",arr_max_index)
    # print("-----------------")
    arr_max_index.to_excel(max_column_values, index=False, header=False)  # Do not create a new row column index

#Find the row with the same elements except for non-zero elements --- call twice + 1  
def find_rows_no_zero(path):
    df = pd.read_excel(path, header=None)
    # Here, replace 1 in the original score matrix with 0
    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            if df.iloc[i, j] == 0:
                df.iloc[i, j] = df.iloc[i, j - 1]  #！！！！----- Replace the value of 0 by twice j+1 and j-1
    df.to_excel(path, index=False, header=False)  # Do not create a new row column index

def get_protein_result(path):
    max_rowsValue=pd.read_excel(path,header=None)
    max_rowsValue=np.array(max_rowsValue)
    protein_file=pd.read_excel(path,header=None)
    protein_file=np.array(protein_file)
    column_index_rowsMax_Result=[]#Save the prediction of in each row - protein results


    for i in range(len(max_rowsValue)):#990
        # print(protein_file[max_rowsValue[i]].flatten())
        column_index_rowsMax_Result.append((protein_file[max_rowsValue[i]]).flatten())

    print(column_index_rowsMax_Result)
    # print(list(column_index_rowsMax_Resultcolumn))
    # print(list(column_index_rowsMax_Resultcolumn)[0])
    pd.DataFrame(np.array(column_index_rowsMax_Result)).to_excel(column_index_rowsMax_Resultcolumn,header=None,columns=None)


if __name__== '__main__':
    row_max_index(scores_path)
    # find_rows_no_zero(path)
    # get_protein_result(to_path)









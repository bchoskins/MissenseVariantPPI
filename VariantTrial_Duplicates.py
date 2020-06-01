#PPI Prediction with a Siamese RCNN 
import pandas as pd 
import sys
import os
import tensorflow as tf
from tensorflow import keras as keras
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras.layers import Dense, Activation, Dropout, Embedding, LSTM, Bidirectional, BatchNormalization, add
from tensorflow.keras.layers import Flatten, Reshapequit
from tensorflow.keras.layers import Concatenate, concatenate, subtract, multiply
from tensorflow.keras.layers import Conv1D
from tensorflow.keras.layers import MaxPooling1D, AveragePooling1D, GlobalAveragePooling1D
from tensorflow.keras.optimizers import Adam,  RMSprop
#from tensorflow.keras.backend import tensorflow_backend as KTF
import numpy as np
from tqdm import tqdm
from tensorflow.keras.layers import Input, LSTM, GRU
from numpy import linalg as LA
import scipy

variantProtein = pd.read_csv('/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/duplicatedTestSeqs5.csv')

def ppiPipeline(data):
    embDict = {}
    #embDict = {'A': 1, 'G':2, 'V':3, 'I':4, 'L':5, 'F':6, 'P':7, 'Y':8, 'M':9, 'T':10, 'S':11, 'H':12, 'N':13, 'Q':14, 'W':15, 'R':16, 'K':17, 'D':18, 'E':19, 'C':20, 'U':21}
    embDict = {'A': 1, 'C':2, 'D':3, 'E':4, 'F':5, 'G':6, 'H':7, 'I':8, 'K':9, 'L':10, 'M':11, 'N':12, 'P':13, 'Q':14, 'R':15, 'S':16, 'T':17, 'U':18, 'V':19, 'W':20, 'Y':21} #'X':21,
    #loop through seq_index to get the sequences, split them, embed them with above dict, one hot encode the values, store into vector matrix 
    def aminoEncoding(data):
      encode_list = []
      for row in variantProtein['variantSeq_matrix'].values:
        row_encode = []
        for code in row:
          print(code)
          row_encode.append(embDict.get(code, 0))
        encode_list.append(np.array(row_encode))
    
      return encode_list
    
    #proteins = pd.read_csv('/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/AllProteins.csv', index_col=0)
    Protein_Encode = aminoEncoding(variantProtein)
    
    #######By this stage there is just encoded proteins (1-22) each in their own array in one large array########
    
    #############ALREADY PADDED SO?????#############
    from tensorflow.keras.preprocessing.sequence import pad_sequences
    #a little more than one standard deviation above the third quartile lengths
    #2500 instead of 1000 since coming from R missense variant work
    max_length = 1000
    protein_pad = pad_sequences(Protein_Encode, maxlen=max_length, padding='post', truncating='post')
    protein_pad.shape
    
    from tensorflow.keras.utils import to_categorical
    padVariant = to_categorical(protein_pad)
    padVariant.shape
    #(249,1000,22)
    type(padVariant)
    
    #np.save("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/variantProteinDuplicatesEmbedded4.npy", padVariant)
    
    #variant_good = np.load('/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/variantProteinDuplicatesEmbedded4.npy')
    variant_good = padVariant
    variant_good.shape
    #(249, 1000, 22)
    
    A_pad = np.array([variant_good[np.where(variantProtein['symbol'][1] == variantProtein['symbol'])]])
    B_pad = np.array([variant_good[np.where(variantProtein['symbol'][len(variantProtein)-1] == variantProtein['symbol'])]])
    
    
    A = np.array(A_pad)
    B = np.array(B_pad)
    
    type(A)
    #<class 'numpy.ndarray'>
    A.shape
    #(1, 248, 1000, 22)
    
    A = A[0,:,:,:]
    B = B[0,:,:,:]
    
    type(A)
    #<class 'numpy.ndarray'>
    A.shape
    #(248, 1000, 22)
    
    #np.save("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/DuplicatesA_pad4.npy", A)
    #np.save("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/DuplicatesB_pad4.npy", B)
    
    ####Post processing######
    
    #run on UI-GPU or UI-GPU-HM
    #A = np.load("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/DuplicatesA_pad.npy")
    #B = np.load("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/DuplicatesB_pad.npy")
    A.shape
    #duplicate B 247 times
    B = np.repeat(B, 134, axis=0)
    B.shape
    
    model3 = keras.models.load_model("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/SeqModel_3.h5")
    
    # evaluate loaded model on test data
    
    model3.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    
    score = model3.predict([A, B], verbose=0)
    
    listScore = score.tolist()
    
    flat_score = []
    for sublist in listScore:
        for item in sublist:
            flat_score.append(item)
    
    import csv
    with open("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/predScoresDuplicated_5.csv", 'w', newline='') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow(flat_score)
    
    seqs_score = variantProtein.drop(variantProtein.index[len(variantProtein)-1])
    
    seqs_score['Score'] = flat_score
    ######HERE used excel to attach predScoresDuplicated_1.csv to duplicatedTestSeqs.csv to make duplicatedTestSeqs_Scored.csv
    
    #seqs_score = pd.read_csv('/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/duplicatedTestSeqs_Scored.csv')
    
    probands = seqs_score[seqs_score['Proband'] == 1]
    #probands = variantProtein[variantProtein['Proband'] == 1]
    proband_scores = probands['Score']
    
    siblings = seqs_score[seqs_score['Sibling'] == 1]
    #siblings = variantProtein[variantProtein['Sibling'] == 1]
    sibling_scores = siblings['Score']
    
    # Wilcoxon signed-rank test
    from numpy.random import seed
    from numpy.random import randn
    from scipy.stats import wilcoxon
    from scipy.stats import ranksums
    # seed the random number generator
    seed(1)
    # compare samples
    stat_ps, p_ps = ranksums(proband_scores, sibling_scores)
    return(stat_ps, p_ps)

test = ppiPipeline(variantProtein)
print('Wilcoxon_proband_vs_sibling = %.3f, p = %.3f' % (test[0], test[1])) 

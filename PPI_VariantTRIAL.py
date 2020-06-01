#PPI Prediction with a Siamese RCNN 
import pandas as pd 
import sys
import keras
from keras.models import Sequential, Model
from keras.layers import Dense, Activation, Dropout, Embedding, LSTM, Bidirectional, BatchNormalization, merge, add
from keras.layers.core import Flatten, Reshape
from keras.layers.merge import Concatenate, concatenate, subtract, multiply
from keras.layers.convolutional import Conv1D
from keras.layers.pooling import MaxPooling1D, AveragePooling1D, GlobalAveragePooling1D
from keras.optimizers import Adam,  RMSprop
import os
import tensorflow as tf
import keras.backend.tensorflow_backend as KTF
import numpy as np
from tqdm import tqdm
from keras.layers import Input, CuDNNGRU, GRU, CuDNNLSTM
from numpy import linalg as LA
import scipy

#import os
#os.system('top')

#data = pd.read_csv('/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/UniprotModelDataWithNegatives.csv')
#proteins = pd.read_csv('/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/AllProteins.csv', index_col=0)
#protein_good = np.load('/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/JustProteinEmbeddings.npy')

variantProtein = pd.read_csv('/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/testSeqs.csv')

#creating dictionary for embedding
#####EDITED TO MATCH RDATA MAPPING WITH 'X'###########
embDict = {}
#embDict = {'A': 1, 'G':2, 'V':3, 'I':4, 'L':5, 'F':6, 'P':7, 'Y':8, 'M':9, 'T':10, 'S':11, 'H':12, 'N':13, 'Q':14, 'W':15, 'R':16, 'K':17, 'D':18, 'E':19, 'C':20, 'U':21}
embDict = {'A': 1, 'C':2, 'D':3, 'E':4, 'F':5, 'G':6, 'H':7, 'I':8, 'K':9, 'L':10, 'M':11, 'N':12, 'P':13, 'Q':14, 'R':15, 'S':16, 'T':17, 'U':18, 'V':19, 'W':20, 'X':21, 'Y':22}
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
from keras.preprocessing.sequence import pad_sequences
#a little more than one standard deviation above the third quartile lengths
#2500 instead of 1000 since coming from R missense variant work
max_length = 2500
protein_pad = pad_sequences(Protein_Encode, maxlen=max_length, padding='post', truncating='post')
protein_pad.shape

from keras.utils import to_categorical
padVariant = to_categorical(protein_pad)
padVariant.shape
#(249,2500,23)
type(padVariant)


from keras.utils import to_categorical
variant_good = to_categorical(padVariant)
variant_good.shape
type(variant_good)

np.save("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/variantProteinEmbeddings.npy", variant_good)

variant_good = np.load('/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/variantProteinEmbeddings.npy')

variant_good.shape
#(249, 2500, 23, 2)

###figure out if order is important here#####
binary = pd.read_csv('/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/UniprotModelDataBinary.csv',index_col=0)

variant_A_pad = np.array([variant_good[np.where("ACOT7" == variantProtein['symbol'])]])
variant_B_pad = np.array([variant_good[np.where("CHD8" == variantProtein['symbol'])]])

#seq_A_pad = np.array([protein_good[np.hstack(np.where(x == proteins['Protein_A']))] for x in binary['Protein_A']])
#seq_B_pad = np.array([protein_good[np.hstack(np.where(x == proteins['Protein_A']))] for x in binary['Protein_B']])

variant_A = np.array(variant_A_pad)
variant_B = np.array(variant_B_pad)

type(variant_A)
#<class 'numpy.ndarray'>
variant_A.shape
#(1, 248, 2500, 23, 2)

#seq_A = np.array(seq_A_pad)
#seq_B = np.array(seq_B_pad)

variant_A = variant_A[0,:,:,:,0]
variant_B = variant_B[0,:,:,:,0]

type(variant_A)
#<class 'numpy.ndarray'>
variant_A.shape
#(248, 2500, 23)

#seq_A = seq_A[:,0,:,:]
#seq_B = seq_B[:,0,:,:]

np.save("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/VariantA_pad.npy", variant_A)
np.save("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/VariantB_pad.npy", variant_B)

####Post processing######

#run on UI-GPU or UI-GPU-HM

#seq_A = np.load("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/SeqA_pad.npy")
#seq_B = np.load("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/SeqB_pad.npy")

var_A = np.load("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/VariantA_pad.npy")
var_B = np.load("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/VariantB_pad.npy")

model3 = keras.models.load_model("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/SeqModel_3.h5")

# evaluate loaded model on test data
model3.compile(loss='mse', optimizer='adam', metrics=['accuracy'])

score = model3.evaluate(var_A, var_B, verbose=0)

#print("%s: %.2f%%" % (model3.metrics_names[1], score[1]*100))


def build_model():
    seq_input1 = Input(shape=(2500, 23), name='seq1')
    seq_input2 = Input(shape=(2500, 23), name='seq2')
    l1=Conv1D(hidden_dim, 3)
    r1=Bidirectional(CuDNNGRU(hidden_dim, return_sequences=True))
    l2=Conv1D(hidden_dim, 3)
    r2=Bidirectional(CuDNNGRU(hidden_dim, return_sequences=True))
    l3=Conv1D(hidden_dim, 3)
    r3=Bidirectional(CuDNNGRU(hidden_dim, return_sequences=True))
    l4=Conv1D(hidden_dim, 3)
    r4=Bidirectional(CuDNNGRU(hidden_dim, return_sequences=True))
    l5=Conv1D(hidden_dim, 3)
    r5=Bidirectional(CuDNNGRU(hidden_dim, return_sequences=True))
    l6=Conv1D(hidden_dim, 3)
    s1=MaxPooling1D(3)(l1(seq_input1))
    s1=concatenate([r1(s1), s1])
    s1=MaxPooling1D(3)(l2(s1))
    s1=concatenate([r2(s1), s1])
    s1=MaxPooling1D(3)(l3(s1))
    s1=concatenate([r3(s1), s1])
    s1=MaxPooling1D(3)(l4(s1))
    s1=concatenate([r4(s1), s1])
    s1=MaxPooling1D(3)(l5(s1))
    s1=concatenate([r5(s1), s1])
    s1=l6(s1)
    s1=GlobalAveragePooling1D()(s1)
    s2=MaxPooling1D(3)(l1(seq_input2))
    s2=concatenate([r1(s2), s2])
    s2=MaxPooling1D(3)(l2(s2))
    s2=concatenate([r2(s2), s2])
    s2=MaxPooling1D(3)(l3(s2))
    s2=concatenate([r3(s2), s2])
    s2=MaxPooling1D(3)(l4(s2))
    s2=concatenate([r4(s2), s2])
    s2=MaxPooling1D(3)(l5(s2))
    s2=concatenate([r5(s2), s2])
    s2=l6(s2)
    s2=GlobalAveragePooling1D()(s2)
    merge_text = multiply([s1, s2])
    x = Dense(100, activation='linear')(merge_text)
    x = keras.layers.LeakyReLU(alpha=0.3)(x)
    x = Dense(int((hidden_dim+7)/2), activation='linear')(x)
    x = keras.layers.LeakyReLU(alpha=0.3)(x)
    main_output = Dense(1, activation='sigmoid')(x)
    merge_model = Model(inputs=[seq_input1, seq_input2], outputs=[main_output])
    return merge_model


hidden_dim = 25
n_epochs=50
batch_size = 256
adam = Adam(lr=0.0001, amsgrad=True, epsilon=1e-6) 

#class label train, validate, test split
class_lbl =  binary.iloc[:, 9].as_matrix()

pos = np.where(class_lbl == 1)
neg = np.where(class_lbl == 0)

import random
random.seed(4)
np.random.shuffle(pos[0])
np.random.shuffle(neg[0])

#train, test indeces
pos = np.array((pos[0][0:125000], pos[0][125000:]))
neg = np.array((neg[0][0:520000], neg[0][520000:]))

train = list(pos[0]) + list(neg[0])
test =  list(pos[1]) + list(neg[1])

np.random.shuffle(train)
np.random.shuffle(test)

train_test = []
train_test.append((train, test))

es = keras.callbacks.EarlyStopping(monitor='val_loss', min_delta=0,patience=5,verbose=0, mode='auto')

cp = keras.callbacks.ModelCheckpoint(filepath="/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/VariantModel_Test.h5",
        verbose=1, save_best_only=True)

from tensorflow.python.client import device_lib
print(device_lib.list_local_devices())

for train, test in train_test:
  merge_model = None
  merge_model = build_model()
  adam = Adam(lr=0.0001, amsgrad=True, epsilon=1e-6)
  merge_model.compile(optimizer=adam, loss='mse', metrics=['accuracy'])
  history = merge_model.fit([variant_A[train], variant_B[train]], class_lbl[train], callbacks=[es,cp], batch_size=batch_size, epochs=n_epochs, validation_split=0.2) 

hist_df = pd.DataFrame(history.history) 
hist_df.to_csv("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/variantHistory1.csv")
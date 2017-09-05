import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from collections import defaultdict
import math
import time

class Heuristic_Search():
    def __init__(self,row, col, preprocess, preprocess_selection, cleanup_selection, khard, score_function):
        self.row = row
        self.col = col
        self.tHardmax = None
        self.preprocess = preprocess
        self.preprocess_selection = preprocess_selection
        self.cleanup_selection = cleanup_selection
        self.score_function = score_function
        self.A = None  #going to be updated
        self.AO = None #original
        self.PM = None #Poles Meter list
        self.MP = None #Meter Poles list
        self.PU = []   #Poles Used for it
        self.Tscore = None #Score is Total number of uncovered nodes
        self.tHard = None  #t-hard to cover
        self.score = None  #score of thard for each pole
        self.khard = khard  # maximum number of t-hard for scoring
        self.supmat = None # support matrix to compute score of each pole


    def read_input(self, data):
        rowindex = np.asarray(data[0], dtype=int)
        colindex = np.asarray(data[1], dtype=int)
        L = max(len(rowindex), len(colindex))
        self.supmat = np.ones((self.row, 1))
        val = np.ones(L) # each value of the matrix is one if in the data else zero
        self.A = csr_matrix((val,(rowindex, colindex)), shape=(self.row, self.col)).toarray()
        self.AO = np.copy(self.A)
        self.PoleMeterList(rowindex, colindex)
        self.MeterPoleList(rowindex, colindex)

    def Search(self):
        iter = 0
        preprocess_selection = self.preprocess_selection
        cleanup_selection = self.cleanup_selection
        self.Hard_To_Cover()
        self.Tscore_func() # Total number of uncovered nodes
        self.initial_score_k() # initial k score for all meters
        print('{0:3s} {1:10s} {2:10s}'.format("iter", "Pol.Old", "Pol.New"))
        #print "  Iter  ", "  Meter_Cover_Old  ", "  Meter_Cover_Curr  "
        while True:
            if (self.preprocess!=0 and (iter / float(preprocess_selection)).is_integer() and (iter / float(preprocess_selection))>0):
                print 'preprocessing', iter
                self.Preprocess()
                #print 'preprocessing end'

            [max_value, max_index] = self.compute_max_score_index()
            if(max_value>0.0):self.update(max_index) # update P and M
            else:
                Lold = len(self.PU)
                self.clean_up_new()
                Lnew = len(self.PU)
                print('{0:3d} {1:5d} {2:9d}'.format(iter, Lold, Lnew))
                break;
            iter = iter + 1
            if((iter/float(cleanup_selection)).is_integer()):
                Lold = len(self.PU)
                self.clean_up_new()
                Lnew = len(self.PU)
                print('{0:3d} {1:5d} {2:9d}'.format(iter, Lold, Lnew))

    def PoleMeterList(self, rowindex, colindex):
        PM = defaultdict(list)
        for i, j in zip(rowindex, colindex):
            PM[j].append(i)
        self.PM = PM

    def MeterPoleList(self, rowindex, colindex):
        MP = defaultdict(list)
        for i, j in zip(rowindex, colindex):
            MP[i].append(j)
        self.MP = MP

    def Preprocess(self):

        col_ones = np.ones((self.col, 1), dtype=int)
        row_ones = np.ones((self.row, 1), dtype=int)


        RSR = False #removing the singleton rows
        RNDC= False #removing non -dominent columns
        RDR = False #removing dominent rows

        while(True):
            ###################################################################
            ################## Removing the Singleton Rows ####################
            ###################################################################
            B = np.reshape(self.A.dot(col_ones), self.row).astype(int)
            meterindexs = np.asarray(np.where(B == 1))
            meterindexs = meterindexs.flatten()
            RSR = False
            if (np.size(meterindexs) != 0):
                for meterind in meterindexs:
                    for poleind in self.MP[meterind]:
                        if(self.A[meterind, poleind]== 1):
                            self.Tscore = self.Tscore - np.sum(self.A[self.PM[poleind],: ], axis= 0).T
                            self.A[self.PM[poleind],: ] = 0 # we want to make zeros for sll those polls that contain it
                            self.Tscore[poleind] =0   # update the Tscores
                            self.PU.append(poleind)   # add pole to the list
                            ############
                            for j in self.PM[poleind]:
                                if self.tHard[j] != 0.0:
                                    for k in self.MP[j]:
                                        if self.score[k ,int(self.tHard[j])]>0:
                                            self.score[k, int(self.tHard[j])] = self.score[k, int(self.tHard[j])] - 1
                            self.score[poleind, :] = 0
                            self.tHard[self.PM[poleind]] = 0
                            RSR = True

            ###############################################################
            ####################Removing Non-Dominent Column###################
            ###############################################################
            RNDC = False
            p = np.reshape(self.A.T.dot(row_ones), self.col).astype(int)
            B = self.A.T.dot(self.A)
            for i in self.PM:
                if p[i]>0:
                    for j in self.PM[i]:
                        breakj=0
                        for k in self.MP[j]:
                            if i != k and B[i, k] == p[i] and p[k] > p[i]:
                                p[i] = 0
                                B[i, :] = 0
                                B[:, i] = 0
                                self.A[self.PM[i], i] = 0
                                self.Tscore[i] = 0

                                for l in self.PM[i]:
                                    if self.tHard[l] != 0.0:
                                        for m in self.MP[l]:
                                            if self.score[m ,int(self.tHard[l])]>0:
                                                self.score[m, int(self.tHard[l])] = self.score[m, int(self.tHard[l])] - 1
                                for l in self.PM[i]:
                                    if self.tHard[l] > 0: self.tHard[l] = self.tHard[l] - 1

                                for l in self.PM[i]:
                                    if self.tHard[l] != 0.0:
                                        for m in self.MP[l]:
                                            self.score[m, int(self.tHard[l])] = self.score[m, int(self.tHard[l])] + 1

                                self.score[i, :] = 0
                                breakj=1
                                RNDC = True
                                break
                        if breakj==1: break

            ###############################################################
            ####################Removing Dominent Row###################
            ###############################################################
            RDR = False
            p = np.reshape(self.A.dot(col_ones), self.row).astype(int)
            B = self.A.dot(self.A.T)

            for i  in self.MP:
                if p[i] > 0:
                    for j in self.MP[i]:
                        for k in self.PM[j]:
                            if i!=k and B[i, k]==p[i] and p[k]>p[i] >0:
                                p[k]=0
                                B[k, :]=0
                                B[:, k]=0
                                self.Tscore = self.Tscore - self.A[k, :].T
                                self.A[k, :]=0
                                if(self.tHard[k]!=0):
                                    for l in self.MP[k]:
                                        if self.score[l, int(self.tHard[k])]>0:
                                            self.score[l, int(self.tHard[k])] = self.score[l, int(self.tHard[k])] - 1
                                self.tHard[k] = 0
                                RDR = True

            #print RSR, RNDC, RDR

            if(RSR==False and RNDC==False and RDR==False): break


    def compute_max_score_index(self):
        if(self.score_function=="greedy_score"):
            [max_value, max_index] = self.greedy_score()
            return max_value, max_index

        if(self.score_function=="modified_greedy_score"):
            [max_value, max_index] = self.modified_greedy_score()
            return max_value, max_index

    def update(self, max_score_pole):
        Tscore = self.Tscore[:]
        score  = self.score[:]
        for j in self.PM[max_score_pole]:
            Tscore = Tscore - self.A[j, :].T

            if self.tHard[j]!=0.0 :
                for k in self.MP[j]:
                    if(score[k, int(self.tHard[j])]!=0):
                        score[k, int(self.tHard[j])] = score[k, int(self.tHard[j])] - 1
            self.A[j, :] = 0.0
            self.tHard[j]= 0.0

        self.Tscore[:] = Tscore[:]
        self.PU.append(max_score_pole)
        score[max_score_pole, :]= 0.0
        self.score[:] = score[:]

    def clean_up_new(self):
        AO = self.AO
        col_ones = np.zeros((self.col, 1), dtype=int)
        col_ones[self.PU, 0] = 1

        metercover = np.reshape(self.AO.dot(col_ones), self.row).astype(int)
        metercoverindex = np.where(metercover>0)
        PURem =[]

        for i in self.PU:
            metercover = metercover - AO[:, i].astype(int)
            Ind = np.where(metercover[metercoverindex]<1)
            if(np.size(Ind)==0):
                PURem.append(i)
            else :
                metercover = metercover + AO[:, i].astype(int)
        for i in PURem:
            self.PU.remove(i)

    def get_input(self):
        return self.A

    def greedy_score(self):
        Tscore = self.Tscore
        max_index = np.argmax(Tscore, axis=0)
        max_value = Tscore[max_index]
        return max_value, max_index

    def modified_greedy_score(self):
        S =self.score_k(self.khard)
        max_index = np.argmax(S, axis=0)
        max_value = S[max_index]
        return max_value, max_index


    def Hard_To_Cover(self):
        Ones = np.ones((self.col, 1), dtype=int)
        L = np.reshape(self.A.dot(Ones), self.row).astype(int)
        self.tHardmax = max(L) +1
        self.tHard = L - np.min(L[1:]) + 1
        self.tHard[0] = 0


    def Tscore_func(self):
        A = self.A
        supmat = np.ones((self.row, 1))
        self.Tscore = np.reshape(A.T.dot(supmat), self.col)

    def initial_score_k(self):
        score = np.zeros((self.col, self.tHardmax))
        for i in range(1, self.col):
            score[i, self.tHard[self.PM[i]]] = score[i, self.tHard[self.PM[i]]] + \
                                                   self.Lfunc(self.tHard[self.PM[i]])[0, self.tHard[self.PM[i]]]
        self.score = score


    def score_k(self, k):
        if k < self.tHardmax:
            score = self.score
            #initial
            S = np.ones(self.col, dtype=float)
            S[0] = 0 # we are not using value at 0 so keep that value as zero
            Sold = np.ones(self.col)
            m = min(self.tHard[1:])

            ### fixing a bug ########
            if(m==0): m=1

            for j in range(m, k+1):
                Sold[:] = S[:]
                for i in range(1,  self.col):
                    S[i] = S[i] * math.pow(score[i, j], 1.00/float(j)) # j has to be greater than zero
                if np.sum(S)==0.0:
                    S[:] = Sold[:]
            if np.sum(S)!=0.0: S = np.multiply(S, self.Tscore)
            else: S = np.multiply(1.00, self.Tscore)
            return S

        else:
            print "error k is more than hardmax"



    def Lfunc(self, set):
        L= np.zeros((1, max(set)+1))
        for j in set:
            L[0, j] = L[0, j]+ 1
        return L


def Set_Covering_Algorithm(data, preprocess, preprocess_selection, cleanup_selection, scoring_type, khard):
    t0 = time.clock()  # initial time stamping
    data = pd.read_csv(data, header=None, sep='\s+')
    row = max(data[0]) + 1
    col = max(data[1]) + 1
    gs = Heuristic_Search(row, col, preprocess, preprocess_selection, cleanup_selection, khard, scoring_type)  # greedy score and modified greedy score
    gs.read_input(data)  # converting data into a sparse matrix
    gs.Search()  # search is performed using greedy and modified greedy sore
    t1 = time.clock()

    print "Time Taken in seconds ", t1 - t0
    print "Number_Of_Poles_used", len(gs.PU)

    ##########################################
    #################TEST#####################
    ##########################################
    row_ones = np.zeros(row, dtype=int)
    for i in gs.PU:
        for j in gs.PM[i]:
            row_ones[j] = 1

    print "Check Number Of Meters Covered", np.sum(row_ones, axis=0)


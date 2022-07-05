from __future__ import print_function
import sys

import cplex
from cplex.exceptions import CplexError
import random
import time
from math import log, exp, fabs
from scipy import stats
from numpy import array, std, var
# import Queue
import numpy
import os

class Coviddata_MD():
    def __init__(self):
        self.group = 8
        self.period = 30
        compartment = 11

        self.ch = 1   # ward cost
        self.cu = 4   # icu ward cost
        self.cx = 12531   # hospitalization cost
        self.cy = 49127   # icu cost
        self.cd = 10   # death loss
        self.cts = 62   # test service cost
        self.ct = 62   # test kit cost
        self.ci = 1 # infection cost

        self.parameter = {'betaD': [0.2238],
                            'alphaD': [0.25],
                            'w': [0.5237, 0.5314, 0.5298, 0.4974, 0.4506, 0.4436, 0.3743, 0.3767],
                            'lambda_a': [115],
                            'lambda_b': [-0.2466],
                            't0': [77],
                            'policy0_63': [1],
                            'policy63_84': [1.1657],
                            'policy84_99': [1.2334],
                            'eta': [0.8808, 0.95, 0.9240, 0.9383, 0.9327, 0.9520, 0.9636, 0.9886],
                            'tau': [0.28087, 0.25877, 0.21773, 0.20165, 0.20956, 0.20843, 0.19064, 0.19771],
                            'a9_64': [0.0528, 0.0918, 0.0791, 0.0789, 0.0388, 0.0535, 0.0048, 0.0171],
                            'b9_64' : [5394, 11250, 11320, 12220, 5571, 8214, 586.6, 2841],
                            'a64_99' : [0.1996, 0.1140, 0.1101, 0.1111, 0.0985, 0.1057, 0.0986, 0.1347],
                            'b64_99' : [19320, 12270, 14010, 15450, 12890, 14200, 14720, 19700],
                            'piD': [0.00393, 0.04265, 0.06982, 0.09805, 0.129, 0.212, 0.4525, 0.943],
                            'ricuD': [0.008, 0.0399, 0.0600, 0.0857, 0.1587, 0.2161, 0.3363, 0.3507],
                            'dqD': [0.0032, 0.0058, 0.0101, 0.0175, 0.0345, 0.0671, 0.1121, 0.1708],
                            'dhD': [0.0185, 0.0154, 0.0336, 0.0886, 0.1227, 0.1496, 0.1964, 0.2366],
                            'dicuD': [0.0050, 0.0042, 0.0092, 0.0241, 0.0335, 0.0408, 0.0536, 0.0645],
                            'gammaD': [0.1052],
                            'lambdahD': [0.0985],
                            'lambdaicuD': [0.1477],
                            'betaO': [2.23],
                            'alphaO': [4/3],
                            'piO': [0.4800, 0.8000, 0.8000, 0.8000, 0.8000, 0.8000, 0.7360, 0.7440],
                            'ricuO': [0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23],
                            'dqO': [0.7469, 0.7469, 0.7469, 0.7469, 0.7469, 0.7469, 0.7469, 0.7469],
                            'dhO': [0.01, 0.01, 0.01, 0.01, 0.05, 0.14, 0.14, 0.14],
                            'dicuO': [0.21, 0.21, 0.21, 0.21, 0.21, 0.21, 0.21, 0.21],
                            'gammaO': [7.1106],
                            'lambdahO': [1.13],
                            'lambdaicuO': [1]}

        self.population = 6045680
        group_prop_population = [25, 13, 14, 13, 14, 12, 7, 4]
        self.group_population = [self.population * group_prop_population[i] / 100 for i in range(self.group)]

        S0 = [1270072,481980,576174,481147,727908,510327,339509,198787]
        E0 = [1162.80000000000,1456.35000000000,1847.75000000000,1814.50000000000,
              1207.45000000000,918.650000000000,218.500000000000,217.550000000000]
        I0 = [452, 916, 938, 623, 615, 474, 361, 90]
        IM0 = [1202.03720310646, 2157.07080705323, 1805.81130925095, 1077.09248576521,
               1067.05282699620, .029526235742, 358.170336501901, 9.64118155893537]
        IS0 = [1.44803425457376, 6.49691729323308, 13.0797773654917, 16.5651710540339,
               18.8132457228630, 17.1742797118848, 8.65883887801696, 8.25675558509870]
        H0 = [10.3124487004104, 51.5622435020520, 103.124487004104, 136.353488372093,
              167.290834473324, 167.290834473324, 100.832831737346, 100.832831737346]
        U0 = [1.97975376196990, 9.89876880984952, 19.7975376196990, 26.1767441860465,
              32.1160054719562, 32.1160054719562, 19.3575923392613, 19.3575923392613]
        DH0 = [10, 48, 120.428571428571, 309.142857142857,
               851.571428571429, 1683.28571428571, 2486.28571428571, 4314.42857142857]
        D0 = [20.69926920611367, 173.00310016634484, 194.8343354976536, 1142.5938896060973,
              3471.90759228442, 6385.222548732182, 10816.167815812161, 7891.229242380346]
        R0 = [self.group_population[j] - S0[j] - E0[j] - I0[j] - IM0[j] - IS0[j] - U0[j] - H0[j] - D0[j]
              for j in range(self.group)]

        self.S = [1221909.83603031, 463439.364496720, 554075.397653914, 463801.289194988,
                704095.100923903, 493887.265367804, 330257.405627254, 193335.810556363]
        self.E = [4394.45484417625, 1691.29554638348, 2015.95732101117, 1583.96898978713,
                2177.70582345691, 1503.75225529999, 848.073650631724, 499.661491065978]
        self.I = [2132.18885677411, 846.990806140762, 1054.52503363789, 850.487557691953,
                1180.30660802910, 845.196594728111, 528.697965563804, 366.840870229350]
        self.IM = [1004.62563529715, 360.006730910099, 371.175697140658, 272.109943174622,
                 375.576375794728, 245.841127849355, 105.960757995562, 23.4804399656852]
        self.IS = [1.50493519792026, 7.45226362802630, 13.1479097586580, 13.6842912978167,
                 25.7739038244388, 29.7659912411751, 34.1538987876377, 49.2739328920495]
        self.H = [12.4457765939109, 55.0459842985154, 96.5822945833096, 100.848029560380,
                176.169377338243, 192.933320949182, 197.280505867159, 245.313844738088]
        self.U = [0.0497034581802129, 1.12871766842250, 3.01796155089700, 4.58952916625200,
                15.7647048702726, 24.7111976809538, 43.9536269199907, 57.7828054414295]
        self.DH = [12.0330127793271, 55.9204165594445, 152.193999624179, 402.116978451194,
                 1047.54902126482, 1942.08506998150, 2803.91257238655, 4720.09411272132]
        self.D = [23.8291172357537, 192.523822160612, 272.041231125302, 1330.82634519198,
                  3941.34579531014, 7268.61987389383, 12273.7074959894, 10159.4564449962]
        self.R = [self.group_population[j] - self.S[j] - self.E[j] - self.I[j] - self.IM[j] - self.IS[j]
                  - self.U[j] - self.H[j] - self.D[j] for j in range(self.group)]

        self.compartment_population = [self.S, self.E, self.I, self.IM, self.IS, self.H, self.U,
                                       [0 for j in range(self.group)],
                                       [0 for j in range(self.group)]]
        self.compartment_population = [[self.compartment_population[n][j] for j in range(self.group)]
                            for n in range(len(self.compartment_population))]



        lam = [1 / (1 + self.parameter['lambda_a'][0] * exp(self.parameter['lambda_b'][0] * (t+100-self.parameter['t0'][0])))
               for t in range(self.period)]
        
        self.w = [self.parameter['w'][j] * self.parameter['policy84_99'][0] for j in range(self.group)]
        betaD = self.parameter['betaD'][0]
        betaO = self.parameter['betaD'][0] * self.parameter['betaO'][0]
        
        rho = self.parameter['a64_99']
        self.bhat = self.parameter['b64_99']
        
        alphaD = self.parameter['alphaD'][0]
        alphaO = self.parameter['alphaD'][0] * self.parameter['alphaO'][0]
        
        gammaD = self.parameter['gammaD'][0]
        gammaO = self.parameter['gammaD'][0] * self.parameter['gammaO'][0]
        
        piD = self.parameter['piD']
        piO = [self.parameter['piD'][j] * self.parameter['piO'][j] for j in range(self.group)]
        
        dqD = [self.parameter['dqD'][j] for j in range(self.group)]
        dqO = [self.parameter['dqD'][j] * self.parameter['dqO'][j] for j in range(self.group)]
        dicuD = [self.parameter['dicuD'][j] for j in range(self.group)]
        dicuO = [self.parameter['dicuD'][j] * self.parameter['dicuO'][j] for j in range(self.group)]
        dhD = [self.parameter['dhD'][j] for j in range(self.group)]
        dhO = [self.parameter['dhD'][j] * self.parameter['dhO'][j] for j in range(self.group)]
        lambdahD = [self.parameter['lambdahD'][0] for j in range(self.group)]
        lambdaicuD = [self.parameter['lambdaicuD'][0] for j in range(self.group)]
        lambdahO = [self.parameter['lambdahD'][0] * self.parameter['lambdahO'][0] for j in range(self.group)]
        lambdaicuO = [self.parameter['lambdaicuD'][0] * self.parameter['lambdaicuO'][0] for j in range(self.group)]
        self.pm = [[(1-lam[t]) * (1-piD[j]) + lam[t] * (1-piO[j]) for t in range(self.period)] for j in range(self.group)]
        self.pq = [[(1-lam[t]) * piD[j] + lam[t] * piO[j] for t in range(self.period)] for j in range(self.group)]
        
        self.beta = [(1 - lam[t]) * betaD + lam[t] * betaO for t in range(self.period)]
        q_EtoI = [[(1 - lam[t]) * alphaD + lam[t] * alphaO for t in range(self.period)] for j in range(self.group)]
        q_ItoR = [[(1-lam[t]) * gammaD * (1-piD[j]) + lam[t] * gammaO * (1-piO[j])
                   for t in range(self.period)] for j in range(self.group)]
        q_ItoD = [[(1-lam[t]) * dqD[j] * piD[j] + lam[t] * dqO[j] * piO[j] for t in range(self.period)] for j in range(self.group)]
        q_MtoR = [[(1-lam[t]) * gammaD + lam[t] * gammaO for t in range(self.period)] for j in range(self.group)]
        q_QtoD = [[(1-lam[t]) * dqD[j] + lam[t] * dqO[j]
                   for t in range(self.period)] for j in range(self.group)]
        q_HtoD = [[(1-lam[t]) * dhD[j] * lambdahD[j] + lam[t] * dhO[j] * lambdahO[j]
                   for t in range(self.period)] for j in range(self.group)]
        q_HtoR = [[(1-lam[t]) * (1-dhD[j]) * lambdahD[j] + lam[t] * (1-dhO[j]) * lambdahO[j]
                   for t in range(self.period)] for j in range(self.group)]
        
        self.q = {('E', 'I'): q_EtoI, ('I', 'R'): q_ItoR,  ('I', 'D'): q_ItoD,  ('M', 'R'): q_MtoR,
            ('Q', 'D'): q_QtoD,  ('H', 'D'): q_HtoD,  ('H', 'R'): q_HtoR}


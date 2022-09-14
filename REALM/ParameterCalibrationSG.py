from __future__ import print_function
import random
from math import exp
import REALM
import xlwt, xlrd
import itertools

user = REALM.REALM()

group = 10
period = 120
user.set_group(group=group)
user.set_time(time=period)
compartment = ['S', 'E', 'I', 'C', 'VS', 'VE', 'VI', 'VC', 'BS', 'BE', 'BI', 'BC','M','VM','BM']
compart_num = len(compartment)
user.set_all_compartment(name=compartment)
total_pop = [263846,277179,397083,845064,1134500,925569,720661,609424,332222,181590]




def readfile(filename:str, sheetname:str, row:list, column:list):
    wb = xlrd.open_workbook(filename=filename)
    sheet = wb.sheet_by_name(sheetname)
    dat = []
    for i in range(row[0], row[1]):
        cells = sheet.row_values(i)
        # print(cells)
        dat.append([float(cells[j]) for j in range(column[0], column[1])])
    return dat

def model(population, alpha, tau, beta, w, policy, x1, x2, Chat):
    for i in range(compart_num):
        user.set_population(compartment=compartment[i], population=population[i])

    user.set_transition_compartment(n='E', m='I', prob=[[alpha for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='VE', m='VI', prob=[[alpha for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='BE', m='BI', prob=[[alpha for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='I', m='C', prob=[[tau*p[j][0] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='VI', m='VC', prob=[[tau*p[j][1] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='BI', m='BC', prob=[[tau*p[j][2] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='I', m='M', prob=[[tau*(1-p[j][0]) for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='VI', m='VM', prob=[[tau*(1-p[j][1]) for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='BI', m='BM', prob=[[tau*(1-p[j][2]) for t in range(period)] for j in range(group)])
    for j in range(group):
        for t in range(period):
            betaSE = [[0 for k in range(group)] for n in range(len(population))]
            for n in range(len(population)):
                if compartment[n] in ['I']:
                    for k in range(group):
                        betaSE[n][k] = beta[0] * w[j][k] * policy[0][j][t] / total_pop[k]
                elif compartment[n] in ['VI']:
                    for k in range(group):
                        betaSE[n][k] = beta[1] * w[j][k] * policy[1][j][t] / total_pop[k]
                elif compartment[n] in ['BI']:
                    for k in range(group):
                        betaSE[n][k] = beta[2] * w[j][k] * policy[2][j][t] / total_pop[k]
            user.set_transition_compartment_dp(flow=('S', 'E', j, t), beta=betaSE)

    for j in range(group):
        for t in range(period):
            betaSE = [[0 for k in range(group)] for n in range(len(population))]
            for n in range(len(population)):
                if compartment[n] in ['I']:
                    for k in range(group):
                        betaSE[n][k] = beta[0] * w[j][k] * policy[0][j][t] / total_pop[k]
                elif compartment[n] in ['VI']:
                    for k in range(group):
                        betaSE[n][k] = beta[1] * w[j][k] * policy[1][j][t] / total_pop[k]
                elif compartment[n] in ['BI']:
                    for k in range(group):
                        betaSE[n][k] = beta[2] * w[j][k] * policy[2][j][t] / total_pop[k]
            user.set_transition_compartment_dp(flow=('VS', 'VE', j, t), beta=betaSE)

    for j in range(group):
        for t in range(period):
            betaSE = [[0 for k in range(group)] for n in range(len(population))]
            for n in range(len(population)):
                if compartment[n] in ['I']:
                    for k in range(group):
                        betaSE[n][k] = beta[0] * w[j][k] * policy[0][j][t] / total_pop[k]
                elif compartment[n] in ['VI']:
                    for k in range(group):
                        betaSE[n][k] = beta[1] * w[j][k] * policy[1][j][t] / total_pop[k]
                elif compartment[n] in ['BI']:
                    for k in range(group):
                        betaSE[n][k] = beta[2] * w[j][k] * policy[2][j][t] / total_pop[k]
            user.set_transition_compartment_dp(flow=('BS', 'BE', j, t), beta=betaSE)

    x = [[[[0 for t in range(len(x1))] for j in range(group)] for m in range(compart_num)] for n in range(compart_num)]
    x[0][4] = [[x1[t][j] for t in range(len(x1))] for j in range(group)]
    x[4][8] = [[x2[t][j] for t in range(len(x2))] for j in range(group)]
    user.x = x

    Xsim_exp = user.Prediction(opt=['D'])
    Cbar = Xsim_exp[compartment.index('C')]
    VCbar = Xsim_exp[compartment.index('VC')]
    BCbar = Xsim_exp[compartment.index('BC')]
    obj = sum([((Cbar[j][t] + VCbar[j][t] + BCbar[j][t]) / Chat[t][j] - 1)**2
               for j in range(group) for t in range(30)]) / (group * 30)
    return obj, Xsim_exp


if __name__ == '__main__':
    bestobj = 1e+20
    bestpopulation = []
    bestalpha = 0
    besttau = 0
    bestbeta = []
    bestw = []

    x1p = readfile(filename='Vaccination and WWT All Numbers and Charts.xlsx',
                  sheetname='AGE (N)', row=[217-6,248+30+90], column=[1, 11])
    x2p = readfile(filename='Vaccination and WWT All Numbers and Charts.xlsx',
                  sheetname='AGE (N)', row=[217-6,248+30+90], column=[13, 23])
    Cnewp = readfile(filename='Vaccination and WWT All Numbers and Charts.xlsx',
                  sheetname='INFECTION (N)', row=[216-6,247+30+90], column=[7, 17])
    cuminfection = readfile(filename='Vaccination and WWT All Numbers and Charts.xlsx',
                  sheetname='INFECTION (N)', row=[216,247+30+90], column=[31, 41])
    Chat = [[cuminfection[t][j]*total_pop[j] for j in range(group)] for t in range(len(cuminfection))]
    # print(x1p)
    # print(Cnewp)
    xx1 = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 11039.864062, 0.0, 20000.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20000.0, 20000.0, 20000.0, 8960.135938, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9598.321527, 10401.678473, 9598.321527, 10401.678473, 5524.718847, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 14120.020346, 5879.979654, 14120.020346, 5879.979654, 9756.911792, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
    xx2 = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4073.60268, 14475.281153, 20000.0, 20000.0, 20000.0, 20000.0, 20000.0, 0.0, 0.0, 0.0, 0.0, 20000.0, 0.0, 0.0, 0.0, 0.0],
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4363.108554, 10243.088208, 20000.0, 20000.0, 10401.678473, 9598.321527, 10401.678473, 9598.321527, 10401.678473, 5524.718847, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
[20000.0, 20000.0, 5879.979654, 14120.020346, 5879.979654, 14120.020346, 5879.979654, 9756.911792, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
    # x1 = [[xx1[j][t-30] if t >= 30 and t<60 else sum([x1p[t+tau][j] for tau in range(7)]) / 7 for j in range(group)] for t in range(period)]
    # x2 = [[xx2[j][t-30] if t >= 30 and t<60 else sum([x2p[t + tau][j] for tau in range(7)]) / 7 for j in range(group)] for t in range(period)]
    x1 = [[sum([x1p[t+tau][j] for tau in range(7)]) / 7 for j in range(group)] for t in range(period)]
    x2 = [[sum([x2p[t + tau][j] for tau in range(7)]) / 7 for j in range(group)] for t in range(period)]
    Cnew = [[sum([Cnewp[t + tau][j] for tau in range(7)]) / 7 for j in range(group)] for t in range(period)]

    LOS = [[1.9544, 2.3698, 2.9677, 4.1774, 4.1405, 4.9579, 5.6474, 5.9455, 6.5570, 6.9782],
           [1.6000, 2.3333, 2.8619, 3.3555, 3.3958, 4.7778, 5.5690, 6.2801, 6.5108, 6.5317],
           [1.6000, 3.5000, 3.0460, 3.2578, 3.3589, 3.8148, 4.3255, 5.0031, 5.2380, 5.5851]]
    lambda1 = [[1/LOS[i][j] for j in range(len(LOS[i]))] for i in range(len(LOS))]
    # alpha_base = 0.3
    # tau_base = 2
    # beta_base = 2
    # w_base = 0.5
    policy_nonvacc = [53.7 if t < 18 else 56.48 for t in range(period)]
    policy_vacc = [41.67 if t < 18 or t > 42 else 44.44 for t in range(period)]
    policy1 = [[53.7/policy_nonvacc[t] for t in range(period)],
              [41.67/policy_vacc[t] for t in range(period)],
              [41.67/policy_vacc[t] for t in range(period)]]
    # social = [1 if t < 15 or t > 45 else 1.5 for t in range(period)]
    social = [[1 + 0.5 * exp(-(t - 32) ** 2 / (2 * 10 ** 2)) for t in range(period)],
              [1 + 0.5 * exp(-(t - 32) ** 2 / (2 * 10 ** 2)) for t in range(period)],
              [1 + 0.5 * exp(-(t - 32) ** 2 / (2 * 10 ** 2)) for t in range(period)],
              [1 + 0.5 * exp(-(t - 32) ** 2 / (2 * 10 ** 2)) for t in range(period)],
              [1 + 0.5 * exp(-(t - 32) ** 2 / (2 * 10 ** 2)) for t in range(period)],
              [1 + 0.5 * exp(-(t - 32) ** 2 / (2 * 10 ** 2)) for t in range(period)],
              [1 + 0.5 * exp(-(t - 32) ** 2 / (2 * 10 ** 2)) for t in range(period)],
              [1 + 0.5 * exp(-(t - 32) ** 2 / (2 * 10 ** 2)) for t in range(period)],
              [1 + 0.5 * exp(-(t - 32) ** 2 / (2 * 10 ** 2)) for t in range(period)],
              [1 + 0.5 * exp(-(t - 32) ** 2 / (2 * 10 ** 2)) for t in range(period)]]
    print(social)
    policy = [[[policy1[0][t] * social[j][t] for t in range(period)] for j in range(group)]]\
              + [[[policy1[1][t] * social[j][t] for t in range(period)] for j in range(group)]]\
              + [[[policy1[2][t] * social[j][t] for t in range(period)] for j in range(group)]]
    # step = 10
    # alpha1 = [round(alpha_base * 0.5 + alpha_base * i / (step-1), 4) for i in range(step)]
    # tau1 = [round(tau_base / (1 + i/2),4) for i in range(step)]
    # beta1 = [round(beta_base * 0.5 + beta_base * i / (step-1), 4) for i in range(step)]
    # w1 = [round(w_base * 0.5 + w_base * i / (step-1), 4) for i in range(step)]
    # c = itertools.product(alpha1, tau1, beta1, beta1, beta1, w1, w1, w1, w1, w1, w1, w1, w1, w1, w1)

    prop = [[0.973724571,0.026275429,0],
            [0.596813125,0.403114823,7.20515E-05],
            [0.070297175,0.568649197,0.361053628],
            [0.081802152,0.30782985,0.610367998],
            [0.098757311,0.317268847,0.583973842],
            [0.080111849,0.257662487,0.662225664],
            [0.0643932,0.233770645,0.701836155],
            [0.060599159,0.171857068,0.767543774],
            [0.070802964,0.178169344,0.751027692],
            [0.115432656,0.25122073,0.633346614]]
    death = [[0.0000158, 0.00000 , 0.00000],
                                [0, 0.00000 , 0.00000],
                                [0, 0.00000 , 0.00000],
                                [0, 0.00000 , 0.00000],
                                [0, 0.00000 , 0.00001],
                                [0.0008306, 0.00013 , 0.00003],
                                [0.000760173, 0.00090 , 0.00003 ],
                                [0.002559568, 0.00307 , 0.00027 ],
                                [0.007675998, 0.00689 , 0.00090 ],
                                [0.0478664, 0.02508 , 0.00460 ]]
    p = [[0.033241357, 0.00168, 0.00000],
                                [0.00853711, 0.00520, 0.22222],
                                [0.013102269, 0.01576, 0.00215],
                                [0.004261432, 0.00319, 0.00213],
                                [0.005579845, 0.00438, 0.00351],
                                [0.007184482, 0.00583, 0.00340],
                                [0.02191819, 0.01557, 0.00787],
                                [0.052791546, 0.03678, 0.01996],
                                [0.139227102, 0.10054, 0.05496],
                                [0.424860886, 0.33714, 0.17464]]
    
    vacprop0 = [0, 0.000869474, 0.874585414, 0.839532864, 0.827357426,
                0.861650509, 0.874687266, 0.864749337, 0.843740029, 0.708822072]
    vacprop3 = [0, 0, 0.052956183, 0.177700151, 0.207553989, 0.426879033,
                0.558551108, 0.652059978, 0.634747247, 0.445178699]
    vacprop2 = [vacprop0[i] - vacprop3[i] for i in range(len(vacprop0))]
    vacprop1 = [1 - vacprop0[i] for i in range(len(vacprop0))]

    # label = 1
    # id = 0
    # alpha_lower = [0.15]
    # alpha_upper = [0.175]
    #
    #
    # # tau_lower = [round(tau_base / 7, 4)]
    # # tau_upper = [round(tau_base / 1, 4)]
    # beta_lower = [1.5, 1.3, 1]
    # beta_upper = [1.62, 1.42, 1.2]
    #
    #
    # w_lower = [0.25 for i in range(group)]
    # w_upper = [0.4 for i in range(group)]
    #
    # step = 1000
    # alpharan = [random.uniform(alpha_lower[0] + (alpha_upper[0] - alpha_lower[0]) / step * j,
    #                            alpha_lower[0] + (alpha_upper[0] - alpha_lower[0]) / step * (j + 1)) for j in range(step)]
    # random.shuffle(alpharan)
    #
    # # tauran = [random.uniform(tau_lower[0] + (tau_upper[0] - tau_lower[0]) / step * j,
    # #                            tau_lower[0] + (tau_upper[0] - tau_lower[0]) / step * (j + 1)) for j in range(step)]
    # # random.shuffle(tauran)
    # tauran = [0.8 for j in range(step)]
    # betaran = [[] for i in range(3)]
    # for i in range(3):
    #     betaran[i] = [random.uniform(beta_lower[i] + (beta_upper[i] - beta_lower[i]) / step * j,
    #                                  beta_lower[i] + (beta_upper[i] - beta_lower[i]) / step * (j+1)) for j in range(step)]
    #     random.shuffle(betaran[i])
    #
    # wran = [[] for i in range(group)]
    # for i in range(group):
    #     wran[i] = [0.7*random.uniform(w_lower[i] + (w_upper[i] - w_lower[i]) / step * j,
    #                                  w_lower[i] + (w_upper[i] - w_lower[i]) / step * (j+1)) for j in range(step)]
    #     random.shuffle(wran[i])
    # temp = [[betaran[i][j] for i in range(3)] for j in range(step)]
    # for j in range(step):
    #     temp[j].sort(reverse=True)
    # cc = [(alpharan[j], tauran[j], temp[j], [wran[i][j] for i in range(group)]) for j in range(step)]

    pp = 0
    currentobj = []

    aa = [1.1, 1.21, 1.05, 0.9, 1.1,
          1.05, 1.05, 1.1, 1.12, 1.1]
    # ww = [0.26969523185668265, 0.261818774553978, 0.24929235891468596, 0.2686408089510693, 0.1808214804854385,
    #        0.17778078625394927, 0.17879090776911108, 0.18018611154356318, 0.18646473940562833, 0.18762363041635252]
    ww = [[0.198408689,0.095931398,0.075377296,0.142649748,0.253510478,0.108750269,0.073748669,0.035325003,0.012413567,0.003884884],
            [0.042348141,0.493675788,0.109501385,0.04957647,0.144000922,0.098354883,0.035741385,0.019017155,0.005522952,0.002260919],
            [0.008397092,0.058551198,0.630588168,0.07494171,0.076106627,0.10390807,0.033625307,0.008642874,0.003740834,0.001498121],
            [0.021373504,0.015515033,0.120991862,0.387622447,0.195210558,0.147935256,0.091932507,0.013800924,0.003935628,0.001682282],
            [0.039496562,0.06440112,0.091816649,0.160731038,0.325826066,0.187764259,0.095721487,0.0258699,0.006377763,0.001995156],
            [0.019274355,0.046058771,0.156463795,0.141086954,0.213844453,0.276874658,0.11243945,0.022746854,0.008147558,0.003063151],
            [0.02611464,0.054701314,0.140714499,0.153019644,0.173112167,0.194484994,0.201579573,0.041879498,0.010127788,0.004265883],
            [0.046704657,0.060300482,0.093816922,0.111580447,0.186902292,0.14634986,0.144322831,0.170207258,0.031640805,0.008174445],
            [0.04306802,0.081068029,0.139702188,0.063191141,0.127455505,0.156192544,0.117513345,0.117726399,0.112221974,0.041860854],
            [0.046512359,0.068146347,0.17516909,0.061971207,0.120615413,0.15547829,0.129778183,0.080576935,0.099652382,0.062099794]
            ]
    alpha = 0.33
    tau = 0.8
    beta = [1.580550580804279*0.68, 1.3003004259503794*0.68, 1.0232515583165258*0.68]
    cc = [[alpha, tau, beta, ww, policy]]
    for (alpha, tau, beta, w, policy) in cc:

        I0 = [Cnew[0][i] * prop[i][0] / tau *1 if i not in [1,2,3,4,5] else Cnew[0][i] * prop[i][0] / tau *6.5 for i in range(group)]
        VI0 = [Cnew[0][i] * prop[i][1] / tau*1 if i not in [1,2,3,4,5] else Cnew[0][i] * prop[i][1] / tau *6.5 for i in range(group)]
        BI0 = [Cnew[0][i] * prop[i][2] / tau*1 if i not in [1,2,3,4,5] else Cnew[0][i] * prop[i][2] / tau *6.5 for i in range(group)]

        E0 = [I0[i] / alpha * (1-tau) *2 for i in range(group)]
        VE0 = [VI0[i] / alpha * (1-tau) *2 for i in range(group)]
        BE0 = [BI0[i] / alpha * (1-tau) *2 for i in range(group)]

        C = [Cnew[0][i] * prop[i][0] * p[i][0] * (1-lambda1[0][i]) / lambda1[0][i]
              + Cnew[0][i] * prop[i][0] * (1-p[i][0]) * (1-1/5)/(1/5) for i in range(group)]
        VC = [Cnew[0][i] * prop[i][1] * p[i][1] * (1-lambda1[1][i]) / lambda1[1][i]
              + Cnew[0][i] * prop[i][1] * (1-p[i][1]) * (1-1/5)/(1/5) for i in range(group)]
        BC = [Cnew[0][i] * prop[i][0] * p[i][2] * (1-lambda1[2][i]) / lambda1[2][i]
              + Cnew[0][i] * prop[i][0] * (1-p[i][2]) * (1-1/5)/(1/5) for i in range(group)]

        C0 = [C[i] / (C[i] + VC[i] + BC[i]) * Chat[0][i] for i in range(group)]
        VC0 = [VC[i] / (C[i] + VC[i] + BC[i]) * Chat[0][i] for i in range(group)]
        BC0 = [BC[i] / (C[i] + VC[i] + BC[i]) * Chat[0][i] for i in range(group)]

        # i = 0
        # print(cuminfection[i], total_pop[i], I0[i], E0[i])
        # sa = [0.5, 0.65, 0.5, 0.4, 0.65, 0.5, 0.4, 0.4, 0.4, 0.4]
        # va = [0.55, 0.55, 0.55, 0.55, 0.65, 0.55, 0.55, 0.5, 0.5, 0.5]
        # ba = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7, 0.75, 0.75]
        S0 = [max(0, vacprop1[i] * total_pop[i] - I0[i] - E0[i] - cuminfection[0][i]*total_pop[i] / 2 ) * 1
               for i in range(group)]
        VS0 = [max(0, vacprop2[i] * total_pop[i] - VI0[i] - VE0[i] - cuminfection[0][i]*total_pop[i] / 3 ) * 1
               for i in range(group)]
        BS0 = [max(0, vacprop3[i] * total_pop[i] - BI0[i] - BE0[i] - cuminfection[0][i]*total_pop[i] / 6 ) * 1
               for i in range(group)]
        print([int(S0[i]) for i in range(group)])
        print([int(VS0[i]) for i in range(group)])
        print([int(BS0[i]) for i in range(group)])
        print([int(E0[i]) for i in range(group)])
        print([int(VE0[i]) for i in range(group)])
        print([int(BE0[i]) for i in range(group)])
        print([int(I0[i]) for i in range(group)])
        print([int(VI0[i]) for i in range(group)])
        print([int(BI0[i]) for i in range(group)])
        print('=====')
        print(sum(S0) + sum(VS0) + sum(BS0))
        print(sum(E0) + sum(VE0) + sum(BE0))
        print(sum(I0) + sum(VI0) + sum(BI0))
        print(sum(C0) + sum(VC0) + sum(BC0))
        print(sum([cuminfection[0][i]*total_pop[i] for i in range(group)]))

        population = [S0, E0, I0, C0, VS0, VE0, VI0, VC0, BS0, BE0, BI0, BC0, [0 for i in range(group)], [0 for i in range(group)], [0 for i in range(group)]]
        obj, Xsim = model(population, alpha, tau, beta, w, policy, x1, x2, Chat)
        currentobj.append([obj, [alpha, tau, beta, w]])
        pp += 1
        if obj < bestobj:
            bestobj = obj
            bestpopulation = population
            bestalpha = alpha
            besttau = tau
            bestbeta = beta
            bestw = w
            beatX = Xsim
            print('============================')
            print('OBJ =', bestobj)
            print('Population = ', bestpopulation)
            print('alpha = ', bestalpha)
            print('tau = ', besttau)
            print('beta = ', bestbeta)
            print('w = ', bestw)
            Mbar = Xsim[compartment.index('M')]
            VMbar = Xsim[compartment.index('VM')]
            BMbar = Xsim[compartment.index('BM')]
            Cbar = Xsim[compartment.index('C')]
            VCbar = Xsim[compartment.index('VC')]
            BCbar = Xsim[compartment.index('BC')]
            for j in range(group):
                print([int(Cbar[j][t] + VCbar[j][t] + BCbar[j][t]) for t in range(period)])
                print([int(Cnew[t][j]) for t in range(len(Cnew))])
                print()
            print()
        else:
            print('========', pp - 1, currentobj[-1][0])
        #     parameter[id] += (parameter_upper[id] - parameter_lower[id]) / 20
        # else:
        #
        #     print(pp, obj)
        #     if parameter[id] >= parameter_upper[id]:
        #         parameter = bestbeta + bestw + [bestalpha] + [besttau]
        #         id += 1
        #     else:
        #         if random.random() < 0.5:
        #
        #             id += 1
        #         else:
        #             parameter[id] += (parameter_upper[id] - parameter_lower[id]) / 20
        #     print(id, parameter)
        # if id >= len(parameter):
        #     break
    currentobj.sort()
    print('--------------------')
    for i in range(min(10,len(currentobj))):
        print(currentobj[i])
    print('--------------------')
    for j in range(group):
        print([int(Cbar[j][t] + VCbar[j][t] + BCbar[j][t] + Mbar[j][t] + VMbar[j][t] + BMbar[j][t]) for t in range(period)])
    print()
    for j in range(group):
        print([int(Cnew[t][j]) for t in range(period)])
    print()

    print('OBJ =', bestobj)
    for j in range(group):
        print([int(Cbar[j][30]-Cbar[j][0]), int(VCbar[j][30]-VCbar[j][0]),  int(BCbar[j][30]-BCbar[j][0])])
    print()
    t = 30
    for j in range(group):
        print([int(Cbar[j][t] + Mbar[j][t])-int(Cbar[j][0] + Mbar[j][0]),
               int(VCbar[j][t] + VMbar[j][t])-int(VCbar[j][0] + VMbar[j][0]),
               int(BCbar[j][t] + BMbar[j][t])-int(BCbar[j][0] + BMbar[j][0])])
    print('I')
    print([int(Xsim[compartment.index('I')][j][t])for j in range(group)])
    print([int(Xsim[compartment.index('VI')][j][t])for j in range(group)])
    print([int(Xsim[compartment.index('BI')][j][t])for j in range(group)])
    print('S')
    print([int(Xsim[compartment.index('S')][j][t])for j in range(group)])
    print([int(Xsim[compartment.index('VS')][j][t])for j in range(group)])
    print([int(Xsim[compartment.index('BS')][j][t])for j in range(group)])
    print('E')
    print([int(Xsim[compartment.index('E')][j][t])for j in range(group)])
    print([int(Xsim[compartment.index('VE')][j][t])for j in range(group)])
    print([int(Xsim[compartment.index('BE')][j][t])for j in range(group)])
    print()
    for j in range(group):
        print([Cbar[j][t] + VCbar[j][t] + BCbar[j][t] for t in range(period)])

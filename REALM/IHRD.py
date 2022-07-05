from __future__ import print_function
import sys

import cplex
from cplex.exceptions import CplexError
import random
import time
from math import log, exp, fabs
from scipy import stats
from numpy import array, std, var
import numpy
# import Datasetnew
from itertools import product
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import REALM
import DatasetMD

data = DatasetMD.Coviddata_MD()
period = 30
# period_num = 4
# period_len = 7
totalperiod = 100
group = data.group
# Various_lambdah = [2, 1.5, 1, 0.9, 0.7, 0.6, 0.5, 0.4]
# Various_lambdaicu = [2, 1.5, 1, 0.9, 0.7, 0.6, 0.5, 0.4]
Various_lambdah = [1 for j in range(group)]
Various_lambdaicu = [1 for j in range(group)]
lam = [1 / (1 + data.parameter['lambda_a'][0] * exp(data.parameter['lambda_b'][0] * (t+100-data.parameter['t0'][0])))
       for t in range(period)]

w = [data.parameter['w'][j] * data.parameter['policy84_99'][0] for j in range(group)]
betaD = data.parameter['betaD'][0]
betaO = data.parameter['betaD'][0] * data.parameter['betaO'][0]

rho = data.parameter['a64_99']
bhat = data.parameter['b64_99']

alphaD = data.parameter['alphaD'][0]
alphaO = data.parameter['alphaD'][0] * data.parameter['alphaO'][0]

gammaD = data.parameter['gammaD'][0]
gammaO = data.parameter['gammaD'][0] * data.parameter['gammaO'][0]

piD = data.parameter['piD']
piO = [data.parameter['piD'][j] * data.parameter['piO'][j] for j in range(group)]

dqD = [data.parameter['dqD'][j] for j in range(group)]
dqO = [data.parameter['dqD'][j] * data.parameter['dqO'][j] for j in range(group)]
dicuD = [data.parameter['dicuD'][j] for j in range(group)]
dicuO = [data.parameter['dicuD'][j] * data.parameter['dicuO'][j] for j in range(group)]
dhD = [data.parameter['dhD'][j] for j in range(group)]
dhO = [data.parameter['dhD'][j] * data.parameter['dhO'][j] for j in range(group)]
lambdahD = [data.parameter['lambdahD'][0] * Various_lambdah[j] for j in range(group)]
lambdaicuD = [data.parameter['lambdaicuD'][0] * Various_lambdaicu[j] for j in range(group)]
lambdahO = [data.parameter['lambdahD'][0] * data.parameter['lambdahO'][0]
            * Various_lambdah[j] for j in range(group)]
lambdaicuO = [data.parameter['lambdaicuD'][0] * data.parameter['lambdaicuO'][0]
              * Various_lambdaicu[j] for j in range(group)]
# pm = [[(1-lam[t]) * (1-piD[j]) + lam[t] * (1-piO[j]) for t in range(period)] for j in range(group-4, group-1)]
# pq = [[(1-lam[t]) * piD[j] + lam[t] * piO[j] for t in range(period)] for j in range(group-4, group-1)]

beta = [(1 - lam[t]) * betaD + lam[t] * betaO for t in range(period)]
q_EtoI = [[(1 - lam[t]) * alphaD + lam[t] * alphaO for t in range(period)] for j in range(group-4, group-1)]
q_ItoR = [[(1-lam[t]) * gammaD * (1-piD[j]) + lam[t] * gammaO * (1-piO[j])
           for t in range(period)] for j in range(group-4, group-1)]
q_ItoD = [[(1-lam[t]) * dqD[j] * piD[j] + lam[t] * dqO[j] * piO[j] for t in range(period)] for j in range(group-4, group-1)]
q_MtoR = [[(1-lam[t]) * gammaD + lam[t] * gammaO for t in range(period)] for j in range(group-4, group-1)]
q_QtoD = [[(1-lam[t]) * dqD[j] + lam[t] * dqO[j]
           for t in range(period)] for j in range(group-4, group-1)]
q_HtoD = [[(1-lam[t]) * dhD[j] * lambdahD[j] + lam[t] * dhO[j] * lambdahO[j]
           for t in range(period)] for j in range(group-4, group-1)]
q_HtoR = [[(1-lam[t]) * (1-dhD[j]) * lambdahD[j] + lam[t] * (1-dhO[j]) * lambdahO[j]
           for t in range(period)] for j in range(group-4, group-1)]
q_UtoD = [[(1-lam[t]) * dicuD[j] * lambdaicuD[j] + lam[t] * dicuO[j] * lambdaicuO[j]
           for t in range(period)] for j in range(group-4, group-1)]
q_UtoR = [[(1-lam[t]) * (1-dicuD[j]) * lambdaicuD[j] + lam[t] * (1-dicuO[j]) * lambdaicuO[j]
           for t in range(period)] for j in range(group-4, group-1)]


population = [[10000000 for j in range(3)]]
population.extend([[data.compartment_population[n][j] for j in range(group - 4, group - 1)] for n in range(4, 9)])
total_pop = [data.group_population for j in range(group - 4, group - 1)]
S = 6

ch = data.ch
cu = data.cu
cx = data.cx
cy = data.cy
cd = data.cd
ci = data.ci
ct = data.ct
cts = data.cts
group = 3
Qdata = [[25.77390382,90.26724523,140.3707237,190.2348299,240.4576826,291.5636321,343.8679376,397.622341,453.0489892,510.3547757,569.7354735,631.3780617,695.4622599,762.1616337,831.6444112,904.0740515,979.609576,1058.405667,1140.542348,1225.402826,1313.147707,1403.848118,1497.553782,1594.300475,1694.111549,1796.997527,1902.955424,2011.967988,2124.002967,2239.012406],
[29.76599124,43.43546201,56.44566214,69.52845559,82.41593275,95.09976973,107.6448553,120.1382058,132.6663719,145.3082102,158.133145,171.2016363,184.5662564,198.2727902,212.3611807,226.8662878,241.8184744,257.2440457,273.1655647,289.6020682,306.5692,324.0793634,342.1425882,360.765774,379.9525566,399.7031828,420.0143987,440.8793313,462.2873683,484.2240389],
[34.15389879,50.10179382,64.63023214,78.74902846,92.27944039,105.2759354,117.8539514,130.1451814,142.2762144,154.361279,166.5004478,178.7803396,191.2756516,204.0508056,217.1614364,230.6556422,244.5749927,258.955319,273.8273158,289.2169855,305.1459505,321.6317583,338.6890213,356.3286545,374.5577312,393.3793452,412.7924829,432.7918773,453.3678488,474.5061345]
]
def model(totalcost):
    '''=============== Step 1. Define model and model size ====================='''
    user = REALM.REALM()

    # user.set_types(types = types) # Expectation or Robustness
    user.set_group(group=group)
    user.set_time(time=period)

    H0 = 1000
    ICU0 = 150

    '''=============== Step 2. Define name of compartment and population ====================='''
    compartment = ['I','Q','H','U','R','D']
    compart_num = len(compartment)
    user.set_all_compartment(name=compartment)
    # print(population)
    for i in range(compart_num):
        print(population[i])
    for i in range(compart_num):
        user.set_population(compartment=compartment[i], population=population[i])
    print('step.2')
    '''=============== Step 3. Define decision variable  ====================='''
    user.set_flow_variable(n='Q', m = 'H', xupper=[[H0 for t in range(period)] for j in range(group)],
                              xlower=[[0 for t in range(period)] for j in range(group)])
    user.set_flow_variable(n='Q', m = 'U', xupper=[[ICU0 for t in range(period)] for j in range(group)],
                              xlower=[[0 for t in range(period)] for j in range(group)])
    user.set_flow_variable(n='I', m = 'Q', xupper=[[Qdata[j][t+1]-Qdata[j][t] if t <period-1 else 0
                                                    for t in range(period)] for j in range(group)],
                              xlower=[[Qdata[j][t+1]-Qdata[j][t] if t <period-1 else 0
                                                    for t in range(period)] for j in range(group)])

    user.set_custom_variable(name='Hmax', xlower=0, xupper=0)
    user.set_custom_variable(name='Umax', xlower=0, xupper=0)

    print('step.3')
    '''====== Step 4.1 Define decision-independent transition between compartment and compartment  ========='''
    user.set_transition_compartment(n='Q', m='D', prob=[[q_QtoD[j][t] for t in range(period)] for j in range(group)], opt = None)
    user.set_transition_compartment(n='H', m='D', prob=[[q_HtoD[j][t] for t in range(period)] for j in range(group)], opt = None)
    user.set_transition_compartment(n='U', m='D', prob=[[q_UtoD[j][t] for t in range(period)] for j in range(group)], opt = None)
    user.set_transition_compartment(n='H', m='R', prob=[[q_HtoR[j][t] for t in range(period)] for j in range(group)], opt = None)
    user.set_transition_compartment(n='U', m='R', prob=[[q_UtoR[j][t] for t in range(period)] for j in range(group)], opt = None)

    print('step.4.2')

    # '''=============== Step 5. Prediction  ====================='''
    # print('step.5p Predicition')
    # user.set_transition_compartment_self()
    # for i in user.compartment_name:
    #     user.set_transition_group(compartment=i, prob=[[[1 if k == j else 0 for t in range(period)]
    #                                                     for k in range(group)] for j in range(group)])
    # B=100
    # user.x = [[[[B if (n,m) == (0,1) else 0 for t in range(period)]
    #                     for j in range(group)] for m in range(compart_num)] for n in range(compart_num)]
    # X = user.DeterPrediction()
    #
    # # X1 = user.StochaPrediction()
    # # for i in range(compart_num):
    # #     print(compartment[i])
    # #     for j in range(group):
    # #         print(X[i][j])
    #         # print(X1[i][j])
    #     # print()
    #     # print([sum([X[i][j][t] for j in range(group)]) for t in range(period)])
    #     # print([sum([X1[i][j][t] for j in range(group)]) for t in range(period)])
    # # print()
    #
    # user.x = [[[[0 if (n,m) == (0,1) else 0 for t in range(period)]
    #                     for j in range(group)] for m in range(compart_num)] for n in range(compart_num)]
    # X1 = user.DeterPrediction()
    #
    # x = [20, 0, 30, 50]
    # user.x = [[[[x[j] if (n,m) == (0,1) else 0 for t in range(period)]
    #                     for j in range(group)] for m in range(compart_num)] for n in range(compart_num)]
    # X3 = user.DeterPrediction()
    # X2,p,q = user.DeterPrediction(types='step')
    #
    # for i in range(compart_num):
    #     print(compartment[i])
    #     for j in range(group):
    #         print([X3[i][j][t] for t in range(period)])
    #         print([X2[i][j][t] for t in range(period)])
    #         print()

    '''=============== Step 5. Define other constraints  ====================='''
    print('step.5')
    # capacity constraints
    for t in range(1,period):
        user.custom_lp(fvar=[('H', j, t) for j in range(group)],
                       fcoef=[1 for j in range(group)],
                       dvar=['Hmax'],
                       dcoef=[-1],
                       sense='L',
                       rhs=H0,
                       name='Hospitalization.' + str(t))
        user.custom_lp(fvar=[('U', j, t) for j in range(group)],
                       fcoef=[1 for j in range(group)],
                       dvar=['Umax'],
                       dcoef=[-1],
                       sense='L',
                       rhs=ICU0,
                       name='ICU.' + str(t))


    
    #########################################
    # for t in range(period):
    #     user.custom_lp(dvar= ['u.'+str(j)+'.'+str(s) for j in range(group) for s in range(S)],
    #                    dcoef= [Test_capacity[t] * Test_probability[s] for j in range(group) for s in range(S)],
    #                    sense='L',
    #                    rhs=Test_capacity[t],
    #                    name='Test_capacity.' + str(t))
    #########################################
    user.custom_lp(fvar=[('Q', 'H', j, t) for j in range(group) for t in range(period-1)] \
                    + [('Q', 'U', j, t) for j in range(group) for t in range(period-1)],
                   fcoef=[cx/Various_lambdah[j] for j in range(group) for t in range(period-1)]
                         + [cy/Various_lambdaicu[j] for j in range(group) for t in range(period-1)],
                   sense='L',
                   rhs=totalcost,
                   name='totalcost')
    #########################################
    # user.custom_lp(fvar=[('H', j, t) for j in range(group) for t in range(period-1)] \
    #                 + [('U', j, t) for j in range(group) for t in range(period-1)],
    #                fcoef=[cx * lambdahD[j] for j in range(group) for t in range(period-1)]
    #                      + [cy * lambdaicuD[j] for j in range(group) for t in range(period-1)],
    #                dvar=['u.' + str(j) + '.' + str(s) for j in range(group) for s in range(S)],
    #                dcoef=[ct * Test_capacity[t] * Test_probability[s] * (period-1) for j in range(group) for s in range(S)],
    #                sense='L',
    #                rhs=totalcost,
    #                name='totalcost')
    #########################################
    # user.custom_lp(dvar=['Hmax','Umax'],
    #                dcoef=[ch, cu],
    #                sense='L',
    #                rhs=totalcost,
    #                name='totalcost')
    #########################################
    # user.custom_lp(fvar=[('I', j, t) for j in range(group) for t in range(period)],
    #                fcoef=[ci for j in range(group) for t in range(period)],
    #                sense='L',
    #                rhs=totalcost[0],
    #                name='Infection')
    # user.custom_lp(fvar=[('D', j, period-1) for j in range(group)],
    #                fcoef=[ci for j in range(group)],
    #                sense='L',
    #                rhs=totalcost[1],
    #                name='Deaths')
    '''=============== Step 6. Define objective  ====================='''
    print('step.6')
    user.set_objectivetype(sense='min')
    # user.set_initial_solution(compartment1='S',compartment2='SS',val=
    # [[50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    # [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    # [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    # [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])

    for j in range(group):
        user.set_objective(state=('D', j, period - 1), value=cd)


    '''=============== Step 7. Solve and output  ====================='''
    user.set_solver(solver='gurobi')
    user.set_approximation(opt='SO')
    user.set_log_stream_SO(label=1, file='SEIHR/log_SEIHR_exp_'+str(totalcost))
    sol = user.solve(label='Expectation', ct=3600)
    (status, obj, Xopt, xopt, Xc) = sol
    print('obj=', obj)
    Xsolution_exp = user.get_solution_compartment(compartment=compartment)
    xsolution_exp = user.get_solution_flow_variable(compartmentfrom=['Q','Q'],compartmentto=['H','U'])
    Xcustom_exp = user.get_solution_custom_variable(name=['Hmax', 'Umax'])
    user.Solution_print(X=Xsolution_exp,x=xsolution_exp,xc=None)
    #
    # '''======================== Simulation =================================='''
    # sample = 1000
    # # randomq = user.Sample_generation(n=['S','SS'], m=['E','EE'], lb=[1,1], ub=[1.01,1.01], size=sample)
    # randomq = [[[0.7 + s * 0.6 / sample if (ii, i) in [(2,5),(3,5)] else 1 for i in range(compart_num)]
    #                 for ii in range(compart_num)] for s in range(sample)]
    # # C = [[sum([Test_capacity[t] * Test_probability[s] * Xcustom['u.' + str(j) + '.' + str(s)] for s in range(S)])
    # #       for t in range(period)] for j in range(group)]
    # # for j in range(group):
    # #     for t in range(period):
    # #         betaIH = [[1 if compartment[n] == 'I' and k == j else 0
    # #                    for k in range(group)] for n in range(len(population))]
    # #         user.set_transition_compartment_dp(flow=('I', 'M', j, t), alpha=C[j][t]*pm[j][t],
    # #                                            betadeno=betaIH, alphadeno=bhat[j], ub=attractive*pm[j][t])
    # #         user.set_transition_compartment_dp(flow=('I', 'Q', j, t), alpha=C[j][t]*pq[j][t],
    # #                                            betadeno=betaIH, alphadeno=bhat[j], ub=attractive*pq[j][t])
    # #         user.set_transition_compartment_dp(flow=('I', 'D', j, t), alpha=(1 - C[j][t] * pm[j][t] - C[j][t] * pq[j][t]) * q_ItoD[j][t],
    # #                                            betadeno=betaIH, alphadeno=bhat[j])
    # #         user.set_transition_compartment_dp(flow=('I', 'R', j, t), alpha=(1 - C[j][t] * pm[j][t] - C[j][t] * pq[j][t]) * q_ItoR[j][t],
    # #                                            betadeno=betaIH, alphadeno=bhat[j])
    # x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)] for n in range(compart_num)]
    # x[0][1] = [[xopt[0][1][j][t] for t in range(period)] for j in range(group)]
    # x[1][2] = [[xopt[1][2][j][t] for t in range(period)] for j in range(group)]
    # x[1][3] = [[xopt[1][3][j][t] for t in range(period)] for j in range(group)]
    # user.x = x
    # Xsim_exp = [0 for s in range(sample)]
    # for s in range(sample):
    #     Xsim_exp[s] = user.Prediction(opt=['D', None, randomq[s]])
    # # user.set_transition_compartment_dp(flow=('I', 'M'), opt='delete')
    # # user.set_transition_compartment_dp(flow=('I', 'Q'), opt='delete')
    # # user.set_transition_compartment_dp(flow=('I', 'R'), opt='delete')
    # # user.set_transition_compartment_dp(flow=('I', 'D'), opt='delete')
    # # user.set_transition_compartment_dp(flow=('I', 'I'), opt='delete')
    # # print(user.compartment_info[2]['tocpm'][7]['q_idp'])
    # # print(user.compartment_info[2]['tocpm'][8]['q_idp'])
    # # print([(Xsim_exp[0][2][0][t+1], Xsim_exp[0][2][0][t] * (1-q_ItoR[0][t]-q_ItoD[0][t]), Xsim_exp[0][1][0][t] * q_EtoI[0][t],
    # #         min(Xsim_exp[0][2][0][t] * 9000 / (Xsim_exp[0][2][0][t]+bhat[j]), attractive*Xsim_exp[0][2][0][t]))
    # #        for t in range(period-1)])
    # # print([sum([Xsim_exp[0][4][j][t] for j in range(group)]) for t in range(period)])
    # # print([sum([Xsim_exp[0][5][j][t] for j in range(group)]) for t in range(period)])
    # '''======================== save =================================='''
    # # lhs = [['totalcost',
    # #         [sum([xsolution[('I','M')][j][t] * cts + xsolution[('I','Q')][j][t] * cts \
    # #               + xsolution[('Q', 'H')][j][t] * cx + xsolution[('Q', 'U')][j][t] * cy
    # #              for j in range(group) for t in range(period)]) + Xcustom['Hmax'] * ch + Xcustom['Umax'] * cu ]]
    # #        ]
    # lhs = [['totalcost',
    #         [sum([xsolution_exp[('Q', 'H')][j][t] * cx + xsolution_exp[('Q', 'U')][j][t] * cy
    #               for j in range(group) for t in range(period)]) + Xcustom_exp['Hmax'] * ch + Xcustom_exp['Umax'] * cu ]
    #         ]]
    # for t in range(period):
    #     lhs.append(['Hospitalization.' + str(t),[sum([Xsim_exp[s][2][j][t] for j in range(group)])
    #                              for s in range(sample)]])
    # for t in range(period):
    #     lhs.append(['ICU.' + str(t), [sum([Xsim_exp[s][3][j][t] for j in range(group)])
    #                              for s in range(sample)]])
    # lhs.append(['objective', [sum([Xsim_exp[s][-1][j][t] * cd for j in range(group)])
    #                           for s in range(sample)]])
    # user.save_result(filename='SEIHR/resultIHR_exp_'+str(totalcost),
    #                  compartment=Xsim_exp, dvar=xopt, custom=Xcustom_exp, lhs=lhs)
    #
    # '''======================== Simulation =================================='''
    # # sample = 1000
    # # # randomq = user.Sample_generation(n=['S','SS'], m=['E','EE'], lb=[1,1], ub=[1.01,1.01], size=sample)
    # # randomq = [[[1 + s * 0.1 / sample if (ii, i) in [(2, 8),(4, 8),(5, 8),(6, 8)] else 1 for i in range(compart_num)]
    # #             for ii in range(compart_num)] for s in range(sample)]
    # # # C = [[sum([Test_capacity[t] * Test_probability[s] * Xcustom['u.' + str(j) + '.' + str(s)] for s in range(S)])
    # # #       for t in range(period)] for j in range(group)]
    # # # for j in range(group):
    # # #     for t in range(period):
    # # #         betaIH = [[1 if compartment[n] == 'I' and k == j else 0
    # # #                    for k in range(group)] for n in range(len(population))]
    # # #         user.set_transition_compartment_dp(flow=('I', 'M', j, t), alpha=C[j][t]*pm[j][t],
    # # #                                            betadeno=betaIH, alphadeno=bhat[j], ub=attractive*pm[j][t])
    # # #         user.set_transition_compartment_dp(flow=('I', 'Q', j, t), alpha=C[j][t]*pq[j][t],
    # # #                                            betadeno=betaIH, alphadeno=bhat[j], ub=attractive*pq[j][t])
    # # #         user.set_transition_compartment_dp(flow=('I', 'D', j, t), alpha=(1 - C[j][t] * pm[j][t] - C[j][t] * pq[j][t]) * q_ItoD[j][t],
    # # #                                            betadeno=betaIH, alphadeno=bhat[j])
    # # #         user.set_transition_compartment_dp(flow=('I', 'R', j, t), alpha=(1 - C[j][t] * pm[j][t] - C[j][t] * pq[j][t]) * q_ItoR[j][t],
    # # #                                            betadeno=betaIH, alphadeno=bhat[j])
    # # x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)] for n in range(compart_num)]
    # # x[2][3] = [[xopt[2][3][j][t] for t in range(period)] for j in range(group)]
    # # x[2][4] = [[xopt[2][4][j][t] for t in range(period)] for j in range(group)]
    # # x[4][5] = [[xopt[4][5][j][t] for t in range(period)] for j in range(group)]
    # # x[4][6] = [[xopt[4][6][j][t] for t in range(period)] for j in range(group)]
    # # user.x = x
    # # Xsim_exp = [0 for s in range(sample)]
    # # for s in range(sample):
    # #     Xsim_exp[s] = user.Prediction(opt=['D', None, randomq[s]])
    # # '''======================== save =================================='''
    # # # lhs = [['totalcost',
    # # #         [sum([xsolution[('I', 'M')][j][t] * cts + xsolution[('I', 'Q')][j][t] * cts \
    # # #               + xsolution[('Q', 'H')][j][t] * cx + xsolution[('Q', 'U')][j][t] * cy
    # # #               for j in range(group) for t in range(period)]) + Xcustom['Hmax'] * ch + Xcustom['Umax'] * cu]]
    # # #        ]
    # # lhs = [['totalcost',
    # #         [sum([xsolution_exp[('Q', 'H')][j][t] * cx + xsolution_exp[('Q', 'U')][j][t] * cy
    # #               for j in range(group) for t in range(period)]) + Xcustom_exp['Hmax'] * ch + Xcustom_exp['Umax'] * cu
    # #          + sum([Xcustom_exp['u.' + str(j) + '.' + str(s)] * ct * Test_capacity[t] * Test_probability[s] * (period-1)
    # #          for j in range(group) for s in range(S)]) ]]
    # #        ]
    # # for t in range(period):
    # #     lhs.append(['Hospitalization.' + str(t), [sum([Xsim_exp[s][5][j][t] for j in range(group)])
    # #                                               for s in range(sample)]])
    # # for t in range(period):
    # #     lhs.append(['ICU.' + str(t), [sum([Xsim_exp[s][6][j][t] for j in range(group)])
    # #                                   for s in range(sample)]])
    # # lhs.append(['objective', [sum([Xsim_exp[s][2][j][t] * ci + Xsim_exp[s][-1][j][t] * cd for j in range(group)])
    # #                           for s in range(sample)]])
    # # user.save_result(filename='SEIHR/resultSEIHR_death_exp_' + str(totalcost),
    # #                  compartment=Xsim_exp, dvar=xopt, custom=Xcustom_exp, lhs=lhs)







    '''====================== robust model =========================='''
    print('Begin to solve robust model')
    # theta = [1+0.01*(1+i) for i in range(20)]
    theta = [1.3]
    # Xoptr = [[] for i in range(len(theta))]
    # Xcr = [[] for i in range(len(theta))]
    # xoptr = [[] for i in range(len(theta))]
    # objr = [0 for i in range(len(theta))]
    # Xsim_ro = [[[] for s in range(sample)] for ii in range(len(theta))]
    #
    # # xlabel = [('I','M'),('I','Q'),('Q','H'),('Q','U')]
    # # for i in range(len(xlabel)):
    # #     user.set_initial_solution(compartment1=xlabel[i][0],compartment2=xlabel[i][1],val=xsolution[xlabel[i]])
    # for i in range(len(theta)):
    #     user.set_solver(solver='cplex')
    #     user.set_log_stream_SO(label=1, file='SEIHR/log_SEIHR_ro'+'_'+str(totalcost)+'_'+str(round(theta[i],2)))
    #     solr = user.default_lp(label='Robustness', ct=3600, target=obj*theta[i])
    #     if solr != 0:
    #         (status, objr[i], Xoptr[i], xoptr[i], Xcr[i]) = solr
    #         Xsolution = user.get_solution_compartment(compartment=compartment)
    #         xsolution = user.get_solution_flow_variable(compartmentfrom=['I','I','Q','Q'],compartmentto=['M','Q','H','U'])
    #         Xcustom = user.get_solution_custom_variable(name=['Hmax', 'Umax']
    #                                                       + ['u.' + str(j) + '.' + str(s) for s in range(S) for j in range(group)]
    #                                                       + ['phi1.' + str(j) + '.' + str(t) + '.'+str(s)
    #                                                          for s in range(S) for j in range(group) for t in range(period)]
    #                                                       + ['phi2.' + str(j) + '.' + str(t) + '.'+str(s)
    #                                                          for s in range(S) for j in range(group) for t in range(period)]
    #                                                       + ['z.' + str(j) + '.' + str(t) + '.' + str(s)
    #                                                          for s in range(S) for j in range(group) for t in range(period)])
    #         user.Solution_print(X=Xsolution,x=xsolution,xc=None)
    #     '''=============== Step 8. Simulation  ====================='''
    #     # user.set_time(time=totalperiod)
    #     # C = [[sum([Test_capacity[t] * Test_probability[s] * Xcustom['u.' + str(j) + '.' + str(s)] for s in range(S)])
    #     #       for t in range(period)] for j in range(group)]
    #     # for j in range(group):
    #     #     for t in range(period):
    #     #         betaIH = [[1 if compartment[n] == 'I' and k == j else 0
    #     #                    for k in range(group)] for n in range(len(population))]
    #     #         user.set_transition_compartment_dp(flow=('I', 'M', j, t), alpha=C[j][t] * pm[j][t],
    #     #                                            betadeno=betaIH, alphadeno=bhat[j], ub=attractive * pm[j][t])
    #     #         user.set_transition_compartment_dp(flow=('I', 'Q', j, t), alpha=C[j][t] * pq[j][t],
    #     #                                            betadeno=betaIH, alphadeno=bhat[j], ub=attractive * pq[j][t])
    #     #         user.set_transition_compartment_dp(flow=('I', 'D', j, t), alpha=(1 - C[j][t] * pm[j][t] - C[j][t] * pq[j][t]) * q_ItoD[j][t],
    #     #                                            betadeno=betaIH, alphadeno=bhat[j])
    #     #         user.set_transition_compartment_dp(flow=('I', 'R', j, t), alpha=(1 - C[j][t] * pm[j][t] - C[j][t] * pq[j][t]) * q_ItoR[j][t],
    #     #                                            betadeno=betaIH, alphadeno=bhat[j])
    #     #
    #     randomq = [[[1 + s * 0.1 / sample if (ii, i) in [(0, 1)] else 1 for i in range(compart_num)]
    #                 for ii in range(compart_num)] for s in range(sample)]
    #     x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)]
    #          for n in range(compart_num)]
    #     x[2][3] = [[xoptr[i][2][3][j][t] for t in range(period)] for j in range(group)]
    #     x[2][4] = [[xoptr[i][2][4][j][t] for t in range(period)] for j in range(group)]
    #     x[4][5] = [[xoptr[i][4][5][j][t] for t in range(period)] for j in range(group)]
    #     x[4][6] = [[xoptr[i][4][6][j][t] for t in range(period)] for j in range(group)]
    #     user.x = x
    #     print(user.x)
    #     Xsim_ro[i] = [0 for s in range(sample)]
    #     for s in range(sample):
    #         Xsim_ro[i][s] = user.Prediction(opt=['D', None, randomq[s]])
    #     # user.set_transition_compartment_dp(flow=('I', 'M'), opt='delete')
    #     # user.set_transition_compartment_dp(flow=('I', 'Q'), opt='delete')
    #     # user.set_transition_compartment_dp(flow=('I', 'R'), opt='delete')
    #     # user.set_transition_compartment_dp(flow=('I', 'D'), opt='delete')
    #     # user.set_transition_compartment_dp(flow=('I', 'I'), opt='delete')
    #
    #     '''======================== save =================================='''
    #     # lhs = [['totalcost',
    #     #         [sum([xsolution[('I', 'M')][j][t] * cts + xsolution[('I', 'Q')][j][t] * cts \
    #     #               + xsolution[('Q', 'H')][j][t] * cx + xsolution[('Q', 'U')][j][t] * cy
    #     #               for j in range(group) for t in range(period)]) + Xcustom['Hmax'] * ch + Xcustom['Umax'] * cu]]
    #     #        ]
    #     lhs = [['totalcost',
    #             [sum([xsolution[('Q', 'H')][j][t] * cx + xsolution[('Q', 'U')][j][t] * cy
    #                   for j in range(group) for t in range(period)]) + Xcustom['Hmax'] * ch + Xcustom['Umax'] * cu
    #              + sum(
    #                 [Xcustom['u.' + str(j) + '.' + str(s)] * ct * Test_capacity[t] * Test_probability[s] * (period - 1)
    #                  for j in range(group) for s in range(S)])]]
    #            ]
    #     for t in range(period):
    #         lhs.append(['Hospitalization.' + str(t), [sum([Xsim_ro[i][s][5][j][t] for j in range(group)])
    #                                                   for s in range(sample)]])
    #     for t in range(period):
    #         lhs.append(['ICU.' + str(t), [sum([Xsim_ro[i][s][6][j][t] for j in range(group)])
    #                                       for s in range(sample)]])
    #     lhs.append(['objective', [sum([Xsim_ro[i][s][2][j][t] * ci + Xsim_ro[i][s][-1][j][t] * cd
    #                                    for j in range(group)]) for s in range(sample)]])
    #
    #     user.save_result(filename='SEIHR/resultSEIHR_ro'+'_'+str(totalcost)+'_'+str(round(theta[i],2)),
    #                      compartment=Xsim_ro[i], dvar=xoptr[i], custom=Xcustom, lhs=lhs)
    #     '''=============== Step 8. Simulation  ====================='''
    #     randomq = [[[1 + s * 0.1 / sample if (ii, i) in [(2, 8),(4, 8),(5, 8),(6, 8)] else 1 for i in range(compart_num)]
    #             for ii in range(compart_num)] for s in range(sample)]
    #     x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)]
    #          for n in range(compart_num)]
    #     x[2][3] = [[xoptr[i][2][3][j][t] for t in range(period)] for j in range(group)]
    #     x[2][4] = [[xoptr[i][2][4][j][t] for t in range(period)] for j in range(group)]
    #     x[4][5] = [[xoptr[i][4][5][j][t] for t in range(period)] for j in range(group)]
    #     x[4][6] = [[xoptr[i][4][6][j][t] for t in range(period)] for j in range(group)]
    #     user.x = x
    #     print(user.x)
    #     Xsim_ro[i] = [0 for s in range(sample)]
    #     for s in range(sample):
    #         Xsim_ro[i][s] = user.Prediction(opt=['D', None, randomq[s]])
    #     '''======================== save =================================='''
    #     # lhs = [['totalcost',
    #     #         [sum([xsolution[('I', 'M')][j][t] * cts + xsolution[('I', 'Q')][j][t] * cts \
    #     #               + xsolution[('Q', 'H')][j][t] * cx + xsolution[('Q', 'U')][j][t] * cy
    #     #               for j in range(group) for t in range(period)]) + Xcustom['Hmax'] * ch + Xcustom['Umax'] * cu]]
    #     #        ]
    #     lhs = [['totalcost',
    #             [sum([xsolution[('Q', 'H')][j][t] * cx + xsolution[('Q', 'U')][j][t] * cy
    #                   for j in range(group) for t in range(period)]) + Xcustom['Hmax'] * ch + Xcustom['Umax'] * cu
    #              + sum(
    #                 [Xcustom['u.' + str(j) + '.' + str(s)] * ct * Test_capacity[t] * Test_probability[s] * (period - 1)
    #                  for j in range(group) for s in range(S)])]]
    #            ]
    #     for t in range(period):
    #         lhs.append(['Hospitalization.' + str(t), [sum([Xsim_ro[i][s][5][j][t] for j in range(group)])
    #                                                   for s in range(sample)]])
    #     for t in range(period):
    #         lhs.append(['ICU.' + str(t), [sum([Xsim_ro[i][s][6][j][t] for j in range(group)])
    #                                       for s in range(sample)]])
    #     lhs.append(['objective', [sum([Xsim_ro[i][s][2][j][t] * ci + Xsim_ro[i][s][-1][j][t] * cd
    #                                    for j in range(group)]) for s in range(sample)]])
    #
    #     user.save_result(filename='SEIHR/resultSEIHR_death_ro'+'_'+str(totalcost)+'_'+str(round(theta[i],2)),
    #                      compartment=Xsim_ro[i], dvar=xoptr[i], custom=Xcustom, lhs=lhs)


    '''====================== robust model with mutliple objective =========================='''
    theta = [1+0.02*(1+i) for i in range(10)]
    Xoptr2 = [[] for i in range(len(theta))]
    xoptr2 = [[] for i in range(len(theta))]
    Xcr2 = [[] for i in range(len(theta))]
    objr2 = [0 for i in range(len(theta))]
    #
    # Xsim_ro2 = [[[] for s in range(sample)] for ii in range(len(theta))]

    for i in range(len(theta)):
        '''=============== Step 5. Define other constraints  ====================='''
        print('step.5')
        user.set_time(time=period)
        # user.custom_lp(option='clear')
        # user.custom_lp(fvar=[('D', j, period-1) for j in range(group)],
        #                fcoef=[1 for j in range(group)],
        #                sense='L',
        #                rhs=sum([Xsolution_exp['D'][j][period-1] for j in range(group)])*theta[i],
        #                name='Death_number')
        user.set_solver(solver='gurobi')
        user.set_log_stream_SO(label=1,
                               file='SEIHR/log_SEIHR_multiro' + '_' + str(totalcost) + '_' + str(round(theta[i], 2)))
        solr2 = user.solve(label='Robustness', ct=3600, target=obj * theta[i])
        if solr2 != 0:
            (status, objr2[i], Xoptr2[i], xoptr2[i], Xcr2[i]) = solr2
            Xsolution = user.get_solution_compartment(compartment=compartment)
            xsolution = user.get_solution_flow_variable(compartmentfrom=['Q', 'Q'],
                                                        compartmentto=['H', 'U'])
            Xcustom = user.get_solution_custom_variable(name=['Hmax', 'Umax'])
            user.Solution_print(X=Xsolution, x=xsolution, xc=None)
        # '''=============== Step 8. Simulation  ====================='''
        # randomq = [[[0.7 + s * 0.6 / sample if (ii, i) in [(2,5),(3,5)] else 1 for i in range(compart_num)]
        #             for ii in range(compart_num)] for s in range(sample)]
        # x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)]
        #      for n in range(compart_num)]
        # x[0][1] = [[xoptr2[i][0][1][j][t] for t in range(period)] for j in range(group)]
        # x[1][2] = [[xoptr2[i][1][2][j][t] for t in range(period)] for j in range(group)]
        # x[1][3] = [[xoptr2[i][1][3][j][t] for t in range(period)] for j in range(group)]
        # user.x = x
        # print(user.x)
        # Xsim_ro2[i] = [0 for s in range(sample)]
        # for s in range(sample):
        #     Xsim_ro2[i][s] = user.Prediction(opt=['D', None, randomq[s]])
        #
        # '''======================== save =================================='''
        # # lhs = [['totalcost',
        # #         [sum([xsolution[('I', 'M')][j][t] * cts + xsolution[('I', 'Q')][j][t] * cts \
        # #               + xsolution[('Q', 'H')][j][t] * cx + xsolution[('Q', 'U')][j][t] * cy
        # #               for j in range(group) for t in range(period)]) + Xcustom['Hmax'] * ch + Xcustom['Umax'] * cu]]
        # #        ]
        # lhs = [['totalcost',
        #     [sum([xsolution_exp[('Q', 'H')][j][t] * cx + xsolution_exp[('Q', 'U')][j][t] * cy
        #           for j in range(group) for t in range(period)]) + Xcustom_exp['Hmax'] * ch + Xcustom_exp['Umax'] * cu ]
        #     ]]
        # for t in range(period):
        #     lhs.append(['Hospitalization.' + str(t), [sum([Xsim_ro2[i][s][2][j][t] for j in range(group)])
        #                                               for s in range(sample)]])
        # for t in range(period):
        #     lhs.append(['ICU.' + str(t), [sum([Xsim_ro2[i][s][3][j][t] for j in range(group)])
        #                                   for s in range(sample)]])
        # lhs.append(['objective', [sum([Xsim_ro2[i][s][-1][j][t] * cd
        #                                for j in range(group)]) for s in range(sample)]])
        #
        # user.save_result(filename='SEIHR/resultIHR_multiro' + '_' + str(totalcost) + '_' + str(round(theta[i], 2)),
        #                  compartment=Xsim_ro2[i], dvar=xoptr2[i], custom=Xcustom, lhs=lhs)
        # '''=============== Step 8. Simulation  ====================='''
        # # randomq = [[[1 + s * 0.1 / sample if (ii, i) in [(2, 8), (4, 8), (5, 8), (6, 8)] else 1 for i in
        # #              range(compart_num)]
        # #             for ii in range(compart_num)] for s in range(sample)]
        # # x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)]
        # #      for n in range(compart_num)]
        # # x[2][3] = [[xoptr2[i][2][3][j][t] for t in range(period)] for j in range(group)]
        # # x[2][4] = [[xoptr2[i][2][4][j][t] for t in range(period)] for j in range(group)]
        # # x[4][5] = [[xoptr2[i][4][5][j][t] for t in range(period)] for j in range(group)]
        # # x[4][6] = [[xoptr2[i][4][6][j][t] for t in range(period)] for j in range(group)]
        # # user.x = x
        # # print(user.x)
        # # Xsim_ro2[i] = [0 for s in range(sample)]
        # # for s in range(sample):
        # #     Xsim_ro2[i][s] = user.Prediction(opt=['D', None, randomq[s]])
        # # '''======================== save =================================='''
        # # # lhs = [['totalcost',
        # # #         [sum([xsolution[('I', 'M')][j][t] * cts + xsolution[('I', 'Q')][j][t] * cts \
        # # #               + xsolution[('Q', 'H')][j][t] * cx + xsolution[('Q', 'U')][j][t] * cy
        # # #               for j in range(group) for t in range(period)]) + Xcustom['Hmax'] * ch + Xcustom['Umax'] * cu]]
        # # #        ]
        # # lhs = [['totalcost',
        # #         [sum([xsolution[('Q', 'H')][j][t] * cx + xsolution[('Q', 'U')][j][t] * cy
        # #               for j in range(group) for t in range(period)]) + Xcustom['Hmax'] * ch + Xcustom['Umax'] * cu
        # #          + sum(
        # #             [Xcustom['u.' + str(j) + '.' + str(s)] * ct * Test_capacity[t] * Test_probability[s] * (
        # #             period - 1)
        # #              for j in range(group) for s in range(S)])]]
        # #        ]
        # # for t in range(period):
        # #     lhs.append(['Hospitalization.' + str(t), [sum([Xsim_ro2[i][s][5][j][t] for j in range(group)])
        # #                                               for s in range(sample)]])
        # # for t in range(period):
        # #     lhs.append(['ICU.' + str(t), [sum([Xsim_ro2[i][s][6][j][t] for j in range(group)])
        # #                                   for s in range(sample)]])
        # # lhs.append(['objective', [sum([Xsim_ro2[i][s][2][j][t] * ci + Xsim_ro2[i][s][-1][j][t] * cd
        # #                                for j in range(group)]) for s in range(sample)]])
        # # user.save_result(
        # #     filename='SEIHR/resultSEIHR_death_multiro' + '_' + str(totalcost) + '_' + str(round(theta[i], 2)),
        # #     compartment=Xsim_ro2[i], dvar=xoptr2[i], custom=Xcustom, lhs=lhs)

if __name__ == '__main__':
    for i in [6]:
        model(totalcost=(0 + i) * 1e+7)
    # for i in [500]:
    #     model(totalcost=i)
    # for i in [[550000,5500], [550000,5000], [500000,5000], [500000,4500]]:
    #     model(totalcost=i)
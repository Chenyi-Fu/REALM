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
totalperiod = 100
group = data.group
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
Qdata = data.Qdata
def model(totalcost):
    '''=============== Step 1. Define model and model size ====================='''
    user = REALM.REALM()
    user.set_group(group=group)
    user.set_time(time=period)
    H0 = 1000
    ICU0 = 150
    '''=============== Step 2. Define name of compartment and population ====================='''
    compartment = ['I','Q','H','U','R','D']
    compart_num = len(compartment)
    user.set_all_compartment(name=compartment)
    for i in range(compart_num):
        print(population[i])
    for i in range(compart_num):
        user.set_population(compartment=compartment[i], population=population[i])
    '''=============== Step 3. Define decision variable  ====================='''
    user.set_flow_variable(n='Q', m = 'H', xupper=[[H0 for t in range(period)] for j in range(group)],
                              xlower=[[0 for t in range(period)] for j in range(group)])
    user.set_flow_variable(n='Q', m = 'U', xupper=[[ICU0 for t in range(period)] for j in range(group)],
                              xlower=[[0 for t in range(period)] for j in range(group)])
    user.set_flow_variable(n='I', m = 'Q', xupper=[[Qdata[j][t+1]-Qdata[j][t] if t <period-1 else 0
                                                    for t in range(period)] for j in range(group)],
                              xlower=[[Qdata[j][t+1]-Qdata[j][t] if t <period-1 else 0
                                                    for t in range(period)] for j in range(group)])

    '''====== Step 4.1 Define decision-independent transition between compartment and compartment  ========='''
    user.set_transition_compartment(n='Q', m='D', prob=[[q_QtoD[j][t] for t in range(period)] for j in range(group)], opt = None)
    user.set_transition_compartment(n='H', m='D', prob=[[q_HtoD[j][t] for t in range(period)] for j in range(group)], opt = None)
    user.set_transition_compartment(n='U', m='D', prob=[[q_UtoD[j][t] for t in range(period)] for j in range(group)], opt = None)
    user.set_transition_compartment(n='H', m='R', prob=[[q_HtoR[j][t] for t in range(period)] for j in range(group)], opt = None)
    user.set_transition_compartment(n='U', m='R', prob=[[q_UtoR[j][t] for t in range(period)] for j in range(group)], opt = None)

    '''=============== Step 5. Define other constraints  ====================='''
    for t in range(1,period):
        user.custom_lp(fvar=[('H', j, t) for j in range(group)], fcoef=[1 for j in range(group)],
                       sense='L', rhs=H0, name='Hospitalization.' + str(t))
        user.custom_lp(fvar=[('U', j, t) for j in range(group)], fcoef=[1 for j in range(group)],
                       sense='L', rhs=ICU0, name='ICU.' + str(t))

    user.custom_lp(fvar=[('Q', 'H', j, t) for j in range(group) for t in range(period-1)] \
                    + [('Q', 'U', j, t) for j in range(group) for t in range(period-1)],
                   fcoef=[cx/Various_lambdah[j] for j in range(group) for t in range(period-1)]
                         + [cy/Various_lambdaicu[j] for j in range(group) for t in range(period-1)],
                   sense='L',
                   rhs=totalcost,
                   name='totalcost')
    '''=============== Step 6. Define objective  ====================='''
    user.set_objectivetype(sense='min')
    for j in range(group):
        user.set_objective(state=('D', j, period - 1), value=cd)

    '''=============== Step 7. Solve and output  ====================='''
    user.set_solver(solver='gurobi')
    user.set_approximation(opt='SO')
    user.set_log_stream_SO(label=1, file='SEIHR/log_SEIHR_exp_'+str(totalcost))
    sol = user.solve(label='Expectation', ct=3600)
    (status, obj, Xopt, xopt, Xc) = sol
    Xsolution_exp = user.get_solution_compartment(compartment=compartment)
    xsolution_exp = user.get_solution_flow_variable(compartmentfrom=['Q','Q'],compartmentto=['H','U'])
    user.Solution_print(X=Xsolution_exp,x=xsolution_exp,xc=None)

    '''======================== Simulation =================================='''
    sample = 1000
    randomq = user.Sample_generation(n=['H','U'], m=['D','D'], lb=[0.7,0.7], ub=[1.3,1.3], size=sample)
    x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)] for n in range(compart_num)]
    x[0][1] = [[xopt[0][1][j][t] for t in range(period)] for j in range(group)]
    x[1][2] = [[xopt[1][2][j][t] for t in range(period)] for j in range(group)]
    x[1][3] = [[xopt[1][3][j][t] for t in range(period)] for j in range(group)]
    user.x = x
    Xsim_exp = [0 for s in range(sample)]
    for s in range(sample):
        Xsim_exp[s] = user.Prediction(opt=['S', None, randomq[s]])
    '''======================== save =================================='''
    lhs = [['totalcost',
            [sum([xsolution_exp[('Q', 'H')][j][t] * cx + xsolution_exp[('Q', 'U')][j][t] * cy
                  for j in range(group) for t in range(period)]) + Xcustom_exp['Hmax'] * ch + Xcustom_exp['Umax'] * cu ]
            ]]
    for t in range(period):
        lhs.append(['Hospitalization.' + str(t), [sum([Xsim_exp[s][2][j][t] for j in range(group)]) for s in range(sample)]])
    for t in range(period):
        lhs.append(['ICU.' + str(t), [sum([Xsim_exp[s][3][j][t] for j in range(group)]) for s in range(sample)]])
    lhs.append(['objective', [sum([Xsim_exp[s][-1][j][t] * cd for j in range(group)]) for s in range(sample)]])
    user.save_result(filename='SEIHR/resultIHR_exp_'+str(totalcost),
                     compartment=Xsim_exp, dvar=xopt, custom=Xcustom_exp, lhs=lhs)

    '''====================== robust model =========================='''
    print('Begin to solve robust model')
    theta = [1+0.02*(1+i) for i in range(10)]
    Xoptr2 = [[] for i in range(len(theta))]
    xoptr2 = [[] for i in range(len(theta))]
    Xcr2 = [[] for i in range(len(theta))]
    objr2 = [0 for i in range(len(theta))]

    for i in range(len(theta)):
        '''=============== Step 5. Define other constraints  ====================='''
        user.set_time(time=period)
        user.set_solver(solver='gurobi')
        user.set_log_stream_SO(label=1,
                               file='SEIHR/log_SEIHR_ro' + '_' + str(totalcost) + '_' + str(round(theta[i], 2)))
        solr2 = user.solve(label='Robustness', ct=3600, target=obj * theta[i])
        if solr2 != 0:
            (status, objr2[i], Xoptr2[i], xoptr2[i], Xcr2[i]) = solr2
            Xsolution = user.get_solution_compartment(compartment=compartment)
            xsolution = user.get_solution_flow_variable(compartmentfrom=['Q', 'Q'],
                                                        compartmentto=['H', 'U'])
            user.Solution_print(X=Xsolution, x=xsolution, xc=None)
        '''=============== Step 8. Simulation  ====================='''
        x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)]
             for n in range(compart_num)]
        x[0][1] = [[xoptr2[i][0][1][j][t] for t in range(period)] for j in range(group)]
        x[1][2] = [[xoptr2[i][1][2][j][t] for t in range(period)] for j in range(group)]
        x[1][3] = [[xoptr2[i][1][3][j][t] for t in range(period)] for j in range(group)]
        user.x = x
        Xsim_ro2[i] = [0 for s in range(sample)]
        for s in range(sample):
            Xsim_ro2[i][s] = user.Prediction(opt=['D', None, randomq[s]])

        '''======================== save =================================='''
        lhs = [['totalcost', [sum([xsolution_exp[('Q', 'H')][j][t] * cx + xsolution_exp[('Q', 'U')][j][t] * cy
                  for j in range(group) for t in range(period)]) + Xcustom_exp['Hmax'] * ch + Xcustom_exp['Umax'] * cu]
            ]]
        for t in range(period):
            lhs.append(['Hospitalization.' + str(t), [sum([Xsim_ro2[i][s][2][j][t] for j in range(group)])
                                                      for s in range(sample)]])
        for t in range(period):
            lhs.append(['ICU.' + str(t), [sum([Xsim_ro2[i][s][3][j][t] for j in range(group)])
                                          for s in range(sample)]])
        lhs.append(['objective', [sum([Xsim_ro2[i][s][-1][j][t] * cd
                                       for j in range(group)]) for s in range(sample)]])
        user.save_result(filename='SEIHR/resultIHR_ro' + '_' + str(totalcost) + '_' + str(round(theta[i], 2)),
                         compartment=Xsim_ro2[i], dvar=xoptr2[i], custom=Xcustom, lhs=lhs)
if __name__ == '__main__':
    for i in [6]:
        model(totalcost=(0 + i) * 1e+7)
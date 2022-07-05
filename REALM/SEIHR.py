from __future__ import print_function
import random
from math import exp
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
pm = [[(1-lam[t]) * (1-piD[j]) + lam[t] * (1-piO[j]) for t in range(period)] for j in range(group)]
pq = [[(1-lam[t]) * piD[j] + lam[t] * piO[j] for t in range(period)] for j in range(group)]

beta = [(1 - lam[t]) * betaD + lam[t] * betaO for t in range(period)]
q_EtoI = [[(1 - lam[t]) * alphaD + lam[t] * alphaO for t in range(period)] for j in range(group)]
q_ItoR = [[(1-lam[t]) * gammaD * (1-piD[j]) + lam[t] * gammaO * (1-piO[j])
           for t in range(period)] for j in range(group)]
q_ItoD = [[(1-lam[t]) * dqD[j] * piD[j] + lam[t] * dqO[j] * piO[j] for t in range(period)] for j in range(group)]
q_MtoR = [[(1-lam[t]) * gammaD + lam[t] * gammaO for t in range(period)] for j in range(group)]
q_QtoD = [[(1-lam[t]) * dqD[j] + lam[t] * dqO[j]
           for t in range(period)] for j in range(group)]
q_HtoD = [[(1-lam[t]) * dhD[j] * lambdahD[j] + lam[t] * dhO[j] * lambdahO[j]
           for t in range(period)] for j in range(group)]
q_HtoR = [[(1-lam[t]) * (1-dhD[j]) * lambdahD[j] + lam[t] * (1-dhO[j]) * lambdahO[j]
           for t in range(period)] for j in range(group)]
q_UtoD = [[(1-lam[t]) * dicuD[j] * lambdaicuD[j] + lam[t] * dicuO[j] * lambdaicuO[j]
           for t in range(period)] for j in range(group)]
q_UtoR = [[(1-lam[t]) * (1-dicuD[j]) * lambdaicuD[j] + lam[t] * (1-dicuO[j]) * lambdaicuO[j]
           for t in range(period)] for j in range(group)]
attractive = 0.5
S = 6
Test_capacity = [30000 for t in range(period)]
Test_probability = [0.05*(i+1) for i in range(S)]
H0 = 1500
ICU0 = 200

ch = data.ch
cu = data.cu
cx = data.cx
cy = data.cy
cd = data.cd
ci = data.ci
ct = data.ct
cts = data.cts

def model(totalcost):
    '''=============== Step 1. Define model and model size ====================='''
    user = REALM.REALM()
    user.set_group(group=group)
    user.set_time(time=period)

    '''=============== Step 2. Define name of compartment and population ====================='''
    compartment = ['S','E','I','M','Q','H','U','R','D']
    compart_num = len(compartment)
    user.set_all_compartment(name=compartment)
    population = data.compartment_population
    total_pop = data.group_population
    for i in range(compart_num):
        user.set_population(compartment=compartment[i], population=population[i])
    print('step.2')
    '''=============== Step 3. Define decision variable  ====================='''
    user.set_flow_variable(n='Q', m = 'H', xupper=[[10*H0 for t in range(period)] for j in range(group)],
                              xlower=[[0 for t in range(period)] for j in range(group)])
    user.set_flow_variable(n='Q', m = 'U', xupper=[[10*ICU0 for t in range(period)] for j in range(group)],
                              xlower=[[0 for t in range(period)] for j in range(group)])
    user.set_flow_variable(n='I', m = 'M', xupper=[[Test_capacity[t]*Test_probability[-1]
                                                    for t in range(period)] for j in range(group)],
                              xlower=[[0 for t in range(period)] for j in range(group)])
    user.set_flow_variable(n='I', m = 'Q', xupper=[[Test_capacity[t]*Test_probability[-1]
                                                    for t in range(period)] for j in range(group)],
                              xlower=[[0 for t in range(period)] for j in range(group)])
    for j in range(group):
        for s in range(S):
            user.set_custom_variable(name='u.' + str(j) + '.' + str(s), types='B', xupper=1, xlower=0)
            for t in range(period):
                user.set_custom_variable(name='z.' + str(j) + '.' + str(t) + '.' + str(s), xlower=0)
                user.set_custom_variable(name='phi1.' + str(j) + '.' + str(t) + '.'+str(s), xlower=-1e20)
                user.set_custom_variable(name='phi2.' + str(j) + '.' + str(t) + '.'+str(s), xlower=0)

    print('step.3')
    '''====== Step 4.1 Define decision-independent transition between compartment and compartment  ========='''
    user.set_transition_compartment(n='E', m='I', prob=[[q_EtoI[j][t] for t in range(period)] for j in range(group)], opt = None)
    user.set_transition_compartment(n='I', m='R', prob=[[q_ItoR[j][t] for t in range(period)] for j in range(group)], opt = None)
    user.set_transition_compartment(n='I', m='D', prob=[[q_ItoD[j][t] for t in range(period)] for j in range(group)], opt = None)
    user.set_transition_compartment(n='M', m='R', prob=[[q_MtoR[j][t] for t in range(period)] for j in range(group)], opt = None)
    user.set_transition_compartment(n='Q', m='D', prob=[[q_QtoD[j][t] for t in range(period)] for j in range(group)], opt = None)
    user.set_transition_compartment(n='H', m='D', prob=[[q_HtoD[j][t] for t in range(period)] for j in range(group)], opt = None)
    user.set_transition_compartment(n='U', m='D', prob=[[q_UtoD[j][t] for t in range(period)] for j in range(group)], opt = None)
    user.set_transition_compartment(n='H', m='R', prob=[[q_HtoR[j][t] for t in range(period)] for j in range(group)], opt = None)
    user.set_transition_compartment(n='U', m='R', prob=[[q_UtoR[j][t] for t in range(period)] for j in range(group)], opt = None)

    print('step.4.1')
    '''======= Step 4.2 Define decision-dependent transition between compartment and compartment  ======'''
    for j in range(group):
        for t in range(period):
            betaSE = [[beta[t] * w[j] * w[k] / total_pop[k]
                       if compartment[n] == 'I' else 0 for k in range(group)] for n in range(len(population))]
            gammaSE = [[0 for k in range(group)] for n in range(len(population))]
            user.set_transition_compartment_dp(flow=('S', 'E', j, t),
                                               beta=betaSE, gamma=gammaSE, alpha=0, opt=None)

    '''=============== Step 5. Define other constraints  ====================='''
    for t in range(period):
        for j in range(group):
            user.custom_lp(fvar=[('I', 'M', j, t),  ('I', 'Q', j, t), ('I', j, t)],
                           fcoef=[1, 1, -attractive],
                           sense='L',
                           rhs=0,
                           name='testAllocation.' + str(j) + '.' + str(t),
                           label='Expectation')
            user.custom_lp(fvar=[('I', 'M', j, t),  ('I', j, t)],
                           fcoef=[1, -attractive*((1-lam[t]) * (1-piD[j]) + lam[t] * (1-piO[j]))],
                           sense='L',
                           rhs=0,
                           name='testAllocation1.' + str(j) + '.' + str(t),
                           label='Expectation')
            user.custom_lp(fvar=[('I', 'Q', j, t), ('I', j, t)],
                           fcoef=[1, -attractive * ((1-lam[t]) * piD[j] + lam[t] * piO[j])],
                           sense='L',
                           rhs=0,
                           name='testAllocation2.' + str(j) + '.' + str(t),
                           label='Expectation')

            user.custom_lp(fvar= [('I', 'M', j, t)],
                           fcoef= [-1],
                           dvar=['z.' + str(j) + '.' + str(t) + '.' + str(s) for s in range(S)],
                           dcoef=[(1-lam[t]) * (1-piD[j]) + lam[t] * (1-piO[j]) for s in range(S)],
                           sense='E',
                           rhs=0,
                           name='Positive_ItoM.' + str(j) + '.' + str(t))
            user.custom_lp(fvar= [('I', 'Q', j, t)],
                           fcoef= [-1],
                           dvar=['z.' + str(j) + '.' + str(t) + '.' + str(s) for s in range(S)],
                           dcoef=[(1-lam[t]) * piD[j] + lam[t] * piO[j] for s in range(S)],
                           sense='E',
                           rhs=0,
                           name='Positive_ItoQ.' + str(j) + '.' + str(t))
            for s in range(S):
                user.custom_lp(fvar=[('I', j, t)],
                               fcoef=[1],
                               dvar=['phi1.'+str(j)+'.'+str(t) + '.' + str(s),
                                     'z.' + str(j) + '.' + str(t) + '.' + str(s)],
                               dcoef=[1, 1],
                               sense='E',
                               rhs=Test_capacity[t]*Test_probability[s]-bhat[j],
                               name='SOCequal1.' + str(j) + '.' + str(t) + '.' + str(s),
                               label='Expectation')
                user.custom_lp(fvar=[('I', j, t)],
                               fcoef=[-1],
                               dvar=['phi2.'+str(j)+'.'+str(t) + '.' + str(s),
                                     'z.' + str(j) + '.' + str(t) + '.' + str(s)],
                               dcoef=[1, 1],
                               sense='E',
                               rhs=Test_capacity[t]*Test_probability[s]+bhat[j],
                               name='SOCequal2.' + str(j) + '.' + str(t) + '.' + str(s),
                               label='Expectation')
                user.custom_qp(dvarqp=['phi1.' + str(j) + '.' + str(t) + '.' + str(s),
                                       'phi2.' + str(j) + '.' + str(t) + '.' + str(s)],
                               dcoefqp=[0.25, -0.25],
                               sense='L',
                               rhs= -Test_capacity[t] * Test_probability[s] * bhat[j],
                               name='SOC.' + str(j) + '.' + str(t) + '.' + str(s))

                user.custom_lp(dvar=['z.' + str(j) + '.' + str(t) + '.' + str(s),
                                     'u.' + str(j) + '.' + str(s)],
                               dcoef=[1, -Test_capacity[t] * Test_probability[s]],
                               sense='L',
                               rhs=0,
                               name='SOC_upperbound_1.' + str(j) + '.' + str(t) + '.' + str(s),
                               label='Expectation')
    for j in range(group):
        user.custom_lp(dvar=['u.' + str(j) + '.' + str(s) for s in range(S)],
                       dcoef=[1 for s in range(S)],
                       sense='E',
                       rhs=1,
                       name='Uequal.' + str(j))

    user.custom_lp(fvar=[('Q', 'H', j, t) for j in range(group) for t in range(period-1)] \
                    + [('Q', 'U', j, t) for j in range(group) for t in range(period-1)],
                   fcoef=[cx/Various_lambdah[j] for j in range(group) for t in range(period-1)]
                         + [cy/Various_lambdaicu[j] for j in range(group) for t in range(period-1)],
                   dvar=['u.' + str(j) + '.' + str(s) for j in range(group) for s in range(S)],
                   dcoef=[ct*Test_capacity[t] * Test_probability[s] * (period-1) for j in range(group) for s in range(S)],
                   sense='L',
                   rhs=totalcost,
                   name='totalcost')
    '''=============== Step 6. Define objective  ====================='''
    user.set_objectivetype(sense='min')
    for j in range(group):
        for t in range(period):
            user.set_objective(state=('I', j, t), value=ci)
        user.set_objective(state=('D', j, period - 1), value=cd)
    '''=============== Step 7. Solve and output  ====================='''
    user.set_solver(solver='cplex')
    user.set_approximation(opt='SO')
    user.set_log_stream_SO(label=1, file='SEIHR/log_SEIHR_exp_'+str(totalcost))
    sol = user.solve(label='Expectation', ct=3600)
    (status, obj, Xopt, xopt, Xc) = sol
    print('obj=', obj)
    Xsolution_exp = user.get_solution_compartment(compartment=compartment)
    xsolution_exp = user.get_solution_flow_variable(compartmentfrom=['I','I','Q','Q'],compartmentto=['M','Q','H','U'])
    Xcustom_exp = user.get_solution_custom_variable(name=['u.' + str(j) + '.' + str(s) for s in range(S) for j in range(group)]
                                                  + ['phi1.' + str(j) + '.' + str(t) + '.'+str(s)
                                                     for s in range(S) for j in range(group) for t in range(period)]
                                                  + ['phi2.' + str(j) + '.' + str(t) + '.'+str(s)
                                                     for s in range(S) for j in range(group) for t in range(period)]
                                                  + ['z.' + str(j) + '.' + str(t) + '.' + str(s)
                                                     for s in range(S) for j in range(group) for t in range(period)])
    user.Solution_print(X=Xsolution_exp,x=xsolution_exp,xc=None)

    '''======================== Simulation =================================='''
    sample = 1000
    randomq = user.Sample_generation(n=['S','I'], m=['E','D'], lb=[1,1], ub=[1.01,1.01], size=sample)
    x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)] for n in range(compart_num)]
    x[2][3] = [[xopt[2][3][j][t] for t in range(period)] for j in range(group)]
    x[2][4] = [[xopt[2][4][j][t] for t in range(period)] for j in range(group)]
    x[4][5] = [[xopt[4][5][j][t] for t in range(period)] for j in range(group)]
    x[4][6] = [[xopt[4][6][j][t] for t in range(period)] for j in range(group)]
    user.x = x
    Xsim_exp = [0 for s in range(sample)]
    for s in range(sample):
        Xsim_exp[s] = user.Prediction(opt=['D', None, randomq[s]])
    '''======================== save =================================='''
    lhs = [['totalcost',
            [sum([xsolution_exp[('Q', 'H')][j][t] * cx + xsolution_exp[('Q', 'U')][j][t] * cy
                  for j in range(group) for t in range(period)]) + Xcustom_exp['Hmax'] * ch + Xcustom_exp['Umax'] * cu
             + sum([Xcustom_exp['u.' + str(j) + '.' + str(s)] * ct * Test_capacity[t] * Test_probability[s] * (period-1)
             for j in range(group) for s in range(S)]) ]]
           ]
    for t in range(period):
        lhs.append(['Hospitalization.' + str(t),[sum([Xsim_exp[s][5][j][t] for j in range(group)])
                                 for s in range(sample)]])
    for t in range(period):
        lhs.append(['ICU.' + str(t), [sum([Xsim_exp[s][6][j][t] for j in range(group)])
                                 for s in range(sample)]])
    lhs.append(['objective', [sum([Xsim_exp[s][2][j][t] * ci + Xsim_exp[s][-1][j][t] * cd for j in range(group)])
                              for s in range(sample)]])
    user.save_result(filename='SEIHR/resultSEIHR_exp_'+str(totalcost),
                     compartment=Xsim_exp, dvar=xopt, custom=Xcustom_exp, lhs=lhs)
    '''====================== robust model =========================='''
    print('Begin to solve robust model')
    theta = [1+0.01*(1+i) for i in range(4)]
    Xoptr = [[] for i in range(len(theta))]
    Xcr = [[] for i in range(len(theta))]
    xoptr = [[] for i in range(len(theta))]
    objr = [0 for i in range(len(theta))]
    Xsim_ro = [[[] for s in range(sample)] for ii in range(len(theta))]

    for i in range(len(theta)):
        user.set_solver(solver='cplex')
        user.set_log_stream_SO(label=1, file='SEIHR/log_SEIHR_ro'+'_'+str(totalcost)+'_'+str(round(theta[i],2)))
        solr = user.default_lp(label='Robustness', ct=3600, target=obj*theta[i])
        if solr != 0:
            (status, objr[i], Xoptr[i], xoptr[i], Xcr[i]) = solr
            Xsolution = user.get_solution_compartment(compartment=compartment)
            xsolution = user.get_solution_flow_variable(compartmentfrom=['I','I','Q','Q'],compartmentto=['M','Q','H','U'])
            Xcustom = user.get_solution_custom_variable(name=['Hmax', 'Umax']
                                                          + ['u.' + str(j) + '.' + str(s) for s in range(S) for j in range(group)]
                                                          + ['phi1.' + str(j) + '.' + str(t) + '.'+str(s)
                                                             for s in range(S) for j in range(group) for t in range(period)]
                                                          + ['phi2.' + str(j) + '.' + str(t) + '.'+str(s)
                                                             for s in range(S) for j in range(group) for t in range(period)]
                                                          + ['z.' + str(j) + '.' + str(t) + '.' + str(s)
                                                             for s in range(S) for j in range(group) for t in range(period)])
            user.Solution_print(X=Xsolution,x=xsolution,xc=None)
        '''=============== Step 8. Simulation  ====================='''
        randomq = [[[1 + s * 0.1 / sample if (ii, i) in [(0, 1)] else 1 for i in range(compart_num)]
                    for ii in range(compart_num)] for s in range(sample)]
        x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)]
             for n in range(compart_num)]
        x[2][3] = [[xoptr[i][2][3][j][t] for t in range(period)] for j in range(group)]
        x[2][4] = [[xoptr[i][2][4][j][t] for t in range(period)] for j in range(group)]
        x[4][5] = [[xoptr[i][4][5][j][t] for t in range(period)] for j in range(group)]
        x[4][6] = [[xoptr[i][4][6][j][t] for t in range(period)] for j in range(group)]
        user.x = x
        print(user.x)
        Xsim_ro[i] = [0 for s in range(sample)]
        for s in range(sample):
            Xsim_ro[i][s] = user.Prediction(opt=['D', None, randomq[s]])

        '''======================== save =================================='''
        lhs = [['totalcost',
                [sum([xsolution[('Q', 'H')][j][t] * cx + xsolution[('Q', 'U')][j][t] * cy
                      for j in range(group) for t in range(period)]) + Xcustom['Hmax'] * ch + Xcustom['Umax'] * cu
                 + sum(
                    [Xcustom['u.' + str(j) + '.' + str(s)] * ct * Test_capacity[t] * Test_probability[s] * (period - 1)
                     for j in range(group) for s in range(S)])]]
               ]
        for t in range(period):
            lhs.append(['Hospitalization.' + str(t), [sum([Xsim_ro[i][s][5][j][t] for j in range(group)])
                                                      for s in range(sample)]])
        for t in range(period):
            lhs.append(['ICU.' + str(t), [sum([Xsim_ro[i][s][6][j][t] for j in range(group)])
                                          for s in range(sample)]])
        lhs.append(['objective', [sum([Xsim_ro[i][s][2][j][t] * ci + Xsim_ro[i][s][-1][j][t] * cd
                                       for j in range(group)]) for s in range(sample)]])

        user.save_result(filename='SEIHR/resultSEIHR_ro'+'_'+str(totalcost)+'_'+str(round(theta[i],2)),
                         compartment=Xsim_ro[i], dvar=xoptr[i], custom=Xcustom, lhs=lhs)
        '''=============== Step 8. Simulation  ====================='''
        x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)]
             for n in range(compart_num)]
        x[2][3] = [[xoptr[i][2][3][j][t] for t in range(period)] for j in range(group)]
        x[2][4] = [[xoptr[i][2][4][j][t] for t in range(period)] for j in range(group)]
        x[4][5] = [[xoptr[i][4][5][j][t] for t in range(period)] for j in range(group)]
        x[4][6] = [[xoptr[i][4][6][j][t] for t in range(period)] for j in range(group)]
        user.x = x
        print(user.x)
        Xsim_ro[i] = [0 for s in range(sample)]
        for s in range(sample):
            Xsim_ro[i][s] = user.Prediction(opt=['D', None, randomq[s]])
        '''======================== save =================================='''
        lhs = [['totalcost',
                [sum([xsolution[('Q', 'H')][j][t] * cx + xsolution[('Q', 'U')][j][t] * cy
                      for j in range(group) for t in range(period)]) + Xcustom['Hmax'] * ch + Xcustom['Umax'] * cu
                 + sum(
                    [Xcustom['u.' + str(j) + '.' + str(s)] * ct * Test_capacity[t] * Test_probability[s] * (period - 1)
                     for j in range(group) for s in range(S)])]]
               ]
        for t in range(period):
            lhs.append(['Hospitalization.' + str(t), [sum([Xsim_ro[i][s][5][j][t] for j in range(group)])
                                                      for s in range(sample)]])
        for t in range(period):
            lhs.append(['ICU.' + str(t), [sum([Xsim_ro[i][s][6][j][t] for j in range(group)])
                                          for s in range(sample)]])
        lhs.append(['objective', [sum([Xsim_ro[i][s][2][j][t] * ci + Xsim_ro[i][s][-1][j][t] * cd
                                       for j in range(group)]) for s in range(sample)]])

        user.save_result(filename='SEIHR/resultSEIHR_death_ro'+'_'+str(totalcost)+'_'+str(round(theta[i],2)),
                         compartment=Xsim_ro[i], dvar=xoptr[i], custom=Xcustom, lhs=lhs)
if __name__ == '__main__':
    for i in [10]:
        model(totalcost=(0 + i) * 1e+7)
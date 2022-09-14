from __future__ import print_function
import random
from math import exp
import REALM
import DatasetSG_4layers
from gurobipy import GRB
import time
import multiprocessing as mp



def model(totalcost):
    '''=============== Step 1. Define model and model size ====================='''
    [capacity_V] = totalcost
    user = REALM.REALM()

    # user.set_types(types = types) # Expectation or Robustness
    user.set_group(group=group)
    user.set_time(time=period)

    '''=============== Step 2. Define name of compartment and population ====================='''
    compartment = ['S', 'E', 'I', 'M', 'H', 'R', 'D',
                   'VS', 'VE', 'VI', 'VM', 'VH', 'VR', 'VD',
                   'BS', 'BE', 'BI', 'BM', 'BH', 'BR', 'BD',
                   'BS2', 'BE2', 'BI2', 'BM2', 'BH2', 'BR2', 'BD2']
    compart_num = len(compartment)
    user.set_all_compartment(name=compartment)
    population = data.population
    total_pop = data.total_pop
    # print(population)
    for i in range(compart_num):
        user.set_population(compartment=compartment[i], population=population[i])
    print('step.2')
    '''=============== Step 3. Define decision variable  ====================='''
    # user.set_flow_variable(n='S', m = 'VS', xupper=[[0 for t in range(period)] for j in range(group)],
    #                           xlower=[[0 for t in range(period)] for j in range(group)])
    user.set_flow_variable(n='VS', m = 'BS',
                           xupper=[[capacity_V if 6>=j >0 else 0 for t in range(period)] for j in range(group)],
                              xlower=[[0 for t in range(period)] for j in range(group)])
    # user.set_flow_variable(n='VS', m = 'BS',
    #                        xupper=[[capacity_V for t in range(period)] for j in range(group)],
    #                        xlower=[[0 for t in range(period)] for j in range(group)])
    capacityBS2 = [[0 for t in range(period)] for j in range(group)]
    for j in range(group):
        for t in range(period):
            if t < period/2:
                if j > 5:
                    capacityBS2[j][t] = capacity_V
            else:
                if j > 5:
                    capacityBS2[j][t] = capacity_V

    user.set_flow_variable(n='BS', m='BS2', xupper=capacityBS2, xlower=[[0 for t in range(period)] for j in range(group)])
    # user.set_flow_variable(n='BS', m='BS2',
    # xupper=[[capacity_V if j > 1 else 0 for t in range(period)] for j in range(group)],
    #                        xlower=[[0 for t in range(period)] for j in range(group)])
    user.set_custom_variable(name='Hmax', xlower=0)

    print('step.3')
    '''====== Step 4.1 Define decision-independent transition between compartment and compartment  ========='''
    user.set_transition_compartment(n='E', m='I', prob=[[q_EtoI[j][t] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='VE', m='VI', prob=[[q_VEtoVI[j][t] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='BE', m='BI', prob=[[q_BEtoBI[j][t] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='BE2', m='BI2', prob=[[q_BE2toBI2[j][t] for t in range(period)] for j in range(group)])

    user.set_transition_compartment(n='I', m='M', prob=[[q_ItoM[j][t] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='VI', m='VM', prob=[[q_VItoVM[j][t] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='BI', m='BM', prob=[[q_BItoBM[j][t] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='BI2', m='BM2', prob=[[q_BI2toBM2[j][t] for t in range(period)] for j in range(group)])

    user.set_transition_compartment(n='I', m='H', prob=[[q_ItoH[j][t] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='VI', m='VH', prob=[[q_VItoVH[j][t] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='BI', m='BH', prob=[[q_BItoBH[j][t] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='BI2', m='BH2', prob=[[q_BI2toBH2[j][t] for t in range(period)] for j in range(group)])

    user.set_transition_compartment(n='M', m='R', prob=[[1/5 for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='VM', m='VR', prob=[[1/5 for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='BM', m='BR', prob=[[1/5 for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='BM2', m='BR2', prob=[[1/5 for t in range(period)] for j in range(group)])

    user.set_transition_compartment(n='H', m='R', prob=[[q_HtoR[j][t] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='VH', m='VR', prob=[[q_VHtoVR[j][t] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='BH', m='BR', prob=[[q_BHtoBR[j][t] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='BH2', m='BR2', prob=[[q_BH2toBR2[j][t] for t in range(period)] for j in range(group)])

    user.set_transition_compartment(n='H', m='D', prob=[[q_HtoD[j][t] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='VH', m='VD', prob=[[q_VHtoVD[j][t] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='BH', m='BD', prob=[[q_BHtoBD[j][t] for t in range(period)] for j in range(group)])
    user.set_transition_compartment(n='BH2', m='BD2', prob=[[q_BH2toBD2[j][t] for t in range(period)] for j in range(group)])


    print('step.4.1')
    '''======= Step 4.2 Define decision-dependent transition between compartment and compartment  ======'''
    for j in range(group):
        for t in range(period):
            betaSE = [[0 for k in range(group)] for n in range(len(population))]
            for n in range(len(population)):
                if compartment[n] in ['I']:
                    for k in range(group):
                        betaSE[n][k] = beta[0] * w[j][k] * policy[0] * social[t] / total_pop[k]
                elif compartment[n] in ['VI']:
                    for k in range(group):
                        betaSE[n][k] = beta[0] * w[j][k] * policy[0] * social[t] / total_pop[k]
                elif compartment[n] in ['BI']:
                    for k in range(group):
                        betaSE[n][k] = beta[0] * w[j][k] * policy[0] * social[t] / total_pop[k]
                elif compartment[n] in ['BI2']:
                    for k in range(group):
                        betaSE[n][k] = beta[0] * w[j][k] * policy[0] * social[t] / total_pop[k]
            user.set_transition_compartment_dp(flow=('S', 'E', j, t), beta=betaSE)

    for j in range(group):
        for t in range(period):
            betaSE = [[0 for k in range(group)] for n in range(len(population))]
            for n in range(len(population)):
                if compartment[n] in ['I']:
                    for k in range(group):
                        betaSE[n][k] = beta[1] * w[j][k] * policy[1] * social[t] / total_pop[k]
                elif compartment[n] in ['VI']:
                    for k in range(group):
                        betaSE[n][k] = beta[1] * w[j][k] * policy[1] * social[t] / total_pop[k]
                elif compartment[n] in ['BI']:
                    for k in range(group):
                        betaSE[n][k] = beta[1] * w[j][k] * policy[1] * social[t] / total_pop[k]
                elif compartment[n] in ['BI2']:
                    for k in range(group):
                        betaSE[n][k] = beta[1] * w[j][k] * policy[1] * social[t] / total_pop[k]
            user.set_transition_compartment_dp(flow=('VS', 'VE', j, t), beta=betaSE)

    for j in range(group):
        for t in range(period):
            betaSE = [[0 for k in range(group)] for n in range(len(population))]
            for n in range(len(population)):
                if compartment[n] in ['I']:
                    for k in range(group):
                        betaSE[n][k] = beta[2] * w[j][k] * policy[2] * social[t] / total_pop[k]
                elif compartment[n] in ['VI']:
                    for k in range(group):
                        betaSE[n][k] = beta[2] * w[j][k] * policy[2] * social[t] / total_pop[k]
                elif compartment[n] in ['BI']:
                    for k in range(group):
                        betaSE[n][k] = beta[2] * w[j][k] * policy[2] * social[t] / total_pop[k]
                elif compartment[n] in ['BI2']:
                    for k in range(group):
                        betaSE[n][k] = beta[2] * w[j][k] * policy[2] * social[t] / total_pop[k]
            user.set_transition_compartment_dp(flow=('BS', 'BE', j, t), beta=betaSE)

    for j in range(group):
        for t in range(period):
            betaSE = [[0 for k in range(group)] for n in range(len(population))]
            for n in range(len(population)):
                if compartment[n] in ['I']:
                    for k in range(group):
                        betaSE[n][k] = beta[3] * w[j][k] * policy[3] * social[t] / total_pop[k]
                elif compartment[n] in ['VI']:
                    for k in range(group):
                        betaSE[n][k] = beta[3] * w[j][k] * policy[3] * social[t] / total_pop[k]
                elif compartment[n] in ['BI']:
                    for k in range(group):
                        betaSE[n][k] = beta[3] * w[j][k] * policy[3] * social[t] / total_pop[k]
                elif compartment[n] in ['BI2']:
                    for k in range(group):
                        betaSE[n][k] = beta[3] * w[j][k] * policy[3] * social[t] / total_pop[k]
            user.set_transition_compartment_dp(flow=('BS2', 'BE2', j, t), beta=betaSE)

    '''=============== Step 5. Define other constraints  ====================='''
    print('step.5')
    for t in range(period):
        user.custom_lp(fvar=[('VS', 'BS', j, t) for j in range(group)]\
                        + [('BS', 'BS2', j, t) for j in range(group)],
                       fcoef=[4 for j in range(group)]
                             + [-1 for j in range(group)],
                       sense='L',
                       rhs=0,
                       name='capacity_vacc1.' + str(t))
    # for t in range(period):
    #     user.custom_lp(fvar=[('VS', 'BS', j, t) for j in range(group)],
    #                    fcoef=[1 for j in range(group)],
    #                    sense='L',
    #                    rhs=capacity_V*0.2,
    #                    name='capacity_vacc1.' + str(t))

    for t in range(period):
        user.custom_lp(fvar=[('VS', 'BS', j, t) for j in range(group)]\
                            +[('BS', 'BS2', j, t) for j in range(group)],
                       fcoef=[1 for j in range(group)]+[1 for j in range(group)],
                       sense='L',
                       rhs=capacity_V,
                       name='capacity_vacc2.' + str(t))

    for t in range(period-1):
        # if t != int(period/2)-1:
        # if t not in [6,13,20,27,34,41,48]:
        if t not in [13,27,41]:
        # if t not in [27]:
            for j in range(group):
                user.custom_lp(fvar=[('VS', 'BS', j, t), ('VS', 'BS', j, t+1)],
                               fcoef=[1, -1],
                               sense='E',
                               rhs=0,
                               name='vacc_VS.' + str(j) + '.' + str(t))
                user.custom_lp(fvar=[('BS', 'BS2', j, t), ('BS', 'BS2', j, t+1)],
                               fcoef=[1, -1],
                               sense='E',
                               rhs=0,
                               name='vacc_BS.' + str(j) + '.' + str(t))

    for t in range(period):
        user.custom_lp(fvar=[('H', j, t) for j in range(group)]
                            + [('VH', j, t) for j in range(group)]
                            + [('BH', j, t) for j in range(group)]
                            + [('BH2', j, t) for j in range(group)],
                       fcoef=[1 for j in range(group)] + [1 for j in range(group)]
                             + [1 for j in range(group)] + [1 for j in range(group)],
                       dvar=['Hmax'],
                       dcoef=[-1],
                       sense='L',
                       rhs=0,
                       name='capacity_hosp.' + str(t))
    '''=============== Step 6. Define objective  ====================='''
    print('step.6')
    user.set_objectivetype(model=user, sense='min')

    for j in range(group):
        for t in range(period):
            user.set_objective(model=user, state=('M', j, t), value=data.gdp/365/max(cd))
            user.set_objective(model=user, state=('VM', j, t), value=data.gdp/365/max(cd))
            user.set_objective(model=user, state=('BM', j, t), value=data.gdp/365/max(cd))
            user.set_objective(model=user, state=('BM2', j, t), value=data.gdp/365/max(cd))
            user.set_objective(model=user, state=('H', j, t), value=(data.gdp/365+0)/max(cd))
            user.set_objective(model=user, state=('VH', j, t), value=(data.gdp/365+0)/max(cd))
            user.set_objective(model=user, state=('BH', j, t), value=(data.gdp/365+0)/max(cd))
            user.set_objective(model=user, state=('BH2', j, t), value=(data.gdp/365+0)/max(cd))
        user.set_objective(model=user, state=('D', j, period-1), value=cd[j]/max(cd))
        user.set_objective(model=user, state=('VD', j, period-1), value=cd[j]/max(cd))
        user.set_objective(model=user, state=('BD', j, period-1), value=cd[j]/max(cd))
        user.set_objective(model=user, state=('BD2', j, period-1), value=cd[j]/max(cd))

    user.set_objective(model=user, state='Hmax', value=ch*period/max(cd))

    '''=============== Step 7. Solve and output  ====================='''
    user.set_solver(solver='gurobi')
    user.set_approximation(opt='SO')
    user.set_log_stream_SO(label=1, file='SVEIHR/log_SVEIHR_exp_'+str(totalcost))
    sol = user.solve(label='Expectation', ct=3600)
    (status, obj, Xopt, xopt, Xc) = sol
    print('obj=', obj)

    print(sum([Xopt[compartment.index('M')][j][t] * (data.gdp/365)/max(cd) \
    + Xopt[compartment.index('VM')][j][t] * (data.gdp/365)/max(cd)
    + Xopt[compartment.index('BM')][j][t] * (data.gdp/365)/max(cd)
    + Xopt[compartment.index('BM2')][j][t] * (data.gdp/365)/max(cd)
    + Xopt[compartment.index('H')][j][t] * (data.gdp/365+0)/max(cd) \
    + Xopt[compartment.index('VH')][j][t] * (data.gdp/365+0)/max(cd)
    + Xopt[compartment.index('BH')][j][t] * (data.gdp/365+0)/max(cd)
    + Xopt[compartment.index('BH2')][j][t] * (data.gdp/365+0)/max(cd)
               for j in range(group) for t in range(period)])
    + sum([Xopt[compartment.index('D')][j][period-1] * cd[j]/max(cd) \
    + Xopt[compartment.index('VD')][j][period-1] * cd[j]/max(cd)
    + Xopt[compartment.index('BD')][j][period-1] * cd[j]/max(cd)
    + Xopt[compartment.index('BD2')][j][period-1] * cd[j]/max(cd) for j in range(group)])
    + Xc['Hmax'] * ch*period/max(cd))

    Xsolution_exp = user.get_solution_compartment(compartment=compartment)
    xsolution_exp = user.get_solution_flow_variable(compartmentfrom=['VS','BS'],
                                                    compartmentto=['BS','BS2'])
    Xcustom_exp = user.get_solution_custom_variable(name=['Hmax'])
    print(Xcustom_exp)
    user.Solution_print(Z=Xsolution_exp, x=xsolution_exp, xc=Xcustom_exp)

    '''======================== Simulation =================================='''
    sample = 1000
    size = 10
    randomq = [[[1 + s * 0.1 / (sample-1) if (ii, i) in [(0, 1), (7, 8), (14, 15), (21, 22)]
                 else 1 for i in range(compart_num)] for ii in range(compart_num)] for s in range(sample)]
    x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)] for n in
         range(compart_num)]
    # x[compartment.index('S')][compartment.index('VS')] = [[xopt[compartment.index('S')][compartment.index('VS')][j][t]
    #      for t in range(period)] for j in range(group)]
    x[compartment.index('VS')][compartment.index('BS')] = [[xopt[compartment.index('VS')][compartment.index('BS')][j][t]
         for t in range(period)] for j in range(group)]
    x[compartment.index('BS')][compartment.index('BS2')] = [[xopt[compartment.index('BS')][compartment.index('BS2')][j][t]
         for t in range(period)] for j in range(group)]
    user.x = x

    for i in range(size):
        Xsim_exp = [0 for s in range(int(sample/size))]
        p = 0
        for s in range(int(i*sample/size), int((i+1)*sample/size)):
            print(s)
            Xsim_exp[p] = user.Prediction(opt=['D', None, randomq[s]])
            p += 1
            # Xsim_exp[s] = user.DeterPrediction(None, randomq[s])
        '''======================== save =================================='''
        lhs = []
        ilabel = [compartment.index('H'), compartment.index('VH'), compartment.index('BH'), compartment.index('BH2')]
        for t in range(period):
            lhs.append(['Hosp_Capacity.' + str(t), [sum([Xsim_exp[s][ii][j][t] for ii in ilabel for j in range(group)])
                                                    for s in range(int(sample/size))]])
        ilabel = [compartment.index('D'), compartment.index('VD'), compartment.index('BD'), compartment.index('BD2')]

        lhs.append(['Death.', [sum([Xsim_exp[s][ii][j][-1] for ii in ilabel for j in range(group)])
                               for s in range(int(sample/size))]])
        costdeath = [sum([Xsim_exp[s][ii][j][-1] * cd[j] for ii in ilabel for j in range(group)])
                               for s in range(int(sample/size))]
        lhs.append(['Deathcost.', costdeath])

        ilabel = [compartment.index('M'), compartment.index('VM'), compartment.index('BM'), compartment.index('BM2')]
        costisolation = [sum([Xsim_exp[s][ii][j][t] * data.gdp/365
                                            for ii in ilabel for j in range(group) for t in range(period)])
                               for s in range(int(sample/size))]
        lhs.append(['Isolationcost.', costisolation])

        ilabel = [compartment.index('H'), compartment.index('VH'), compartment.index('BH'), compartment.index('BH2')]
        costhos = [sum([Xsim_exp[s][ii][j][t] * (data.gdp/365+0)
                                            for ii in ilabel for j in range(group) for t in range(period)])
                               for s in range(int(sample/size))]
        lhs.append(['Hospitalcost.', costhos])

        lhs.append(['Wardcost.', [Xcustom_exp['Hmax'] * ch*period for s in range(int(sample/size))]])
        lhs.append(['Totalcost.', [Xcustom_exp['Hmax'] * ch*period + costdeath[s] + costisolation[s] + costhos[s]
                                   for s in range(int(sample/size))]])

        user.save_result(filename='SVEIHR/resultSVEIHR_exp_4layers_'+str(i)+'_'+str(totalcost),
                         compartment=Xsim_exp, dvar=x, custom=Xcustom_exp, lhs=lhs)
    '''=========================================================='''
    for i in range(size):
        Xsim_exp = [0 for s in range(int(sample/size))]
        p = 0
        for s in range(int(i*sample/size), int((i+1)*sample/size)):
            print(s)
            Xsim_exp[p] = user.Prediction(opt=['S'])
            p += 1
            # Xsim_exp[s] = user.DeterPrediction(None, randomq[s])
        '''======================== save =================================='''
        lhs = []
        ilabel = [compartment.index('H'), compartment.index('VH'), compartment.index('BH'), compartment.index('BH2')]
        for t in range(period):
            lhs.append(['Hosp_Capacity.' + str(t), [sum([Xsim_exp[s][ii][j][t] for ii in ilabel for j in range(group)])
                                                    for s in range(int(sample/size))]])
        ilabel = [compartment.index('D'), compartment.index('VD'), compartment.index('BD'), compartment.index('BD2')]

        lhs.append(['Death.', [sum([Xsim_exp[s][ii][j][-1] for ii in ilabel for j in range(group)])
                               for s in range(int(sample/size))]])
        costdeath = [sum([Xsim_exp[s][ii][j][-1] * cd[j] for ii in ilabel for j in range(group)])
                               for s in range(int(sample/size))]
        lhs.append(['Deathcost.', costdeath])

        ilabel = [compartment.index('M'), compartment.index('VM'), compartment.index('BM'), compartment.index('BM2')]
        costisolation = [sum([Xsim_exp[s][ii][j][t] * data.gdp/365
                                            for ii in ilabel for j in range(group) for t in range(period)])
                               for s in range(int(sample/size))]
        lhs.append(['Isolationcost.', costisolation])

        ilabel = [compartment.index('H'), compartment.index('VH'), compartment.index('BH'), compartment.index('BH2')]
        costhos = [sum([Xsim_exp[s][ii][j][t] * (data.gdp/365+0)
                                            for ii in ilabel for j in range(group) for t in range(period)])
                               for s in range(int(sample/size))]
        lhs.append(['Hospitalcost.', costhos])

        lhs.append(['Wardcost.', [Xcustom_exp['Hmax'] * ch*period for s in range(int(sample/size))]])
        lhs.append(['Totalcost.', [Xcustom_exp['Hmax'] * ch*period + costdeath[s] + costisolation[s] + costhos[s]
                                   for s in range(int(sample/size))]])

        user.save_result(filename='SVEIHR/resultSVEIHR_exp_stochastic_4layers_'+str(i)+'_'+str(totalcost),
                         compartment=Xsim_exp, dvar=xopt, custom=Xcustom_exp, lhs=lhs)

    '''======================== Simulation =================================='''
    # sample = 1000
    # randomq = [[[1 + s * 0.1 / sample if (ii, i) in [(4, 6), (11,13), (18, 20)] else 1 for i in range(compart_num)]
    #                 for ii in range(compart_num)] for s in range(sample)]
    # x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)] for n in
    #      range(compart_num)]
    # x[compartment.index('S')][compartment.index('VS')] = [[xopt[compartment.index('S')][compartment.index('VS')][j][t]
    #      for t in range(period)] for j in range(group)]
    # x[compartment.index('VS')][compartment.index('BS')] = [[xopt[compartment.index('VS')][compartment.index('BS')][j][t]
    #      for t in range(period)] for j in range(group)]
    # x[compartment.index('BS')][compartment.index('BS2')] = [[xopt[compartment.index('BS')][compartment.index('BS2')][j][t]
    #      for t in range(period)] for j in range(group)]
    # user.x = x
    # Xsim_exp = [0 for s in range(sample)]
    # for s in range(sample):
    #     print(s)
    #     Xsim_exp[s] = user.Prediction(opt=['D', None, randomq[s]])
    # '''======================== save =================================='''
    # lhs = []
    # ilabel = [compartment.index('H'), compartment.index('VH'), compartment.index('BH'), compartment.index('BH2')]
    # for t in range(period):
    #     lhs.append(['Hosp_Capacity.' + str(t), [sum([Xsim_exp[s][ii][j][t] for ii in ilabel for j in range(group)])
    #                                             for s in range(sample)]])
    # ilabel = [compartment.index('D'), compartment.index('VD'), compartment.index('BD'), compartment.index('BD2')]
    #
    # lhs.append(['Death.', [sum([Xsim_exp[s][ii][j][-1] for ii in ilabel for j in range(group)])
    #                        for s in range(sample)]])
    #
    # user.save_result(filename='SVEIHR/resultSVEIHR_exp_4layers_'+str(totalcost),
    #                  compartment=Xsim_exp, dvar=xopt, custom=Xcustom_exp, lhs=lhs)

    '''====================== robust model =========================='''
    print('Begin to solve robust model')
    # obj = 33.586389523637514
    # Xcustom_exp = {'Hmax':673.2040662324}
    theta = [1.02]
    Xoptr = [[] for i in range(len(theta))]
    Xcr = [[] for i in range(len(theta))]
    xoptr = [[] for i in range(len(theta))]
    objr = [0 for i in range(len(theta))]
    Xsim_ro = [[[] for s in range(sample)] for ii in range(len(theta))]
    for i in range(len(theta)):
        user.set_custom_variable(name='Hmax', xlower=0, xupper=Xcustom_exp['Hmax']*1.0)
        # user.custom_lp(fvar=[('M', j, t) for j in range(group) for t in range(period)]
        #                     + [('VM', j, t) for j in range(group) for t in range(period)]
        #                     + [('BM', j, t) for j in range(group) for t in range(period)]
        #                     + [('BM2', j, t) for j in range(group) for t in range(period)]
        #                     +[('H', j, t) for j in range(group) for t in range(period)]
        #                     + [('VH', j, t) for j in range(group) for t in range(period)]
        #                     + [('BH', j, t) for j in range(group) for t in range(period)]
        #                     + [('BH2', j, t) for j in range(group) for t in range(period)]
        #                     +[('D', j, period-1) for j in range(group)]
        #                     + [('VD', j, period-1) for j in range(group)]
        #                     + [('BD', j, period-1) for j in range(group)]
        #                     + [('BD2', j, period-1) for j in range(group)],
        #                fcoef=[data.gdp/365/max(cd) for j in range(group) for t in range(period)]
        #                      + [data.gdp/365/max(cd) for j in range(group) for t in range(period)]
        #                      + [data.gdp/365/max(cd) for j in range(group) for t in range(period)]
        #                      + [data.gdp/365/max(cd) for j in range(group) for t in range(period)]
        #                      + [data.gdp/365/max(cd) for j in range(group) for t in range(period)]
        #                      + [data.gdp/365/max(cd) for j in range(group) for t in range(period)]
        #                      + [data.gdp/365/max(cd) for j in range(group) for t in range(period)]
        #                      + [data.gdp/365/max(cd) for j in range(group) for t in range(period)]
        #                      + [cd[j]/max(cd) for j in range(group) for t in range(period)]
        #                      + [cd[j]/max(cd) for j in range(group) for t in range(period)]
        #                      + [cd[j]/max(cd) for j in range(group) for t in range(period)]
        #                      + [cd[j]/max(cd) for j in range(group) for t in range(period)],
        #                dvar=['Hmax'],
        #                dcoef=[ch*period/max(cd)],
        #                sense='L',
        #                rhs=obj,
        #                name='capacity_obj.')
        # user.set_initial_solution(compartment1='BS',compartment2='BS2',val=xopt[compartment.index('BS')][compartment.index('BS2')])
        user.set_solver(solver='gurobi')
        user.set_approximation(opt='SO')
        user.set_log_stream_SO(label=1,
                               file='SVEIHR/log_SVEIHR_ro_4layers'+'_'+str(totalcost)+'_'+str(round(theta[i],2)))
        solr = user.solve(label='Robustness', ct=300, target=obj*theta[i])
        if solr != 0:
            (status, objr[i], Xoptr[i], xoptr[i], Xcr[i]) = solr
            Xsolution = user.get_solution_compartment(compartment=compartment)
            xsolution = user.get_solution_flow_variable(compartmentfrom=['VS','BS'],
                                                    compartmentto=['BS','BS2'])
            Xcustom_ro = user.get_solution_custom_variable(name=['Hmax'])
            user.Solution_print(Z=Xsolution, x=xsolution, xc=Xcustom_ro)
        '''=============== Step 8. Simulation  ====================='''
        # randomq = [[[0.9 + s * 0.2 / sample if (ii, i) in [(0, 1), (7, 8), (14, 15), (21, 22)]
        #              else 1 for i in range(compart_num)] for ii in range(compart_num)] for s in range(sample)]
        x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)] for n in
             range(compart_num)]
        # x[compartment.index('S')][compartment.index('VS')] = \
        #     [[xoptr[i][compartment.index('S')][compartment.index('VS')][j][t]
        #       for t in range(period)] for j in range(group)]
        x[compartment.index('VS')][compartment.index('BS')] = \
            [[xoptr[i][compartment.index('VS')][compartment.index('BS')][j][t]
              for t in range(period)] for j in range(group)]
        x[compartment.index('BS')][compartment.index('BS2')] = \
            [[xoptr[i][compartment.index('BS')][compartment.index('BS2')][j][t]
              for t in range(period)] for j in range(group)]
        user.x = x
        for j in range(size):
            Xsim_ro[i] = [0 for s in range(int(sample/size))]
            p = 0
            for s in range(j * int(sample/size), (j+1) * int(sample/size)):
                print(s)
                Xsim_ro[i][p] = user.Prediction(opt=['D', None, randomq[s]])
                p += 1
            '''======================== save =================================='''
            lhs = []
            ilabel = [compartment.index('H'), compartment.index('VH'), compartment.index('BH'), compartment.index('BH2')]
            for t in range(period):
                lhs.append(['Hosp_Capacity.' + str(t), [sum([Xsim_ro[i][s][ii][j][t]
                                                             for ii in ilabel for j in range(group)])
                             for s in range(int(sample/size))]])
            ilabel = [compartment.index('D'), compartment.index('VD'), compartment.index('BD'), compartment.index('BD2')]

            lhs.append(['Death.', [sum([Xsim_ro[i][s][ii][j][-1] for ii in ilabel for j in range(group)])
                         for s in range(int(sample/size))]])
            costdeath = [sum([Xsim_ro[i][s][ii][j][-1] * cd[j] for ii in ilabel for j in range(group)])
                                       for s in range(int(sample / size))]
            lhs.append(['Deathcost.', costdeath])

            ilabel = [compartment.index('M'), compartment.index('VM'), compartment.index('BM'),
                      compartment.index('BM2')]
            costisolation = [sum([Xsim_ro[i][s][ii][j][t] * data.gdp / 365
                                                for ii in ilabel for j in range(group) for t in range(period)])
                                           for s in range(int(sample / size))]
            lhs.append(['Isolationcost.', costisolation])

            ilabel = [compartment.index('H'), compartment.index('VH'), compartment.index('BH'),
                      compartment.index('BH2')]
            costhos = [sum([Xsim_ro[i][s][ii][j][t] * (data.gdp / 365 +0)
                                               for ii in ilabel for j in range(group) for t in range(period)])
                                          for s in range(int(sample / size))]
            lhs.append(['Hospitalcost.', costhos])

            lhs.append(['Wardcost.', [Xcustom_ro['Hmax'] * ch * period for s in range(int(sample / size))]])
            lhs.append(['Totalcost.', [Xcustom_ro['Hmax'] * ch * period + costdeath[s] + costisolation[s] + costhos[s]
                                       for s in range(int(sample / size))]])

            # user.save_result(filename='SVEIHR/resultSVEIHR_ro_4layers_' + str(totalcost)+'_'+str(round(theta[i],2)),
            #                  compartment=Xsim_ro[i], dvar=xoptr[i], custom=Xcustom_ro, lhs=lhs)

            user.save_result(filename='SVEIHR/resultSVEIHR_ro_4layers_' + str(j) + '_' + str(totalcost) + '_'
                                      + str(round(theta[i],2)),
                             compartment=Xsim_ro[i], dvar=xoptr[i], custom=Xcustom_ro, lhs=lhs)
        '''=========================================================='''
        for j in range(size):
            Xsim_ro[i] = [0 for s in range(int(sample/size))]
            p = 0
            for s in range(j * int(sample/size), (j+1) * int(sample/size)):
                print(s)
                Xsim_ro[i][p] = user.Prediction(opt=['S'])
                p += 1
            '''======================== save =================================='''
            lhs = []
            ilabel = [compartment.index('H'), compartment.index('VH'), compartment.index('BH'), compartment.index('BH2')]
            for t in range(period):
                lhs.append(['Hosp_Capacity.' + str(t), [sum([Xsim_ro[i][s][ii][j][t]
                                                             for ii in ilabel for j in range(group)])
                             for s in range(int(sample/size))]])
            ilabel = [compartment.index('D'), compartment.index('VD'), compartment.index('BD'), compartment.index('BD2')]

            lhs.append(['Death.', [sum([Xsim_ro[i][s][ii][j][-1] for ii in ilabel for j in range(group)])
                         for s in range(int(sample/size))]])
            costdeath = [sum([Xsim_ro[i][s][ii][j][-1] * cd[j] for ii in ilabel for j in range(group)])
                                       for s in range(int(sample / size))]
            lhs.append(['Deathcost.', costdeath])

            ilabel = [compartment.index('M'), compartment.index('VM'), compartment.index('BM'),
                      compartment.index('BM2')]
            costisolation = [sum([Xsim_ro[i][s][ii][j][t] * data.gdp / 365
                                                for ii in ilabel for j in range(group) for t in range(period)])
                                           for s in range(int(sample / size))]
            lhs.append(['Isolationcost.', costisolation])

            ilabel = [compartment.index('H'), compartment.index('VH'), compartment.index('BH'),
                      compartment.index('BH2')]
            costhos = [sum([Xsim_ro[i][s][ii][j][t] * (data.gdp / 365 +0)
                                               for ii in ilabel for j in range(group) for t in range(period)])
                                          for s in range(int(sample / size))]
            lhs.append(['Hospitalcost.', costhos])

            lhs.append(['Wardcost.', [Xcustom_ro['Hmax'] * ch * period for s in range(int(sample / size))]])
            lhs.append(['Totalcost.', [Xcustom_ro['Hmax'] * ch * period + costdeath[s] + costisolation[s] + costhos[s]
                                       for s in range(int(sample / size))]])

            # user.save_result(filename='SVEIHR/resultSVEIHR_ro_4layers_' + str(totalcost)+'_'+str(round(theta[i],2)),
            #                  compartment=Xsim_ro[i], dvar=xoptr[i], custom=Xcustom_ro, lhs=lhs)

            user.save_result(filename='SVEIHR/resultSVEIHR_ro_stochastic_4layers_' + str(j) + '_' + str(totalcost) + '_'
                                      + str(round(theta[i],2)),
                             compartment=Xsim_ro[i], dvar=xoptr[i], custom=Xcustom_ro, lhs=lhs)
        # '''=============== Step 8. Simulation  ====================='''
        # sample = 1000
        # randomq = [[[1 + s * 0.1 / sample if (ii, i) in [(4, 6), (11, 13), (18, 20)] else 1 for i in range(compart_num)]
        #             for ii in range(compart_num)] for s in range(sample)]
        # x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)] for n in
        #      range(compart_num)]
        # x[0][7] = [[xoptr[i][0][7][j][t] for t in range(period)] for j in range(group)]
        # x[7][14] = [[xoptr[i][7][14][j][t] for t in range(period)] for j in range(group)]
        # user.x = x
        #
        # Xsim_ro[i] = [0 for s in range(sample)]
        # for s in range(sample):
        #     Xsim_ro[i][s] = user.Prediction(opt=['D', None, randomq[s]])
        #
        # '''======================== save =================================='''
        # lhs = []
        # for t in range(period):
        #     lhs.append(['Hosp_Capacity.' + str(t),
        #                 [sum([Xsim_ro[i][s][4][j][t] + Xsim_ro[i][s][11][j][t] + Xsim_ro[i][s][18][j][t] for j in
        #                       range(group)]) for s in range(sample)]])
        #
        # lhs.append(['Death.',
        #             [sum([Xsim_ro[i][s][6][j][-1] + Xsim_ro[i][s][13][j][-1] + Xsim_ro[i][s][20][j][-1]
        #                   for j in range(group)]) for s in range(sample)]])
        #
        # # lhs.append(['Hosp_budget.' + str(t),
        # #             [ch * sum([Xsim_ro[s][4][j][t] + Xsim_ro[s][11][j][t] + Xsim_ro[s][18][j][t] for j in range(group)
        # #                  for t in range(period)]) for s in range(sample)]])
        # #
        # # lhs.append(['objective', [sum([Xsim_ro[s][2][j][t] * cd[j] for j in range(group) for t in range(period)])
        # #                           + Xcustom_ro['Hmax'] * sum(cd) / group * 1.5 for s in range(sample)]])
        # user.save_result(filename='SVEIHR/resultSVEIHR_death_ro_' + str(totalcost)+'_'+str(round(theta[i],2)),
        #                  compartment=Xsim_ro[i], dvar=xoptr[i], custom=Xcustom_ro, lhs=lhs)


if __name__ == '__main__':
    data = DatasetSG_4layers.Coviddata_SG()
    period = data.period
    group = data.group

    beta = data.parameter['beta']
    policy = data.parameter['a']
    social = data.parameter['s']
    w = data.parameter['w']

    q_EtoI = [[data.parameter['alpha'] for t in range(period)] for j in range(group)]
    q_VEtoVI = [[data.parameter['alpha'] for t in range(period)] for j in range(group)]
    q_BEtoBI = [[data.parameter['alpha'] for t in range(period)] for j in range(group)]
    q_BE2toBI2 = [[data.parameter['alpha'] for t in range(period)] for j in range(group)]

    q_ItoM = [[data.parameter['tau'] * (1 - data.parameter['p'][j][0]) for t in range(period)] for j in range(group)]
    q_VItoVM = [[data.parameter['tau'] * (1 - data.parameter['p'][j][1]) for t in range(period)] for j in range(group)]
    q_BItoBM = [[data.parameter['tau'] * (1 - data.parameter['p'][j][2]) for t in range(period)] for j in range(group)]
    q_BI2toBM2 = [[data.parameter['tau'] * (1 - data.parameter['p'][j][3]) for t in range(period)] for j in
                  range(group)]

    q_ItoH = [[data.parameter['tau'] * data.parameter['p'][j][0] for t in range(period)] for j in range(group)]
    q_VItoVH = [[data.parameter['tau'] * data.parameter['p'][j][1] for t in range(period)] for j in range(group)]
    q_BItoBH = [[data.parameter['tau'] * data.parameter['p'][j][2] for t in range(period)] for j in range(group)]
    q_BI2toBH2 = [[data.parameter['tau'] * data.parameter['p'][j][3] for t in range(period)] for j in range(group)]

    q_HtoR = [[data.parameter['lambda'][0][j] * (1 - data.parameter['d'][j][0]) for t in range(period)] for j in
              range(group)]
    q_VHtoVR = [[data.parameter['lambda'][1][j] * (1 - data.parameter['d'][j][1]) for t in range(period)] for j in
                range(group)]
    q_BHtoBR = [[data.parameter['lambda'][2][j] * (1 - data.parameter['d'][j][2]) for t in range(period)] for j in
                range(group)]
    q_BH2toBR2 = [[data.parameter['lambda'][3][j] * (1 - data.parameter['d'][j][3]) for t in range(period)] for j in
                  range(group)]

    q_HtoD = [[data.parameter['lambda'][0][j] * data.parameter['d'][j][0] for t in range(period)] for j in range(group)]
    q_VHtoVD = [[data.parameter['lambda'][1][j] * data.parameter['d'][j][1] for t in range(period)] for j in
                range(group)]
    q_BHtoBD = [[data.parameter['lambda'][2][j] * data.parameter['d'][j][2] for t in range(period)] for j in
                range(group)]
    q_BH2toBD2 = [[data.parameter['lambda'][3][j] * data.parameter['d'][j][3] for t in range(period)] for j in
                  range(group)]

    ch = data.cost_hosp
    cd = data.cost_death
    crd = data.cost_risk_death
    cv = data.cost_vacc

    for i in [5000]:
        model(totalcost=[i])
    # for i in [20000]:,6000,7000,8000,9000,4000,3000,2000,1000
    #     for j in [1000]:
    #         if (i,j) not in [(10000,700), (10000,800), (10000,900)]:
    #             model(totalcost=[i,j])
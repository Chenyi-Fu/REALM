from __future__ import print_function
import random
from math import exp
import REALM
import DatasetMD

data = DatasetMD.Coviddata_MD()
period = 30
# period_num = 4
# period_len = 7
totalperiod = 100
group = data.group
Various_lambdah = [2, 1.5, 1, 0.9, 0.7, 0.6, 0.5, 0.4]
Various_lambdaicu = [2, 1.5, 1, 0.9, 0.7, 0.6, 0.5, 0.4]
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

print([pm[j][-1] for j in range(group)])
print([pq[j][-1] for j in range(group)])
# for t in range(period):
#     print([q_ItoD[j][t] for j in range(group)])
#     print([q_QtoD[j][t] for j in range(group)])
#     print([q_HtoD[j][t] for j in range(group)])
#     print([q_UtoD[j][t] for j in range(group)])
#     print()
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

    # user.set_types(types = types) # Expectation or Robustness
    user.set_group(group=group)
    user.set_time(time=period)

    '''=============== Step 2. Define name of compartment and population ====================='''
    compartment = ['S','E','I','M','Q','H','U','R','D']
    compart_num = len(compartment)
    user.set_all_compartment(name=compartment)
    population = data.compartment_population
    total_pop = data.group_population
    # print(population)
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

    user.set_custom_variable(name='Hmax', xlower=0, xupper=0)
    user.set_custom_variable(name='Umax', xlower=0, xupper=0)
    for j in range(group):
        for s in range(S):
            user.set_custom_variable(name='u.' + str(j) + '.' + str(s), types='B', xupper=1, xlower=0)
            # if s == S-1:
            #     user.set_custom_variable(name='u.' + str(j) + '.' + str(s), xupper=1, xlower=1)
            # else:
            #     user.set_custom_variable(name='u.' + str(j) + '.' + str(s), xupper=0, xlower=0)
            for t in range(period):
                # user.set_custom_variable(name='C.'+str(j)+'.'+str(t), xupper=None, xlower=0)
                user.set_custom_variable(name='z.' + str(j) + '.' + str(t) + '.' + str(s), xlower=0)
                user.set_custom_variable(name='phi1.' + str(j) + '.' + str(t) + '.'+str(s), xlower=-1e20)
                user.set_custom_variable(name='phi2.' + str(j) + '.' + str(t) + '.'+str(s), xlower=0)

    print('step.3')
    '''====== Step 4.1 Define decision-independent transition between compartment and compartment  ========='''
    # q_StoE = [[beta[t]*w[j] * sum([w[k] * I[k][t] / total_pop[k] for k in range(group)]) for t in range(period)] for j in range(group)]
    # user.set_transition_compartment(n='S', m='E', prob=[[q_StoE[j][t] for t in range(period)] for j in range(group)], opt = None)
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
    # for t in range(period):
    #     for j in range(group):
    #         user.custom_lp(fvar=[('I', 'M', j, t)],
    #                        fcoef=[1],
    #                        sense='E',
    #                        rhs=ItoM[j][t],
    #                        name='priorIM.' + str(j) + '.' + str(t))
    #         user.custom_lp(fvar=[('I', 'Q', j, t)],
    #                        fcoef=[1],
    #                        sense='E',
    #                        rhs=ItoQ[j][t],
    #                        name='priorIQ.' + str(j) + '.' + str(t))
    # for t in range(period):
    #     for j in range(group):
    #         user.custom_lp(fvar=[('Q', 'U', j, t)],
    #                        fcoef=[1],
    #                        sense='E',
    #                        rhs=QtoU[j][t],
    #                        name='priorQU.' + str(j) + '.' + str(t))
    #         user.custom_lp(fvar=[('Q', 'H', j, t)],
    #                        fcoef=[1],
    #                        sense='E',
    #                        rhs=QtoH[j][t],
    #                        name='priorQH.' + str(j) + '.' + str(t))
    #
    # us = [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    # [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
    # [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    # [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    # [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    # [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    # [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
    # [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]
    #
    # for s in range(S):
    #     for j in range(group):
    #         user.custom_lp(dvar=['u.' + str(j) + '.' + str(s)],
    #                        dcoef=[1],
    #                        sense='E',
    #                        rhs=us[j][s],
    #                        name='prioru.' + str(j) + '.' + str(s))

    # zs = [[[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    # [169.3957860857372, 166.23511774031383, 176.09538928196866, 186.31553236024098, 197.13621185955796, 208.83027199702337, 221.5171914050908, 235.24919393169128, 250.05066412981174, 265.9327591254017, 282.89931574555317, 300.9480210527813, 320.0767481802033, 340.2702659114708, 361.50282763516594, 383.7395418564804, 406.93257613836977, 431.0119780096601, 455.8855997622922, 482.41895141997264, 510.33587600790594, 539.2395951364263, 564.5452365372108, 589.3658707322429, 614.151898435647, 638.4551657725876, 0.00010311289671765416, 710.8012105257757, 712.5490249452253, 239.59127689193264],
    # [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    # [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    # [236.0587108387026, 227.46902209868912, 238.3700940337656, 251.40555027905324, 265.8817979103364, 281.7519275337484, 299.0334741095686, 317.7545090547754, 337.9399730217103, 359.60805731162793, 382.7705377068421, 407.4322663475403, 431.92888853290333, 455.39599958885424, 479.81539322615896, 505.0066527182441, 48.35395635749987, 7.925326086902475e-06, 5.098544323984111e-06, 1.899818057222968e-06, 1.1846213802371285e-06, 1.0436262504703972e-06, 1.0704852542440199e-06, 1.1957808449208398e-06, 1.4698679685588835e-06, 2.134142653389168e-06, 4.9453486271149485e-06, 891.1482214429215, 885.329828384805, 323.3573049696869],
    # [168.5286984298582, 163.49643342694506, 170.776313688477, 179.62670413539016, 189.47549823533188, 200.2156992528679, 211.8162731591134, 224.26787973354226, 237.56410679752497, 251.69481993485334, 266.64475529189735, 282.39283358318045, 298.9168204981283, 316.18367795669064, 334.1515615101918, 352.77152427003244, 371.98555557772363, 200.97152620093425, 92.84938195171041, 2.419964484110387e-06, 1.0066068452057811e-06, 7.046977844729862e-07, 6.875110302555954e-07, 7.027610722120372e-07, 7.447134260088197e-07, 9.254521841252383e-07, 2.1647946204547135e-06, 660.86621318772, 652.3696755599965, 226.54476724790396],
    # [104.01218538652392, 101.22234534800742, 105.46168949784317, 111.01880226754373, 117.40032489928873, 124.44442652797545, 132.08629331682056, 140.30359570289693, 149.09091602077217, 158.44779940219289, 168.37378086064004, 178.86604821409003, 189.92129946864432, 201.52990330073123, 213.67646783367277, 226.34064012277813, 239.49585579520746, 253.10554476798993, 267.12297369709995, 281.9117517591101, 297.4065050490295, 238.36016502249402, 205.4701443657331, 170.99836531693495, 135.37078918386442, 98.98610233400518, 62.24468357911046, 454.5525687405128, 448.36036165198203, 146.37372960975108],
    # [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    # ],
    #     [[426.43375782729345, 420.72154669062195, 447.28193905933404, 473.6675500557293, 501.22716216244623, 530.9187186040792, 563.1216208754279, 597.9863784415849, 635.5775752429251, 669.9791164214979, 706.2803382130048, 744.2395928378063, 783.6969851519455, 824.5023856800753, 866.495865760428, 909.5067247561326, 953.3475389149177, 997.8005519701474, 1042.62302595491, 1089.258362072006, 1137.0511424525414, 1185.2102248251404, 1233.2004280155545, 1281.882880059118, 1329.8763348271536, 1376.4803725092127, 1421.1502386820302, 1463.344028015742, 1504.5383858252612, 627.2886946006864],
    # [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    # [210.90251862449986, 203.45506697539716, 213.7806723925254, 225.60080252834734, 238.53539557225596, 252.6543700294038, 268.0126549221158, 284.6453291151359, 302.5740652079969, 321.81069234589074, 342.3598810674445, 364.21928620776305, 387.3867563701769, 411.8444964686211, 437.5619619095743, 464.49748791759004, 492.59378855779244, 521.7669679607366, 551.9062832804651, 584.0423159226771, 617.8590549041331, 652.882784514221, 688.766726396113, 726.1965335058451, 764.1765205082123, 3.3673964261912445e-05, 2.56636303372268e-05, 938.9510029515695, 926.9485794258635, 346.359394115164],
    # [170.09474749626315, 163.15584514213435, 170.8172536949888, 180.06040094800704, 190.34721658054315, 201.63449458379785, 213.93059651497524, 227.25263816588816, 241.61573045243054, 257.03006592038037, 273.50109971965463, 291.029101498063, 309.61468226457066, 329.24639198364497, 349.90287752093883, 371.55419075593886, 394.1582949756191, 417.65239359941586, 441.95263386478223, 467.87607614804716, 2.4805582016698973e-05, 8.387468658170726e-06, 5.748458895652365e-06, 4.895378878891343e-06, 4.8211096005532484e-06, 5.712056522513977e-06, 1.0142216351430595e-05, 764.8387505985146, 750.1063375351197, 266.88200541918553],
    # [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    # [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    # [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    # [73.36456304587409, 73.7131114790902, 77.92969019416527, 83.10943712877702, 88.96284284632203, 95.37060203760318, 102.27495948613692, 109.65306055791052, 117.50148699545412, 125.82687153052858, 134.64038000404548, 143.9542708649587, 153.78167054242027, 164.13207249898295, 175.01060581671493, 186.41792938590538, 198.34903847695657, 210.7902038460428, 223.71793762026326, 237.39426343873947, 251.8309526152461, 266.94304514551584, 282.6341886391397, 299.05901800660143, 316.01906534779835, 333.32800455983937, 350.7980771062083, 368.22357892778325, 385.8471141797807, 123.58255482543102]
    # ]]

    # for s in [1,2]:
    #     for j in range(group):
    #         for t in range(period):
    #             user.custom_lp(dvar=['z.' + str(j) + '.' + str(t) + '.' + str(s)],
    #                            dcoef=[1],
    #                            sense='E',
    #                            rhs=zs[s-1][j][t],
    #                            name='priorz.' + str(j) + '.' + str(t) + '.' + str(s))
    # test constraints
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
                # user.custom_qp(dvarqp=['phi1.' + str(j) + '.' + str(t) + '.' + str(s),
                #                        'phi2.' + str(j) + '.' + str(t) + '.' + str(s),
                #                        'u.' + str(j) + '.' + str(s)],
                #                dcoefqp=[0.25, -0.25, sum(total_pop)],
                #                sense='L',
                #                rhs=sum(total_pop) - Test_capacity[t] * Test_probability[s] * bhat[j],
                #                name='SOC.' + str(j) + '.' + str(t) + '.' + str(s))
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

    # capacity constraints
    # for t in range(1,period):
    #     user.custom_lp(fvar=[('H', j, t) for j in range(group)],
    #                    fcoef=[1 for j in range(group)],
    #                    dvar=['Hmax'],
    #                    dcoef=[-1],
    #                    sense='L',
    #                    rhs=H0,
    #                    name='Hospitalization.' + str(t))
    #     user.custom_lp(fvar=[('U', j, t) for j in range(group)],
    #                    fcoef=[1 for j in range(group)],
    #                    dvar=['Umax'],
    #                    dcoef=[-1],
    #                    sense='L',
    #                    rhs=ICU0,
    #                    name='ICU.' + str(t))

    # test allocation constraints
    for j in range(group):
        user.custom_lp(dvar=['u.' + str(j) + '.' + str(s) for s in range(S)],
                       dcoef=[1 for s in range(S)],
                       sense='E',
                       rhs=1,
                       name='Uequal.' + str(j))

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
                   dvar=['u.' + str(j) + '.' + str(s) for j in range(group) for s in range(S)],
                   dcoef=[ct*Test_capacity[t] * Test_probability[s] * (period-1) for j in range(group) for s in range(S)],
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
        for t in range(period):
            user.set_objective(state=('I', j, t), value=ci)
        user.set_objective(state=('D', j, period - 1), value=cd)

    # for j in range(group):
    #     for t in range(period-1):
    #         user.set_objective(state=('Q', 'H', j, t), value=cx/1000)
    #         user.set_objective(state=('Q', 'U', j, t), value=cy/1000)
    #     for s in range(S):
    #         user.set_objective(state='u.' + str(j) + '.' + str(s),
    #                            value=ct * Test_capacity[t] * Test_probability[s] * (period - 1)/1000)

    '''=============== Step 7. Solve and output  ====================='''
    user.set_solver(solver='cplex')
    user.set_approximation(opt='SO')
    user.set_log_stream_SO(label=1, file='SEIHR/log_SEIHR_exp_'+str(totalcost))
    sol = user.solve(label='Expectation', ct=3600)
    (status, obj, Xopt, xopt, Xc) = sol
    print('obj=', obj)
    Xsolution_exp = user.get_solution_compartment(compartment=compartment)
    xsolution_exp = user.get_solution_flow_variable(compartmentfrom=['I','I','Q','Q'],compartmentto=['M','Q','H','U'])
    Xcustom_exp = user.get_solution_custom_variable(name=['Hmax', 'Umax']
                                                  + ['u.' + str(j) + '.' + str(s) for s in range(S) for j in range(group)]
                                                  + ['phi1.' + str(j) + '.' + str(t) + '.'+str(s)
                                                     for s in range(S) for j in range(group) for t in range(period)]
                                                  + ['phi2.' + str(j) + '.' + str(t) + '.'+str(s)
                                                     for s in range(S) for j in range(group) for t in range(period)]
                                                  + ['z.' + str(j) + '.' + str(t) + '.' + str(s)
                                                     for s in range(S) for j in range(group) for t in range(period)])
    user.Solution_print(X=Xsolution_exp,x=xsolution_exp,xc=None)

    '''======================== Simulation =================================='''
    sample = 1000
    # randomq = user.Sample_generation(n=['S','SS'], m=['E','EE'], lb=[1,1], ub=[1.01,1.01], size=sample)
    randomq = [[[1 + s * 0.1 / sample if (ii, i) in [(0, 1)] else 1 for i in range(compart_num)]
                    for ii in range(compart_num)] for s in range(sample)]
    # C = [[sum([Test_capacity[t] * Test_probability[s] * Xcustom['u.' + str(j) + '.' + str(s)] for s in range(S)])
    #       for t in range(period)] for j in range(group)]
    # for j in range(group):
    #     for t in range(period):
    #         betaIH = [[1 if compartment[n] == 'I' and k == j else 0
    #                    for k in range(group)] for n in range(len(population))]
    #         user.set_transition_compartment_dp(flow=('I', 'M', j, t), alpha=C[j][t]*pm[j][t],
    #                                            betadeno=betaIH, alphadeno=bhat[j], ub=attractive*pm[j][t])
    #         user.set_transition_compartment_dp(flow=('I', 'Q', j, t), alpha=C[j][t]*pq[j][t],
    #                                            betadeno=betaIH, alphadeno=bhat[j], ub=attractive*pq[j][t])
    #         user.set_transition_compartment_dp(flow=('I', 'D', j, t), alpha=(1 - C[j][t] * pm[j][t] - C[j][t] * pq[j][t]) * q_ItoD[j][t],
    #                                            betadeno=betaIH, alphadeno=bhat[j])
    #         user.set_transition_compartment_dp(flow=('I', 'R', j, t), alpha=(1 - C[j][t] * pm[j][t] - C[j][t] * pq[j][t]) * q_ItoR[j][t],
    #                                            betadeno=betaIH, alphadeno=bhat[j])
    x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)] for n in range(compart_num)]
    x[2][3] = [[xopt[2][3][j][t] for t in range(period)] for j in range(group)]
    x[2][4] = [[xopt[2][4][j][t] for t in range(period)] for j in range(group)]
    x[4][5] = [[xopt[4][5][j][t] for t in range(period)] for j in range(group)]
    x[4][6] = [[xopt[4][6][j][t] for t in range(period)] for j in range(group)]
    user.x = x
    Xsim_exp = [0 for s in range(sample)]
    for s in range(sample):
        Xsim_exp[s] = user.Prediction(opt=['D', None, randomq[s]])
    # user.set_transition_compartment_dp(flow=('I', 'M'), opt='delete')
    # user.set_transition_compartment_dp(flow=('I', 'Q'), opt='delete')
    # user.set_transition_compartment_dp(flow=('I', 'R'), opt='delete')
    # user.set_transition_compartment_dp(flow=('I', 'D'), opt='delete')
    # user.set_transition_compartment_dp(flow=('I', 'I'), opt='delete')
    # print(user.compartment_info[2]['tocpm'][7]['q_idp'])
    # print(user.compartment_info[2]['tocpm'][8]['q_idp'])
    # print([(Xsim_exp[0][2][0][t+1], Xsim_exp[0][2][0][t] * (1-q_ItoR[0][t]-q_ItoD[0][t]), Xsim_exp[0][1][0][t] * q_EtoI[0][t],
    #         min(Xsim_exp[0][2][0][t] * 9000 / (Xsim_exp[0][2][0][t]+bhat[j]), attractive*Xsim_exp[0][2][0][t]))
    #        for t in range(period-1)])
    # print([sum([Xsim_exp[0][4][j][t] for j in range(group)]) for t in range(period)])
    # print([sum([Xsim_exp[0][5][j][t] for j in range(group)]) for t in range(period)])
    '''======================== save =================================='''
    # lhs = [['totalcost',
    #         [sum([xsolution[('I','M')][j][t] * cts + xsolution[('I','Q')][j][t] * cts \
    #               + xsolution[('Q', 'H')][j][t] * cx + xsolution[('Q', 'U')][j][t] * cy
    #              for j in range(group) for t in range(period)]) + Xcustom['Hmax'] * ch + Xcustom['Umax'] * cu ]]
    #        ]
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

    '''======================== Simulation =================================='''
    # sample = 1000
    # # randomq = user.Sample_generation(n=['S','SS'], m=['E','EE'], lb=[1,1], ub=[1.01,1.01], size=sample)
    # randomq = [[[1 + s * 0.1 / sample if (ii, i) in [(2, 8),(4, 8),(5, 8),(6, 8)] else 1 for i in range(compart_num)]
    #             for ii in range(compart_num)] for s in range(sample)]
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
    # x[2][3] = [[xopt[2][3][j][t] for t in range(period)] for j in range(group)]
    # x[2][4] = [[xopt[2][4][j][t] for t in range(period)] for j in range(group)]
    # x[4][5] = [[xopt[4][5][j][t] for t in range(period)] for j in range(group)]
    # x[4][6] = [[xopt[4][6][j][t] for t in range(period)] for j in range(group)]
    # user.x = x
    # Xsim_exp = [0 for s in range(sample)]
    # for s in range(sample):
    #     Xsim_exp[s] = user.Prediction(opt=['D', None, randomq[s]])
    # '''======================== save =================================='''
    # # lhs = [['totalcost',
    # #         [sum([xsolution[('I', 'M')][j][t] * cts + xsolution[('I', 'Q')][j][t] * cts \
    # #               + xsolution[('Q', 'H')][j][t] * cx + xsolution[('Q', 'U')][j][t] * cy
    # #               for j in range(group) for t in range(period)]) + Xcustom['Hmax'] * ch + Xcustom['Umax'] * cu]]
    # #        ]
    # lhs = [['totalcost',
    #         [sum([xsolution_exp[('Q', 'H')][j][t] * cx + xsolution_exp[('Q', 'U')][j][t] * cy
    #               for j in range(group) for t in range(period)]) + Xcustom_exp['Hmax'] * ch + Xcustom_exp['Umax'] * cu
    #          + sum([Xcustom_exp['u.' + str(j) + '.' + str(s)] * ct * Test_capacity[t] * Test_probability[s] * (period-1)
    #          for j in range(group) for s in range(S)]) ]]
    #        ]
    # for t in range(period):
    #     lhs.append(['Hospitalization.' + str(t), [sum([Xsim_exp[s][5][j][t] for j in range(group)])
    #                                               for s in range(sample)]])
    # for t in range(period):
    #     lhs.append(['ICU.' + str(t), [sum([Xsim_exp[s][6][j][t] for j in range(group)])
    #                                   for s in range(sample)]])
    # lhs.append(['objective', [sum([Xsim_exp[s][2][j][t] * ci + Xsim_exp[s][-1][j][t] * cd for j in range(group)])
    #                           for s in range(sample)]])
    # user.save_result(filename='SEIHR/resultSEIHR_death_exp_' + str(totalcost),
    #                  compartment=Xsim_exp, dvar=xopt, custom=Xcustom_exp, lhs=lhs)







    '''====================== robust model =========================='''
    print('Begin to solve robust model')
    theta = [1+0.01*(1+i) for i in range(4)]
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
    # theta = [1+0.02*(1+i) for i in range(10)]
    Xoptr2 = [[] for i in range(len(theta))]
    xoptr2 = [[] for i in range(len(theta))]
    Xcr2 = [[] for i in range(len(theta))]
    objr2 = [0 for i in range(len(theta))]

    Xsim_ro2 = [[[] for s in range(sample)] for ii in range(len(theta))]

    for i in range(len(theta)):
        '''=============== Step 5. Define other constraints  ====================='''
        print('step.5')
        user.set_time(time=period)
        # user.custom_lp(option='clear')
        user.custom_lp(fvar=[('I', j, t) for j in range(group) for t in range(period)],
                       fcoef=[1 for j in range(group) for t in range(period)],
                       sense='L',
                       rhs=sum([Xsolution_exp['I'][j][t] for j in range(group) for t in range(period)])*theta[i],
                       name='Infection_number',
                       label='Expectation')
        # user.custom_lp(fvar=[('D', j, period-1) for j in range(group)],
        #                fcoef=[1 for j in range(group)],
        #                sense='L',
        #                rhs=sum([Xsolution_exp['D'][j][period-1] for j in range(group)])*theta[i],
        #                name='Death_number')
        user.set_solver(solver='cplex')
        user.set_log_stream_SO(label=1,
                               file='SEIHR/log_SEIHR_multiro' + '_' + str(totalcost) + '_' + str(round(theta[i], 2)))
        solr2 = user.solve(label='Robustness', ct=3600, target=obj * theta[i])
        if solr2 != 0:
            (status, objr2[i], Xoptr2[i], xoptr2[i], Xcr2[i]) = solr2
            Xsolution = user.get_solution_compartment(compartment=compartment)
            xsolution = user.get_solution_flow_variable(compartmentfrom=['I', 'I', 'Q', 'Q'],
                                                        compartmentto=['M', 'Q', 'H', 'U'])
            Xcustom = user.get_solution_custom_variable(name=['Hmax', 'Umax']
                                                             + ['u.' + str(j) + '.' + str(s) for s in range(S) for j
                                                                in range(group)]
                                                             + ['phi1.' + str(j) + '.' + str(t) + '.' + str(s)
                                                                for s in range(S) for j in range(group) for t in
                                                                range(period)]
                                                             + ['phi2.' + str(j) + '.' + str(t) + '.' + str(s)
                                                                for s in range(S) for j in range(group) for t in
                                                                range(period)]
                                                             + ['z.' + str(j) + '.' + str(t) + '.' + str(s)
                                                                for s in range(S) for j in range(group) for t in
                                                                range(period)])
            user.Solution_print(X=Xsolution, x=xsolution, xc=None)
        '''=============== Step 8. Simulation  ====================='''
        randomq = [[[1 + s * 0.1 / sample if (ii, i) in [(0, 1)] else 1 for i in range(compart_num)]
                    for ii in range(compart_num)] for s in range(sample)]
        x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)]
             for n in range(compart_num)]
        x[2][3] = [[xoptr2[i][2][3][j][t] for t in range(period)] for j in range(group)]
        x[2][4] = [[xoptr2[i][2][4][j][t] for t in range(period)] for j in range(group)]
        x[4][5] = [[xoptr2[i][4][5][j][t] for t in range(period)] for j in range(group)]
        x[4][6] = [[xoptr2[i][4][6][j][t] for t in range(period)] for j in range(group)]
        user.x = x
        print(user.x)
        Xsim_ro2[i] = [0 for s in range(sample)]
        for s in range(sample):
            Xsim_ro2[i][s] = user.Prediction(opt=['D', None, randomq[s]])

        '''======================== save =================================='''
        # lhs = [['totalcost',
        #         [sum([xsolution[('I', 'M')][j][t] * cts + xsolution[('I', 'Q')][j][t] * cts \
        #               + xsolution[('Q', 'H')][j][t] * cx + xsolution[('Q', 'U')][j][t] * cy
        #               for j in range(group) for t in range(period)]) + Xcustom['Hmax'] * ch + Xcustom['Umax'] * cu]]
        #        ]
        lhs = [['totalcost',
                [sum([xsolution[('Q', 'H')][j][t] * cx + xsolution[('Q', 'U')][j][t] * cy
                      for j in range(group) for t in range(period)]) + Xcustom['Hmax'] * ch + Xcustom['Umax'] * cu
                 + sum(
                    [Xcustom['u.' + str(j) + '.' + str(s)] * ct * Test_capacity[t] * Test_probability[s] * (
                    period - 1)
                     for j in range(group) for s in range(S)])]]
               ]
        for t in range(period):
            lhs.append(['Hospitalization.' + str(t), [sum([Xsim_ro2[i][s][5][j][t] for j in range(group)])
                                                      for s in range(sample)]])
        for t in range(period):
            lhs.append(['ICU.' + str(t), [sum([Xsim_ro2[i][s][6][j][t] for j in range(group)])
                                          for s in range(sample)]])
        lhs.append(['objective', [sum([Xsim_ro2[i][s][2][j][t] * ci + Xsim_ro2[i][s][-1][j][t] * cd
                                       for j in range(group)]) for s in range(sample)]])

        user.save_result(filename='SEIHR/resultSEIHR_multiro' + '_' + str(totalcost) + '_' + str(round(theta[i], 2)),
                         compartment=Xsim_ro2[i], dvar=xoptr2[i], custom=Xcustom, lhs=lhs)
        '''=============== Step 8. Simulation  ====================='''
        # randomq = [[[1 + s * 0.1 / sample if (ii, i) in [(2, 8), (4, 8), (5, 8), (6, 8)] else 1 for i in
        #              range(compart_num)]
        #             for ii in range(compart_num)] for s in range(sample)]
        # x = [[[[0 for t in range(period)] for j in range(group)] for m in range(compart_num)]
        #      for n in range(compart_num)]
        # x[2][3] = [[xoptr2[i][2][3][j][t] for t in range(period)] for j in range(group)]
        # x[2][4] = [[xoptr2[i][2][4][j][t] for t in range(period)] for j in range(group)]
        # x[4][5] = [[xoptr2[i][4][5][j][t] for t in range(period)] for j in range(group)]
        # x[4][6] = [[xoptr2[i][4][6][j][t] for t in range(period)] for j in range(group)]
        # user.x = x
        # print(user.x)
        # Xsim_ro2[i] = [0 for s in range(sample)]
        # for s in range(sample):
        #     Xsim_ro2[i][s] = user.Prediction(opt=['D', None, randomq[s]])
        # '''======================== save =================================='''
        # # lhs = [['totalcost',
        # #         [sum([xsolution[('I', 'M')][j][t] * cts + xsolution[('I', 'Q')][j][t] * cts \
        # #               + xsolution[('Q', 'H')][j][t] * cx + xsolution[('Q', 'U')][j][t] * cy
        # #               for j in range(group) for t in range(period)]) + Xcustom['Hmax'] * ch + Xcustom['Umax'] * cu]]
        # #        ]
        # lhs = [['totalcost',
        #         [sum([xsolution[('Q', 'H')][j][t] * cx + xsolution[('Q', 'U')][j][t] * cy
        #               for j in range(group) for t in range(period)]) + Xcustom['Hmax'] * ch + Xcustom['Umax'] * cu
        #          + sum(
        #             [Xcustom['u.' + str(j) + '.' + str(s)] * ct * Test_capacity[t] * Test_probability[s] * (
        #             period - 1)
        #              for j in range(group) for s in range(S)])]]
        #        ]
        # for t in range(period):
        #     lhs.append(['Hospitalization.' + str(t), [sum([Xsim_ro2[i][s][5][j][t] for j in range(group)])
        #                                               for s in range(sample)]])
        # for t in range(period):
        #     lhs.append(['ICU.' + str(t), [sum([Xsim_ro2[i][s][6][j][t] for j in range(group)])
        #                                   for s in range(sample)]])
        # lhs.append(['objective', [sum([Xsim_ro2[i][s][2][j][t] * ci + Xsim_ro2[i][s][-1][j][t] * cd
        #                                for j in range(group)]) for s in range(sample)]])
        # user.save_result(
        #     filename='SEIHR/resultSEIHR_death_multiro' + '_' + str(totalcost) + '_' + str(round(theta[i], 2)),
        #     compartment=Xsim_ro2[i], dvar=xoptr2[i], custom=Xcustom, lhs=lhs)

if __name__ == '__main__':
    for i in [10]:
        model(totalcost=(0 + i) * 1e+7)
    # for i in [500]:
    #     model(totalcost=i)
    # for i in [[550000,5500], [550000,5000], [500000,5000], [500000,4500]]:
    #     model(totalcost=i)
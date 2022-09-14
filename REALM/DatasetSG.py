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

class Coviddata_SG():
    def __init__(self):
        self.group = 10
        self.period = 30
        compartment = 21
        LOS = [[1.9544, 2.3698, 2.9677, 4.1774, 4.1405, 4.9579, 5.6474, 5.9455, 6.5570, 6.9782],
               [1.6000, 2.3333, 2.8619, 3.3555, 3.3958, 4.7778, 5.5690, 6.2801, 6.5108, 6.5317],
               [1.6000, 3.5000, 3.0460, 3.2578, 3.3589, 3.8148, 4.3255, 5.0031, 5.2380, 5.5851]]
        self.cost_hosp = 1000
        self.cost_vacc = 60
        a = [3, 9, 16, 25, 35, 45, 55, 60, 65, 75, 80]
        gdp = 59797
        self.cost_death = [a[i] * gdp for i in range(self.group)]
        self.cost_risk_death = [1.5 * self.cost_death[i] for i in range(self.group)]

        self.parameter = {'beta': [1.0747743949469097, 0.8842042896462581, 0.6958110596552376],
                          'w': [[0.198408689, 0.095931398, 0.075377296, 0.142649748, 0.253510478,
                                 0.108750269, 0.073748669, 0.035325003, 0.012413567, 0.003884884],
                                [0.042348141, 0.493675788, 0.109501385, 0.04957647, 0.144000922,
                                 0.098354883, 0.035741385, 0.019017155, 0.005522952, 0.002260919],
                                [0.008397092, 0.058551198, 0.630588168, 0.07494171, 0.076106627,
                                 0.10390807, 0.033625307, 0.008642874, 0.003740834, 0.001498121],
                                [0.021373504, 0.015515033, 0.120991862, 0.387622447, 0.195210558,
                                 0.147935256, 0.091932507, 0.013800924, 0.003935628, 0.001682282],
                                [0.039496562, 0.06440112, 0.091816649, 0.160731038, 0.325826066,
                                 0.187764259, 0.095721487, 0.0258699, 0.006377763, 0.001995156],
                                [0.019274355, 0.046058771, 0.156463795, 0.141086954, 0.213844453,
                                 0.276874658, 0.11243945, 0.022746854, 0.008147558, 0.003063151],
                                [0.02611464, 0.054701314, 0.140714499, 0.153019644, 0.173112167,
                                 0.194484994, 0.201579573, 0.041879498, 0.010127788, 0.004265883],
                                [0.046704657, 0.060300482, 0.093816922, 0.111580447, 0.186902292,
                                 0.14634986, 0.144322831, 0.170207258, 0.031640805, 0.008174445],
                                [0.04306802, 0.081068029, 0.139702188, 0.063191141, 0.127455505,
                                 0.156192544, 0.117513345, 0.117726399, 0.112221974, 0.041860854],
                                [0.046512359, 0.068146347, 0.17516909, 0.061971207, 0.120615413,
                                 0.15547829, 0.129778183, 0.080576935, 0.099652382, 0.062099794]],
                          'a': [53.7/56.48, 1, 1],
                          's': [1+0.5 *exp(-(t-1) ** 2 / (2 * 10**2)) for t in range(self.period)],
                          'alpha': 0.33,
                          'tau': 0.8,
                          'lambda': [[1/LOS[i][j] for j in range(len(LOS[i]))] for i in range(len(LOS))],
                          # 'lambda_dis': [[0.160383237, 0.157501168, 0.098068235, 0.078906372, 0.073765384, 0.063327621,
                          #                 0.051331983, 0.04276367, 0.035130083, 0.028197539, 0.024380745, 0.026094407,
                          #                 0.021031313, 0.017681882, 0.012307213, 0.012852469, 0.011294594, 0.008256738,
                          #                 0.008490419, 0.007711482, 0.006932544, 0.004907306, 0.003816794, 0.0049852,
                          #                 0.003583113, 0.002726281, 0.002570494, 0.002025238, 0.002570494, 0.02641],
                          #                [0.128832808, 0.142208202, 0.124542587, 0.10977918, 0.090536278, 0.080946372,
                          #                 0.055962145, 0.044164038, 0.040504732, 0.028264984, 0.022523659, 0.018044164,
                          #                 0.015205047, 0.011987382, 0.010157729, 0.007949527, 0.007507886, 0.005488959,
                          #                 0.006498423, 0.005362776, 0.004542587, 0.003028391, 0.002776025, 0.002208202,
                          #                 0.003028391, 0.001577287, 0.001892744, 0.001451104, 0.001829653, 0.02120],
                          #                [0.136098231, 0.176283703, 0.141421947, 0.111196978, 0.09067491, 0.082002404,
                          #                 0.057616349, 0.044049459, 0.038897476, 0.024729521, 0.01760261, 0.015971149,
                          #                 0.013652756, 0.008414906, 0.006955178, 0.005581315, 0.004808518, 0.003949854,
                          #                 0.003606388, 0.002833591, 0.002575992, 0.001202129, 0.001373862, 0.001116263,
                          #                 0.001717328, 0.000429332, 0.000601065, 0.000858664, 0.000515198, 0.00326]],
                          'lambda_dis': [
                              [[0.4099, 0.9172, 0.9715, 0.9843, 0.9893, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                           0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                           0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                           0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],
                                            [0.2926, 0.3859, 0.4307, 0.4569, 0.4740, 0.4862, 0.4952, 0.5022, 0.5077,
                                             0.5122, 0.5160, 0.5192, 0.5219, 0.5243, 0.5264, 0.5281, 0.5293, 0.5313,
                                             0.5328, 0.5330, 0.5346, 0.5355, 0.5357, 0.5368, 0.5380, 0.5393, 0.5396,
                                             0.5398, 0.5405, 0.5420, 0.5442, 0.5468, 0.5535, 0.5671, 0.6003, 0.6866, 1.0000],
                                            [0.2258, 0.3290, 0.3820, 0.4139, 0.4350, 0.4500, 0.4612, 0.4699, 0.4768,
                                             0.4825, 0.4872, 0.4911, 0.4945, 0.4974, 0.5000, 0.5019, 0.5045, 0.5063,
                                             0.5078, 0.5090, 0.5103, 0.5115, 0.5128, 0.5133, 0.5151, 0.5152, 0.5165,
                                             0.5180, 0.5186, 0.5193, 0.5225, 0.5262, 0.5332, 0.5498, 0.5855, 0.6759, 1.0000],
                                            [0.1811, 0.2820, 0.3371, 0.3711, 0.3940, 0.4105, 0.4228, 0.4325, 0.4402,
                                             0.4465, 0.4517, 0.4561, 0.4599, 0.4632, 0.4660, 0.4685, 0.4707, 0.4728,
                                             0.4751, 0.4764, 0.4780, 0.4798, 0.4802, 0.4813, 0.4830, 0.4839, 0.4844,
                                             0.4860, 0.4873, 0.4890, 0.4919, 0.4974, 0.5061, 0.5248, 0.5649, 0.6621, 1.0000],
                                            [0.1545, 0.2357, 0.2816, 0.3109, 0.3310, 0.3457, 0.3569, 0.3657, 0.3728,
                                             0.3786, 0.3835, 0.3877, 0.3912, 0.3943, 0.3971, 0.3995, 0.4016, 0.4036,
                                             0.4053, 0.4068, 0.4083, 0.4093, 0.4113, 0.4122, 0.4136, 0.4144, 0.4153,
                                             0.4173, 0.4196, 0.4227, 0.4275, 0.4343, 0.4485, 0.4728, 0.5217, 0.6330, 1.0000],
                                            [0.1488, 0.1836, 0.2025, 0.2147, 0.2233, 0.2296, 0.2346, 0.2386, 0.2419,
                                             0.2446, 0.2469, 0.2489, 0.2507, 0.2522, 0.2536, 0.2549, 0.2561, 0.2572,
                                             0.2583, 0.2593, 0.2605, 0.2616, 0.2630, 0.2645, 0.2663, 0.2686, 0.2716,
                                             0.2755, 0.2808, 0.2882, 0.2988, 0.3140, 0.3365, 0.3734, 0.4398, 0.5759, 1.0000],
                                            [0.1255, 0.1740, 0.2019, 0.2203, 0.2333, 0.2430, 0.2506, 0.2567, 0.2616,
                                             0.2658, 0.2693, 0.2723, 0.2749, 0.2772, 0.2793, 0.2811, 0.2828, 0.2843,
                                             0.2857, 0.2871, 0.2884, 0.2897, 0.2911, 0.2926, 0.2944, 0.2965, 0.2991,
                                             0.3025, 0.3070, 0.3133, 0.3225, 0.3358, 0.3570, 0.3922, 0.4547, 0.5868, 1.0000],
                                            [0.1192, 0.1557, 0.1765, 0.1902, 0.2000, 0.2074, 0.2131, 0.2178, 0.2217,
                                             0.2249, 0.2276, 0.2300, 0.2321, 0.2340, 0.2357, 0.2372, 0.2386, 0.2400,
                                             0.2413, 0.2426, 0.2440, 0.2454, 0.2470, 0.2489, 0.2511, 0.2538, 0.2572,
                                             0.2616, 0.2674, 0.2754, 0.2865, 0.3025, 0.3266, 0.3650, 0.4322, 0.5714, 1.0000],
                                            [0.0888, 0.1312, 0.1575, 0.1754, 0.1885, 0.1985, 0.2064, 0.2128, 0.2180,
                                             0.2225, 0.2263, 0.2295, 0.2324, 0.2350, 0.2373, 0.2393, 0.2412, 0.2430,
                                             0.2447, 0.2463, 0.2480, 0.2497, 0.2515, 0.2535, 0.2559, 0.2586, 0.2621,
                                             0.2665, 0.2723, 0.2802, 0.2911, 0.3069, 0.3307, 0.3682, 0.4348, 0.5736, 1.0000],
                                            [0.0533, 0.0979, 0.1299, 0.1533, 0.1711, 0.1849, 0.1960, 0.2051, 0.2126,
                                             0.2190, 0.2244, 0.2291, 0.2332, 0.2369, 0.2401, 0.2431, 0.2457, 0.2482,
                                             0.2505, 0.2527, 0.2548, 0.2569, 0.2591, 0.2614, 0.2639, 0.2669, 0.2704,
                                             0.2748, 0.2806, 0.2883, 0.2990, 0.3144, 0.3376, 0.3751, 0.4403, 0.5771, 1.0000]],
                                            [[0.4922, 0.9197, 0.9690, 0.9821, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                              0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                              0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                              0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                              0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                              0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                              0.0000, 0.0000, 0.0000],
                                            [0.3215, 0.4368, 0.4901, 0.5205, 0.5400, 0.5536, 0.5636, 0.5713, 0.5775, 0.5826,
                                             0.5873, 0.5922, 0.5991, 0.6110, 0.6383, 0.7128, 1.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000],
                                            [0.2415, 0.3423, 0.3931, 0.4234, 0.4434, 0.4577, 0.4683, 0.4765, 0.4831, 0.4885,
                                             0.4930, 0.4970, 0.5006, 0.5042, 0.5082, 0.5139, 0.5234, 0.5416, 0.5788, 0.6718,
                                             1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000],
                                            [0.1867, 0.2800, 0.3304, 0.3616, 0.3826, 0.3978, 0.4092, 0.4181, 0.4253, 0.4311,
                                             0.4360, 0.4402, 0.4438, 0.4470, 0.4499, 0.4528, 0.4558, 0.4590, 0.4637, 0.4704,
                                             0.4819, 0.5035, 0.5480, 0.6507, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000],
                                            [0.1834, 0.2756, 0.3256, 0.3567, 0.3777, 0.3928, 0.4042, 0.4131, 0.4203, 0.4262,
                                             0.4311, 0.4353, 0.4389, 0.4421, 0.4450, 0.4479, 0.4508, 0.4543, 0.4593, 0.4663,
                                             0.4775, 0.4997, 0.5440, 0.6477, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000],
                                            [0.1310, 0.1787, 0.2059, 0.2236, 0.2362, 0.2456, 0.2529, 0.2587, 0.2635, 0.2675,
                                             0.2708, 0.2737, 0.2762, 0.2784, 0.2804, 0.2821, 0.2837, 0.2851, 0.2864, 0.2876,
                                             0.2887, 0.2898, 0.2908, 0.2919, 0.2929, 0.2941, 0.2954, 0.2967, 0.2987, 0.3012,
                                             0.3047, 0.3092, 0.3147, 0.3240, 0.3376, 0.3590, 0.3938, 0.4557, 0.5871, 1.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000],
                                            [0.0940, 0.1412, 0.1704, 0.1902, 0.2046, 0.2155, 0.2241, 0.2310, 0.2367, 0.2415,
                                             0.2455, 0.2490, 0.2520, 0.2547, 0.2571, 0.2592, 0.2611, 0.2628, 0.2643, 0.2658,
                                             0.2671, 0.2683, 0.2695, 0.2706, 0.2717, 0.2728, 0.2739, 0.2750, 0.2763, 0.2779,
                                             0.2797, 0.2817, 0.2846, 0.2884, 0.2931, 0.2991, 0.3097, 0.3241, 0.3465, 0.3828,
                                             0.4467, 0.5807, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000],
                                            [0.0791, 0.1195, 0.1451, 0.1629, 0.1759, 0.1860, 0.1939, 0.2004, 0.2057, 0.2102,
                                             0.2141, 0.2174, 0.2203, 0.2229, 0.2251, 0.2272, 0.2290, 0.2306, 0.2321, 0.2335,
                                             0.2348, 0.2359, 0.2370, 0.2380, 0.2390, 0.2399, 0.2408, 0.2416, 0.2425, 0.2433,
                                             0.2442, 0.2452, 0.2462, 0.2475, 0.2488, 0.2502, 0.2522, 0.2553, 0.2577, 0.2629,
                                             0.2675, 0.2759, 0.2871, 0.3029, 0.3270, 0.3647, 0.4320, 0.5714, 1.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000],
                                            [0.0737, 0.1129, 0.1380, 0.1556, 0.1686, 0.1786, 0.1866, 0.1931, 0.1985, 0.2030,
                                             0.2069, 0.2103, 0.2132, 0.2158, 0.2181, 0.2201, 0.2220, 0.2237, 0.2252, 0.2266,
                                             0.2279, 0.2290, 0.2301, 0.2312, 0.2321, 0.2330, 0.2339, 0.2347, 0.2355, 0.2363,
                                             0.2371, 0.2379, 0.2388, 0.2397, 0.2408, 0.2417, 0.2430, 0.2447, 0.2468, 0.2495,
                                             0.2532, 0.2580, 0.2632, 0.2713, 0.2825, 0.2991, 0.3229, 0.3620, 0.4303, 0.5696,
                                             1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000],
                                            [0.0699, 0.1103, 0.1367, 0.1553, 0.1691, 0.1798, 0.1883, 0.1953, 0.2010, 0.2059,
                                             0.2100, 0.2136, 0.2167, 0.2195, 0.2220, 0.2242, 0.2261, 0.2279, 0.2295, 0.2310,
                                             0.2324, 0.2336, 0.2348, 0.2359, 0.2369, 0.2379, 0.2388, 0.2397, 0.2406, 0.2414,
                                             0.2423, 0.2432, 0.2442, 0.2453, 0.2466, 0.2479, 0.2494, 0.2514, 0.2539, 0.2574,
                                             0.2614, 0.2667, 0.2752, 0.2861, 0.3020, 0.3262, 0.3649, 0.4321, 0.5697, 1.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000]],
                                            [[0.4922, 0.9197, 0.9690, 0.9821, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                              0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                              0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                              0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                              0.0000],
                                            [0.0571, 0.2408, 0.3889, 0.4821, 0.5415, 0.5816, 0.6101, 0.6313, 0.6476, 0.6605,
                                             0.6713, 0.6811, 0.6922, 0.7125, 0.7666, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000],
                                            [0.2480, 0.3059, 0.3348, 0.3523, 0.3642, 0.3727, 0.3792, 0.3843, 0.3884, 0.3918,
                                             0.3946, 0.3970, 0.3992, 0.4011, 0.4028, 0.4045, 0.4062, 0.4082, 0.4107, 0.4142,
                                             0.4196, 0.4272, 0.4413, 0.4668, 0.5170, 0.6288, 1.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000],
                                            [0.2018, 0.2902, 0.3370, 0.3657, 0.3850, 0.3989, 0.4094, 0.4176, 0.4242, 0.4296,
                                             0.4341, 0.4379, 0.4412, 0.4442, 0.4469, 0.4495, 0.4523, 0.4555, 0.4601, 0.4675,
                                             0.4789, 0.5008, 0.5445, 0.6498, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000],
                                            [0.1991, 0.2782, 0.3201, 0.3460, 0.3635, 0.3762, 0.3858, 0.3934, 0.3994, 0.4044,
                                             0.4085, 0.4121, 0.4152, 0.4180, 0.4205, 0.4229, 0.4253, 0.4281, 0.4317, 0.4365,
                                             0.4445, 0.4572, 0.4813, 0.5290, 0.6379, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000],
                                            [0.1669, 0.2365, 0.2748, 0.2991, 0.3158, 0.3280, 0.3374, 0.3447, 0.3507, 0.3556,
                                             0.3598, 0.3633, 0.3663, 0.3690, 0.3714, 0.3736, 0.3756, 0.3775, 0.3794, 0.3815,
                                             0.3838, 0.3866, 0.3905, 0.3968, 0.4072, 0.4222, 0.4497, 0.5030, 0.6201, 1.0000,
                                             0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000],
                                            [0.1398, 0.2007, 0.2354, 0.2579, 0.2737, 0.2854, 0.2944, 0.3016, 0.3074, 0.3122,
                                             0.3163, 0.3198, 0.3228, 0.3254, 0.3278, 0.3299, 0.3317, 0.3335, 0.3351, 0.3367,
                                             0.3382, 0.3398, 0.3416, 0.3439, 0.3464, 0.3496, 0.3555, 0.3619, 0.3734, 0.3918,
                                             0.4229, 0.4796, 0.6042, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                             0.0000],
                                            [0.1230, 0.1685, 0.1946, 0.2118, 0.2241, 0.2333, 0.2405, 0.2462, 0.2509, 0.2548,
                                             0.2582, 0.2610, 0.2635, 0.2657, 0.2677, 0.2694, 0.2709, 0.2723, 0.2736, 0.2748,
                                             0.2760, 0.2770, 0.2780, 0.2790, 0.2801, 0.2812, 0.2823, 0.2837, 0.2854, 0.2872,
                                             0.2901, 0.2938, 0.2986, 0.3046, 0.3133, 0.3283, 0.3503, 0.3853, 0.4488, 0.5829,
                                             1.0000],
                                            [0.1070, 0.1555, 0.1845, 0.2040, 0.2179, 0.2284, 0.2367, 0.2433, 0.2487, 0.2532,
                                             0.2571, 0.2604, 0.2633, 0.2658, 0.2680, 0.2700, 0.2718, 0.2734, 0.2749, 0.2763,
                                             0.2776, 0.2788, 0.2799, 0.2810, 0.2822, 0.2833, 0.2846, 0.2860, 0.2878, 0.2897,
                                             0.2924, 0.2958, 0.3002, 0.3078, 0.3168, 0.3306, 0.3524, 0.3874, 0.4502, 0.5850,
                                             1.0000],
                                            [0.0913, 0.1396, 0.1697, 0.1902, 0.2051, 0.2164, 0.2253, 0.2325, 0.2384, 0.2434,
                                             0.2476, 0.2512, 0.2543, 0.2571, 0.2596, 0.2618, 0.2637, 0.2655, 0.2672, 0.2687,
                                             0.2701, 0.2714, 0.2727, 0.2739, 0.2751, 0.2764, 0.2778, 0.2794, 0.2812, 0.2835,
                                             0.2863, 0.2898, 0.2951, 0.3024, 0.3112, 0.3254, 0.3479, 0.3839, 0.4474, 0.5817,
                                             1.0000]]],
                          'p': [[0.033241357, 0.00168, 0.00000],
                                [0.00853711, 0.00520, 0.22222],
                                [0.013102269, 0.01576, 0.00215],
                                [0.004261432, 0.00319, 0.00213],
                                [0.005579845, 0.00438, 0.00351],
                                [0.007184482, 0.00583, 0.00340],
                                [0.02191819, 0.01557, 0.00787],
                                [0.052791546, 0.03678, 0.01996],
                                [0.139227102, 0.10054, 0.05496],
                                [0.424860886, 0.33714, 0.17464]],
                          'd': [[0.0000158, 0.00000 , 0.00000],
                                [0, 0.00000 , 0.00000],
                                [0, 0.00000 , 0.00000],
                                [0, 0.00000 , 0.00000],
                                [0, 0.00000 , 0.00001],
                                [0.0008306, 0.00013 , 0.00003],
                                [0.000760173, 0.00090 , 0.00003 ],
                                [0.002559568, 0.00307 , 0.00027 ],
                                [0.007675998, 0.00689 , 0.00090 ],
                                [0.0478664, 0.02508 , 0.00460 ]]}

        self.total_pop = [263846, 277179, 397083, 845064, 1134500, 925569, 720661, 609424, 332222, 181590]

        self.S = [256567, 238046, 40879, 118367, 171311, 107924, 73267, 69477, 44957, 48432]
        self.VS = [50, 28734, 295386, 329905, 420608, 246702, 157756, 90576, 49342, 35878]
        self.BS = [0, 1, 45961, 355049, 484256, 521849, 449325, 418567, 221787, 88744]

        self.E = [840, 954, 132, 340, 513, 317, 216, 205, 136, 147]
        self.VE = [0, 72, 957, 972, 1290, 742, 472, 271, 151, 110]
        self.BE = [0, 0, 141, 979, 1399, 1510, 1313, 1233, 669, 267]

        self.I = [309, 358, 49, 125, 189, 117, 79, 75, 50, 54]
        self.VI = [0, 20, 354, 365, 484, 278, 176, 100, 56, 40]
        self.BI = [0, 0, 51, 353, 506, 552, 482, 453, 246, 98]

        self.H = [0 for i in range(self.group)]
        self.VH = [0 for i in range(self.group)]
        self.BH = [0 for i in range(self.group)]

        self.M = [0 for i in range(self.group)]
        self.VM = [0 for i in range(self.group)]
        self.BM = [0 for i in range(self.group)]

        self.R = [0 for i in range(self.group)]
        self.VR = [0 for i in range(self.group)]
        self.BR = [0 for i in range(self.group)]

        self.D = [0 for i in range(self.group)]
        self.VD = [0 for i in range(self.group)]
        self.BD = [0 for i in range(self.group)]

        self.population = [self.S,  self.E,  self.I,  self.M,  self.H,  self.R,  self.D,
                           self.VS, self.VE, self.VI, self.VM, self.VH, self.VR, self.VD,
                           self.BS, self.BE, self.BI, self.BM, self.BH, self.BR, self.BD]
        self.populationH = [self.S,  self.E,  self.I,  self.M] + [self.H for t in range(self.period)] + [self.R,  self.D,
                           self.VS, self.VE, self.VI, self.VM] + [self.VH for t in range(self.period)] + [self.VR, self.VD,
                           self.BS, self.BE, self.BI, self.BM] + [self.BH for t in range(self.period)] + [self.BR, self.BD]
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 

def find_closest_time(true, lower_target_percent, upper_target_percent, sh_ratios):
    """ Purpose: Look for value closest to a specified target
        Args :  lower_target :  Lower limit of target i.e. 0.9*target
                upper_target : Upper limit of target i.e. 1.1*target
                sh_ratios : list of spitzer ratios
                dt = time step
        Returns : spitzer ratio closest to target and time step
        
    """
    lower_target = true * lower_target_percent
    upper_target = true * upper_target_percent
    print('Lower', lower_target, 'Upper', upper_target)
    array = np.asarray(sh_ratios)
    array = array[~np.isnan(array)]
    #i = (np.abs(array- lower_target)).argmin()
    #j = (np.abs(array- upper_target)).argmin()
    i = np.where(np.isclose(array, lower_target, 1e-3))[0]
    j = np.where(np.isclose(array, upper_target, 1e-3))[0]
    print('Lower Value: ', sh_ratios[i])
    #print('indices : ', i)
    print('Upper Value: ', sh_ratios[j])
    #print('indices : ', j)
    if len(i) ==0 and len(j) ==0:
    #    find_closest_time(true, lower_target_percent*1.1, upper_target_percent*1.1, sh_ratios)
        return np.nan
    if len(j) == 0:
        print('j ==0')
        return np.max(i)
    if len(i) == 0:
        print('i ==0')
        return np.max(j)
    k = np.max(np.array(np.concatenate((i,j)))) #max as we have flipped the values.
    print('Percentage is :', upper_target_percent)
    print('Percentage is :', lower_target_percent)
    return k
#createDicts('/home/abetharan/HYDRO_KINETIC_COUPLING/ES_Investigation/Tables_Of_Q_Qsh_ratio/3_DP_f_0_0_q_qsh_ratio.csv')
#df_1 = pd.DataFrame(data = all_sh_ratios)
all_sh_ratios= pd.read_csv('/home/abetharan/HYDRO_KINETIC_COUPLING/ES_Investigation/Tables_Of_Q_Qsh_ratio/local_init_f_1_q_qsh_ratio.csv', header=0)
keys = list(all_sh_ratios.keys())
#converged_ratio, time = convergence_times(correct_q_ratio[4], all_sh_ratios['02'].tolist(), dt = 0.01)
dict_of_times = {}
dt = np.array([0.1, 0.1, 0.1, 0.1, 0.01, 0.01, 0.01, 0.01])
columns = ['0.0075', '0.02', '0.04', '0.075', '0.2', '0.4', '0.75', '2']
i = 7
#sh_ratio, time = find_convergence_times(all_sh_ratios[keys[i]].tolist(), dt[i])
klambda = [0.0075, 0.02, 0.04, 0.075, 0.2, 0.4, 0.75, 2]
#f_1 = 0
#f_0_convergence_steps = np.array([262, 385 ,353 ,263 ,1711 ,1646 ,1417 ,128]) #to 2dp
#f_0_convergence_times = np.array([262 * 0.1,385 * 0.1,353 * 0.1,263 * 0.1,1711 * 0.01,1646 * 0.01,1417 * 0.01,128 * 0.01]) #to 2dp
#f_0_first_correct = np.array([0,0,130,76,277,136,72,6])
f_0_convergence_values = [0.995,0.987,0.963,0.902,0.69,0.48,0.31,0.13]
#f_1 = local
f_0_first_correct = np.array([0, 0, 0 ,276,566,296,389,152])
f_0_convergence_steps = np.array([274, 314, 47*5, 375, 2514, 2407, 2048, 1084])


f_0_convergence_times_001 = []#
f_0_convergence_times_005 = []

for i in ([0.01, 0.05]):
    for k in range(len(keys)):
        print('key', k, 'converged value is ', f_0_convergence_values[k])
        nan_sh_ratios = np.array(all_sh_ratios[keys[k]].tolist())[f_0_first_correct[k]:f_0_convergence_steps[k]]
        klambda_sh_ratios = nan_sh_ratios[~np.isnan(nan_sh_ratios)]
        if i == 0.01:
            f_0_convergence_times_001.append((find_closest_time(f_0_convergence_values[k], 1-i, 1+i, klambda_sh_ratios) + f_0_first_correct[k])*dt[k]) 
        if i == 0.05:
            f_0_convergence_times_005.append((find_closest_time(f_0_convergence_values[k], 1-i, 1+i, klambda_sh_ratios) + f_0_first_correct[k]) * dt[k]) 

plt.rcParams.update({'font.size': 40})
plt.rcParams.update({'xtick.major.size':10})
plt.rcParams.update({'xtick.minor.size':6})
plt.rcParams.update({'ytick.major.size':10})
plt.rcParams.update({'ytick.minor.size':6})
plt.rcParams.update({'lines.markeredgewidth':6})
plt.rcParams.update({'lines.markersize':10})
plt.rcParams.update({'legend.fontsize':28})
plt.rcParams.update({'legend.frameon': False})
plt.rcParams.update({'lines.linewidth' : 3})
fig = plt.figure()
fig.set_size_inches(15.5,12.5,forward=True)
ax = fig.add_subplot(111)
ax.set_autoscale_on(True)
#ax.set_ylim(0.002, 0.006)
ax.minorticks_on()
ax.tick_params(which='major', length=10, width=2,color='gray', direction='inout', bottom = True, top= True, left= True, right = True)
ax.tick_params(which='minor', length=5, width=2, color='gray', direction='in', bottom = True, top= True, left= True, right = True)
print('Convergence time for Errors 1%', f_0_convergence_times_001)
print('Convergence time for Errors 5%', f_0_convergence_times_005)

#f_0_
#f_0_first_correct = np.array([0,0,130,76,275,135,72,6])#
#f_0_convergence_steps = np.array([262, 385 ,353 ,263 ,1711 ,1646 ,1417 ,128]) #to 2dp
#percent_1 = np.multiply(np.array([121, 129 , 27, nan, 2466, 2406, 1892, 1015]), dt) #klambda = 0.04 has max error of 0.35
#percent_5 = np.multiply(np.array([50, 55, nan, nan, 1160, 1705, 1635, 1074]), dt)  #klambda = 0.75 has max error of 2.1%, klambda = 2 4.15% error
#percent_1 = f_0_convergence_times_001 
#percent_5 = f_0_convergence_times_005

f_0_convergence_steps = np.array([274, 314, 47*5, 375, 2514, 2407, 2048, 1084])
percent_1 = np.multiply(np.array([118, 126 , 28, 276, 2466, 2397, 1892, 1016]), dt) #klambda = 0.04 has max error of 0.35
percent_5 = np.multiply(np.array([50, 55, 15*5, 276, 1160, 1705, 1637, 1074]), dt)  #klambda = 0.75 has max error of 2.1%, klambda = 2 4.15% error
plt.plot(klambda, np.multiply(f_0_convergence_steps,dt), 'kv-', label = 'Converged')
plt.plot(klambda, percent_1, 'rx-', label = '1% error')
plt.plot(klambda, percent_5, 'go-', label = '5% error')
ax.set_ylabel(r'$Time/\tau_{ei}$')
ax.set_xlabel(r'$k\lambda$')
ax.set_xscale('log')
#ax.set_yscale('log')
# ax.axis('equal')
ax.legend(loc=0, frameon=None)
plt.tight_layout()
plt.savefig('/home/abetharan/Documents/MY_CONFERENCES_POSTER_PRESENTATIONS/HPL_CHRISTMAS_MEETING_2019/q_qsh_convergence_f1_local_initpng', papertyper = 'a0')
plt.show()

# klambda = [0.0075, 0.02, 0.04, 0.075, 0.2, 0.4, 0.75, 2]
# f_0_convergence_times = np.array([262, 385, 353, 263, 1070, 1722, 1860, 238])
# dt = np.array([0.1,0.1,0.1,0.1,0.01,0.01,0.01,0.01])
# f_0_error = np.array([86, 82, 0, 0, 876, 543, 1160, 126])
# f_0_error_in_convergence = f_0_convergence_times - f_0_error
# plt.errorbar(klambda, f_0_convergence_times * dt, yerr = f_0_error_in_convergence * dt)
# plt.show()
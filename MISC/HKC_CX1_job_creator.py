from string import Template
import os
import sys
sys.path.insert(0, "/home/abetharan/HYDRO_KINETIC_COUPLING/MISC/HKC_CX1_job_creator.py")
import yaml
from distutils.dir_util import copy_tree
def templateCX1YML(write_path_, base_dir_, input_path_, run_name_, cycles_, K_NT_, K_DT_, F_NT_, F_DT_):
    
    file_path = '/home/abetharan/HYDRO_KINETIC_COUPLING/MISC/CX1_input.yml'
    with open(file_path) as f:
        list_doc = yaml.load(f)

    list_doc['BASE_DIR'] = base_dir_

    list_doc['RUN_NAME'] = run_name_

    list_doc['CYCLES'] = cycles_
    
    list_doc['F_INIT_PATH'] = input_path_
    
    list_doc['K_DT'] = K_DT_

    list_doc['K_NT'] = K_NT_

    list_doc['F_INITIAL_DT'] = F_DT_

    list_doc['F_STEPS'] = F_NT_

    dump_path = write_path_ + '/INPUT.yml'

    with open(dump_path, "w") as f:
        yaml.dump(list_doc, f)

def load_template(filename):

    template_file = open(filename)

    template_string = ''

    for lines in template_file:
        template_string = template_string + lines

    batch_template = Template(template_string)

    return batch_template


def templateCX1PBS(write_path_, 
                   file_name_,
                   input_path,
                   run_path,
                   nump='16',
                   time_h='24',
                   time_m='00',
                   mem='60GB',
                   ):


    templ = load_template('/home/abetharan/HYDRO_KINETIC_COUPLING/MISC/Hydro_kinetic_coupling_job_template')

    mins = int(time_h)*60 + int(time_m) - 20

    batch = templ.safe_substitute(num_proc=nump,
                                  job_mem=mem,
                                  job_time_h=time_h,
                                  job_time_m=time_m,
                                  INPUT_PATH=input_path,
                                  RUN_PATH=run_path,
                                  MINS=str(mins)+'m')

    file = open(write_path_ + '/' + file_name_,'w')
    file.write(batch)

def creator(main_job_name, job_names, cycles_, K_NT_, K_DT_, F_NT_, F_DT_, np_ = [16], mem_ = ['60GB']):
    j = 0

    for i in job_names:
        file_name = i + '_PBS'
        main_writePath = '/home/abetharan/HYDRO_KINETIC_COUPLING/MISC/' + main_job_name
        writePath = main_writePath + '/' +  'INPUT_FOR_JOB_' + i 
        yml_input_path = '/rds/general/user/aa10516/home/' + main_job_name + '/' +  'INPUT_FOR_JOB_' + i + '/INPUT.yml'
        init_input_path = '/rds/general/user/aa10516/home/' + main_job_name + '/' + 'INPUT_FOR_JOB_' + i + '/init/'
        make_init_input_path = '/home/abetharan/HYDRO_KINETIC_COUPLING/MISC/' + main_job_name + '/' +  'INPUT_FOR_JOB_' + i + '/init/'
        
        if not os.path.exists(main_writePath):
            os.makedirs(main_writePath)

        if not os.path.exists(writePath):
            os.makedirs(writePath)

        if not os.path.exists(make_init_input_path):
            copy_tree('/home/abetharan/ELH1/init_data/', make_init_input_path)
                     
        if not os.path.exists(main_writePath + '/submit_batch.sh'):
            import shutil
            shutil.copy('/home/abetharan/HYDRO_KINETIC_COUPLING/MISC/submit_batch.sh', main_writePath + '/submit_batch.sh')

        
        base_dir = '/rds/general/user/aa10516/home/' + main_job_name 
        run_path = base_dir + '/' + i + '/'
        cycles = cycles_[j]
        K_NT = K_NT_[j]
        K_DT = K_DT_[j]
        F_NT = F_NT_[j]
        F_DT = F_DT_[j]
        np = np_[j]
        mem = mem_[j]
        
        if len(cycles_) != 1:
            j+=1    
        
        templateCX1PBS(main_writePath, file_name, yml_input_path, run_path, nump = np)
        templateCX1YML(writePath, base_dir, init_input_path, i, cycles, K_NT, K_DT, F_NT, F_DT)




K_time_step = []
f_time_step = []
names = []
f_dt = []
k_dt = [] 
cycles = []
nump = []
mem = []
k_dt_ = 0.01
k_NT = 842
f_dt_ = 1e-8
f_NT = 1590
import numpy as np
steps = [1, 1.5, 2, 2.5]
import math
for i in steps:
    for j in steps:
        k_step = math.ceil(k_NT * i)
        f_step = math.ceil(f_NT * j)

        K_time_step.append(k_step)
        f_time_step.append(f_step)            
        names.append(str(i) + "_KNT_" + str(j) + "_KFNT")
        nump.append(32)
        mem.append('60GB')
        cycles.append(7)
        k_dt.append(k_dt_)
        f_dt.append(f_dt_)


name = 'EPPERLEIN_SHORT_COUPLED_02'
creator(name, names, cycles, K_time_step, k_dt, f_time_step, f_dt,nump, mem)
#creator('test_restart', 'restart_test', [50], [1000], [0.01], [1000], [1e-7])

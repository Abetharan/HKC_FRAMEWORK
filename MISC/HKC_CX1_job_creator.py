from string import Template
import os
import yaml

def templateCX1YML(write_path_, base_dir_, run_name_, cycles_, K_NT_, K_DT_, F_NT_, F_DT_):
    
    file_path = os.getcwd() + '/CX1_input.yml'
    with open(file_path) as f:
        list_doc = yaml.load(f)

    list_doc['BASE_DIR'] = base_dir_

    list_doc['RUN_NAME'] = run_name_

    list_doc['CYCLES'] = cycles_

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


    templ = load_template(os.getcwd() + '/Hydro_kinetic_coupling_job_template')

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

def creator(main_job_name, job_names, cycles_, K_NT_, K_DT_, F_NT_, F_DT_, np_, mem_):
    j = 0

    for i in job_names:
        file_name = i + '_PBS'
        main_writePath = os.getcwd() + '/' + main_job_name
        writePath = main_writePath + '/' + i 
        input_path = '$HOME/' + i
        
        if not os.path.exists(main_writePath):
            os.makedirs(main_writePath)
        if not os.path.exists(writePath):
            
        
        base_dir = '$HOME/' + i
        run_path = base_dir + '/' + i
        cycles = cycles_[j]
        K_NT = K_NT_[j]
        K_DT = K_DT_[j]
        F_NT = F_NT_[j]
        F_DT = F_DT_[j]
        np = np_[j]
        mem = mem_[j]
        
        if len(cycles_) != 1:
            j+=1    
        
        templateCX1PBS(main_writePath, file_name, input_path, run_path, nump = np, mem = mem)
        templateCX1YML(writePath, base_dir, input_path, cycles, K_NT, K_DT, F_NT, F_DT)





creator('kappa', ['kappa_1', 'kappa_2'], [50], [1], [1],[1] ,[1],[1],[1])
#----
import numpy as np
from scipy import constants
k_b = constants.value("Boltzmann constant")
m_e = constants.value("electron mass")
m_p = constants.value("proton mass")
q_e = constants.value("elementary charge")
epsilon0 = 8.854188E-12    # Vacuum dielectric constant
planck_h = constants.value("Planck constant")
bohr_radi = constants.value("Bohr radius")
#-----------------------------------------------------------------------

def extract_power(x):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return b    
#-----------------------------------------------------------------------


def impact_inputs(ne,Te,Z,Ar, Bz):
    '''
        ne, Te, Z, Bz, Ar are the IMPACT reference density [cm**-3], temperature [eV], ionisation (Z)
        magnetic field [T], relative atomic mass (Ar)
        
        returns a dictionary containing the normalisation parameters for the IMPACT reference material
    '''
    # Reference density [cm**-3], temperature [eV], ionisation and B-field [T]
    #ne = 1.0e+20#10.0e19#0.5*9.08e21
    #Te = 30.0
    #Z = 2.0 #ls
    #Bz = 3.75 #3.75
    #Ar = 40.0
    
    dict = {}
    dict['ne'] = ne
    dict['Te'] = Te
    dict['Z'] = Z
    dict['Ar'] = Ar
    # Convert ne to 10**21 cm**-3
    ne = ne / 1.0e21
    ni = ne / Z

    # Calculate Coulomb np.logarithm.
    
    #Lamda bar de brogile length for a thermal electron
    lambda_bar = (2.76e-10) / np.sqrt(Te) 
    #90 scatter cross scattering
    b0 = (1.44e-9) * (Z / Te)

    if (lambda_bar > b0):
        log_lambda = np.log(Te) - np.log(ne)/2.0 - 0.16
    else: 
        log_lambda = (3.0/2.0)*np.log(Te) - np.log(Z) - np.log(ne)/2.0 - 1.8


    # Reference values
    #   v0 : thermal velocity [m/s],      t0 : collision time [s].
    #   nu0 : collision frequency [s**-1], l0 : mean free path [m])
    v0 = (0.5931e6) * (Te**0.5)
    t0 = (2.588e-16) * (Te**(1.5)) / (Z*Z*ni*log_lambda)
    nu0 = 1.0 / (t0)
    l0 = (1.535e-10) * (Te*Te) / (Z*Z*ni*log_lambda)
    # IMPACT inputs 
    wpe_by_nu_ei = 0.4618 * (Te**(1.5)) / (Z * np.sqrt(ne) * log_lambda)
    # for tau_B instead of tau_ei as in IMPACTA - note that ln_lam will be
    # different too
    #wpe_by_nu_ei = 0.4618 * (Te**(3/2)) / (Z * sqrt(ne) * log_lambda) * (3*sqrt(pi)/4)
    c_by_v0 = 1.0 / ( (1.9784e-3) * np.sqrt(Te) )
    prof_Bz_ave = (1.756e11) * Bz * t0

    # Display
    
    print('\nINPUT QUANTITIES')
    print('Density \t', ne * 1e21,'[cm**3]')

    print('Temperature \t', Te, '[eV]')
    print('Ionisation \t',Z, '[ ]')
    print('Bz \t\t', Bz, '[T]')

    print('\n IMPACT VARIABLES')
    print('log_lambda \t', log_lambda)

    print('wpe_by_nu_ei \t', wpe_by_nu_ei)
 
    print('c / v0 \t\t', c_by_v0)
    print('prof_Bz_ave \t', prof_Bz_ave)


    #clear arry
    print('\n\nPLASMA REFERENCE VARIABLES')
    print('Reference thermal velocity \t %1.5e [m/s]' % (v0))
    print('Reference collision time \t %1.5e [s]' % (t0))
    print('Reference collision frequency \t%1.5e [s**-1]' % (nu0))
    print('Reference mfp \t\t \t %1.5e [m]\n' % (l0))

    dict['vte'] = v0
    dict['tau_ei'] = t0
    dict['nu_ei'] = nu0
    dict['lambda_mfp'] = l0
    dict['c_over_vte'] = c_by_v0
    dict['log_lambda'] = log_lambda
    dict['wpe_over_nu_ei'] = wpe_by_nu_ei 
    dict['Bz_norm'] = prof_Bz_ave
    dict["ne"] = ne
    dict["ni"] = ni
    dict["Te"] = 2 * Te
    # --- get transport coeffs
    # if T in Kelvin...
    #kappa*gradT = (1/(msec))*(J/K) * (K/m) = kappa*k_b *grad (T/Kelvin)  
    # W/m^2= [kappa*gradT] = (1/(m*sec))*J/m = J/(m^2 sec) 
    
    
    return dict

#-----------------------------------------------------------------------
def calc_wo(Te, lamda, ne, Z, log_lamda):
    """ 
    Purpose: Calculates Wo in IMPACT units
    Args:
        Te = temperature in eVs.
        lamda = Wavelength in meters.
        ne = electron density in cm^3.
        Z = Ionisation.
        log_lamda = Coulomb log. 

    Returns:
        laser frequency in Impac Norms.
    """
    
    wo = (0.6 * Te ** 3/2) / ((lamda/1E-6) * (ne/1E21) * Z**2 * log_lamda)

    return(wo)


#-----------------------------------------------------------------------
def calc_norms(var,normal_dict,sample =0.0,forced_power=[]):
    '''
        norm_const, ylab = calc_norms(var,sample,forced_power)
        Obtains the normalising constant to convert variable into units specified in ylab.
        
        This function will automatically quote the variable to one significant figure by default.
        To override this set forced_power=0
        
    '''
    c_fmt = '%3.2f'
    
    v_te = normal_dict['vte']
    tau_ei = normal_dict['tau_ei']
    nu_ei = normal_dict['nu_ei']
    lambda_mfp = normal_dict['lambda_mfp']
    n0,T0 = normal_dict['ne'], normal_dict['Te']
    #-----------------

    if var== 'Cx' or var=='Cy':
        
        norm_const = v_te*1e-3
        title = r'$' + var[0] + '_' + var[1] + r'$ [ $  kms^{-1}$ ]'
        
    elif var=='n':
        power = extract_power(n0)
        mod = r'$ 10^{' +str(power)  + '} $'

        norm_const = n0*(10**-power)
        title = r'$n_e$ [' + mod + r'$\si{cm^{-3}}$'+ r']'

    elif var=='Te':
        norm_const = 2.0*T0
        title = r'$T_e$ [ $eV$ ]'
        c_fmt = '%3.0f'
    elif var=='Ui':
        norm_const = (2.0/3.0)*(2.0*T0)
        title = r'$T_i$ [ $eV$ ]'
        c_fmt = '%3.0f'
    
    elif var=='Bz':
        norm_const = (m_e/(q_e*tau_ei))
        
        power = extract_power(sample*norm_const)
        norm_const*= (10**-power)
        ##print ' sample = ', sample, 'power = ', power

        var_name = '$B_z$'
        if power==0:
            mod = ''            
            units = r'[$si{T}$]'
            title = r'$B_z$ [$\si{T}$ ]'

        else:
            mod = r'$ 10^{' +str(power)  + '} $'
        
            units = r'[' + mod + '$\si{T}$]'
            title = r'$B_z$ [' + mod + '$\si{T}$ ]'
        c_fmt = '%1.1f'
        
    elif var=='wt':
        if len(forced_power)!=0:
            power = forced_power[0]
        else:
            power = extract_power(sample)
        
        if power!= 0:
            mod = r'$ 10^{' +str(power)  + '} $'
        else:
            mod = r''            
        ##mod = r'$ 10^{' +str(power)  + '} $'
        
        norm_const = 1.0*(10**(-power))
        
        title = mod + r'$\omega \tau_{ei}$'
        c_fmt = '%1.1f'
        
    elif var[0]=='E':
        #print 'DOING E-field - min sample = ', sample
        norm_const = (lambda_mfp/(tau_ei**2))*(m_e/q_e)
        power = extract_power(norm_const*sample)
        #print ' power extracted = ', power
        c_fmt = '%1.1f'
        mod = r'$ 10^{' +str(power)  + '}$'
        norm_const = norm_const*(10**-power)
        if power==0:
            mod=''           
        title = r'$' + var[0] + '_' + var[1] + r'$ [ ' + mod + '$V/m$ ]'
    
    elif var[0] =='q':
        
        c_fmt = '%1.1f'
        power = extract_power(sample)
        norm_const = 1.0*(10**-power)
        mod = r'$ 10^{' +str(power)  + '} $'
        if power==0:
            mod=''
        #c_fmt = '%1.0e'
        #norm_const = m_e*(v_te**3)*n_04
        title = r'$' + var[0] + '_' + var[1] + r'$ [ ' + mod + '$q_0$ ]'# + r' [ $Jms^{-1}$ ]'

    elif var[0] =='j':
        c_fmt = '%1.1f'
        power = extract_power(sample)
        norm_const = 1.0*(10**-power)
        mod = r'$ 10^{' +str(power)  + '} $'
        if power==0:
            mod=''
        #norm_const = 1.0#q_e*n0*v_te
        
        title = r'$' + var[0] + '_' + var[1] + r'$ [ ' +mod + '$j_0$ ]'#r' [ $C/m^2 s$ ]'    
    elif var =='U':
        
        c_fmt = '%1.1f'
        power = extract_power(sample)
        norm_const = 1.0*(10**-power)
        mod = r'$ 10^{' +str(power)  + '} $'
        #c_fmt = '%1.0e'
        #norm_const = m_e*(v_te**3)*n_04
        title = r'$' + var[0] + r'$ [ ' + mod + '$m_e v_n^2 n_0$ ]'# + r' [ $Jms^{-1}$ ]'

    return norm_const, title, c_fmt

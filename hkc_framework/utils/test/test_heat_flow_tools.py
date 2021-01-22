import pytest
import numpy as np 
import os 
import sys
from scipy import constants
mp = constants.value("proton mass")
BOLTZMANN_CONSTANT = constants.value("Boltzmann constant")
ELECTRON_MASS = constants.value("electron mass")
PROTON_MASS = constants.value("proton mass")
ELEMENTARY_CHARGE = constants.value("elementary charge")
VACUUM_PERMITTIVITY = 8.854188E-12    # Vacuum dielectric constant
PLANCK_CONSTANT = constants.value("Planck constant")
BOHR_RADIUS = constants.value("Bohr radius")
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../../../')
from hkc_framework.utils.heat_flow_coupling_tools import HeatFlowCouplingTools


class TestHFCT():

    def test_spitzer_harm(self):
        true_sh = np.array([5.75922290e+15, 7.37788102e+15, 7.96809003e+15, 8.59090381e+15,
       9.75644723e+15, 1.17071165e+16, 1.45091910e+16, 1.80604684e+16,
       2.18283952e+16, 2.53601545e+16, 2.86952761e+16, 3.17744431e+16,
       3.45701054e+16, 3.70855052e+16, 3.93610384e+16, 4.14955795e+16,
       4.36013905e+16, 4.57347134e+16, 4.78701648e+16, 4.99643907e+16,
       5.20029743e+16, 5.39898943e+16, 5.59476581e+16, 5.79269115e+16,
       5.99404220e+16, 6.19555714e+16, 6.39387694e+16, 6.58594512e+16,
       6.77078241e+16, 6.94877844e+16, 7.12162756e+16, 7.29286163e+16,
       7.46842842e+16, 7.65222609e+16, 7.84311979e+16, 8.04106311e+16,
       8.24702721e+16, 8.46264818e+16, 8.69028960e+16, 8.93359184e+16,
       9.18275763e+16, 9.44008917e+16, 9.72095153e+16, 1.00327566e+17,
       1.03832478e+17, 1.07816241e+17, 1.12396120e+17, 1.17739322e+17,
       1.23633523e+17, 1.29957811e+17, 1.36961601e+17, 1.44780680e+17,
       1.53568456e+17, 1.63507389e+17, 1.74814592e+17, 1.87741725e+17,
       2.03201972e+17, 2.08664083e+17, 2.14450860e+17, 2.20533822e+17,
       2.26929054e+17, 2.33652129e+17, 2.40720810e+17, 2.48154740e+17,
       2.55975485e+17, 2.64206636e+17, 2.72873913e+17, 2.82004610e+17,
       2.91630328e+17, 3.01783631e+17, 3.12500450e+17, 3.23819028e+17,
       3.35782252e+17, 3.48436915e+17, 3.61832348e+17, 3.76021693e+17,
       3.91064749e+17, 4.07022579e+17, 4.23962260e+17, 4.41952934e+17,
       4.61058188e+17, 4.82362639e+17, 4.83697533e+17, 4.85459940e+17,
       4.87232916e+17, 4.89028229e+17, 4.90811422e+17, 4.92628445e+17,
       4.94445188e+17, 4.96261227e+17, 4.98111569e+17, 4.99961837e+17,
       5.01835354e+17, 5.03697057e+17, 5.05593713e+17, 5.07478649e+17,
       5.09398866e+17, 5.11319475e+17, 5.13252105e+17, 5.15196881e+17,
       5.17153863e+17, 5.19123212e+17, 5.21117189e+17, 5.23099844e+17,
       5.25119133e+17, 5.27139476e+17, 5.29172738e+17, 5.31219033e+17,
       5.33278483e+17, 5.35351202e+17, 5.37437302e+17, 5.39549366e+17,
       5.41650610e+17, 5.43790106e+17, 5.45931462e+17, 5.48086818e+17,
       5.50256301e+17, 5.52440042e+17, 5.54638166e+17, 5.56850794e+17,
       5.59090807e+17, 5.61320584e+17, 5.63590341e+17, 5.65849999e+17,
       5.68150029e+17, 5.70453016e+17, 5.72771448e+17, 5.75105473e+17,
       5.77455226e+17, 5.79833894e+17, 5.82202997e+17, 5.84600804e+17,
       5.87028062e+17, 5.89459176e+17, 5.91893688e+17, 5.94370876e+17,
       5.96852271e+17, 5.99350734e+17, 6.01866422e+17, 6.04399492e+17,
       6.06963522e+17, 6.09518950e+17, 6.12105141e+17, 6.14722891e+17,
       6.17332322e+17, 6.19986707e+17, 6.22646608e+17, 6.25325234e+17,
       6.28009058e+17, 6.30738824e+17, 6.33474618e+17, 6.36243678e+17,
       6.39005289e+17, 6.41799962e+17, 6.44614622e+17, 6.47463403e+17,
       6.50305209e+17, 6.53194995e+17, 6.56078058e+17, 6.59009632e+17,
       6.61948872e+17, 6.64909480e+17, 6.67891667e+17, 6.70895650e+17,
       6.73921632e+17, 6.76969835e+17, 6.80040473e+17, 6.83133761e+17,
       6.86249921e+17, 6.89403643e+17, 6.92552355e+17, 6.95753055e+17,
       6.98949064e+17, 7.02197686e+17, 7.05441950e+17, 7.08739450e+17,
       7.12047667e+17, 7.15380929e+17, 7.18724671e+17, 7.22122953e+17,
       7.25532844e+17, 7.28968803e+17, 7.32431101e+17, 7.35920015e+17,
       7.39435810e+17, 7.42993853e+17, 7.46549840e+17, 7.50148018e+17,
       7.53789433e+17, 7.57429446e+17, 7.61112651e+17, 7.64840110e+17,
       7.68566853e+17, 7.72353216e+17, 7.76139348e+17, 7.79985893e+17,
       7.83848219e+17, 7.87725881e+17, 7.91665193e+17, 7.95621216e+17,
       7.99609088e+17, 8.03613583e+17, 8.07681430e+17, 8.11767241e+17,
       8.15886413e+17, 8.20039336e+17, 8.24226386e+17, 8.28447957e+17,
       8.32704443e+17, 8.36996244e+17, 8.41323769e+17, 8.45687424e+17,
       8.50087631e+17, 8.54541009e+17, 8.59000162e+17, 8.63512651e+17,
       8.68063457e+17, 8.72653020e+17, 8.77298205e+17, 8.81951050e+17,
       8.86659735e+17, 8.91425606e+17, 8.96200368e+17, 9.01032568e+17,
       9.05923593e+17, 9.10824811e+17, 9.15785145e+17, 9.20806010e+17,
       9.25838452e+17, 9.30948615e+17, 9.36071327e+17, 9.41273068e+17,
       9.46488382e+17, 9.51767048e+17, 9.57110596e+17, 9.62469297e+17,
       9.67910464e+17, 9.73367913e+17, 9.78909288e+17, 9.84468110e+17,
       9.90112365e+17, 9.95775289e+17, 1.00152520e+18, 1.00729505e+18,
       1.01315350e+18, 1.01903322e+18, 1.02498566e+18, 1.03101255e+18,
       1.03706280e+18, 1.04320594e+18, 1.04937384e+18, 1.05563650e+18,
       1.06192552e+18, 1.06831111e+18, 1.07472471e+18, 1.08121891e+18,
       1.08779565e+18, 1.09440295e+18, 1.10109379e+18, 1.10787020e+18,
       1.11467989e+18, 1.12157625e+18, 1.12856139e+18, 1.13558274e+18,
       1.14269409e+18, 1.14987929e+18, 1.15713949e+18, 1.16449424e+18,
       1.17189055e+18, 1.17938288e+18, 1.18695503e+18, 1.19460829e+18,
       1.20234394e+18, 1.21016333e+18, 1.21806781e+18, 1.22605879e+18,
       1.23413767e+18, 1.24230594e+18, 1.25056507e+18, 1.25891662e+18,
       1.26736214e+18, 1.27588430e+18, 1.28454044e+18, 1.29327762e+18,
       1.30209637e+18, 1.31105446e+18, 1.32009879e+18, 1.32923003e+18,
       1.33850636e+18, 1.34785534e+18, 1.35735349e+18, 1.36692822e+18,
       1.37663707e+18, 1.38648344e+18, 1.39641269e+18, 1.40648279e+18,
       1.41667789e+18, 1.42700035e+18, 1.43747217e+18, 1.44803861e+18,
       1.45875834e+18, 1.46959617e+18, 1.48061230e+18, 1.49175329e+18,
       1.50304041e+18, 1.51447667e+18, 1.52604541e+18, 1.53780766e+18,
       1.54969031e+18, 1.56177317e+18, 1.57398317e+18, 1.58640047e+18,
       1.59895208e+18, 1.61169862e+18, 1.62464541e+18, 1.63773813e+18,
       1.65103792e+18, 1.66453050e+18, 1.67822033e+18, 1.69211196e+18,
       1.70619010e+18, 1.72051819e+18, 1.73504416e+18, 1.74977173e+18,
       1.76476496e+18, 1.77995218e+18, 1.79541642e+18, 1.81108644e+18,
       1.82702556e+18, 1.84322165e+18, 1.85970149e+18, 1.87641333e+18,
       1.89340129e+18, 1.91071294e+18, 1.92829899e+18, 1.94618594e+18,
       1.96436202e+18, 1.98289446e+18, 2.00171534e+18, 2.02091089e+18,
       2.04041415e+18, 2.06029188e+18, 2.08055644e+18, 2.10116073e+18,
       2.12217265e+18, 2.14356610e+18, 2.16541157e+18, 2.18766621e+18,
       2.21034196e+18, 2.23351079e+18, 2.25711245e+18, 2.28123786e+18,
       2.30582867e+18, 2.33095737e+18, 2.35662451e+18, 2.38284928e+18,
       2.40965186e+18, 2.43705354e+18, 2.46505732e+18, 2.49374261e+18,
       2.52306152e+18, 2.55311281e+18, 2.58385227e+18, 2.61536203e+18,
       2.64765649e+18, 2.68076876e+18, 2.71473399e+18, 2.74958960e+18,
       2.78535690e+18, 2.82213101e+18, 2.85988895e+18, 2.89872867e+18,
       2.93870429e+18, 2.97982050e+18, 3.02218509e+18, 3.06584722e+18,
       3.11085884e+18, 3.15734347e+18, 3.20531952e+18, 3.25493656e+18,
       3.30622946e+18, 3.35934600e+18, 3.41438673e+18, 3.47147673e+18,
       3.53075367e+18, 3.59236967e+18, 3.65647847e+18, 3.72330875e+18,
       3.79301947e+18, 3.86588297e+18, 3.94215627e+18, 4.02214438e+18,
       4.10619757e+18, 4.19472173e+18, 4.28818001e+18, 4.38716590e+18,
       4.49231978e+18, 4.60449454e+18, 4.72472298e+18, 4.85434843e+18,
       4.99518449e+18, 5.14984860e+18, 5.32272707e+18, 5.48886945e+18,
       5.63181212e+18, 5.75039073e+18, 5.84394376e+18, 5.91168487e+18,
       5.95321151e+18, 5.96864796e+18, 5.95867324e+18, 5.92443077e+18,
       5.86751281e+18, 5.78983371e+18, 5.69351074e+18, 5.58083597e+18,
       5.45412325e+18, 5.31569031e+18, 5.16781381e+18, 5.01261006e+18,
       4.85210297e+18, 4.68813592e+18, 4.52239482e+18, 4.35638115e+18,
       4.19141321e+18, 4.02866102e+18, 3.86909964e+18, 3.71356006e+18,
       3.56273022e+18, 3.41713916e+18, 3.27722530e+18, 3.14328107e+18,
       3.01551354e+18, 2.89404053e+18, 2.77890387e+18, 2.67007480e+18,
       2.56746695e+18, 2.47095142e+18, 2.38035643e+18, 2.29548331e+18,
       2.21610583e+18, 2.14198331e+18, 2.07285953e+18, 2.00847337e+18,
       1.94855726e+18, 1.89284558e+18, 1.84107337e+18, 1.79297940e+18,
       1.74831269e+18, 1.70682784e+18, 1.66828623e+18, 1.61931123e+18,
       1.57624471e+18, 1.55214994e+18, 1.52985054e+18, 1.50919427e+18,
       1.49003438e+18, 1.47223312e+18, 1.45565802e+18, 1.44018690e+18,
       1.42570339e+18, 1.41209809e+18, 1.39926941e+18, 1.38712151e+18,
       1.37556573e+18, 1.36451847e+18, 1.35390388e+18, 1.34364850e+18,
       1.33368649e+18, 1.32395670e+18, 1.31440117e+18, 1.30496730e+18,
       1.29560675e+18, 1.28627432e+18, 1.27692924e+18, 1.26753250e+18,
       1.25805092e+18, 1.24845138e+18, 1.23870567e+18, 1.22878715e+18,
       1.21867226e+18, 1.20833939e+18, 1.19776898e+18, 1.18694408e+18,
       1.17584922e+18, 1.16447117e+18, 1.15279864e+18, 1.14082175e+18,
       1.12853225e+18, 1.11592319e+18, 1.10299008e+18, 1.08972911e+18,
       1.07613864e+18, 1.06221670e+18, 1.04796520e+18, 1.03338527e+18,
       1.01848012e+18, 1.00325434e+18, 9.87713447e+17, 9.71863920e+17,
       9.55714176e+17, 9.39272737e+17, 9.22550193e+17, 9.05557331e+17,
       8.88306480e+17, 8.70810563e+17, 8.53084105e+17, 8.35141806e+17,
       8.16999611e+17, 7.98674443e+17, 7.80183446e+17, 7.61545137e+17,
       7.42778401e+17, 7.23902693e+17, 7.04938252e+17, 6.85905575e+17,
       6.66826110e+17, 6.47720975e+17, 6.28612153e+17, 6.09521860e+17,
       5.90472144e+17, 5.71485620e+17, 5.52584551e+17, 5.33791361e+17,
       5.15128355e+17, 4.96617725e+17, 4.78281166e+17, 4.60140215e+17,
       4.42215850e+17, 4.24528768e+17, 4.07098843e+17, 3.89945422e+17,
       3.73087147e+17, 3.56541810e+17, 3.40326320e+17, 3.24456700e+17,
       3.08948017e+17, 2.93814148e+17, 2.79068014e+17, 2.64721301e+17,
       2.50784476e+17, 2.37266807e+17, 2.24176232e+17, 2.11519405e+17,
       1.99301665e+17, 1.87526925e+17, 1.76197768e+17, 1.65315462e+17,
       1.54879798e+17, 1.44889304e+17, 1.35341133e+17, 1.26231163e+17,
       1.17553999e+17, 1.09303065e+17, 1.01470653e+17, 9.40480161e+16,
       8.70254963e+16, 8.03926614e+16, 7.41384649e+16, 6.82514661e+16,
       6.27201336e+16, 5.75331391e+16, 5.26798013e+16, 4.81506386e+16,
       4.39379539e+16, 4.00368976e+16, 3.64491578e+16, 3.30352114e+16,
       2.91130987e+16, 2.59815849e+16, 2.30823721e+16, 2.03420353e+16,
       1.56575853e+16, 1.13854566e+16, 8.08944835e+15, 5.77694062e+15,
       4.19897579e+15, 3.09933533e+15, 2.36144868e+15, 2.08011806e+15,
       1.86535649e+15, 1.43735680e+15, 1.05911602e+15, 8.40847901e+14,
       7.17794248e+14, 5.90993955e+14, 4.70522425e+14, 3.68505810e+14,
       2.92159757e+14, 2.29411710e+14, 1.78291737e+14, 1.50278094e+14,
       1.29679407e+14, 1.10111155e+14, 1.06752948e+14, 9.97678992e+13,
       9.81352424e+13, 8.86190879e+13, 7.62182245e+13])
        hfct_obj = HeatFlowCouplingTools()
        fluid_Te = np.loadtxt(os.path.join(myPath,'test/gdhohlraum_xmic_5ps_TekeV_interp'), skiprows=1)[:, 1] * 1E3 *(ELEMENTARY_CHARGE/BOLTZMANN_CONSTANT)
        fluid_Z = np.loadtxt(os.path.join(myPath,'test/gdhohlraum_xmic_Z_interp'), skiprows=1)[:, 1] 
        fluid_x_grid = np.loadtxt(os.path.join(myPath,'test/gdhohlraum_xmic_5ps_separatedsnbWcm2'), skiprows=1)[:, 0]*1e-6
        fluid_ne = np.loadtxt(os.path.join(myPath,'test/gdhohlraum_xmic_ne1e20cm3_interp'), skiprows=1)[:, 1] * (1e20 * 1e6)
        fluid_x_centered_grid = np.array([(fluid_x_grid[i] + fluid_x_grid[i+1])/2 for i in range(len(fluid_x_grid) - 1)])
        hfct_obj.electron_temperature = fluid_Te
        hfct_obj.electron_number_density = fluid_ne
        hfct_obj.zbar = fluid_Z
        hfct_obj.cell_wall_coord = fluid_x_grid
        hfct_obj.cell_centered_coord = fluid_x_centered_grid
        # hfct_obj.mass = fluid_mass
        hfct_obj.lambda_ei(hfct_obj.electron_temperature * (BOLTZMANN_CONSTANT/ELEMENTARY_CHARGE), 
                            hfct_obj.electron_number_density,
                            hfct_obj.zbar)
        hfct_obj.spitzerHarmHeatFlow()
        rel_err = 100*(abs(hfct_obj.spitzer_harm_heat[1:-1] - true_sh) / true_sh)
        assert np.nanmax(rel_err) < 1e-6 


    def test_snb(self):
        true_snb_heat_flow = np.array([0, 2.63172238e+15, 5.10400192e+15, 7.49830626e+15, 9.87855348e+15,
        1.22955265e+16, 1.47819239e+16, 1.73506431e+16, 1.99942102e+16,
        2.26666720e+16, 2.53321496e+16, 2.79968812e+16, 3.06487864e+16,
        3.32808999e+16, 3.58914581e+16, 3.84837994e+16, 4.10658680e+16,
        4.36259163e+16, 4.61672491e+16, 4.87198626e+16, 5.12866885e+16,
        5.38714453e+16, 5.64790410e+16, 5.91156267e+16, 6.17882677e+16,
        6.44748504e+16, 6.71875875e+16, 6.99479802e+16, 7.27606073e+16,
        7.56321642e+16, 7.85717364e+16, 8.15909481e+16, 8.47040215e+16,
        8.79050755e+16, 9.12144598e+16, 9.46630811e+16, 9.82675830e+16,
        1.02046500e+17, 1.06020532e+17, 1.10212749e+17, 1.14648637e+17,
        1.19330099e+17, 1.24293646e+17, 1.29598346e+17, 1.35289454e+17,
        1.41415725e+17, 1.48029500e+17, 1.55187446e+17, 1.62949387e+17,
        1.71337728e+17, 1.80469585e+17, 1.90413295e+17, 2.01271846e+17,
        2.13164975e+17, 2.26232579e+17, 2.40639163e+17, 2.56580204e+17,
        2.74346440e+17, 2.80233885e+17, 2.86324027e+17, 2.92624734e+17,
        2.99144627e+17, 3.05892830e+17, 3.12879018e+17, 3.20113450e+17,
        3.27607012e+17, 3.35371256e+17, 3.43418460e+17, 3.51761676e+17,
        3.60414813e+17, 3.69392684e+17, 3.78711093e+17, 3.88386925e+17,
        3.98438261e+17, 4.08884487e+17, 4.19746404e+17, 4.31046401e+17,
        4.42808656e+17, 4.55059316e+17, 4.67826811e+17, 4.81142198e+17,
        4.95039636e+17, 5.09725539e+17, 5.10918965e+17, 5.12118524e+17,
        5.13322400e+17, 5.14530663e+17, 5.15743182e+17, 5.16960177e+17,
        5.18181520e+17, 5.19407229e+17, 5.20637478e+17, 5.21872138e+17,
        5.23111331e+17, 5.24354927e+17, 5.25603150e+17, 5.26855820e+17,
        5.28113162e+17, 5.29375046e+17, 5.30641545e+17, 5.31912683e+17,
        5.33188481e+17, 5.34468963e+17, 5.35754205e+17, 5.37044074e+17,
        5.38338800e+17, 5.39638252e+17, 5.40942504e+17, 5.42251582e+17,
        5.43565508e+17, 5.44884307e+17, 5.46208004e+17, 5.47536675e+17,
        5.48870188e+17, 5.50208776e+17, 5.51552309e+17, 5.52900862e+17,
        5.54254461e+17, 5.55613131e+17, 5.56976898e+17, 5.58345787e+17,
        5.59719878e+17, 5.61099038e+17, 5.62483503e+17, 5.63873089e+17,
        5.65268033e+17, 5.66668205e+17, 5.68073682e+17, 5.69484493e+17,
        5.70900663e+17, 5.72322275e+17, 5.73749194e+17, 5.75181608e+17,
        5.76619545e+17, 5.78062928e+17, 5.79511782e+17, 5.80966296e+17,
        5.82426340e+17, 5.83891995e+17, 5.85363289e+17, 5.86840250e+17,
        5.88322964e+17, 5.89811297e+17, 5.91305437e+17, 5.92805419e+17,
        5.94311107e+17, 5.95822749e+17, 5.97340213e+17, 5.98863582e+17,
        6.00392831e+17, 6.01928156e+17, 6.03469425e+17, 6.05016778e+17,
        6.06570081e+17, 6.08129529e+17, 6.09695100e+17, 6.11266881e+17,
        6.12844738e+17, 6.14428924e+17, 6.16019251e+17, 6.17615972e+17,
        6.19218955e+17, 6.20828287e+17, 6.22444001e+17, 6.24066130e+17,
        6.25694709e+17, 6.27329771e+17, 6.28971350e+17, 6.30619480e+17,
        6.32274196e+17, 6.33935589e+17, 6.35603526e+17, 6.37278263e+17,
        6.38959615e+17, 6.40647838e+17, 6.42342746e+17, 6.44044598e+17,
        6.45753264e+17, 6.47468834e+17, 6.49191288e+17, 6.50920831e+17,
        6.52657335e+17, 6.54400891e+17, 6.56151537e+17, 6.57909309e+17,
        6.59674246e+17, 6.61446443e+17, 6.63225770e+17, 6.65012430e+17,
        6.66806467e+17, 6.68607749e+17, 6.70416483e+17, 6.72232711e+17,
        6.74056302e+17, 6.75887522e+17, 6.77726185e+17, 6.79572558e+17,
        6.81426512e+17, 6.83288088e+17, 6.85157496e+17, 6.87034609e+17,
        6.88919524e+17, 6.90812225e+17, 6.92712927e+17, 6.94621502e+17,
        6.96538048e+17, 6.98462608e+17, 7.00395225e+17, 7.02335942e+17,
        7.04284804e+17, 7.06241854e+17, 7.08207135e+17, 7.10180694e+17,
        7.12162573e+17, 7.14152876e+17, 7.16151476e+17, 7.18158588e+17,
        7.20174201e+17, 7.22198363e+17, 7.24231176e+17, 7.26272515e+17,
        7.28322597e+17, 7.30381470e+17, 7.32449010e+17, 7.34525433e+17,
        7.36610790e+17, 7.38704957e+17, 7.40808150e+17, 7.42920422e+17,
        7.45041648e+17, 7.47172104e+17, 7.49311613e+17, 7.51460451e+17,
        7.53618440e+17, 7.55785800e+17, 7.57962583e+17, 7.60148668e+17,
        7.62344332e+17, 7.64549399e+17, 7.66764146e+17, 7.68988399e+17,
        7.71222434e+17, 7.73466079e+17, 7.75719610e+17, 7.77982854e+17,
        7.80256088e+17, 7.82539141e+17, 7.84832232e+17, 7.87135418e+17,
        7.89448581e+17, 7.91771999e+17, 7.94105502e+17, 7.96449367e+17,
        7.98803424e+17, 8.01167952e+17, 8.03542781e+17, 8.05928133e+17,
        8.08324064e+17, 8.10730461e+17, 8.13147544e+17, 8.15575373e+17,
        8.18013833e+17, 8.20463146e+17, 8.22923370e+17, 8.25394392e+17,
        8.27876434e+17, 8.30369499e+17, 8.32873642e+17, 8.35388977e+17,
        8.37915391e+17, 8.40453105e+17, 8.43002124e+17, 8.45562502e+17,
        8.48134297e+17, 8.50717566e+17, 8.53312365e+17, 8.55918750e+17,
        8.58536779e+17, 8.61166507e+17, 8.63807991e+17, 8.66461289e+17,
        8.69126456e+17, 8.71803493e+17, 8.74492621e+17, 8.77193735e+17,
        8.79906888e+17, 8.82632300e+17, 8.85369865e+17, 8.88119637e+17,
        8.90881833e+17, 8.93656294e+17, 8.96443289e+17, 8.99242659e+17,
        9.02054616e+17, 9.04879216e+17, 9.07716354e+17, 9.10566238e+17,
        9.13428871e+17, 9.16304304e+17, 9.19192642e+17, 9.22093777e+17,
        9.25007914e+17, 9.27935001e+17, 9.30875242e+17, 9.33828532e+17,
        9.36794967e+17, 9.39774594e+17, 9.42767406e+17, 9.45773601e+17,
        9.48793020e+17, 9.51825907e+17, 9.54872101e+17, 9.57931842e+17,
        9.61004968e+17, 9.64091664e+17, 9.67191967e+17, 9.70305761e+17,
        9.73433224e+17, 9.76574338e+17, 9.79729131e+17, 9.82897627e+17,
        9.86079800e+17, 9.89275817e+17, 9.92485553e+17, 9.95709021e+17,
        9.98946375e+17, 1.00219744e+18, 1.00546240e+18, 1.00874108e+18,
        1.01203361e+18, 1.01533995e+18, 1.01866013e+18, 1.02199399e+18,
        1.02534161e+18, 1.02870305e+18, 1.03207816e+18, 1.03546694e+18,
        1.03886930e+18, 1.04228533e+18, 1.04571480e+18, 1.04915781e+18,
        1.05261413e+18, 1.05608381e+18, 1.05956675e+18, 1.06306273e+18,
        1.06657178e+18, 1.07009368e+18, 1.07362843e+18, 1.07717577e+18,
        1.08073553e+18, 1.08430768e+18, 1.08789186e+18, 1.09148802e+18,
        1.09509581e+18, 1.09871508e+18, 1.10234555e+18, 1.10598694e+18,
        1.10963894e+18, 1.11330123e+18, 1.11697342e+18, 1.12065524e+18,
        1.12434614e+18, 1.12804581e+18, 1.13175365e+18, 1.13546923e+18,
        1.13919197e+18, 1.14292128e+18, 1.14665651e+18, 1.15039698e+18,
        1.15414190e+18, 1.15789057e+18, 1.16164200e+18, 1.16539535e+18,
        1.16914963e+18, 1.17290366e+18, 1.17665637e+18, 1.18040647e+18,
        1.18415260e+18, 1.18789337e+18, 1.19162711e+18, 1.19535222e+18,
        1.19906678e+18, 1.20276889e+18, 1.20645639e+18, 1.21012696e+18,
        1.21377812e+18, 1.21740717e+18, 1.22101116e+18, 1.22458697e+18,
        1.22813110e+18, 1.23163988e+18, 1.23510922e+18, 1.23853472e+18,
        1.24191157e+18, 1.24523453e+18, 1.24849787e+18, 1.25169534e+18,
        1.25482002e+18, 1.25786435e+18, 1.26081990e+18, 1.26367733e+18,
        1.26642617e+18, 1.26905462e+18, 1.27154974e+18, 1.27386689e+18,
        1.27597413e+18, 1.27785804e+18, 1.27951304e+18, 1.28093715e+18,
        1.28213204e+18, 1.28308674e+18, 1.28379051e+18, 1.28423577e+18,
        1.28441794e+18, 1.28433502e+18, 1.28398728e+18, 1.28337685e+18,
        1.28250740e+18, 1.28138376e+18, 1.28001172e+18, 1.27839775e+18,
        1.27654882e+18, 1.27447228e+18, 1.27217577e+18, 1.26966709e+18,
        1.26695420e+18, 1.26404520e+18, 1.26094824e+18, 1.25767157e+18,
        1.25422352e+18, 1.25061249e+18, 1.24684694e+18, 1.24293540e+18,
        1.23888645e+18, 1.23470874e+18, 1.23041094e+18, 1.22600174e+18,
        1.22148985e+18, 1.21688396e+18, 1.21219272e+18, 1.20742473e+18,
        1.20258848e+18, 1.19769238e+18, 1.19274469e+18, 1.18775353e+18,
        1.18272681e+18, 1.17767226e+18, 1.17259737e+18, 1.16750938e+18,
        1.16241522e+18, 1.15732138e+18, 1.15223352e+18, 1.14715266e+18,
        1.14202253e+18, 1.13697118e+18, 1.13199077e+18, 1.12707963e+18,
        1.12223682e+18, 1.11746139e+18, 1.11275221e+18, 1.10810783e+18,
        1.10352651e+18, 1.09900617e+18, 1.09454439e+18, 1.09013845e+18,
        1.08578530e+18, 1.08148156e+18, 1.07722357e+18, 1.07300731e+18,
        1.06882850e+18, 1.06468254e+18, 1.06056454e+18, 1.05646931e+18,
        1.05239141e+18, 1.04832508e+18, 1.04426434e+18, 1.04020290e+18,
        1.03613427e+18, 1.03205168e+18, 1.02794814e+18, 1.02381644e+18,
        1.01964917e+18, 1.01543871e+18, 1.01117725e+18, 1.00685682e+18,
        1.00246906e+18, 9.98005448e+17, 9.93457325e+17, 9.88815899e+17,
        9.84072280e+17, 9.79217503e+17, 9.74242558e+17, 9.69138413e+17,
        9.63896041e+17, 9.58506447e+17, 9.52960711e+17, 9.47250004e+17,
        9.41365626e+17, 9.35299044e+17, 9.29041919e+17, 9.22586141e+17,
        9.15923872e+17, 9.09047572e+17, 9.01950045e+17, 8.94624467e+17,
        8.87064433e+17, 8.79263985e+17, 8.71217665e+17, 8.62920537e+17,
        8.54368234e+17, 8.45557000e+17, 8.36483715e+17, 8.27145948e+17,
        8.17541982e+17, 8.07670854e+17, 7.97532390e+17, 7.87127238e+17,
        7.76456900e+17, 7.65523756e+17, 7.54331098e+17, 7.42883157e+17,
        7.31185110e+17, 7.19243116e+17, 7.07064319e+17, 6.94656861e+17,
        6.82029893e+17, 6.69193577e+17, 6.56159070e+17, 6.42938530e+17,
        6.29545092e+17, 6.15992861e+17, 6.02296869e+17, 5.88473054e+17,
        5.74538220e+17, 5.60509990e+17, 5.46406751e+17, 5.32247606e+17,
        5.18052306e+17, 5.03841172e+17, 4.89635033e+17, 4.75455130e+17,
        4.61323035e+17, 4.47260552e+17, 4.33289615e+17, 4.19432182e+17,
        4.05710120e+17, 3.92145070e+17, 3.78758333e+17, 3.65570728e+17,
        3.52602438e+17, 3.39872880e+17, 3.27400539e+17, 3.15202828e+17,
        3.03295938e+17, 2.91694699e+17, 2.80412466e+17, 2.69461011e+17,
        2.58850464e+17, 2.48589281e+17, 2.38684256e+17, 2.29140594e+17,
        2.19962042e+17, 2.11151075e+17, 2.02709155e+17, 1.94637042e+17,
        1.86935135e+17, 1.79603799e+17, 1.72642312e+17, 1.66016593e+17,
        1.57198831e+17, 1.50029417e+17, 1.43323942e+17, 1.37034232e+17,
        1.23444669e+17, 1.11105346e+17, 1.00199045e+17, 9.06758772e+16,
        8.22858243e+16, 7.47539206e+16, 6.79224040e+16, 6.18169087e+16,
        5.60946901e+16, 5.05047680e+16, 4.52627120e+16, 4.04588879e+16,
        3.60307461e+16, 3.18890931e+16, 2.80431363e+16, 2.44302985e+16,
        2.11315170e+16, 1.81260303e+16, 1.54000164e+16, 1.29496040e+16,
        1.07453474e+16, 8.76595927e+15, 7.00500446e+15, 5.42062451e+15,
        3.98988996e+15, 2.65667160e+15, 1.36940104e+15, 0])

        hfct_obj = HeatFlowCouplingTools()
        fluid_Te = np.loadtxt(os.path.join(myPath,'test/gdhohlraum_xmic_5ps_TekeV_interp'), skiprows=1)[:, 1] * 1E3 *(ELEMENTARY_CHARGE/BOLTZMANN_CONSTANT)
        fluid_Z = np.loadtxt(os.path.join(myPath,'test/gdhohlraum_xmic_Z_interp'), skiprows=1)[:, 1] 
        fluid_x_grid = np.loadtxt(os.path.join(myPath,'test/gdhohlraum_xmic_5ps_separatedsnbWcm2'), skiprows=1)[:, 0]*1e-6
        fluid_ne = np.loadtxt(os.path.join(myPath,'test/gdhohlraum_xmic_ne1e20cm3_interp'), skiprows=1)[:, 1] * (1e20 * 1e6)
        fluid_x_centered_grid = np.array([(fluid_x_grid[i] + fluid_x_grid[i+1])/2 for i in range(len(fluid_x_grid) - 1)])
        hfct_obj.electron_temperature = fluid_Te
        hfct_obj.electron_number_density = fluid_ne
        hfct_obj.zbar = fluid_Z
        hfct_obj.cell_wall_coord = fluid_x_grid
        hfct_obj.cell_centered_coord = fluid_x_centered_grid
        # hfct_obj.mass = fluid_mass
        hfct_obj.lambda_ei(hfct_obj.electron_temperature * (BOLTZMANN_CONSTANT/ELEMENTARY_CHARGE), 
                            hfct_obj.electron_number_density,
                            hfct_obj.zbar)
        hfct_obj.spitzerHarmHeatFlow()
        hfct_obj.snb = True
        hfct_obj.snb_heat_flow(20, 12, 2)
        rel_err = 100*(abs(hfct_obj.q_snb - true_snb_heat_flow) / true_snb_heat_flow)
        assert np.nanmax(rel_err) < 1e-6

# Size of variable arrays:
sizeAlgebraic = 45
sizeStates = 11
sizeConstants = 64
from math import *
from numpy import *

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (millisecond)"
    legend_algebraic[3] = "SOVFThick in component sarcomere_geometry (dimensionless)"
    legend_algebraic[4] = "SOVFThin in component sarcomere_geometry (dimensionless)"
    legend_algebraic[0] = "sovr_ze in component sarcomere_geometry (micrometre)"
    legend_algebraic[1] = "sovr_cle in component sarcomere_geometry (micrometre)"
    legend_algebraic[2] = "len_sovr in component sarcomere_geometry (micrometre)"
    legend_constants[0] = "SLmax in component normalised_active_and_passive_force (micrometre)"
    legend_constants[1] = "SLmin in component normalised_active_and_passive_force (micrometre)"
    legend_constants[2] = "len_thin in component model_parameters (micrometre)"
    legend_constants[3] = "len_thick in component model_parameters (micrometre)"
    legend_constants[4] = "len_hbare in component model_parameters (micrometre)"
    legend_states[0] = "SL in component normalised_active_and_passive_force (micrometre)"
    legend_states[1] = "TRPNCaL in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)"
    legend_states[2] = "TRPNCaH in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)"
    legend_algebraic[42] = "dTRPNCaL in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)"
    legend_algebraic[43] = "dTRPNCaH in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)"
    legend_algebraic[10] = "kn_pT in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)"
    legend_algebraic[14] = "kp_nT in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)"
    legend_constants[55] = "konT in component Ca_binding_to_troponin_to_thin_filament_regulation (second_order_rate_constant)"
    legend_constants[56] = "koffLT in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)"
    legend_constants[57] = "koffHT in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)"
    legend_constants[5] = "Qkon in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)"
    legend_constants[6] = "Qkoff in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)"
    legend_constants[7] = "Qkn_p in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)"
    legend_constants[8] = "Qkp_n in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)"
    legend_constants[9] = "kon in component Ca_binding_to_troponin_to_thin_filament_regulation (second_order_rate_constant)"
    legend_constants[10] = "koffL in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)"
    legend_constants[11] = "koffH in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)"
    legend_constants[12] = "perm50 in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)"
    legend_constants[13] = "nperm in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)"
    legend_constants[14] = "kn_p in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)"
    legend_constants[15] = "kp_n in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)"
    legend_constants[16] = "koffmod in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)"
    legend_algebraic[6] = "Tropreg in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)"
    legend_algebraic[8] = "permtot in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)"
    legend_algebraic[12] = "inprmt in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)"
    legend_constants[17] = "TmpC in component model_parameters (celsius)"
    legend_algebraic[40] = "Cai in component equation_for_simulated_calcium_transient (micromolar)"
    legend_constants[58] = "fappT in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)"
    legend_algebraic[17] = "gappT in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)"
    legend_algebraic[20] = "hfT in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)"
    legend_algebraic[21] = "hbT in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)"
    legend_algebraic[23] = "gxbT in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)"
    legend_constants[18] = "fapp in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)"
    legend_constants[19] = "gapp in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)"
    legend_constants[20] = "hf in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)"
    legend_constants[21] = "hb in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)"
    legend_constants[22] = "gxb in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)"
    legend_constants[23] = "gslmod in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)"
    legend_algebraic[18] = "hfmd in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)"
    legend_algebraic[19] = "hbmd in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)"
    legend_constants[24] = "hfmdc in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)"
    legend_constants[25] = "hbmdc in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)"
    legend_constants[26] = "sigmap in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)"
    legend_constants[27] = "sigman in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)"
    legend_constants[28] = "xbmodsp in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)"
    legend_constants[29] = "Qfapp in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)"
    legend_constants[30] = "Qgapp in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)"
    legend_constants[31] = "Qhf in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)"
    legend_constants[32] = "Qhb in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)"
    legend_constants[33] = "Qgxb in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)"
    legend_algebraic[22] = "gxbmd in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)"
    legend_algebraic[16] = "gapslmd in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)"
    legend_constants[34] = "x_0 in component model_parameters (micrometre)"
    legend_states[3] = "xXBpostr in component mean_strain_of_strongly_bound_states (micrometre)"
    legend_states[4] = "xXBprer in component mean_strain_of_strongly_bound_states (micrometre)"
    legend_states[5] = "XBpostr in component regulation_and_crossbridge_cycling_state_equations (dimensionless)"
    legend_states[6] = "XBprer in component regulation_and_crossbridge_cycling_state_equations (dimensionless)"
    legend_algebraic[24] = "dXBpostr in component regulation_and_crossbridge_cycling_state_equations (first_order_rate_constant)"
    legend_algebraic[26] = "dXBprer in component regulation_and_crossbridge_cycling_state_equations (first_order_rate_constant)"
    legend_states[7] = "N_NoXB in component regulation_and_crossbridge_cycling_state_equations (dimensionless)"
    legend_states[8] = "P_NoXB in component regulation_and_crossbridge_cycling_state_equations (dimensionless)"
    legend_algebraic[25] = "P in component regulation_and_crossbridge_cycling_state_equations (dimensionless)"
    legend_states[9] = "N in component regulation_and_crossbridge_cycling_state_equations (dimensionless)"
    legend_algebraic[31] = "dxXBpostr in component mean_strain_of_strongly_bound_states (micrometre_per_millisecond)"
    legend_algebraic[30] = "dxXBprer in component mean_strain_of_strongly_bound_states (micrometre_per_millisecond)"
    legend_constants[35] = "xPsi in component mean_strain_of_strongly_bound_states (dimensionless)"
    legend_algebraic[27] = "dutyprer in component mean_strain_of_strongly_bound_states (dimensionless)"
    legend_algebraic[28] = "dutypostr in component mean_strain_of_strongly_bound_states (dimensionless)"
    legend_algebraic[29] = "dSL in component normalised_active_and_passive_force (micrometre_per_millisecond)"
    legend_constants[61] = "SSXBpostr in component normalised_active_and_passive_force (dimensionless)"
    legend_constants[59] = "SSXBprer in component normalised_active_and_passive_force (dimensionless)"
    legend_constants[36] = "kxb in component normalised_active_and_passive_force (millinewton_per_millimetre2)"
    legend_constants[62] = "Fnordv in component normalised_active_and_passive_force (millinewton_micrometre_per_millimetre2)"
    legend_algebraic[5] = "force in component normalised_active_and_passive_force (millinewton_micrometre_per_millimetre2)"
    legend_algebraic[7] = "active in component normalised_active_and_passive_force (unit_normalised_force)"
    legend_algebraic[13] = "ppforce in component normalised_active_and_passive_force (unit_normalised_force)"
    legend_algebraic[9] = "ppforce_t in component normalised_active_and_passive_force (unit_normalised_force)"
    legend_algebraic[11] = "ppforce_c in component normalised_active_and_passive_force (unit_normalised_force)"
    legend_constants[63] = "preload in component normalised_active_and_passive_force (unit_normalised_force)"
    legend_algebraic[15] = "afterload in component normalised_active_and_passive_force (unit_normalised_force)"
    legend_states[10] = "intf in component normalised_active_and_passive_force (unit_normalised_force_millisecond)"
    legend_constants[37] = "SL_c in component normalised_active_and_passive_force (micrometre)"
    legend_constants[38] = "SLrest in component normalised_active_and_passive_force (micrometre)"
    legend_constants[39] = "SLset in component normalised_active_and_passive_force (micrometre)"
    legend_constants[40] = "PCon_t in component normalised_active_and_passive_force (unit_normalised_force)"
    legend_constants[41] = "PExp_t in component normalised_active_and_passive_force (per_micrometre)"
    legend_constants[42] = "PCon_c in component normalised_active_and_passive_force (unit_normalised_force)"
    legend_constants[43] = "PExp_c in component normalised_active_and_passive_force (per_micrometre)"
    legend_constants[44] = "massf in component normalised_active_and_passive_force (unit_normalised_force_millisecond2_per_micrometre)"
    legend_constants[45] = "visc in component normalised_active_and_passive_force (unit_normalised_force_millisecond_per_micrometre)"
    legend_constants[46] = "KSE in component normalised_active_and_passive_force (unit_normalised_force_per_micrometre)"
    legend_constants[47] = "SEon in component normalised_active_and_passive_force (dimensionless)"
    legend_algebraic[32] = "FrSBXB in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (dimensionless)"
    legend_algebraic[33] = "dFrSBXB in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (first_order_rate_constant)"
    legend_algebraic[35] = "dsovr_ze in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (micrometre_per_millisecond)"
    legend_algebraic[36] = "dsovr_cle in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (micrometre_per_millisecond)"
    legend_algebraic[37] = "dlen_sovr in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (micrometre_per_millisecond)"
    legend_algebraic[39] = "dSOVFThick in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (first_order_rate_constant)"
    legend_algebraic[38] = "dSOVFThin in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (first_order_rate_constant)"
    legend_constants[48] = "kxb in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (millinewton_per_millimetre2)"
    legend_algebraic[41] = "dforce in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (millinewton_micrometre_per_millimetre2_per_millisecond)"
    legend_constants[49] = "Trop_conc in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (micromolar)"
    legend_algebraic[34] = "TropTot in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (micromolar)"
    legend_algebraic[44] = "dTropTot in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (micromolar_per_millisecond)"
    legend_constants[60] = "beta in component equation_for_simulated_calcium_transient (dimensionless)"
    legend_constants[50] = "tau1 in component equation_for_simulated_calcium_transient (millisecond)"
    legend_constants[51] = "tau2 in component equation_for_simulated_calcium_transient (millisecond)"
    legend_constants[52] = "start_time in component equation_for_simulated_calcium_transient (millisecond)"
    legend_constants[53] = "Ca_amplitude in component equation_for_simulated_calcium_transient (micromolar)"
    legend_constants[54] = "Ca_diastolic in component equation_for_simulated_calcium_transient (micromolar)"
    legend_rates[1] = "d/dt TRPNCaL in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)"
    legend_rates[2] = "d/dt TRPNCaH in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)"
    legend_rates[7] = "d/dt N_NoXB in component regulation_and_crossbridge_cycling_state_equations (dimensionless)"
    legend_rates[8] = "d/dt P_NoXB in component regulation_and_crossbridge_cycling_state_equations (dimensionless)"
    legend_rates[9] = "d/dt N in component regulation_and_crossbridge_cycling_state_equations (dimensionless)"
    legend_rates[6] = "d/dt XBprer in component regulation_and_crossbridge_cycling_state_equations (dimensionless)"
    legend_rates[5] = "d/dt XBpostr in component regulation_and_crossbridge_cycling_state_equations (dimensionless)"
    legend_rates[4] = "d/dt xXBprer in component mean_strain_of_strongly_bound_states (micrometre)"
    legend_rates[3] = "d/dt xXBpostr in component mean_strain_of_strongly_bound_states (micrometre)"
    legend_rates[0] = "d/dt SL in component normalised_active_and_passive_force (micrometre)"
    legend_rates[10] = "d/dt intf in component normalised_active_and_passive_force (unit_normalised_force_millisecond)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    constants[0] = 2.4
    constants[1] = 1.4
    constants[2] = 1.2
    constants[3] = 1.65
    constants[4] = 0.1
    states[0] = 1.89999811516093
    states[1] = 0.0147730085063734
    states[2] = 0.13066096561522
    constants[5] = 1.5
    constants[6] = 1.3
    constants[7] = 1.6
    constants[8] = 1.6
    constants[9] = 0.05
    constants[10] = 0.25
    constants[11] = 0.025
    constants[12] = 0.5
    constants[13] = 15
    constants[14] = 0.5
    constants[15] = 0.05
    constants[16] = 1
    constants[17] = 24
    constants[18] = 0.5
    constants[19] = 0.07
    constants[20] = 2
    constants[21] = 0.4
    constants[22] = 0.07
    constants[23] = 6
    constants[24] = 5
    constants[25] = 0
    constants[26] = 8
    constants[27] = 1
    constants[28] = 1
    constants[29] = 6.25
    constants[30] = 2.5
    constants[31] = 6.25
    constants[32] = 6.25
    constants[33] = 6.25
    constants[34] = 0.007
    states[3] = 0.00700005394873882
    states[4] = 3.41212828972468e-8
    states[5] = 1.81017564383744e-6
    states[6] = 3.0494964880038e-7
    states[7] = 0.999999959256274
    states[8] = 4.07437173988636e-8
    states[9] = 0.999997834540066
    constants[35] = 2
    constants[36] = 120
    states[10] = -4.5113452510363e-6
    constants[37] = 2.25
    constants[38] = 1.85
    constants[39] = 1.9
    constants[40] = 0.002
    constants[41] = 10
    constants[42] = 0.02
    constants[43] = 70
    constants[44] = 50
    constants[45] = 3
    constants[46] = 1
    constants[47] = 1
    constants[48] = 120
    constants[49] = 70
    constants[50] = 20
    constants[51] = 110
    constants[52] = 5
    constants[53] = 1.45
    constants[54] = 0.09
    constants[55] = constants[9]*constants[5]**(constants[17]-37.0000)/10.0000
    constants[56] = constants[10]*constants[16]*constants[6]**(constants[17]-37.0000)/10.0000
    constants[57] = constants[11]*constants[16]*constants[6]**(constants[17]-37.0000)/10.0000
    constants[58] = constants[18]*constants[28]*constants[29]**(constants[17]-37.0000)/10.0000
    constants[59] = (constants[21]*constants[18]+constants[22]*constants[18])/(constants[18]*constants[20]+constants[22]*constants[20]+constants[22]*constants[19]+constants[21]*constants[18]+constants[21]*constants[19]+constants[22]*constants[18])
    constants[60] = constants[50]/constants[51]**-1.00000/(constants[50]/constants[51]-1.00000)-constants[50]/constants[51]**-1.00000/(1.00000-constants[51]/constants[50])
    constants[61] = constants[18]*constants[20]/(constants[18]*constants[20]+constants[22]*constants[20]+constants[22]*constants[19]+constants[21]*constants[18]+constants[21]*constants[19]+constants[22]*constants[18])
    constants[62] = constants[36]*constants[34]*constants[61]
    constants[63] = fabs(constants[39]-constants[38])/(constants[39]-constants[38])*constants[40]*(exp(constants[41]*fabs(constants[39]-constants[38]))-1.00000)
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    algebraic[0] = custom_piecewise([less(constants[3]/2.00000 , states[0]/2.00000), constants[3]/2.00000 , True, states[0]/2.00000])
    algebraic[1] = custom_piecewise([greater(states[0]/2.00000-(states[0]-constants[2]) , constants[4]/2.00000), states[0]/2.00000-(states[0]-constants[2]) , True, constants[4]/2.00000])
    algebraic[2] = algebraic[0]-algebraic[1]
    algebraic[4] = algebraic[2]/constants[2]
    algebraic[6] = (1.00000-algebraic[4])*states[1]+algebraic[4]*states[2]
    algebraic[8] = fabs(1.00000/(1.00000+constants[12]/algebraic[6]**constants[13]))**(1.0/2)
    algebraic[10] = constants[14]*algebraic[8]*constants[7]**(constants[17]-37.0000)/10.0000
    algebraic[12] = custom_piecewise([less(1.00000/algebraic[8] , 100.000), 1.00000/algebraic[8] , True, 100.000])
    algebraic[14] = constants[15]*algebraic[12]*constants[8]**(constants[17]-37.0000)/10.0000
    rates[7] = algebraic[14]*states[8]-algebraic[10]*states[7]
    rates[8] = algebraic[10]*states[7]-algebraic[14]*states[8]
    algebraic[3] = algebraic[2]*2.00000/(constants[3]-constants[4])
    algebraic[5] = constants[36]*algebraic[3]*(states[3]*states[5]+states[4]*states[6])
    algebraic[7] = 1.00000*algebraic[5]/constants[62]
    algebraic[9] = (states[0]-constants[38])/fabs(states[0]-constants[38])*constants[40]*(exp(constants[41]*fabs(states[0]-constants[38]))-1.00000)
    algebraic[11] = custom_piecewise([greater(states[0] , constants[37]), constants[42]*(exp(constants[43]*fabs(states[0]-constants[37]))-1.00000) , True, 0.00000])
    algebraic[13] = algebraic[9]+algebraic[11]
    algebraic[15] = custom_piecewise([equal(constants[47] , 1.00000), constants[46]*(constants[39]-states[0]) , True, 0.00000])
    rates[10] = constants[63]+algebraic[15]-(algebraic[13]+algebraic[7])
    algebraic[18] = exp(-states[4]/fabs(states[4])*constants[24]*states[4]/constants[34]**2.00000)
    algebraic[20] = constants[20]*algebraic[18]*constants[28]*constants[31]**(constants[17]-37.0000)/10.0000
    algebraic[19] = exp((states[3]-constants[34])/fabs(states[3]-constants[34])*constants[25]*(states[3]-constants[34])/constants[34]**2.00000)
    algebraic[21] = constants[21]*algebraic[19]*constants[28]*constants[32]**(constants[17]-37.0000)/10.0000
    algebraic[22] = custom_piecewise([less(states[3] , constants[34]), exp(constants[26]*(constants[34]-states[3])/constants[34]**2.00000) , True, exp(constants[27]*(states[3]-constants[34])/constants[34]**2.00000)])
    algebraic[23] = constants[22]*algebraic[22]*constants[28]*constants[33]**(constants[17]-37.0000)/10.0000
    algebraic[24] = algebraic[20]*states[6]-(algebraic[21]*states[5]+algebraic[23]*states[5])
    rates[5] = algebraic[24]
    algebraic[25] = 1.00000-states[9]-states[6]-states[5]
    rates[9] = algebraic[14]*algebraic[25]-algebraic[10]*states[9]
    algebraic[16] = 1.00000+(1.00000-algebraic[3])*constants[23]
    algebraic[17] = constants[19]*algebraic[16]*constants[28]*constants[30]**(constants[17]-37.0000)/10.0000
    algebraic[26] = constants[58]*algebraic[25]+algebraic[21]*states[5]-(algebraic[17]*states[6]+algebraic[20]*states[6])
    rates[6] = algebraic[26]
    algebraic[29] = custom_piecewise([less_equal(states[0] , constants[0]) & greater(states[0] , constants[1]), (states[10]+(constants[39]-states[0])*constants[45])/constants[44] , True, 0.00000])
    rates[0] = algebraic[29]
    algebraic[27] = (algebraic[21]*constants[58]+algebraic[23]*constants[58])/(constants[58]*algebraic[20]+algebraic[23]*algebraic[20]+algebraic[23]*algebraic[17]+algebraic[21]*constants[58]+algebraic[21]*algebraic[17]+algebraic[23]*constants[58])
    algebraic[30] = algebraic[29]/2.00000+constants[35]/algebraic[27]*(constants[58]*-states[4]+algebraic[21]*(states[3]-constants[34]+states[4]))
    rates[4] = algebraic[30]
    algebraic[28] = constants[58]*algebraic[20]/(constants[58]*algebraic[20]+algebraic[23]*algebraic[20]+algebraic[23]*algebraic[17]+algebraic[21]*constants[58]+algebraic[21]*algebraic[17]+algebraic[23]*constants[58])
    algebraic[31] = algebraic[29]/2.00000+constants[35]/algebraic[28]*algebraic[20]*(states[4]+constants[34]-states[3])
    rates[3] = algebraic[31]
    algebraic[40] = custom_piecewise([greater(voi , constants[52]), (constants[53]-constants[54])/constants[60]*(exp(-(voi-constants[52])/constants[50])-exp(-(voi-constants[52])/constants[51]))+constants[54] , True, constants[54]])
    algebraic[42] = constants[55]*algebraic[40]*(1.00000-states[1])-constants[56]*states[1]
    rates[1] = algebraic[42]
    algebraic[43] = constants[55]*algebraic[40]*(1.00000-states[2])-constants[57]*states[2]
    rates[2] = algebraic[43]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([less(constants[3]/2.00000 , states[0]/2.00000), constants[3]/2.00000 , True, states[0]/2.00000])
    algebraic[1] = custom_piecewise([greater(states[0]/2.00000-(states[0]-constants[2]) , constants[4]/2.00000), states[0]/2.00000-(states[0]-constants[2]) , True, constants[4]/2.00000])
    algebraic[2] = algebraic[0]-algebraic[1]
    algebraic[4] = algebraic[2]/constants[2]
    algebraic[6] = (1.00000-algebraic[4])*states[1]+algebraic[4]*states[2]
    algebraic[8] = fabs(1.00000/(1.00000+constants[12]/algebraic[6]**constants[13]))**(1.0/2)
    algebraic[10] = constants[14]*algebraic[8]*constants[7]**(constants[17]-37.0000)/10.0000
    algebraic[12] = custom_piecewise([less(1.00000/algebraic[8] , 100.000), 1.00000/algebraic[8] , True, 100.000])
    algebraic[14] = constants[15]*algebraic[12]*constants[8]**(constants[17]-37.0000)/10.0000
    algebraic[3] = algebraic[2]*2.00000/(constants[3]-constants[4])
    algebraic[5] = constants[36]*algebraic[3]*(states[3]*states[5]+states[4]*states[6])
    algebraic[7] = 1.00000*algebraic[5]/constants[62]
    algebraic[9] = (states[0]-constants[38])/fabs(states[0]-constants[38])*constants[40]*(exp(constants[41]*fabs(states[0]-constants[38]))-1.00000)
    algebraic[11] = custom_piecewise([greater(states[0] , constants[37]), constants[42]*(exp(constants[43]*fabs(states[0]-constants[37]))-1.00000) , True, 0.00000])
    algebraic[13] = algebraic[9]+algebraic[11]
    algebraic[15] = custom_piecewise([equal(constants[47] , 1.00000), constants[46]*(constants[39]-states[0]) , True, 0.00000])
    algebraic[18] = exp(-states[4]/fabs(states[4])*constants[24]*states[4]/constants[34]**2.00000)
    algebraic[20] = constants[20]*algebraic[18]*constants[28]*constants[31]**(constants[17]-37.0000)/10.0000
    algebraic[19] = exp((states[3]-constants[34])/fabs(states[3]-constants[34])*constants[25]*(states[3]-constants[34])/constants[34]**2.00000)
    algebraic[21] = constants[21]*algebraic[19]*constants[28]*constants[32]**(constants[17]-37.0000)/10.0000
    algebraic[22] = custom_piecewise([less(states[3] , constants[34]), exp(constants[26]*(constants[34]-states[3])/constants[34]**2.00000) , True, exp(constants[27]*(states[3]-constants[34])/constants[34]**2.00000)])
    algebraic[23] = constants[22]*algebraic[22]*constants[28]*constants[33]**(constants[17]-37.0000)/10.0000
    algebraic[24] = algebraic[20]*states[6]-(algebraic[21]*states[5]+algebraic[23]*states[5])
    algebraic[25] = 1.00000-states[9]-states[6]-states[5]
    algebraic[16] = 1.00000+(1.00000-algebraic[3])*constants[23]
    algebraic[17] = constants[19]*algebraic[16]*constants[28]*constants[30]**(constants[17]-37.0000)/10.0000
    algebraic[26] = constants[58]*algebraic[25]+algebraic[21]*states[5]-(algebraic[17]*states[6]+algebraic[20]*states[6])
    algebraic[29] = custom_piecewise([less_equal(states[0] , constants[0]) & greater(states[0] , constants[1]), (states[10]+(constants[39]-states[0])*constants[45])/constants[44] , True, 0.00000])
    algebraic[27] = (algebraic[21]*constants[58]+algebraic[23]*constants[58])/(constants[58]*algebraic[20]+algebraic[23]*algebraic[20]+algebraic[23]*algebraic[17]+algebraic[21]*constants[58]+algebraic[21]*algebraic[17]+algebraic[23]*constants[58])
    algebraic[30] = algebraic[29]/2.00000+constants[35]/algebraic[27]*(constants[58]*-states[4]+algebraic[21]*(states[3]-constants[34]+states[4]))
    algebraic[28] = constants[58]*algebraic[20]/(constants[58]*algebraic[20]+algebraic[23]*algebraic[20]+algebraic[23]*algebraic[17]+algebraic[21]*constants[58]+algebraic[21]*algebraic[17]+algebraic[23]*constants[58])
    algebraic[31] = algebraic[29]/2.00000+constants[35]/algebraic[28]*algebraic[20]*(states[4]+constants[34]-states[3])
    algebraic[40] = custom_piecewise([greater(voi , constants[52]), (constants[53]-constants[54])/constants[60]*(exp(-(voi-constants[52])/constants[50])-exp(-(voi-constants[52])/constants[51]))+constants[54] , True, constants[54]])
    algebraic[42] = constants[55]*algebraic[40]*(1.00000-states[1])-constants[56]*states[1]
    algebraic[43] = constants[55]*algebraic[40]*(1.00000-states[2])-constants[57]*states[2]
    algebraic[32] = (states[5]+states[6])/(constants[61]+constants[59])
    algebraic[33] = (algebraic[24]+algebraic[26])/(constants[61]+constants[59])
    algebraic[34] = constants[49]*((1.00000-algebraic[4])*states[1]+algebraic[4]*(algebraic[32]*states[2]+(1.00000-algebraic[32])*states[1]))
    algebraic[35] = custom_piecewise([less(states[0] , constants[3]), -0.500000*algebraic[29] , True, 0.00000])
    algebraic[36] = custom_piecewise([greater(2.00000*constants[2]-states[0] , constants[4]), -0.500000*algebraic[29] , True, 0.00000])
    algebraic[37] = algebraic[35]-algebraic[36]
    algebraic[38] = algebraic[37]/constants[2]
    algebraic[39] = 2.00000*algebraic[37]/(constants[3]-constants[4])
    algebraic[41] = constants[48]*algebraic[39]*(states[3]*states[5]+states[4]*states[6])+constants[48]*algebraic[3]*(algebraic[31]*states[5]+states[3]*algebraic[24]+algebraic[30]*states[6]+states[4]*algebraic[26])
    algebraic[44] = constants[49]*(-algebraic[38]*states[1]+(1.00000-algebraic[4])*algebraic[42]+algebraic[38]*(algebraic[32]*states[2]+(1.00000-algebraic[32])*states[1])+algebraic[4]*(algebraic[33]*states[2]+algebraic[32]*algebraic[43]+(1.00000-algebraic[32])*algebraic[42]-algebraic[33]*states[1]))
    return algebraic

def custom_piecewise(cases):
    """Compute result of a piecewise function"""
    return select(cases[0::2],cases[1::2])

def solve_model():
    """Solve model with ODE solver"""
    from scipy.integrate import ode
    # Initialise constants and state variables
    (init_states, constants) = initConsts()

    # Set timespan to solve over
    voi = linspace(0, 10, 500)

    # Construct ODE object to solve
    r = ode(computeRates)
    r.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=1)
    r.set_initial_value(init_states, voi[0])
    r.set_f_params(constants)

    # Solve model
    states = array([[0.0] * len(voi)] * sizeStates)
    states[:,0] = init_states
    for (i,t) in enumerate(voi[1:]):
        if r.successful():
            r.integrate(t)
            states[:,i+1] = r.y
        else:
            break

    # Compute algebraic variables
    algebraic = computeAlgebraic(constants, states, voi)
    return (voi, states, algebraic)

def plot_model(voi, states, algebraic):
    """Plot variables against variable of integration"""
    import pylab
    (legend_states, legend_algebraic, legend_voi, legend_constants) = createLegends()
    pylab.figure(1)
    pylab.plot(voi,vstack((states,algebraic)).T)
    pylab.xlabel(legend_voi)
    pylab.legend(legend_states + legend_algebraic, loc='best')
    pylab.show()

if __name__ == "__main__":
    (voi, states, algebraic) = solve_model()
    plot_model(voi, states, algebraic)


function [VOI, STATES, ALGEBRAIC, CONSTANTS] = mainFunction()
    % This is the "main function".  In Matlab, things work best if you rename this function to match the filename.
   [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel();
end

function [algebraicVariableCount] = getAlgebraicVariableCount() 
    % Used later when setting a global variable with the number of algebraic variables.
    % Note: This is not the "main method".  
    algebraicVariableCount =45;
end
% There are a total of 11 entries in each of the rate and state variable arrays.
% There are a total of 64 entries in the constant variable array.
%

function [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel()
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  algebraicVariableCount = getAlgebraicVariableCount();
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;

    % Set timespan to solve over 
    tspan = [0, 10];

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 1);

    % Solve model with ODE solver
    [VOI, STATES] = ode15s(@(VOI, STATES)computeRates(VOI, STATES, CONSTANTS), tspan, INIT_STATES, options);

    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI);

    % Plot state variables against variable of integration
    [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends();
    figure();
    plot(VOI, STATES);
    xlabel(LEGEND_VOI);
    l = legend(LEGEND_STATES);
    set(l,'Interpreter','none');
end

function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_VOI = ''; LEGEND_CONSTANTS = '';
    LEGEND_VOI = strpad('time in component environment (millisecond)');
    LEGEND_ALGEBRAIC(:,4) = strpad('SOVFThick in component sarcomere_geometry (dimensionless)');
    LEGEND_ALGEBRAIC(:,5) = strpad('SOVFThin in component sarcomere_geometry (dimensionless)');
    LEGEND_ALGEBRAIC(:,1) = strpad('sovr_ze in component sarcomere_geometry (micrometre)');
    LEGEND_ALGEBRAIC(:,2) = strpad('sovr_cle in component sarcomere_geometry (micrometre)');
    LEGEND_ALGEBRAIC(:,3) = strpad('len_sovr in component sarcomere_geometry (micrometre)');
    LEGEND_CONSTANTS(:,1) = strpad('SLmax in component normalised_active_and_passive_force (micrometre)');
    LEGEND_CONSTANTS(:,2) = strpad('SLmin in component normalised_active_and_passive_force (micrometre)');
    LEGEND_CONSTANTS(:,3) = strpad('len_thin in component model_parameters (micrometre)');
    LEGEND_CONSTANTS(:,4) = strpad('len_thick in component model_parameters (micrometre)');
    LEGEND_CONSTANTS(:,5) = strpad('len_hbare in component model_parameters (micrometre)');
    LEGEND_STATES(:,1) = strpad('SL in component normalised_active_and_passive_force (micrometre)');
    LEGEND_STATES(:,2) = strpad('TRPNCaL in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)');
    LEGEND_STATES(:,3) = strpad('TRPNCaH in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)');
    LEGEND_ALGEBRAIC(:,43) = strpad('dTRPNCaL in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)');
    LEGEND_ALGEBRAIC(:,44) = strpad('dTRPNCaH in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)');
    LEGEND_ALGEBRAIC(:,11) = strpad('kn_pT in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)');
    LEGEND_ALGEBRAIC(:,15) = strpad('kp_nT in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,56) = strpad('konT in component Ca_binding_to_troponin_to_thin_filament_regulation (second_order_rate_constant)');
    LEGEND_CONSTANTS(:,57) = strpad('koffLT in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,58) = strpad('koffHT in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,6) = strpad('Qkon in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)');
    LEGEND_CONSTANTS(:,7) = strpad('Qkoff in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)');
    LEGEND_CONSTANTS(:,8) = strpad('Qkn_p in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)');
    LEGEND_CONSTANTS(:,9) = strpad('Qkp_n in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)');
    LEGEND_CONSTANTS(:,10) = strpad('kon in component Ca_binding_to_troponin_to_thin_filament_regulation (second_order_rate_constant)');
    LEGEND_CONSTANTS(:,11) = strpad('koffL in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,12) = strpad('koffH in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,13) = strpad('perm50 in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)');
    LEGEND_CONSTANTS(:,14) = strpad('nperm in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)');
    LEGEND_CONSTANTS(:,15) = strpad('kn_p in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,16) = strpad('kp_n in component Ca_binding_to_troponin_to_thin_filament_regulation (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,17) = strpad('koffmod in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)');
    LEGEND_ALGEBRAIC(:,7) = strpad('Tropreg in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)');
    LEGEND_ALGEBRAIC(:,9) = strpad('permtot in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)');
    LEGEND_ALGEBRAIC(:,13) = strpad('inprmt in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)');
    LEGEND_CONSTANTS(:,18) = strpad('TmpC in component model_parameters (celsius)');
    LEGEND_ALGEBRAIC(:,41) = strpad('Cai in component equation_for_simulated_calcium_transient (micromolar)');
    LEGEND_CONSTANTS(:,59) = strpad('fappT in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)');
    LEGEND_ALGEBRAIC(:,18) = strpad('gappT in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)');
    LEGEND_ALGEBRAIC(:,21) = strpad('hfT in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)');
    LEGEND_ALGEBRAIC(:,22) = strpad('hbT in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)');
    LEGEND_ALGEBRAIC(:,24) = strpad('gxbT in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,19) = strpad('fapp in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,20) = strpad('gapp in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,21) = strpad('hf in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,22) = strpad('hb in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,23) = strpad('gxb in component thin_filament_regulation_and_crossbridge_cycling_rates (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,24) = strpad('gslmod in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)');
    LEGEND_ALGEBRAIC(:,19) = strpad('hfmd in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)');
    LEGEND_ALGEBRAIC(:,20) = strpad('hbmd in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)');
    LEGEND_CONSTANTS(:,25) = strpad('hfmdc in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)');
    LEGEND_CONSTANTS(:,26) = strpad('hbmdc in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)');
    LEGEND_CONSTANTS(:,27) = strpad('sigmap in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)');
    LEGEND_CONSTANTS(:,28) = strpad('sigman in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)');
    LEGEND_CONSTANTS(:,29) = strpad('xbmodsp in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)');
    LEGEND_CONSTANTS(:,30) = strpad('Qfapp in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)');
    LEGEND_CONSTANTS(:,31) = strpad('Qgapp in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)');
    LEGEND_CONSTANTS(:,32) = strpad('Qhf in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)');
    LEGEND_CONSTANTS(:,33) = strpad('Qhb in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)');
    LEGEND_CONSTANTS(:,34) = strpad('Qgxb in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)');
    LEGEND_ALGEBRAIC(:,23) = strpad('gxbmd in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)');
    LEGEND_ALGEBRAIC(:,17) = strpad('gapslmd in component thin_filament_regulation_and_crossbridge_cycling_rates (dimensionless)');
    LEGEND_CONSTANTS(:,35) = strpad('x_0 in component model_parameters (micrometre)');
    LEGEND_STATES(:,4) = strpad('xXBpostr in component mean_strain_of_strongly_bound_states (micrometre)');
    LEGEND_STATES(:,5) = strpad('xXBprer in component mean_strain_of_strongly_bound_states (micrometre)');
    LEGEND_STATES(:,6) = strpad('XBpostr in component regulation_and_crossbridge_cycling_state_equations (dimensionless)');
    LEGEND_STATES(:,7) = strpad('XBprer in component regulation_and_crossbridge_cycling_state_equations (dimensionless)');
    LEGEND_ALGEBRAIC(:,25) = strpad('dXBpostr in component regulation_and_crossbridge_cycling_state_equations (first_order_rate_constant)');
    LEGEND_ALGEBRAIC(:,27) = strpad('dXBprer in component regulation_and_crossbridge_cycling_state_equations (first_order_rate_constant)');
    LEGEND_STATES(:,8) = strpad('N_NoXB in component regulation_and_crossbridge_cycling_state_equations (dimensionless)');
    LEGEND_STATES(:,9) = strpad('P_NoXB in component regulation_and_crossbridge_cycling_state_equations (dimensionless)');
    LEGEND_ALGEBRAIC(:,26) = strpad('P in component regulation_and_crossbridge_cycling_state_equations (dimensionless)');
    LEGEND_STATES(:,10) = strpad('N in component regulation_and_crossbridge_cycling_state_equations (dimensionless)');
    LEGEND_ALGEBRAIC(:,32) = strpad('dxXBpostr in component mean_strain_of_strongly_bound_states (micrometre_per_millisecond)');
    LEGEND_ALGEBRAIC(:,31) = strpad('dxXBprer in component mean_strain_of_strongly_bound_states (micrometre_per_millisecond)');
    LEGEND_CONSTANTS(:,36) = strpad('xPsi in component mean_strain_of_strongly_bound_states (dimensionless)');
    LEGEND_ALGEBRAIC(:,28) = strpad('dutyprer in component mean_strain_of_strongly_bound_states (dimensionless)');
    LEGEND_ALGEBRAIC(:,29) = strpad('dutypostr in component mean_strain_of_strongly_bound_states (dimensionless)');
    LEGEND_ALGEBRAIC(:,30) = strpad('dSL in component normalised_active_and_passive_force (micrometre_per_millisecond)');
    LEGEND_CONSTANTS(:,62) = strpad('SSXBpostr in component normalised_active_and_passive_force (dimensionless)');
    LEGEND_CONSTANTS(:,60) = strpad('SSXBprer in component normalised_active_and_passive_force (dimensionless)');
    LEGEND_CONSTANTS(:,37) = strpad('kxb in component normalised_active_and_passive_force (millinewton_per_millimetre2)');
    LEGEND_CONSTANTS(:,63) = strpad('Fnordv in component normalised_active_and_passive_force (millinewton_micrometre_per_millimetre2)');
    LEGEND_ALGEBRAIC(:,6) = strpad('force in component normalised_active_and_passive_force (millinewton_micrometre_per_millimetre2)');
    LEGEND_ALGEBRAIC(:,8) = strpad('active in component normalised_active_and_passive_force (unit_normalised_force)');
    LEGEND_ALGEBRAIC(:,14) = strpad('ppforce in component normalised_active_and_passive_force (unit_normalised_force)');
    LEGEND_ALGEBRAIC(:,10) = strpad('ppforce_t in component normalised_active_and_passive_force (unit_normalised_force)');
    LEGEND_ALGEBRAIC(:,12) = strpad('ppforce_c in component normalised_active_and_passive_force (unit_normalised_force)');
    LEGEND_CONSTANTS(:,64) = strpad('preload in component normalised_active_and_passive_force (unit_normalised_force)');
    LEGEND_ALGEBRAIC(:,16) = strpad('afterload in component normalised_active_and_passive_force (unit_normalised_force)');
    LEGEND_STATES(:,11) = strpad('intf in component normalised_active_and_passive_force (unit_normalised_force_millisecond)');
    LEGEND_CONSTANTS(:,38) = strpad('SL_c in component normalised_active_and_passive_force (micrometre)');
    LEGEND_CONSTANTS(:,39) = strpad('SLrest in component normalised_active_and_passive_force (micrometre)');
    LEGEND_CONSTANTS(:,40) = strpad('SLset in component normalised_active_and_passive_force (micrometre)');
    LEGEND_CONSTANTS(:,41) = strpad('PCon_t in component normalised_active_and_passive_force (unit_normalised_force)');
    LEGEND_CONSTANTS(:,42) = strpad('PExp_t in component normalised_active_and_passive_force (per_micrometre)');
    LEGEND_CONSTANTS(:,43) = strpad('PCon_c in component normalised_active_and_passive_force (unit_normalised_force)');
    LEGEND_CONSTANTS(:,44) = strpad('PExp_c in component normalised_active_and_passive_force (per_micrometre)');
    LEGEND_CONSTANTS(:,45) = strpad('massf in component normalised_active_and_passive_force (unit_normalised_force_millisecond2_per_micrometre)');
    LEGEND_CONSTANTS(:,46) = strpad('visc in component normalised_active_and_passive_force (unit_normalised_force_millisecond_per_micrometre)');
    LEGEND_CONSTANTS(:,47) = strpad('KSE in component normalised_active_and_passive_force (unit_normalised_force_per_micrometre)');
    LEGEND_CONSTANTS(:,48) = strpad('SEon in component normalised_active_and_passive_force (dimensionless)');
    LEGEND_ALGEBRAIC(:,33) = strpad('FrSBXB in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (dimensionless)');
    LEGEND_ALGEBRAIC(:,34) = strpad('dFrSBXB in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (first_order_rate_constant)');
    LEGEND_ALGEBRAIC(:,36) = strpad('dsovr_ze in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (micrometre_per_millisecond)');
    LEGEND_ALGEBRAIC(:,37) = strpad('dsovr_cle in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (micrometre_per_millisecond)');
    LEGEND_ALGEBRAIC(:,38) = strpad('dlen_sovr in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (micrometre_per_millisecond)');
    LEGEND_ALGEBRAIC(:,40) = strpad('dSOVFThick in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (first_order_rate_constant)');
    LEGEND_ALGEBRAIC(:,39) = strpad('dSOVFThin in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (first_order_rate_constant)');
    LEGEND_CONSTANTS(:,49) = strpad('kxb in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (millinewton_per_millimetre2)');
    LEGEND_ALGEBRAIC(:,42) = strpad('dforce in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (millinewton_micrometre_per_millimetre2_per_millisecond)');
    LEGEND_CONSTANTS(:,50) = strpad('Trop_conc in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (micromolar)');
    LEGEND_ALGEBRAIC(:,35) = strpad('TropTot in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (micromolar)');
    LEGEND_ALGEBRAIC(:,45) = strpad('dTropTot in component calculation_of_micromolar_per_millisecondes_of_Ca_for_apparent_Ca_binding (micromolar_per_millisecond)');
    LEGEND_CONSTANTS(:,61) = strpad('beta in component equation_for_simulated_calcium_transient (dimensionless)');
    LEGEND_CONSTANTS(:,51) = strpad('tau1 in component equation_for_simulated_calcium_transient (millisecond)');
    LEGEND_CONSTANTS(:,52) = strpad('tau2 in component equation_for_simulated_calcium_transient (millisecond)');
    LEGEND_CONSTANTS(:,53) = strpad('start_time in component equation_for_simulated_calcium_transient (millisecond)');
    LEGEND_CONSTANTS(:,54) = strpad('Ca_amplitude in component equation_for_simulated_calcium_transient (micromolar)');
    LEGEND_CONSTANTS(:,55) = strpad('Ca_diastolic in component equation_for_simulated_calcium_transient (micromolar)');
    LEGEND_RATES(:,2) = strpad('d/dt TRPNCaL in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)');
    LEGEND_RATES(:,3) = strpad('d/dt TRPNCaH in component Ca_binding_to_troponin_to_thin_filament_regulation (dimensionless)');
    LEGEND_RATES(:,8) = strpad('d/dt N_NoXB in component regulation_and_crossbridge_cycling_state_equations (dimensionless)');
    LEGEND_RATES(:,9) = strpad('d/dt P_NoXB in component regulation_and_crossbridge_cycling_state_equations (dimensionless)');
    LEGEND_RATES(:,10) = strpad('d/dt N in component regulation_and_crossbridge_cycling_state_equations (dimensionless)');
    LEGEND_RATES(:,7) = strpad('d/dt XBprer in component regulation_and_crossbridge_cycling_state_equations (dimensionless)');
    LEGEND_RATES(:,6) = strpad('d/dt XBpostr in component regulation_and_crossbridge_cycling_state_equations (dimensionless)');
    LEGEND_RATES(:,5) = strpad('d/dt xXBprer in component mean_strain_of_strongly_bound_states (micrometre)');
    LEGEND_RATES(:,4) = strpad('d/dt xXBpostr in component mean_strain_of_strongly_bound_states (micrometre)');
    LEGEND_RATES(:,1) = strpad('d/dt SL in component normalised_active_and_passive_force (micrometre)');
    LEGEND_RATES(:,11) = strpad('d/dt intf in component normalised_active_and_passive_force (unit_normalised_force_millisecond)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end

function [STATES, CONSTANTS] = initConsts()
    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];
    CONSTANTS(:,1) = 2.4;
    CONSTANTS(:,2) = 1.4;
    CONSTANTS(:,3) = 1.2;
    CONSTANTS(:,4) = 1.65;
    CONSTANTS(:,5) = 0.1;
    STATES(:,1) = 1.89999811516093;
    STATES(:,2) = 0.0147730085063734;
    STATES(:,3) = 0.13066096561522;
    CONSTANTS(:,6) = 1.5;
    CONSTANTS(:,7) = 1.3;
    CONSTANTS(:,8) = 1.6;
    CONSTANTS(:,9) = 1.6;
    CONSTANTS(:,10) = 0.05;
    CONSTANTS(:,11) = 0.25;
    CONSTANTS(:,12) = 0.025;
    CONSTANTS(:,13) = 0.5;
    CONSTANTS(:,14) = 15;
    CONSTANTS(:,15) = 0.5;
    CONSTANTS(:,16) = 0.05;
    CONSTANTS(:,17) = 1;
    CONSTANTS(:,18) = 24;
    CONSTANTS(:,19) = 0.5;
    CONSTANTS(:,20) = 0.07;
    CONSTANTS(:,21) = 2;
    CONSTANTS(:,22) = 0.4;
    CONSTANTS(:,23) = 0.07;
    CONSTANTS(:,24) = 6;
    CONSTANTS(:,25) = 5;
    CONSTANTS(:,26) = 0;
    CONSTANTS(:,27) = 8;
    CONSTANTS(:,28) = 1;
    CONSTANTS(:,29) = 1;
    CONSTANTS(:,30) = 6.25;
    CONSTANTS(:,31) = 2.5;
    CONSTANTS(:,32) = 6.25;
    CONSTANTS(:,33) = 6.25;
    CONSTANTS(:,34) = 6.25;
    CONSTANTS(:,35) = 0.007;
    STATES(:,4) = 0.00700005394873882;
    STATES(:,5) = 3.41212828972468e-8;
    STATES(:,6) = 1.81017564383744e-6;
    STATES(:,7) = 3.0494964880038e-7;
    STATES(:,8) = 0.999999959256274;
    STATES(:,9) = 4.07437173988636e-8;
    STATES(:,10) = 0.999997834540066;
    CONSTANTS(:,36) = 2;
    CONSTANTS(:,37) = 120;
    STATES(:,11) = -4.5113452510363e-6;
    CONSTANTS(:,38) = 2.25;
    CONSTANTS(:,39) = 1.85;
    CONSTANTS(:,40) = 1.9;
    CONSTANTS(:,41) = 0.002;
    CONSTANTS(:,42) = 10;
    CONSTANTS(:,43) = 0.02;
    CONSTANTS(:,44) = 70;
    CONSTANTS(:,45) = 50;
    CONSTANTS(:,46) = 3;
    CONSTANTS(:,47) = 1;
    CONSTANTS(:,48) = 1;
    CONSTANTS(:,49) = 120;
    CONSTANTS(:,50) = 70;
    CONSTANTS(:,51) = 20;
    CONSTANTS(:,52) = 110;
    CONSTANTS(:,53) = 5;
    CONSTANTS(:,54) = 1.45;
    CONSTANTS(:,55) = 0.09;
    CONSTANTS(:,56) =  CONSTANTS(:,10).*CONSTANTS(:,6) .^ (CONSTANTS(:,18) - 37.0000)./10.0000;
    CONSTANTS(:,57) =  CONSTANTS(:,11).*CONSTANTS(:,17).*CONSTANTS(:,7) .^ (CONSTANTS(:,18) - 37.0000)./10.0000;
    CONSTANTS(:,58) =  CONSTANTS(:,12).*CONSTANTS(:,17).*CONSTANTS(:,7) .^ (CONSTANTS(:,18) - 37.0000)./10.0000;
    CONSTANTS(:,59) =  CONSTANTS(:,19).*CONSTANTS(:,29).*CONSTANTS(:,30) .^ (CONSTANTS(:,18) - 37.0000)./10.0000;
    CONSTANTS(:,60) = ( CONSTANTS(:,22).*CONSTANTS(:,19)+ CONSTANTS(:,23).*CONSTANTS(:,19))./( CONSTANTS(:,19).*CONSTANTS(:,21)+ CONSTANTS(:,23).*CONSTANTS(:,21)+ CONSTANTS(:,23).*CONSTANTS(:,20)+ CONSTANTS(:,22).*CONSTANTS(:,19)+ CONSTANTS(:,22).*CONSTANTS(:,20)+ CONSTANTS(:,23).*CONSTANTS(:,19));
    CONSTANTS(:,61) = CONSTANTS(:,51)./CONSTANTS(:,52) .^  - 1.00000./(CONSTANTS(:,51)./CONSTANTS(:,52) - 1.00000) - CONSTANTS(:,51)./CONSTANTS(:,52) .^  - 1.00000./(1.00000 - CONSTANTS(:,52)./CONSTANTS(:,51));
    CONSTANTS(:,62) =  CONSTANTS(:,19).*CONSTANTS(:,21)./( CONSTANTS(:,19).*CONSTANTS(:,21)+ CONSTANTS(:,23).*CONSTANTS(:,21)+ CONSTANTS(:,23).*CONSTANTS(:,20)+ CONSTANTS(:,22).*CONSTANTS(:,19)+ CONSTANTS(:,22).*CONSTANTS(:,20)+ CONSTANTS(:,23).*CONSTANTS(:,19));
    CONSTANTS(:,63) =  CONSTANTS(:,37).*CONSTANTS(:,35).*CONSTANTS(:,62);
    CONSTANTS(:,64) =  abs(CONSTANTS(:,40) - CONSTANTS(:,39))./(CONSTANTS(:,40) - CONSTANTS(:,39)).*CONSTANTS(:,41).*(exp( CONSTANTS(:,42).*abs(CONSTANTS(:,40) - CONSTANTS(:,39))) - 1.00000);
    if (isempty(STATES)), warning('Initial values for states not set');, end
end

function [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS)
    global algebraicVariableCount;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        ALGEBRAIC = zeros(1, algebraicVariableCount);
    else
        statesRowCount = statesSize(1);
        ALGEBRAIC = zeros(statesRowCount, algebraicVariableCount);
        RATES = zeros(statesRowCount, statesColumnCount);
    end
    ALGEBRAIC(:,1) = piecewise({CONSTANTS(:,4)./2.00000<STATES(:,1)./2.00000, CONSTANTS(:,4)./2.00000 }, STATES(:,1)./2.00000);
    ALGEBRAIC(:,2) = piecewise({STATES(:,1)./2.00000 - (STATES(:,1) - CONSTANTS(:,3))>CONSTANTS(:,5)./2.00000, STATES(:,1)./2.00000 - (STATES(:,1) - CONSTANTS(:,3)) }, CONSTANTS(:,5)./2.00000);
    ALGEBRAIC(:,3) = ALGEBRAIC(:,1) - ALGEBRAIC(:,2);
    ALGEBRAIC(:,5) = ALGEBRAIC(:,3)./CONSTANTS(:,3);
    ALGEBRAIC(:,7) =  (1.00000 - ALGEBRAIC(:,5)).*STATES(:,2)+ ALGEBRAIC(:,5).*STATES(:,3);
    ALGEBRAIC(:,9) =  abs(1.00000./(1.00000+CONSTANTS(:,13)./ALGEBRAIC(:,7) .^ CONSTANTS(:,14))) .^ (1.0 / 2);
    ALGEBRAIC(:,11) =  CONSTANTS(:,15).*ALGEBRAIC(:,9).*CONSTANTS(:,8) .^ (CONSTANTS(:,18) - 37.0000)./10.0000;
    ALGEBRAIC(:,13) = piecewise({1.00000./ALGEBRAIC(:,9)<100.000, 1.00000./ALGEBRAIC(:,9) }, 100.000);
    ALGEBRAIC(:,15) =  CONSTANTS(:,16).*ALGEBRAIC(:,13).*CONSTANTS(:,9) .^ (CONSTANTS(:,18) - 37.0000)./10.0000;
    RATES(:,8) =  ALGEBRAIC(:,15).*STATES(:,9) -  ALGEBRAIC(:,11).*STATES(:,8);
    RATES(:,9) =  ALGEBRAIC(:,11).*STATES(:,8) -  ALGEBRAIC(:,15).*STATES(:,9);
    ALGEBRAIC(:,4) =  ALGEBRAIC(:,3).*2.00000./(CONSTANTS(:,4) - CONSTANTS(:,5));
    ALGEBRAIC(:,6) =  CONSTANTS(:,37).*ALGEBRAIC(:,4).*( STATES(:,4).*STATES(:,6)+ STATES(:,5).*STATES(:,7));
    ALGEBRAIC(:,8) =  1.00000.*ALGEBRAIC(:,6)./CONSTANTS(:,63);
    ALGEBRAIC(:,10) =  (STATES(:,1) - CONSTANTS(:,39))./abs(STATES(:,1) - CONSTANTS(:,39)).*CONSTANTS(:,41).*(exp( CONSTANTS(:,42).*abs(STATES(:,1) - CONSTANTS(:,39))) - 1.00000);
    ALGEBRAIC(:,12) = piecewise({STATES(:,1)>CONSTANTS(:,38),  CONSTANTS(:,43).*(exp( CONSTANTS(:,44).*abs(STATES(:,1) - CONSTANTS(:,38))) - 1.00000) }, 0.00000);
    ALGEBRAIC(:,14) = ALGEBRAIC(:,10)+ALGEBRAIC(:,12);
    ALGEBRAIC(:,16) = piecewise({CONSTANTS(:,48)==1.00000,  CONSTANTS(:,47).*(CONSTANTS(:,40) - STATES(:,1)) }, 0.00000);
    RATES(:,11) = CONSTANTS(:,64)+ALGEBRAIC(:,16) - (ALGEBRAIC(:,14)+ALGEBRAIC(:,8));
    ALGEBRAIC(:,19) = exp(  - STATES(:,5)./abs(STATES(:,5)).*CONSTANTS(:,25).*STATES(:,5)./CONSTANTS(:,35) .^ 2.00000);
    ALGEBRAIC(:,21) =  CONSTANTS(:,21).*ALGEBRAIC(:,19).*CONSTANTS(:,29).*CONSTANTS(:,32) .^ (CONSTANTS(:,18) - 37.0000)./10.0000;
    ALGEBRAIC(:,20) = exp( (STATES(:,4) - CONSTANTS(:,35))./abs(STATES(:,4) - CONSTANTS(:,35)).*CONSTANTS(:,26).*(STATES(:,4) - CONSTANTS(:,35))./CONSTANTS(:,35) .^ 2.00000);
    ALGEBRAIC(:,22) =  CONSTANTS(:,22).*ALGEBRAIC(:,20).*CONSTANTS(:,29).*CONSTANTS(:,33) .^ (CONSTANTS(:,18) - 37.0000)./10.0000;
    ALGEBRAIC(:,23) = piecewise({STATES(:,4)<CONSTANTS(:,35), exp( CONSTANTS(:,27).*(CONSTANTS(:,35) - STATES(:,4))./CONSTANTS(:,35) .^ 2.00000) }, exp( CONSTANTS(:,28).*(STATES(:,4) - CONSTANTS(:,35))./CONSTANTS(:,35) .^ 2.00000));
    ALGEBRAIC(:,24) =  CONSTANTS(:,23).*ALGEBRAIC(:,23).*CONSTANTS(:,29).*CONSTANTS(:,34) .^ (CONSTANTS(:,18) - 37.0000)./10.0000;
    ALGEBRAIC(:,25) =  ALGEBRAIC(:,21).*STATES(:,7) - ( ALGEBRAIC(:,22).*STATES(:,6)+ ALGEBRAIC(:,24).*STATES(:,6));
    RATES(:,6) = ALGEBRAIC(:,25);
    ALGEBRAIC(:,26) = 1.00000 - STATES(:,10) - STATES(:,7) - STATES(:,6);
    RATES(:,10) =  ALGEBRAIC(:,15).*ALGEBRAIC(:,26) -  ALGEBRAIC(:,11).*STATES(:,10);
    ALGEBRAIC(:,17) = 1.00000+ (1.00000 - ALGEBRAIC(:,4)).*CONSTANTS(:,24);
    ALGEBRAIC(:,18) =  CONSTANTS(:,20).*ALGEBRAIC(:,17).*CONSTANTS(:,29).*CONSTANTS(:,31) .^ (CONSTANTS(:,18) - 37.0000)./10.0000;
    ALGEBRAIC(:,27) =  CONSTANTS(:,59).*ALGEBRAIC(:,26)+ ALGEBRAIC(:,22).*STATES(:,6) - ( ALGEBRAIC(:,18).*STATES(:,7)+ ALGEBRAIC(:,21).*STATES(:,7));
    RATES(:,7) = ALGEBRAIC(:,27);
    ALGEBRAIC(:,30) = piecewise({STATES(:,1)<=CONSTANTS(:,1)&STATES(:,1)>CONSTANTS(:,2), (STATES(:,11)+ (CONSTANTS(:,40) - STATES(:,1)).*CONSTANTS(:,46))./CONSTANTS(:,45) }, 0.00000);
    RATES(:,1) = ALGEBRAIC(:,30);
    ALGEBRAIC(:,28) = ( ALGEBRAIC(:,22).*CONSTANTS(:,59)+ ALGEBRAIC(:,24).*CONSTANTS(:,59))./( CONSTANTS(:,59).*ALGEBRAIC(:,21)+ ALGEBRAIC(:,24).*ALGEBRAIC(:,21)+ ALGEBRAIC(:,24).*ALGEBRAIC(:,18)+ ALGEBRAIC(:,22).*CONSTANTS(:,59)+ ALGEBRAIC(:,22).*ALGEBRAIC(:,18)+ ALGEBRAIC(:,24).*CONSTANTS(:,59));
    ALGEBRAIC(:,31) = ALGEBRAIC(:,30)./2.00000+ CONSTANTS(:,36)./ALGEBRAIC(:,28).*( CONSTANTS(:,59).* - STATES(:,5)+ ALGEBRAIC(:,22).*(STATES(:,4) - CONSTANTS(:,35)+STATES(:,5)));
    RATES(:,5) = ALGEBRAIC(:,31);
    ALGEBRAIC(:,29) =  CONSTANTS(:,59).*ALGEBRAIC(:,21)./( CONSTANTS(:,59).*ALGEBRAIC(:,21)+ ALGEBRAIC(:,24).*ALGEBRAIC(:,21)+ ALGEBRAIC(:,24).*ALGEBRAIC(:,18)+ ALGEBRAIC(:,22).*CONSTANTS(:,59)+ ALGEBRAIC(:,22).*ALGEBRAIC(:,18)+ ALGEBRAIC(:,24).*CONSTANTS(:,59));
    ALGEBRAIC(:,32) = ALGEBRAIC(:,30)./2.00000+ CONSTANTS(:,36)./ALGEBRAIC(:,29).*ALGEBRAIC(:,21).*(STATES(:,5)+CONSTANTS(:,35) - STATES(:,4));
    RATES(:,4) = ALGEBRAIC(:,32);
    ALGEBRAIC(:,41) = piecewise({VOI>CONSTANTS(:,53),  (CONSTANTS(:,54) - CONSTANTS(:,55))./CONSTANTS(:,61).*(exp( - (VOI - CONSTANTS(:,53))./CONSTANTS(:,51)) - exp( - (VOI - CONSTANTS(:,53))./CONSTANTS(:,52)))+CONSTANTS(:,55) }, CONSTANTS(:,55));
    ALGEBRAIC(:,43) =  CONSTANTS(:,56).*ALGEBRAIC(:,41).*(1.00000 - STATES(:,2)) -  CONSTANTS(:,57).*STATES(:,2);
    RATES(:,2) = ALGEBRAIC(:,43);
    ALGEBRAIC(:,44) =  CONSTANTS(:,56).*ALGEBRAIC(:,41).*(1.00000 - STATES(:,3)) -  CONSTANTS(:,58).*STATES(:,3);
    RATES(:,3) = ALGEBRAIC(:,44);
   RATES = RATES';
end

% Calculate algebraic variables
function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI)
    ALGEBRAIC(:,1) = piecewise({CONSTANTS(:,4)./2.00000<STATES(:,1)./2.00000, CONSTANTS(:,4)./2.00000 }, STATES(:,1)./2.00000);
    ALGEBRAIC(:,2) = piecewise({STATES(:,1)./2.00000 - (STATES(:,1) - CONSTANTS(:,3))>CONSTANTS(:,5)./2.00000, STATES(:,1)./2.00000 - (STATES(:,1) - CONSTANTS(:,3)) }, CONSTANTS(:,5)./2.00000);
    ALGEBRAIC(:,3) = ALGEBRAIC(:,1) - ALGEBRAIC(:,2);
    ALGEBRAIC(:,5) = ALGEBRAIC(:,3)./CONSTANTS(:,3);
    ALGEBRAIC(:,7) =  (1.00000 - ALGEBRAIC(:,5)).*STATES(:,2)+ ALGEBRAIC(:,5).*STATES(:,3);
    ALGEBRAIC(:,9) =  abs(1.00000./(1.00000+CONSTANTS(:,13)./ALGEBRAIC(:,7) .^ CONSTANTS(:,14))) .^ (1.0 / 2);
    ALGEBRAIC(:,11) =  CONSTANTS(:,15).*ALGEBRAIC(:,9).*CONSTANTS(:,8) .^ (CONSTANTS(:,18) - 37.0000)./10.0000;
    ALGEBRAIC(:,13) = piecewise({1.00000./ALGEBRAIC(:,9)<100.000, 1.00000./ALGEBRAIC(:,9) }, 100.000);
    ALGEBRAIC(:,15) =  CONSTANTS(:,16).*ALGEBRAIC(:,13).*CONSTANTS(:,9) .^ (CONSTANTS(:,18) - 37.0000)./10.0000;
    ALGEBRAIC(:,4) =  ALGEBRAIC(:,3).*2.00000./(CONSTANTS(:,4) - CONSTANTS(:,5));
    ALGEBRAIC(:,6) =  CONSTANTS(:,37).*ALGEBRAIC(:,4).*( STATES(:,4).*STATES(:,6)+ STATES(:,5).*STATES(:,7));
    ALGEBRAIC(:,8) =  1.00000.*ALGEBRAIC(:,6)./CONSTANTS(:,63);
    ALGEBRAIC(:,10) =  (STATES(:,1) - CONSTANTS(:,39))./abs(STATES(:,1) - CONSTANTS(:,39)).*CONSTANTS(:,41).*(exp( CONSTANTS(:,42).*abs(STATES(:,1) - CONSTANTS(:,39))) - 1.00000);
    ALGEBRAIC(:,12) = piecewise({STATES(:,1)>CONSTANTS(:,38),  CONSTANTS(:,43).*(exp( CONSTANTS(:,44).*abs(STATES(:,1) - CONSTANTS(:,38))) - 1.00000) }, 0.00000);
    ALGEBRAIC(:,14) = ALGEBRAIC(:,10)+ALGEBRAIC(:,12);
    ALGEBRAIC(:,16) = piecewise({CONSTANTS(:,48)==1.00000,  CONSTANTS(:,47).*(CONSTANTS(:,40) - STATES(:,1)) }, 0.00000);
    ALGEBRAIC(:,19) = exp(  - STATES(:,5)./abs(STATES(:,5)).*CONSTANTS(:,25).*STATES(:,5)./CONSTANTS(:,35) .^ 2.00000);
    ALGEBRAIC(:,21) =  CONSTANTS(:,21).*ALGEBRAIC(:,19).*CONSTANTS(:,29).*CONSTANTS(:,32) .^ (CONSTANTS(:,18) - 37.0000)./10.0000;
    ALGEBRAIC(:,20) = exp( (STATES(:,4) - CONSTANTS(:,35))./abs(STATES(:,4) - CONSTANTS(:,35)).*CONSTANTS(:,26).*(STATES(:,4) - CONSTANTS(:,35))./CONSTANTS(:,35) .^ 2.00000);
    ALGEBRAIC(:,22) =  CONSTANTS(:,22).*ALGEBRAIC(:,20).*CONSTANTS(:,29).*CONSTANTS(:,33) .^ (CONSTANTS(:,18) - 37.0000)./10.0000;
    ALGEBRAIC(:,23) = piecewise({STATES(:,4)<CONSTANTS(:,35), exp( CONSTANTS(:,27).*(CONSTANTS(:,35) - STATES(:,4))./CONSTANTS(:,35) .^ 2.00000) }, exp( CONSTANTS(:,28).*(STATES(:,4) - CONSTANTS(:,35))./CONSTANTS(:,35) .^ 2.00000));
    ALGEBRAIC(:,24) =  CONSTANTS(:,23).*ALGEBRAIC(:,23).*CONSTANTS(:,29).*CONSTANTS(:,34) .^ (CONSTANTS(:,18) - 37.0000)./10.0000;
    ALGEBRAIC(:,25) =  ALGEBRAIC(:,21).*STATES(:,7) - ( ALGEBRAIC(:,22).*STATES(:,6)+ ALGEBRAIC(:,24).*STATES(:,6));
    ALGEBRAIC(:,26) = 1.00000 - STATES(:,10) - STATES(:,7) - STATES(:,6);
    ALGEBRAIC(:,17) = 1.00000+ (1.00000 - ALGEBRAIC(:,4)).*CONSTANTS(:,24);
    ALGEBRAIC(:,18) =  CONSTANTS(:,20).*ALGEBRAIC(:,17).*CONSTANTS(:,29).*CONSTANTS(:,31) .^ (CONSTANTS(:,18) - 37.0000)./10.0000;
    ALGEBRAIC(:,27) =  CONSTANTS(:,59).*ALGEBRAIC(:,26)+ ALGEBRAIC(:,22).*STATES(:,6) - ( ALGEBRAIC(:,18).*STATES(:,7)+ ALGEBRAIC(:,21).*STATES(:,7));
    ALGEBRAIC(:,30) = piecewise({STATES(:,1)<=CONSTANTS(:,1)&STATES(:,1)>CONSTANTS(:,2), (STATES(:,11)+ (CONSTANTS(:,40) - STATES(:,1)).*CONSTANTS(:,46))./CONSTANTS(:,45) }, 0.00000);
    ALGEBRAIC(:,28) = ( ALGEBRAIC(:,22).*CONSTANTS(:,59)+ ALGEBRAIC(:,24).*CONSTANTS(:,59))./( CONSTANTS(:,59).*ALGEBRAIC(:,21)+ ALGEBRAIC(:,24).*ALGEBRAIC(:,21)+ ALGEBRAIC(:,24).*ALGEBRAIC(:,18)+ ALGEBRAIC(:,22).*CONSTANTS(:,59)+ ALGEBRAIC(:,22).*ALGEBRAIC(:,18)+ ALGEBRAIC(:,24).*CONSTANTS(:,59));
    ALGEBRAIC(:,31) = ALGEBRAIC(:,30)./2.00000+ CONSTANTS(:,36)./ALGEBRAIC(:,28).*( CONSTANTS(:,59).* - STATES(:,5)+ ALGEBRAIC(:,22).*(STATES(:,4) - CONSTANTS(:,35)+STATES(:,5)));
    ALGEBRAIC(:,29) =  CONSTANTS(:,59).*ALGEBRAIC(:,21)./( CONSTANTS(:,59).*ALGEBRAIC(:,21)+ ALGEBRAIC(:,24).*ALGEBRAIC(:,21)+ ALGEBRAIC(:,24).*ALGEBRAIC(:,18)+ ALGEBRAIC(:,22).*CONSTANTS(:,59)+ ALGEBRAIC(:,22).*ALGEBRAIC(:,18)+ ALGEBRAIC(:,24).*CONSTANTS(:,59));
    ALGEBRAIC(:,32) = ALGEBRAIC(:,30)./2.00000+ CONSTANTS(:,36)./ALGEBRAIC(:,29).*ALGEBRAIC(:,21).*(STATES(:,5)+CONSTANTS(:,35) - STATES(:,4));
    ALGEBRAIC(:,41) = piecewise({VOI>CONSTANTS(:,53),  (CONSTANTS(:,54) - CONSTANTS(:,55))./CONSTANTS(:,61).*(exp( - (VOI - CONSTANTS(:,53))./CONSTANTS(:,51)) - exp( - (VOI - CONSTANTS(:,53))./CONSTANTS(:,52)))+CONSTANTS(:,55) }, CONSTANTS(:,55));
    ALGEBRAIC(:,43) =  CONSTANTS(:,56).*ALGEBRAIC(:,41).*(1.00000 - STATES(:,2)) -  CONSTANTS(:,57).*STATES(:,2);
    ALGEBRAIC(:,44) =  CONSTANTS(:,56).*ALGEBRAIC(:,41).*(1.00000 - STATES(:,3)) -  CONSTANTS(:,58).*STATES(:,3);
    ALGEBRAIC(:,33) = (STATES(:,6)+STATES(:,7))./(CONSTANTS(:,62)+CONSTANTS(:,60));
    ALGEBRAIC(:,34) = (ALGEBRAIC(:,25)+ALGEBRAIC(:,27))./(CONSTANTS(:,62)+CONSTANTS(:,60));
    ALGEBRAIC(:,35) =  CONSTANTS(:,50).*( (1.00000 - ALGEBRAIC(:,5)).*STATES(:,2)+ ALGEBRAIC(:,5).*( ALGEBRAIC(:,33).*STATES(:,3)+ (1.00000 - ALGEBRAIC(:,33)).*STATES(:,2)));
    ALGEBRAIC(:,36) = piecewise({STATES(:,1)<CONSTANTS(:,4),   - 0.500000.*ALGEBRAIC(:,30) }, 0.00000);
    ALGEBRAIC(:,37) = piecewise({ 2.00000.*CONSTANTS(:,3) - STATES(:,1)>CONSTANTS(:,5),   - 0.500000.*ALGEBRAIC(:,30) }, 0.00000);
    ALGEBRAIC(:,38) = ALGEBRAIC(:,36) - ALGEBRAIC(:,37);
    ALGEBRAIC(:,39) = ALGEBRAIC(:,38)./CONSTANTS(:,3);
    ALGEBRAIC(:,40) =  2.00000.*ALGEBRAIC(:,38)./(CONSTANTS(:,4) - CONSTANTS(:,5));
    ALGEBRAIC(:,42) =  CONSTANTS(:,49).*ALGEBRAIC(:,40).*( STATES(:,4).*STATES(:,6)+ STATES(:,5).*STATES(:,7))+ CONSTANTS(:,49).*ALGEBRAIC(:,4).*( ALGEBRAIC(:,32).*STATES(:,6)+ STATES(:,4).*ALGEBRAIC(:,25)+ ALGEBRAIC(:,31).*STATES(:,7)+ STATES(:,5).*ALGEBRAIC(:,27));
    ALGEBRAIC(:,45) =  CONSTANTS(:,50).*(  - ALGEBRAIC(:,39).*STATES(:,2)+ (1.00000 - ALGEBRAIC(:,5)).*ALGEBRAIC(:,43)+ ALGEBRAIC(:,39).*( ALGEBRAIC(:,33).*STATES(:,3)+ (1.00000 - ALGEBRAIC(:,33)).*STATES(:,2))+ ALGEBRAIC(:,5).*( ALGEBRAIC(:,34).*STATES(:,3)+ ALGEBRAIC(:,33).*ALGEBRAIC(:,44)+ (1.00000 - ALGEBRAIC(:,33)).*ALGEBRAIC(:,43) -  ALGEBRAIC(:,34).*STATES(:,2)));
end

% Compute result of a piecewise function
function x = piecewise(cases, default)
    set = [0];
    for i = 1:2:length(cases)
        if (length(cases{i+1}) == 1)
            x(cases{i} & ~set,:) = cases{i+1};
        else
            x(cases{i} & ~set,:) = cases{i+1}(cases{i} & ~set);
        end
        set = set | cases{i};
        if(set), break, end
    end
    if (length(default) == 1)
        x(~set,:) = default;
    else
        x(~set,:) = default(~set);
    end
end

% Pad out or shorten strings to a set length
function strout = strpad(strin)
    req_length = 160;
    insize = size(strin,2);
    if insize > req_length
        strout = strin(1:req_length);
    else
        strout = [strin, blanks(req_length - insize)];
    end
end


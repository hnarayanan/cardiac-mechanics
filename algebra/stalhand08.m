clear all;
clc;

# Define system level parameters
global q_1 = 1.0;
global q_2 = 1.0;
global f_1 = 1.3;
global f_2 = 85.5;
global E_1 = 12.5;
global E_2 = 8.8;
global lambda_opt = 1.4;
global v = 5;
global xi = 1/sqrt(2);

function alpha_dot = hai_murphy_98(alpha, t)

  global lambda;
  global beta;

  C_0 = -13.723;
  C_1 = 19.807;
  beta_0 = C_0 + C_1/lambda;

  # Rate constants
  k_1 = beta**4/(beta**4 + beta_0**4);
  k_2 = 0.5;
  k_3 = 0.4;
  k_4 = 0.1;
  k_5 = 0.5;
  k_6 = k_1;
  k_7 = 0.01;

  # Put rate constants into column vectors
  gamma_1 = [k_1; -k_1; 0; 0];
  gamma_2 = [-k_2; k_2 + k_3; -k_3; 0];
  gamma_3 = [0; -k_4; k_4 + k_5; -k_5];
  gamma_4 = [-k_7; 0; -k_6; -k_6 + k_7];

  # Define the eta functions
  lambda_bar = 1.12;
  A_1 =  2.069;
  A_2 = -3.289;
  A_3 = -1.744;
  A_4 =  0.989;
  eta_1 = A_1*(lambda - lambda_bar) + 1;
  eta_2 = A_2*(lambda - lambda_bar) + 1;
  eta_3 = A_3*(lambda - lambda_bar) + 1;
  eta_4 = A_4*(lambda - lambda_bar) + 1;

  alpha_dot = -[eta_1*gamma_1 eta_2*gamma_2 eta_3*gamma_3 eta_4*gamma_4]*alpha;

endfunction

function lambda_a_dot = stalhand_08(lambda_a, t)

  global alpha_steady;
  global lambda_opt;
  global f_1;
  global f_2;
  global E_1;
  global E_2;
  global v;
  global lambda;
  global xi;
  lambda_c = lambda/lambda_a;


  lambda_a_dot = (-alpha_steady(3)*f_1*v*exp(-(lambda_a - lambda_opt)**2/(2*xi**2)) + (2*lambda_a - 2*lambda_opt)*(lambda_c - 1)**2*(E_1*alpha_steady(3)/2 + E_2*alpha_steady(4)/2)*exp(-(lambda_a - lambda_opt)**2/(2*xi**2))/(2*xi**2) + lambda_c*(2*lambda_c - 2)*(E_1*alpha_steady(3)/2 + E_2*alpha_steady(4)/2)*exp(-(lambda_a - lambda_opt)**2/(2*xi**2))/lambda_a)*exp((lambda_a - lambda_opt)**2/(2*xi**2))/(alpha_steady(3)*f_1 + alpha_steady(4)*f_2);

endfunction

function stress = total_stress(alpha, lambda, lambda_a)

  global q_1;
  global q_2;
  global E_1;
  global E_2;
  global lambda_opt;
  global xi;

  lambda_c = lambda/lambda_a;

  # Generated from stalhand08.py
  stress = q_1*(q_1*exp(q_1*(lambda - 1)) - q_2)/q_2 + \
      (2*lambda_c - 2)*(E_1*alpha(3)/2 + E_2*alpha(4)/2) * \
      exp(-(lambda_a - lambda_opt)**2/(2*xi**2))/lambda_a;

endfunction


# Main driver

# Step 1: Solve for chemical state, given lambda and beta
global lambda = 0.75;
global beta = 0.1;

alpha_0 = [0.25; 0.25; 0.25; 0.25];
steps = 1000;
t = linspace(0, 100, steps)';

alpha = lsode("hai_murphy_98", alpha_0, t);
global alpha_steady = alpha(steps, :);

# Step 2: Solve for the active stretch, given the steady-state
#         chemical concentrations above

lambda_a_0 = 1.0;
t = linspace(0, 100, steps)';

lambda_a = lsode("stalhand_08", lambda_a_0, t);

lambda_a_steady = lambda_a(steps);

P = total_stress(alpha_steady, lambda, lambda_a_steady)


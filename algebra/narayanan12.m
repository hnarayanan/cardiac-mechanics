1;

clear all;
clc;

# Define system level parameters
global q_1 = 10.0;
global q_2 = 10.0;
global f_1 = 1.3;
global f_2 = 85.5;
global E_1 = 12.5;
global E_2 = 8.8;
global lambda_opt = 1.4;
global v = 5;
global xi = 1/sqrt(2);

global lambda = 1.5;
global beta = 1.1;

function lambda_a_dot = stalhand_08(lambda_a, t)

  global alphas;
  alpha_steady = alphas(1, :);
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
global alphas = [0.9987702708891195 0.00000589174448707144 0.00017297197566117962 0.0010508653907322332;
0.998762573843706 0.000005929923993356007 0.0001740544110342275 0.00105744182126638;
0.9986643685097072 0.0000071128544475041605 0.00018787560249589286 0.001140643033349413;
0.8208941057402434 0.0475614686013137 0.02486416276029727 0.10668026289814568;
0.031302907149558856 0.3117142373634343 0.12501012182413113 0.5319727336628757;
0.031178422520788545 0.31176826838826843 0.12502404223129404 0.532029266859649];

# Step 2: Solve for the active stretch, given the steady-state
#         chemical concentrations above

steps = 1000;
lambda_a_0 = 1.0;
t = linspace(0, 100, steps)';

lambda_a = lsode("stalhand_08", lambda_a_0, t);
lambda_a_steady = lambda_a(steps)
plot(t, lambda_a)

# Step 3: Compute the total stress, given the steady state
#         chemical state and active stretch

P = total_stress(alphas(1, :), lambda, lambda_a_steady)


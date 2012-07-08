1;

function xdot = hai_murphy_98(x, t)

  lambda = 0.75;
  beta = 1.0;

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

  xdot = -[eta_1*gamma_1 eta_2*gamma_2 eta_3*gamma_3 eta_4*gamma_4]*x;

endfunction

x0 = [0.25; 0.25; 0.25; 0.25];
t = linspace(0, 100, 1000)';

x = lsode("hai_murphy_98", x0, t);

plot(t, x)

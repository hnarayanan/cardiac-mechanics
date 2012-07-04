# Size of variable arrays:
sizeAlgebraic = 4
sizeStates = 4
sizeConstants = 5
from math import *
from numpy import *

lmbda = 1.50
beta = 50.0
C_0 = -13.723
C_1 = 19.807
beta_0 = C_0 + C_1/lmbda

k_1 = beta**4/(beta**4 + beta_0**4)

def createLegends():
    legend_states = [""] * sizeStates
    legend_rates = [""] * sizeStates
    legend_algebraic = [""] * sizeAlgebraic
    legend_voi = ""
    legend_constants = [""] * sizeConstants
    legend_voi = "time in component environment (second)"
    legend_states[0] = "M in component M (dimensionless)"
    legend_states[1] = "AM in component AM (dimensionless)"
    legend_states[2] = "Mp in component Mp (dimensionless)"
    legend_algebraic[0] = "k1 in component model_parameters (first_order_rate_constant)"
    legend_constants[0] = "k2 in component model_parameters (first_order_rate_constant)"
    legend_constants[1] = "k7 in component model_parameters (first_order_rate_constant)"
    legend_states[3] = "AMp in component AMp (dimensionless)"
    legend_constants[2] = "k3 in component model_parameters (first_order_rate_constant)"
    legend_constants[3] = "k4 in component model_parameters (first_order_rate_constant)"
    legend_constants[4] = "k5 in component model_parameters (first_order_rate_constant)"
    legend_algebraic[3] = "k6 in component model_parameters (first_order_rate_constant)"
    legend_algebraic[1] = "phosphorylation in component phosphorylation (dimensionless)"
    legend_algebraic[2] = "stress in component stress (dimensionless)"
    legend_rates[0] = "d/dt M in component M (dimensionless)"
    legend_rates[2] = "d/dt Mp in component Mp (dimensionless)"
    legend_rates[3] = "d/dt AMp in component AMp (dimensionless)"
    legend_rates[1] = "d/dt AM in component AM (dimensionless)"
    return (legend_states, legend_algebraic, legend_voi, legend_constants)

def initConsts():
    constants = [0.0] * sizeConstants; states = [0.0] * sizeStates;
    states[0] = 1.0
    states[1] = 0.0
    states[2] = 0.0
    constants[0] = 0.5
    constants[1] = 0.01
    states[3] = 0.0
    constants[2] = 0.4
    constants[3] = 0.1
    constants[4] = 0.5
    return (states, constants)

def computeRates(voi, states, constants):
    rates = [0.0] * sizeStates; algebraic = [0.0] * sizeAlgebraic
    # algebraic[0] = custom_piecewise([greater_equal(voi , 0.00000) & less(voi , 5.00000), 0.550000 , True, 0.300000])
    algebraic[0] = k_1
    rates[0] = -(algebraic[0]*states[0])+constants[0]*states[2]+constants[1]*states[1]
    rates[2] = constants[3]*states[3]+algebraic[0]*states[0]-(constants[0]+constants[2])*states[2]
    algebraic[3] = algebraic[0]
    rates[3] = constants[2]*states[2]+algebraic[3]*states[1]-(constants[4]+constants[3])*states[3]
    rates[1] = constants[4]*states[3]-(algebraic[3]+constants[1])*states[1]
    return(rates)

def computeAlgebraic(constants, states, voi):
    algebraic = array([[0.0] * len(voi)] * sizeAlgebraic)
    states = array(states)
    voi = array(voi)
    algebraic[0] = custom_piecewise([greater_equal(voi , 0.00000) & less(voi , 5.00000), 0.550000 , True, 0.300000])
    algebraic[3] = algebraic[0]
    algebraic[1] = states[3]+states[2]
    algebraic[2] = states[3]+states[1]
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
    voi = linspace(0, 100, 5000)

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
    pylab.plot(voi,vstack((states,algebraic)).T[:, 0:4])
    pylab.xlabel(legend_voi)
    pylab.legend(legend_states + legend_algebraic, loc='best')
    pylab.show()

if __name__ == "__main__":
    (voi, states, algebraic) = solve_model()
    plot_model(voi, states, algebraic)

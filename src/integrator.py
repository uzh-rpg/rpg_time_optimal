from casadi import MX, Function


def Euler(dynamics):
  dt = MX.sym('dt', 1)
  x = MX.sym('x', dynamics.size1_in(0))
  u = MX.sym('u', dynamics.size1_in(1))
  integrator = Function('integrator',
    [x, u, dt],
    [x + dt * dynamics(x, u)],
    ['x', 'u', 'dt'], ['xn'])
  return integrator

def RungeKutta4(dynamics):
  dt = MX.sym('dt', 1)
  x = MX.sym('x', dynamics.size1_in(0))
  u = MX.sym('u', dynamics.size1_in(1))
  k1 = dynamics(x, u)
  k2 = dynamics(x + dt/2 * k1, u)
  k3 = dynamics(x + dt/2 * k2, u)
  k4 = dynamics(x + dt * k3, u)
  integrator = Function('integrator',
    [x, u, dt],
    [x + dt/6 * (k1 + 2*k2 + 2*k3 + k4)],
    ['x', 'u', 'dt'], ['xn'])
  return integrator
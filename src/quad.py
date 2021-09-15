from casadi import MX, DM, vertcat, mtimes, Function, inv, cross, sqrt, norm_2
import yaml
from quaternion import *


class Quad:
  def __init__(self, filename = ""):
    self.m = 1                                        # mass in [kg]
    self.l = 1                                        # arm length
    self.I = DM([(1, 0, 0), (0, 1, 0), (0, 0, 1)])    # Inertia
    self.I_inv = inv(self.I)                          # Inertia inverse
    self.T_max = 5                                    # max thrust [N]
    self.T_min = 0                                    # min thrust [N]
    self.omega_max = 3                                # max bodyrate [rad/s]
    self.ctau = 0.5                                   # thrust torque coeff.
    self.rampup_dist = 0
    self.T_ramp_start = 5
    self.omega_ramp_start = 3

    self.v_max = None
    self.cd = 0.0

    self.g = 9.801

    if filename:
      self.load(filename)

  def load(self, filename):
    print("Loading track from " + filename)
    with open(filename, 'r') as file:
      quad = yaml.load(file, Loader=yaml.FullLoader)

    if 'mass' in quad:
      self.m = quad['mass']
    else:
      print("No mass specified in " + filename)

    if 'arm_length' in quad:
      self.l = quad['arm_length']
    else:
      print("No arm length specified in " + filename)

    if 'inertia' in quad:
      self.I = DM(quad['inertia'])
      self.I_inv = inv(self.I)
    else:
      print("No inertia specified in " + filename)


    if 'TWR_max' in quad:
      self.T_max = quad['TWR_max'] * 9.81 * self.m / 4
    elif 'thrust_max' in quad:
      self.T_max = quad['thrust_max']
    else:
      print("No max thrust specified in " + filename)

    if 'TWR_min' in quad:
      self.T_min = quad['TWR_min'] * 9.81 * self.m / 4
    elif 'thrust_min' in quad:
      self.T_min = quad['thrust_min']
    else:
      print("No min thrust specified in " + filename)

    if 'omega_max_xy' in quad:
      self.omega_max_xy = quad['omega_max_xy']
    else:
      print("No max omega_xy specified in " + filename)

    if 'omega_max_z' in quad:
      self.omega_max_z = quad['omega_max_z']
    else:
      print("No max omega_z specified in " + filename)

    if 'torque_coeff' in quad:
      self.ctau = quad['torque_coeff']
    else:
      print("No thrust to drag coefficient specified in " + filename)

    if 'v_max' in quad:
      self.v_max = quad['v_max']
      a_max = 4 * self.T_max / self.m
      a_hmax = sqrt(a_max**2 - self.g**2)
      self.cd = a_hmax / self.v_max
    if 'drag_coeff' in quad:
      self.cd = quad['drag_coeff']

    if 'rampup_dist' in quad:
      self.rampup_dist = quad['rampup_dist']
      if 'TWR_ramp_start' in quad and 'omega_ramp_start' in quad:
        self.T_ramp_start = min(quad['TWR_ramp_start'] * 9.81 * self.m / 4, self.T_max)
        self.omega_ramp_start = min(quad['omega_ramp_start'], self.omega_max_xy)
      else:
        print("No TWR_ramp_start or omega_ramp_start specified. Disabling rampup")
        rampup_dist = 0


  def dynamics(self):
    p = MX.sym('p', 3)
    v = MX.sym('v', 3)
    q = MX.sym('q', 4)
    w = MX.sym('w', 3)
    T = MX.sym('thrust', 4)

    x = vertcat(p, v, q, w)
    u = vertcat(T)

    g = DM([0, 0, -self.g])

    x_dot = vertcat(
      v,
      rotate_quat(q, vertcat(0, 0, (T[0]+T[1]+T[2]+T[3])/self.m)) + g - v * self.cd,
      0.5*quat_mult(q, vertcat(0, w)),
      mtimes(self.I_inv, vertcat(
        self.l*(T[0]-T[1]-T[2]+T[3]),
        self.l*(-T[0]-T[1]+T[2]+T[3]),
        self.ctau*(T[0]-T[1]+T[2]-T[3]))
      -cross(w,mtimes(self.I,w)))
    )
    fx = Function('f',  [x, u], [x_dot], ['x', 'u'], ['x_dot'])
    return fx

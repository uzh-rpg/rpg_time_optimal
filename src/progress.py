import casadi as ca

class Progress4:
  def __init__(self, waypoints, options = {}):
    self.wp = waypoints
    self.NW = waypoints.shape[1]
    if 'distance_threshold' in options:
      self.dist = options['distance_threshold']
    else:
      self.dist = 0.3
    self.th = 5/3*self.dist**4

  def step(self):
    p_gate = ca.MX.sym('p_gate', 3)
    p_adj = ca.MX.sym('p_adj', 3)
    p = ca.MX.sym('p', 3)
    v = ca.MX.sym('v', 3)
    q = ca.MX.sym('q', 4)
    w = ca.MX.sym('w', 3)
    x = ca.vertcat(p, v, q, w)

    f_prog = ca.Function('f_prog',
      [x, p_gate, p_adj],
      [1-self.th/(ca.dot((p-p_gate-p_adj), (p-p_gate-p_adj))**2+self.th)])

    mu = ca.MX.sym('mu', self.NW)
    tau = ca.MX.sym('tau', 3, self.NW)
    mu_step = ca.Function('mu_step',
      [x, mu, tau],
      [ca.mtimes(ca.diag(f_prog(ca.repmat(x,1,self.NW),self.wp[:,0:self.NW], tau))
      ,mu)], ['x', 'mu', 'tau'], ['mun'])
    return mu_step

class Progress2:
  def __init__(self, waypoints, options = {}):
    self.wp = waypoints
    self.NW = waypoints.shape[1]
    if 'distance_threshold' in options:
      self.dist = options['distance_threshold']
    else:
      self.dist = 0.05
    self.th = 3*self.dist**2

  def step(self):
    p_gate = ca.MX.sym('p_gate', 3)
    p_adj = ca.MX.sym('p_adj', 3)
    p = ca.MX.sym('p', 3)
    v = ca.MX.sym('v', 3)
    q = ca.MX.sym('q', 4)
    w = ca.MX.sym('w', 3)
    x = ca.vertcat(p, v, q, w)

    f_prog = ca.Function('f_prog',
      [x, p_gate, p_adj],
      [1-self.th/(ca.dot((p-p_gate-p_adj), (p-p_gate-p_adj))+self.th)])

    mu = ca.MX.sym('mu', self.NW)
    tau = ca.MX.sym('tau', 3, self.NW)
    mu_step = ca.Function('mu_step',
      [x, mu, tau],
      [ca.mtimes(ca.diag(f_prog(ca.repmat(x,1,self.NW),self.wp[:,0:self.NW], tau))
      ,mu)], ['x', 'mu', 'tau'], ['mun'])
    return mu_step
  
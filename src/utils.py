import casadi as ca
from trajectory import Trajectory

def progressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    if iteration == total: print()

class CallbackSaver(ca.Callback):
  def __init__(self, folderpath, opts={}):
    ca.Callback.__init__(self)

    self.nx = None
    self.ng = None
    self.NPW = None
    self.opts = opts

    assert type(folderpath)==str
    self.folderpath = folderpath
    self.i = 0
  
  def set_size(self, nx, ng, NPW):
    self.nx = nx
    self.ng = ng
    self.NPW = NPW
    self.construct('CallbackSaver', self.opts)
  
  def set_wp(self, wp):
    self.wp = wp

  def get_n_in(self): return ca.nlpsol_n_out()
  def get_n_out(self): return 1
  def get_name_in(self, i): return ca.nlpsol_out(i)
  def get_name_out(self, i): return "ret"

  def get_sparsity_in(self, i):
    n = ca.nlpsol_out(i)
    if n=='f':
      return ca.Sparsity.scalar()
    elif n in ('x', 'lam_x'):
      return ca.Sparsity.dense(self.nx)
    elif n in ('g', 'lam_g'):
      return ca.Sparsity.dense(self.ng)
    else:
      return ca.Sparsity(0,0)
      
  def eval(self, arg):
    # Create dictionary
    darg = {}
    for (i,s) in enumerate(ca.nlpsol_out()): darg[s] = arg[i]

    X_opt = darg['x'].full().flatten()
    traj = Trajectory(X_opt, NPW=self.NPW, wp=self.wp)

    traj.save(self.folderpath+'/iteration_%05d.csv' % self.i)

    self.i += 1
    return [0]
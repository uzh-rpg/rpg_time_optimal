#!/usr/bin/python3
import matplotlib.pyplot as plt
import argparse
import math
import casadi as ca
from trajectory import Trajectory

class CallbackPlot(ca.Callback):
  def __init__(self,
               pos='xy', vel=None, ori=None, rate=None,
               inputs=None, prog=None, save=None, fig=None, opts={}):
    ca.Callback.__init__(self)

    self.nx = None
    self.ng = None
    self.NPW = None
    self.wp = None
    self.i = 0
    self.opts = opts

    if pos is not None: assert type(pos) == str
    if vel is not None: assert type(vel) == str
    if ori is not None: assert type(ori) == str
    if rate is not None: assert type(rate) == str
    if inputs is not None: assert type(inputs) == str
    if prog is not None: assert type(prog) == str
    if save is not None: assert type(save) == str

    self.pos = pos
    self.vel = vel
    self.ori = ori
    self.rate = rate
    self.inputs = inputs
    self.prog = prog
    self.save = save
    self.n_plots = 0

    for plot in [self.pos, self.vel, self.ori, self.rate, self.inputs, self.prog]:
      self.n_plots += int((plot is not None) & (type(plot)==str))
    assert self.n_plots > 0

    if self.n_plots<=3: 
      self.plot_r = self.n_plots
      self.plot_c = 1
    else:
      self.plot_r = 2
      self.plot_c = int((self.n_plots+1)/2)

    if fig is None:
      self.fig = plt.figure()
    elif type(fig) == int:
      self.fig = plt.figure(fig)
    elif type(fig) == plt.Figure:
      self.fig = plt.figure(fig.number)
    else:
      raise Exception('No valid figure passed!')

    print(self.plot_r)
    print(self.plot_c)
    self.axes = []
    if self.n_plots>1:
      for i in range(self.n_plots):
        self.axes.append(plt.subplot(self.plot_r, self.plot_c, 1+i))
    else:
      self.axes.append(self.fig.gca())
    
    self.fig.show()
  
  def set_size(self, nx, ng, NPW):
    self.nx = nx
    self.ng = ng
    self.NPW = NPW
    self.construct('CallbackPlot', self.opts)
  
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
    
    i_plot = 0
    if self.pos is not None:
      self.axes[i_plot].cla()
      traj.plotPos(self.axes[i_plot], plot_axis=self.pos)
      i_plot += 1

    if self.ori is not None:
      self.axes[i_plot].cla()
      traj.plotOrientation(self.axes[i_plot], plot_axis=self.ori)
      i_plot += 1

    if self.prog is not None:
      self.axes[i_plot].cla()
      traj.plotProgress(self.axes[i_plot], plot_axis=self.prog)
      i_plot += 1

    if self.vel is not None:
      self.axes[i_plot].cla()
      traj.plotVel(self.axes[i_plot], plot_axis=self.vel)
      i_plot += 1

    if self.rate is not None:
      self.axes[i_plot].cla()
      traj.plotOmega(self.axes[i_plot], plot_axis=self.rate)
      i_plot += 1

    if self.inputs is not None:
      self.axes[i_plot].cla()
      traj.plotThrust(self.axes[i_plot], plot_axis=self.inputs)
      i_plot += 1

    if self.save is not None:
      traj.save(self.save+'/iteration_%05d.csv' % self.i)
    else:
      plt.draw()
      self.fig.canvas.start_event_loop(0.0002)

    self.i += 1
    return [0]


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Plots a given trajectory from .csv format.')
  parser.add_argument('filename', metavar='file', type=str,
                      help='filename of the trajectory csv file')
  parser.add_argument('-p', '--position', dest='pos', type=str,
                      help='plot position over specified axes')
  parser.add_argument('-v', '--velocity', dest='vel', type=str,
                      help='plot velocity over specified axes')
  parser.add_argument('-w', '--omega', dest='omega', type=str,
                      help='plot omega over specified axes')
  parser.add_argument('-q', '--orientation', dest='ori', type=str,
                      help='plot orientation quaternion over specified axes')
  parser.add_argument('-u', '--thrust', dest='thrust', type=str,
                      help='plot thrust at rotors or acceleration over specified axes')
  parser.add_argument('-m', '--progress', dest='prog', type=str,
                      help='plot progress variables')
  parser.add_argument('-s', '--sub', dest='subp', metavar=('rows', 'cols'), type=int, nargs=2,
                      help='plot all in one figure arrange by rows*cols plots')
  parser.add_argument('-o', '--output', dest='output', type=str,
                      help='filename to output figure')
  args = parser.parse_args()

  traj = Trajectory(args.filename)
  print("Plotting trajectory with %d waypoints and %d nodes." % (traj.NW, traj.N))
  print("Overall time: %1.3f" % traj.t_total)

  if all(arg==None for arg in [
    args.pos, args.vel, args.omega, args.ori, args.thrust, args.prog]):
    args.subp = [2, 3]
    args.pos = 'xy'
    args.vel = 'xyza'
    args.omega = 'xyz'
    args.ori = 'wxyza'
    args.thrust = 'u'
    args.prog = 'mn'

  n_plots = 0
  if args.pos is not None: n_plots += 1
  if args.vel is not None: n_plots += 1
  if args.omega is not None: n_plots += 1
  if args.ori is not None: n_plots += 1
  if args.thrust is not None: n_plots += 1
  if args.prog is not None: n_plots += 1

  use_subplot = False
  if args.subp is not None:
    subr = args.subp[0]
    subc = args.subp[1]
    use_subplot = True
    if subr * subc < n_plots:
      Warning('Subplot too small')
      subc = math.ceil(n_plots/subr)
      
  if use_subplot:
    fig = plt.figure(0)
    
  axes = []
  for i in range(n_plots):
    if use_subplot:
      axes += [plt.subplot(subr, subc, i+1)]
    else:
      axes += [plt.figure(i)]

  i_fig = 0
  if args.pos is not None:
    traj.plotPos(axes[i_fig], plot_axis=args.pos, wp_style='rx')
    i_fig += 1

  if args.ori is not None:
    traj.plotOrientation(axes[i_fig], plot_axis=args.ori)
    i_fig += 1

  if args.prog is not None:
    traj.plotProgress(axes[i_fig], plot_axis=args.prog)
    i_fig += 1

  if args.vel is not None:
    traj.plotVel(axes[i_fig], plot_axis=args.vel)
    i_fig += 1

  if args.omega is not None:
    traj.plotOmega(axes[i_fig], plot_axis=args.omega)
    i_fig += 1

  if args.thrust is not None:
    traj.plotThrust(axes[i_fig], plot_axis=args.thrust)
    i_fig += 1

  if args.output is not None:
    plt.savefig(fig, args.output)
  else:
    plt.show()

import casadi as ca
import numpy as np
import csv
import matplotlib.pyplot as plt
from quaternion import rotate_quat
from mpl_toolkits.mplot3d import Axes3D
from quad import Quad
import warnings

# Plot and analyse trajectories
class Trajectory:
  def __init__(self, x=None, NPW=None, wp = [], NW = 0):
    assert type(wp) == ca.MX or type(wp) == ca.DM or type(wp) == np.ndarray or type(wp)==list
    if x is not None:
      assert isinstance(x, (ca.MX, ca.DM, np.ndarray, list, str))

    if x is None:
      self.NX = 0
      self.NU = 0
      self.x = []
      self.NPW = 0
      return
    elif isinstance(x, (ca.MX, ca.DM)):
      assert x.shape[1] == 1
      self.x = x.full().flatten()
    elif isinstance(x, (np.ndarray, list)):
      self.x = x
    elif isinstance(x, str):
      self.load(x)
      self.parse()
      return
    else:
      raise Exception('Unknown state vector x passed.')
    
    if NPW is None:
      raise Exception('Number of nodes per waypoint (NPW) must be provided!')

    self.NX = 13
    self.NU = 4
    self.NPW = NPW

    if type(wp) == ca.MX or type(wp) == ca.DM:
      self.wp = wp
      self.NW = self.wp.shape[1] 
    elif NW>0:
      self.NW = NW
    else:
      raise Exception('need valid waypointlist [wp] or number of waypoints [NW]!')

    self.N = self.NW * self.NPW

    self.parse()

  def parse(self):
    n_slice = self.NX + self.NU + 3 * self.NW
    n_start = 1 + self.NX + self.NW

    self.t_total = self.x[0]
    self.t_x = ca.DM(np.linspace(0, self.t_total, self.N+1, True))
    self.t_u = self.t_x[0:self.N]

    idx = np.array([0, *list(range(n_start+self.NU-1, len(self.x), n_slice))])

    self.p = np.array([
      self.x[1+idx],
      self.x[2+idx],
      self.x[3+idx]
    ])

    self.v = np.array([
      self.x[4+idx],
      self.x[5+idx],
      self.x[6+idx]
    ])

    self.q = np.array([
      self.x[7+idx],
      self.x[8+idx],
      self.x[9+idx],
      self.x[10+idx]
    ])

    self.w = np.array([
      self.x[11+idx],
      self.x[12+idx],
      self.x[13+idx]
    ])
    
    self.u = np.array([
      self.x[n_start+0::n_slice],
      self.x[n_start+1::n_slice],
      self.x[n_start+2::n_slice],
      self.x[n_start+3::n_slice]
    ])

    dt = self.t_total / self.N
    self.a_lin = np.zeros((3, self.N+1))
    self.a_rot = np.zeros((3, self.N+1))

    self.a_lin[:,0:-1] = np.diff(self.v) / dt
    self.a_rot[:,0:-1] = np.diff(self.w) / dt

    self.mu = np.zeros((self.NW, self.N+1))
    self.nu = np.zeros((self.NW, self.N))
    self.tau = np.zeros((self.NW, self.N))
    for i in range(self.NW):
      self.mu[i,:] = self.x[1+self.NX+i::n_slice]
      self.nu[i,:] = self.x[n_start+self.NU+self.NX+i::n_slice]
      self.tau[i,:] = self.x[n_start+self.NU+self.NX+self.NW+i::n_slice]

    self.thrust = np.zeros((3,self.N))
    self.dir = np.zeros((3,self.N))
    for i in range(self.N):
      self.thrust[:,i] = np.array(
        rotate_quat(ca.DM(self.q[:,i]), ca.vertcat(0, 0, ca.cumsum(self.u[:,i]))))[:,0]
      self.dir[:,i] = self.thrust[:,i] / np.linalg.norm(self.thrust[:,i])

  def unparse(self):
    n_slice = self.NX + self.NU + 3 * self.NW
    n_start = 1 + self.NX + self.NW

    self.x = np.zeros(n_start + self.N * n_slice)
    idx = np.array([0, *list(range(n_start+self.NU-1, len(self.x)-1, n_slice))])
    
    self.x[0] = self.t_total
    self.x[1+idx] = self.p[0,:]
    self.x[2+idx] = self.p[1,:]
    self.x[3+idx] = self.p[2,:]

    self.x[4+idx] = self.v[0,:]
    self.x[5+idx] = self.v[1,:]
    self.x[6+idx] = self.v[2,:]
    
    self.x[7+idx] = self.q[0,:]
    self.x[8+idx] = self.q[1,:]
    self.x[9+idx] = self.q[2,:]
    self.x[10+idx] = self.q[3,:]

    self.x[11+idx] = self.w[0,:]
    self.x[12+idx] = self.w[1,:]
    self.x[13+idx] = self.w[2,:]

    self.x[n_start+0::n_slice] = self.u[0,:]
    self.x[n_start+1::n_slice] = self.u[1,:]
    self.x[n_start+2::n_slice] = self.u[2,:]
    self.x[n_start+3::n_slice] = self.u[3,:]


  def getAxesHandle(self, fig, plot3d=False):
    kwargs = {}
    if plot3d:
      kwargs = {'projection': '3d'}
    if fig is None:
      return plt.figure(**kwargs).gca()
    elif type(fig) == int or type(fig) == str:
      return plt.figure(fig, **kwargs).gca()
    elif type(fig) == plt.Figure:
      return fig.gca(**kwargs)
    elif type(fig) == plt.Axes or type(fig) == plt.Subplot:
      return fig
    else:
      raise Exception('Provided figure or axis handle is invalid')

  def getDataAxes(self, axes_str, cset="xyz"):
    assert type(cset) == str
    assert type(axes_str) == str
    axes_str = axes_str.lower()
    cset = cset.lower()
    axes = []
    for i in range(len(axes_str)):
      axes += [cset.find(axes_str[i])]
      if axes[-1] < 0:
        raise Exception('Invalid axes %c specified' % axes_str[i])

    return axes

  def plotWaypoints(self, wp, fig=None, plot_axis='xy', style='rx', **kwargs):
    if wp is None: wp = self.wp
    p = []
    o = []

    if type(wp) == list:
      if len(wp) > 0:
        wp = ca.DM(wp)
      else:
        return p
    if type(wp) == ca.DM:
      if wp.shape[0]!=3:
        if wp.shape[1]==3:
          wp = wp.T
        else:
          raise Exception('waypoints have incorect format')
      elif wp.shape[1]<1:
        return p

    data = self.getDataAxes(plot_axis) 
    ax = self.getAxesHandle(fig, plot3d=(len(data)>2))

    if len(data) == 2:
      p = ax.plot(wp[data[0],:], wp[data[1],:], style, **kwargs)
      o = ax.plot(self.p[data[0],0], self.p[data[1],0], '.k')
    elif len(data) == 3:
      p = ax.plot(wp[data[0],:], wp[data[1],:], wp[data[2],:], style, **kwargs)
      o = ax.plot(self.p[data[0],[0]], self.p[data[1],[0]], self.p[data[2],[0]], '.k')
    else:
      raise Exception('No valid data axes specified to plot')
    plt.draw()
    return [o, p]

  def plotPos(self,
              fig=None, title=None, plot_axis='xy',
              wp=None, wp_style='rx', arrow_nth=None, arrow_size=0.5, arrow_args=None, **kwargs):
    if 'color' not in kwargs:
      kwargs['color'] = 'b'

    if arrow_args is None:
      arrow_args = {}
      arrow_args['color'] = kwargs['color']
    
    if not 'width' in arrow_args:
      arrow_args['width'] = 0.01
    arrow_args['zorder'] = 100

    if title is None:
      title = 'Position%'
    else:
      assert type(title) == str
    
    if arrow_nth is None:
      arrow_nth = int(self.N/10)

    title = title.replace('%', '   $t_{N}= %1.3fs$' % self.t_total)
    
    plot_ori = False
    if 'q' in plot_axis:
      plot_ori = True
      plot_axis = plot_axis.replace('q', '')

    data = self.getDataAxes(plot_axis) 
    ax = self.getAxesHandle(fig, plot3d=(len(data)>2))
    p = []
    pdir = []
    if len(data)==2:
      p += ax.plot(self.p[data[0],:], self.p[data[1],:], **kwargs)
      if plot_ori:
        idx = [i for i in range(0, self.N, arrow_nth)]
        px = [self.p[data[0], i] for i in idx]
        py = [self.p[data[1], i] for i in idx]
        dirx = [self.dir[data[0], i] for i in idx]
        diry = [self.dir[data[1], i] for i in idx]
        for i in range(len(idx)):
          pdir += [ax.arrow(px[i], py[i], arrow_size * dirx[i], arrow_size * diry[i], **arrow_args)]
      ax.set_xlabel('$p_%c$ $[m]$' % plot_axis[0])
      ax.set_ylabel('$p_%c$ $[m]$' % plot_axis[1])
      ax.axis('equal')
    elif len(data)==3:
      assert type(fig) == plt.Figure
      p += ax.plot(self.p[data[0],:], self.p[data[1],:], self.p[data[2],:], **kwargs)
      smax = np.max(self.p)
      smin = np.min(self.p)
      if plot_ori:
        idx = [i for i in range(0, self.N+1, arrow_nth)]
        dist = smax-smin
        px = [self.p[data[0], i] for i in idx]
        py = [self.p[data[1], i] for i in idx]
        pz = [self.p[data[2], i] for i in idx]
        dirx = [self.dir[data[0], i] for i in idx]
        diry = [self.dir[data[1], i] for i in idx]
        dirz = [self.dir[data[2], i] for i in idx]
        pdir += [ax.quiver(px, py, pz, dirx, diry, dirz)]
        # for i in range(len(idx)):
        #   pdir += [ax.arrow(px[i], py[i], dirx[i], diry[i], width=0.003*dist,
        #            ec='k', fc='k')]

      ax.set_xlabel('$p_%c$ $[m]$' % plot_axis[0])
      ax.set_ylabel('$p_%c$ $[m]$' % plot_axis[1])
      ax.set_zlabel('$p_%c$ $[m]$' % plot_axis[2])
      ax.set_xlim(smin, smax)
      ax.set_ylim(smin, smax)
      ax.set_zlim(smin, smax)
    else:
      raise Exception('No valid axes specified to plot')

    wpp = self.plotWaypoints(wp, fig, plot_axis, wp_style, ms=5)
    if len(title) > 0: ax.set_title(title)
    plt.draw()
    if len(pdir)>0: p += pdir
    if len(wpp)>0: p += wpp
    return p
  

  def plotVel(self, fig=None, title='Velocity', plot_axis='xyza', **kwargs):
    ax = self.getAxesHandle(fig)
    plot_abs = False
    if 'a' in plot_axis:
      plot_abs = True
      plot_axis = plot_axis.replace('a', '')
    data = self.getDataAxes(plot_axis)
  
    if len(data)<1 and not plot_abs: return ax 
    if len(data)>3:
      warnings.warn('Can only print 3 axes for velocity')

    p = [] 
    for i in range(min(len(data), 3)):
      if not 'label' in kwargs:
        kwargs['label'] = '$v_%c$' % plot_axis[i]
      p += ax.plot(self.t_x, self.v[data[i]], **kwargs)

    if plot_abs:
      if not 'label' in kwargs:
        kwargs['label'] = '$\|v\|$'
      p += ax.plot(self.t_x, np.linalg.norm(self.v, axis=0), **kwargs)

    if len(title) > 0:
      ax.set_title(title)
    ax.legend(loc='right')
    ax.set_xlabel('$t$ $[s]$')
    ax.set_ylabel("$v$ $[m/s]$")
    plt.draw()
    return p
  
  def plotOmega(self, fig=None, title='Bodyrate', plot_axis='xyz', **kwargs):
    ax = self.getAxesHandle(fig)
    data = self.getDataAxes(plot_axis)
  
    if len(data)<1: return fig 
    if len(data)>3:
      warnings.warn('Can only print 3 axes for bodyrate')

    p = [] 
    lgnd = []
    for i in range(min(len(data), 3)):
      p += ax.plot(self.t_x, self.w[data[i]], **kwargs)
      lgnd += ['$\omega_%c$' % plot_axis[i]]

    if len(title) > 0:
      ax.set_title(title)
    ax.legend(lgnd)
    ax.set_xlabel('$t$ $[s]$')
    ax.set_ylabel("$\omega$ $[rad/s]$")
    plt.draw()
    return p
  
  def plotOrientation(self, fig=None, title='Orientation', plot_axis='wxyz', **kwargs):
    ax = self.getAxesHandle(fig)
    plot_abs = False
    if 'a' in plot_axis:
      plot_abs = True
      plot_axis = plot_axis.replace('a', '')
    data = self.getDataAxes(plot_axis, 'wxyz')
  
    if len(data)<1: return fig 
    if len(data)>4:
      warnings.warn('Can only print 4 axes for bodyrate')

    p = [] 
    lgnd = []
    for i in range(len(data)):
      p += ax.plot(self.t_x, self.q[data[i]], **kwargs)
      lgnd += ['$q_%c$' % plot_axis[i]]

    if plot_abs:
      p += ax.plot(self.t_x, np.linalg.norm(self.q, axis=0), **kwargs)
      lgnd += ['$\|q\|$']

    if len(title) > 0:
      ax.set_title(title)
    ax.legend(lgnd)
    ax.set_xlabel('$t$ $[s]$')
    ax.set_ylabel("$q$ $[unit]$")
    plt.draw()
    return p

  def plotThrust(self,
    fig=None, title='Thrust', plot_axis='xyza', astyle='-', ustyle=None, **kwargs):
    if ustyle is None: ustyle = astyle

    ax = self.getAxesHandle(fig)
    plot_abs = False
    if 'a' in plot_axis:
      plot_abs = True
      plot_axis = plot_axis.replace('a', '')
    plot_thrusts = False
    if 'u' in plot_axis:
      plot_thrusts = True
      plot_axis = plot_axis.replace('u', '')
    data = self.getDataAxes(plot_axis, 'xyz')
  
    if len(data)>3:
      warnings.warn('Can only print 4 axes for directional force')

    p = []
    lgnd = []
    if plot_thrusts:
      for i in range(self.NU):
        p += ax.plot(self.t_u, self.u[i,:], ustyle, label='$u_%d$' % i, **kwargs)

    for i in range(min(len(data), 3)):
      p += ax.plot(self.t_u, self.m * self.a[data[i]], astyle, label='$f_%c$' % plot_axis[i], **kwargs)

    if plot_abs:
      if len(data) < 1:
        label = '$\|u\|$'
      else:
        label = '$\|f\|$'
      p += ax.plot(self.t_u, np.sum(self.u, axis=0), astyle, label=label, **kwargs)

    if len(title) > 0:
      ax.set_title(title)
    # ax.legend(lgnd)
    ax.set_xlabel('$t$ $[s]$')
    if plot_thrusts:
      ax.set_ylabel("$u$ $[N]$")
    else:
      ax.set_ylabel("$f$ $[N]$")
    plt.draw()
    return p

  def plotProgress(self,
    fig=None, title='Progress', plot_axis='mn', mstyle='-', nstyle='--', tstyle='.', **kwargs):

    lgnd = []
    ax = self.getAxesHandle(fig)

    plot_mu = 'm' in plot_axis
    plot_nu = 'n' in plot_axis
    plot_tau = 't' in plot_axis

    all_labels = (plot_mu ^ plot_nu ^ plot_tau) and self.NW <=6

    p = []
    for i in range(self.NW):
      prev = None
      if plot_mu:
        prev = ax.plot(self.t_x, self.mu[i,:], mstyle, **kwargs)
        if all_labels:
          lgnd +=['$\mu_%d$' % i]
        elif i==0:
          lgnd +=['$\mu$']
        p += prev
      if plot_nu:
        if prev is None:
          prev = ax.plot(self.t_u, self.nu[i,:], nstyle, **kwargs)
        else:
          prev = ax.plot(self.t_u, self.nu[i,:], nstyle, color=prev[0].get_color(), **kwargs)
        if all_labels:
          lgnd +=['$\\nu_%d $' % i]
        elif i==0:
          lgnd +=['$\\nu $']
        p += prev
      if plot_tau:
        if prev is None:
          prev = ax.plot(self.t_u, self.tau[i,:], tstyle, **kwargs)
        else:
          prev = ax.plot(self.t_u, self.tau[i,:], tstyle, color=prev[0].get_color(), **kwargs)
        if all_labels:
          lgnd +=['$\\tau_%d$' % i]
        elif i==0:
          lgnd +=['$\\tau$']
        p += prev

    if len(title) > 0:
      ax.set_title(title)
    ax.legend(lgnd, loc='right')
    ax.set_xlabel('$t$ $[s]$')
    ax.set_ylabel("progress")
    plt.draw()
    return p

  def save(self, filename, readable=False):
    assert type(filename) == str
    with open(filename, 'w') as csvfile:
      traj_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
      if not readable:
        traj_writer.writerow([self.NW, self.NPW, self.NX, self.NU, self.N])
        traj_writer.writerow(self.x)
        if type(self.wp) == np.ndarray:
          wp = self.wp.T.flatten()
        else:
          wp = self.wp.T.full().flatten()
        traj_writer.writerow(wp)
      if readable:
        labels = ['t', 'p_x', 'p_y', 'p_z',
                  'q_w', 'q_x', 'q_y', 'q_z',
                  'v_x', 'v_y', 'v_z',
                  'w_x', 'w_y', 'w_z',
                  'a_lin_x', 'a_lin_y', 'a_lin_z',
                  'a_rot_x', 'a_rot_y', 'a_rot_z',
                  'u_1', 'u_2', 'u_3', 'u_4']
        for i in range(self.NW):
          labels += ['mu_' + str(i)]
          labels += ['nu_' + str(i)]
          labels += ['tau_' + str(i)]
        
        traj_writer.writerow(labels)
        for i in range(self.N+1):
          # States
          row = [self.t_x[i],
            self.p[0,i], self.p[1,i], self.p[2,i],
            self.q[0,i], self.q[1,i], self.q[2,i], self.q[3,i],
            self.v[0,i], self.v[1,i], self.v[2,i],
            self.w[0,i], self.w[1,i], self.w[2,i],
            self.a_lin[0, i], self.a_lin[1, i], self.a_lin[2, i],
            self.a_rot[0, i], self.a_rot[1, i], self.a_rot[2, i]]

          # Inputs
          if i<self.N:
            for j in range(self.NU):
              row += [self.u[j,i]]
          else:
            row += [0]*self.NU

          # Progress
          for j in range(self.NW):
            row += [self.mu[j,i]]
            if i>0:
              row += [self.nu[j,i-1]]
              row += [self.tau[j,i-1]]
            else:
              row += [0]*2

          traj_writer.writerow(row)

  def load(self, filename):
    assert type(filename) == str
    with open(filename, 'r') as csvfile:
      traj_reader = csv.reader(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
      params = next(traj_reader)
      [self.NW, self.NPW, self.NX, self.NU, self.N] = params[:5]
      if len(params)>5:
        self.m = params[5]
      self.NW = int(self.NW)
      self.NPW = int(self.NPW)
      self.NX = int(self.NX)
      self.NU = int(self.NU)
      self.N = int(self.N)
      self.x = np.array(next(traj_reader))
      wp = np.array(next(traj_reader))
      self.wp = wp.reshape(self.NW, 3).T

      self.parse()
      return self


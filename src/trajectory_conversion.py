from trajectory import Trajectory
from quaternion import *
from pandas import read_csv

def viconToTrajectory(filename, interest=None, dt=0.005, wp=None):
  data = read_csv(filename)
  keys = data.keys()
  t = a = w = None
  if 'ts' in keys:
    t = data['ts'].to_numpy().T
  else:
    t = dt * np.arange(len(data))
  
  p = data[['TX', 'TY', 'TZ']].to_numpy().T
  r = data[['RX', 'RY', 'RZ']].to_numpy().T
  if ' accSmooth[0] (m/s/s)' in keys:
    a = data[[' accSmooth[0] (m/s/s)', ' accSmooth[1] (m/s/s)', ' accSmooth[2] (m/s/s)']].to_numpy().T
  if ' gyroADC[0] (deg/s)' in keys:
    w = data[[' gyroADC[0] (deg/s)', ' gyroADC[1] (deg/s)', ' gyroADC[2] (deg/s)']].to_numpy().T

  if interest is not None:
    p = p[:,interest[0]:interest[1]]
    r = r[:,interest[0]:interest[1]]
    if a is not None: a = a[:,interest[0]:interest[1]]
    if w is not None: w = w[:,interest[0]:interest[1]]
    if t is not None: t = t[interest[0]:interest[1]]

  ## Selecting segment
  if interest is not None:
    idx = ~(np.isnan(p).any(axis=0) | np.isnan(r).any(axis=0))
    t = t[idx]
    p = p[:, idx]
    r = r[:, idx]
    if a is not None: a = a[:, idx]
    if w is not None: w = w[:, idx]

  n = len(t)

  # Polishing Data
  if (p>100.0).any(): p *= 1e-3
  q = np.zeros((4, n))
  for i in range(n):
    q[:,i] = angleAxisToQuaternion(axis=r[:,i]).flatten()

  wind = 11
  k = np.ones(wind, 'd')
  k /= np.linalg.norm(k)
  dt = np.diff(t)

  pf = np.apply_along_axis(lambda m: np.convolve(m, k, mode='same'), axis=1, arr=p)
  v = np.zeros((3, n))
  for i in range(n-1):
    v[:,i] = (pf[:,i+1] - pf[:,i]) / dt[i]
  inv_wind = int(wind-1)
  v[:,0:inv_wind] = np.repeat(v[:,[inv_wind]], inv_wind, axis=1)
  v[:,n-inv_wind:n] = np.repeat(v[:,[-1-inv_wind]], inv_wind, axis=1)


  if a is None:
    a = np.diff(v)
    for i in range(n-1):
      a[:,i] /= dt[i]
    a = np.hstack((a, a[:,[-1]]))

  n -= 1-(n%2)
  traj = Trajectory()
  traj.N = n-1
  traj.NX = 13
  traj.NU = 4
  if wp is not None:
    traj.NW = wp.shape[1]
  else:
    traj.NW = 1
  traj.NPW = int(n/traj.NW)
  if wp is not None:
    traj.wp = wp
  else:
    traj.wp = p[:,[-1]]

  traj.t_x = t[:n] - t[0]
  traj.t_u = t[:n-1] - t[0]
  traj.t_total = traj.t_x[-1] - traj.t_x[0]
  traj.p = p[:,:n]
  traj.v = v[:,:n]
  traj.q = q[:,:n]
  traj.a = a[:,:n]
  if w is not None: traj.w = w[:,:n]
  traj.dir = np.zeros((3, n))
  for i in range(n):
    traj.dir[:,i] = qRotv(traj.q[:,i], np.array([0, 0, 1]))
  traj.u = np.zeros((4,n-1))
  traj.mu = np.zeros((traj.NW, n))
  traj.nu = np.zeros((traj.NW, n))
  traj.tau = np.zeros((traj.NW, n-1))
  return traj

  
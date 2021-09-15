import sys
import os
BASEPATH = os.path.abspath(__file__).split('rpg_time_optimal', 1)[0]+'rpg_time_optimal/'
sys.path += [BASEPATH + 'src']
from trajectory import Trajectory
import matplotlib.pyplot as plt

path = './example'

traj = Trajectory(path + '/result_cpc_format.csv')
print('Total Time: %1.4fs' % traj.t_total)

fig_pos_xy = plt.figure(0, (3,3))
axhxy = fig_pos_xy.gca()
traj.plotPos(fig_pos_xy, '', plot_axis='xyq', arrow_size=3.0, arrow_nth=15, arrow_args={'color': 'k'})
axhxy.set_title('')
axhxy.grid(True)
fig_pos_xy.tight_layout()
fig_pos_xy.savefig(path + '/pos_xy.pdf')

fig_pos_xz = plt.figure(1, (3,3))
axhxz = fig_pos_xz.gca()
traj.plotPos(fig_pos_xz, '', plot_axis='xzq', arrow_size=3.0, arrow_nth=15, arrow_args={'color': 'k'})
axhxz.set_title('')
axhxz.grid(True)
fig_pos_xz.tight_layout()
fig_pos_xz.savefig(path + '/pos_xz.pdf')

fig_vel = plt.figure(2, (6, 2))
axhv = fig_vel.gca()
traj.plotVel(fig_vel, '', label=None)
axhv.grid(True)
fig_vel.legend(['$v_x$', '$v_y$', '$v_z$', '$\|v\|$'], loc='right')
fig_vel.tight_layout()
fig_vel.savefig(path + '/vel.pdf')

fig_prog = plt.figure(3, (6, 2))
axhm = fig_prog.gca()
traj.plotProgress(fig_prog, '')
axhm.grid(True)
fig_prog.tight_layout()
fig_prog.savefig(path + '/prog.pdf')

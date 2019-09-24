import matplotlib.pyplot as plt
from zjh_fermi_daily_geometry import Fermi_daily_geometry


time = '2018-07-25 00:08:00'

fermi_object = Fermi_daily_geometry(time)

mp = fermi_object.plot_position_on_earth(plot_fermi = True)
plt.savefig('/home/laojin/shiyan/plot_SAA.png')



















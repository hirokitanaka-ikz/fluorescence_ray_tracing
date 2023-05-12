import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

temp = 300

data_pi = np.loadtxt(f'YbYLF_pi_{temp}K.txt', skiprows=0)
data_sigma = np.loadtxt(f'YbYLF_sigma_{temp}K.txt', skiprows=0)

wl = data_pi[:,0]
If_pi = data_pi[:,1]
If_sigma = data_sigma[:,1]
abs_xs_pi = data_pi[:,3]
abs_xs_sigma = data_sigma[:,3]

fig = plt.figure(figsize=(4, 4))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
fig.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95, hspace=0.0)


ax1.plot(wl, If_sigma, label='E$\perp$c')
ax1.plot(wl, If_pi, label='E$\parallel$c')
ax1.set_ylabel('$I_f$ (a.u.)', fontsize=10)
ax1.set_ylim(0, 1.1)
ax1.set_xticks([])
ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))

ax2.plot(wl, abs_xs_sigma * 1e20, label='E$\perp$c')
ax2.plot(wl, abs_xs_pi * 1e20, label='E$\parallel$c')
ax2.set_xlabel('Wavelength (nm)', fontsize=10)
ax2.set_ylabel('$\sigma_\mathrm{abs}$ ($10^{-20}$ cm$^2$)', fontsize=10)
ax2.set_ylim(0, 2.9)
ax2.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))

ax1.text(x=930, y=0.95, s=f'T = {temp}K', fontsize=10)
ax1.legend(loc='upper right', frameon=False, fontsize=10)

for ax in [ax1, ax2]:
    ax.set_xlim(925, 1075)
    ax.tick_params(axis='both', labelsize=9)

plt.show()
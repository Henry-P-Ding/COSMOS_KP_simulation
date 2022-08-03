import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
import matplotlib.font_manager as font_manager
import scipy.constants as const
from matplotlib import rcParams
import matplotlib
matplotlib.rc('font',family='cmss10')

a, P, E_range = 1, 1.33, np.arange(0, 40, 0.1)
alpha_a = np.sqrt(E_range) * a
# V0b product equals 1b
ka_values = np.arccos(np.cos(alpha_a) + P * np.sin(alpha_a) / alpha_a)
cos_ka_values = np.cos(alpha_a) + P * np.sin(alpha_a) / alpha_a
band_slices = np.ma.clump_unmasked(np.ma.masked_invalid(ka_values))
low_ka = ka_values[band_slices[0]]
low_ka = np.concatenate([np.flip(np.multiply(-1, low_ka)), low_ka])
high_ka = ka_values[band_slices[1]]
high_ka = np.concatenate([np.multiply(-1, high_ka), np.flip(high_ka)])
low_E = E_range[band_slices[0]]
low_E = np.concatenate([np.flip(low_E), low_E])
high_E = E_range[band_slices[1]]
high_E = np.concatenate([high_E, np.flip(high_E)])

''' for conclusion sections
metal_E = valid_E[valid_E < 5]
print(metal_E.size)
metal_E = np.concatenate([np.flip(metal_E), metal_E])
print(metal_E.size)
metal_ka = valid_ka[valid_E < 5]
metal_ka = np.concatenate([np.flip(np.multiply(-1, metal_ka)), metal_ka])

c_metal_E = valid_E[valid_E < 7]
print(c_metal_E.size)
c_metal_E = np.concatenate([np.flip(c_metal_E)[40:], c_metal_E])
print(c_metal_E.size)
c_metal_ka = valid_ka[valid_E < 7]
c_metal_ka = np.concatenate([np.flip(np.multiply(-1, c_metal_ka))[40:], c_metal_ka])

insulator_E = low_E
insulator_ka = low_ka
'''
plt.figure(figsize=(14, 8))
ax = plt.gca()
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)
# Move left y-axis and bottim x-axis to centre, passing through (0,0)
ax.spines['left'].set_position('center')
# Eliminate upper and right axes
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.set_ylabel(r'$E\;(\mathrm{eV})$', rotation=0, fontsize=18, loc='top', math_fontfamily='cm')
plt.xlabel(r'$k\;(\mathrm{m^{-1}})$', fontsize=18, math_fontfamily='cm')
ax.yaxis.set_label_coords(0.45, 0.95)
ax.legend(frameon=False, fontsize=18)
plt.xticks(fontfamily="serif")
plt.yticks(fontfamily="serif")
#ax.axes.xaxis.set_ticks([])
#ax.axes.yaxis.set_ticks([])
titles = ['Filled Energies for an Insulator', 'Filled Energies for a Conducting Metal', 'Filled Energies for a Ground State Metal', 'Kronig-Penney Dispersion Relation']
ax.set_title(titles[0], pad=20,fontsize=36, font="cmss10")

ev_hbar = const.physical_constants['Planck constant in eV/Hz'][0]
print(f"a: {a}, P: {P}")

plt.plot(low_ka / a, low_E * ev_hbar * ev_hbar / 2 / const.m_e, '-', color='#bfdfb0', linewidth=4.0, label='Energy Bands')
plt.plot(high_ka / a, high_E * ev_hbar * ev_hbar / 2 / const.m_e, '-', color='#bfdfb0', linewidth=4.0)
plt.plot([low_ka[np.argmax(low_E)] / a + 0.5, low_ka[np.argmax(low_E)] / a + 0.5], np.array([np.max(low_E), np.min(high_E)]) * ev_hbar * ev_hbar / 2 / const.m_e, color='black')
plt.plot([low_ka[np.argmax(low_E)] / a + 0.5, low_ka[np.argmax(low_E)] / a + 0.5], np.array([np.max(low_E), np.min(high_E)]) * ev_hbar * ev_hbar / 2 / const.m_e, '_', color='black')
plt.text(low_ka[np.argmax(low_E)] / a + 0.7, (np.max(low_E) + np.min(high_E)) * ev_hbar * ev_hbar / 4 / const.m_e, 'band gap', fontsize=18, fontfamily="cmss10")

#plt.plot(metal_ka / a, metal_E * ev_hbar * ev_hbar / 2 / const.m_e, linewidth=6, color="#5f6353", label='Filled Energies')
#plt.plot(c_metal_ka / a, c_metal_E * ev_hbar * ev_hbar / 2 / const.m_e, linewidth=6, color="#5f6353", label='Filled Energies')
#plt.arrow(0.4, 100, 1.3, 0, width=2, head_width=8, head_length=0.06, edgecolor=None, color='black')
#plt.annotate(r'$\mathbf{\vec{E}}$', (1, 110), fontsize=18, math_fontfamily='cm')
#plt.plot(insulator_ka / a, insulator_E * ev_hbar * ev_hbar / 2 / const.m_e, '-', color='#5f6353', linewidth=6, label='Filled Energies')

print(f"Band Gap: {(np.min(high_E) - np.max(low_E)) * ev_hbar * ev_hbar / 2 / const.m_e} eV")

plt.show()


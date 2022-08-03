import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
import matplotlib.font_manager as font_manager
import scipy.constants as const
from matplotlib import rcParams
import matplotlib
import sys
matplotlib.rc('font', family='cmss10')

EV_HBAR = const.physical_constants['Planck constant in eV/Hz'][0]

a, b_min, b_max, db, V_0, E_min, E_max, dE = 1, 0.05, 0.32, 0.01, 250, 0, 80, 0.1


def calculate_bands(a, b, V_0, E_min, E_max, dE):
    print(f"P = {const.m_e * a * b * V_0 / EV_HBAR / EV_HBAR}")
    E_range = np.arange(E_min, E_max, dE)
    alpha = np.sqrt(E_range)
    beta = np.sqrt(E_range + V_0)
    # V0b product equals 1b
    ka_values = np.arccos(
        np.cos(beta * b) * np.cos(alpha * (a - b)) + (np.square(alpha) + np.square(beta)) / 2 / alpha / beta * np.sin(
            beta * b) * np.sin(alpha * (a - b)))
    band_slices = np.ma.clump_unmasked(np.ma.masked_invalid(ka_values))
    low_ka = ka_values[band_slices[0]]
    low_ka = np.concatenate([np.flip(np.multiply(-1, low_ka)), low_ka])
    high_ka = ka_values[band_slices[1]]
    high_ka = np.concatenate([np.multiply(-1, high_ka), np.flip(high_ka)])
    low_E = E_range[band_slices[0]]
    low_E = np.concatenate([np.flip(low_E), low_E])
    high_E = E_range[band_slices[1]]
    high_E = np.concatenate([high_E, np.flip(high_E)])
    return low_ka, low_E, high_ka, high_E


def setup_dispersion_plot():
    plt.figure(figsize=(14, 8))
    ax = plt.gca()
    for axis in ['top', 'bottom', 'left', 'right']:
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
    # ax.axes.xaxis.set_ticks([])
    # ax.axes.yaxis.set_ticks([])
    titles = ['Filled Energies for an Insulator', 'Filled Energies for a Conducting Metal',
              'Filled Energies for a Ground State Metal', 'Kronig-Penney Dispersion Relation']
    ax.set_title(titles[0], pad=20, fontsize=36, font="cmss10")


def setup_plot_bg_b_plot(b_list, gap_list):
    plt.figure(figsize=(14, 8))
    ax = plt.gca()
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(2)
    # Move left y-axis and bottim x-axis to centre, passing through (0,0)
    # Eliminate upper and right axes
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.set_ylabel(r'$E_{bg}\;(\mathrm{eV})$', rotation=0, fontsize=18, loc='top', math_fontfamily='cm')
    plt.xlabel(r'$b\;(\mathrm{m})$', fontsize=18, math_fontfamily='cm')
    ax.legend(frameon=False, fontsize=18)
    plt.xticks(fontfamily="serif")
    plt.yticks(fontfamily="serif")
    ax.set_title('Band Gap v.s. Well Width', pad=20, fontsize=36, font="cmss10")

    plt.plot(b_list, gap_list, '-', linewidth=4.0, color='#5c7750')
    plt.show()


def plot_bands(low_ka, low_E, high_ka, high_E):
    plt.plot(low_ka / a, low_E * EV_HBAR * EV_HBAR / 2 / const.m_e, '-', color='#bfdfb0', linewidth=4.0,
             label='Energy Bands')
    plt.plot(high_ka / a, high_E * EV_HBAR * EV_HBAR / 2 / const.m_e, '-', color='#bfdfb0', linewidth=4.0)
    plt.plot([low_ka[np.argmax(low_E)] / a + 0.5, low_ka[np.argmax(low_E)] / a + 0.5],
             np.array([np.max(low_E), np.min(high_E)]) * EV_HBAR * EV_HBAR / 2 / const.m_e, color='black')
    plt.plot([low_ka[np.argmax(low_E)] / a + 0.5, low_ka[np.argmax(low_E)] / a + 0.5],
             np.array([np.max(low_E), np.min(high_E)]) * EV_HBAR * EV_HBAR / 2 / const.m_e, '_', color='black')
    plt.text(low_ka[np.argmax(low_E)] / a + 0.7, (np.max(low_E) + np.min(high_E)) * EV_HBAR * EV_HBAR / 4 / const.m_e,
             'band gap', fontsize=18, fontfamily="cmss10")


b_list = np.arange(b_min, b_max, db)
gap_list = []
'''
for b in b_list:
    print(f"b = {b}")
    low_ka, low_E, high_ka, high_E = calculate_bands(a, b, V_0, E_min, E_max, dE)
    bg = (np.min(high_E) - np.max(low_E)) * EV_HBAR * EV_HBAR / 2 / const.m_e
    gap_list.append(bg)
    print(f"Band Gap: {bg} eV")
setup_plot_bg_b_plot(b_list, gap_list)
'''

for d in np.arange(0.2, 0.3, 0.1):
    a = d
    b = d
    print(f"a = {a}, b = {b}")
    low_ka, low_E, high_ka, high_E = calculate_bands(a, b, V_0, E_min, E_max, dE)
    bg = (np.min(high_E) - np.max(low_E)) * EV_HBAR * EV_HBAR / 2 / const.m_e
    gap_list.append(bg)
    print(f"Band Gap: {bg} eV")
setup_plot_bg_b_plot(np.arange(1, 2, 0.1), gap_list)

setup_dispersion_plot()
low_ka, low_E, high_ka, high_E = calculate_bands(a, 0.2, V_0, E_min, E_max, dE)
plot_bands(low_ka, low_E, high_ka, high_E)
plt.show()
#plot_bands(low_ka, low_E, high_ka, high_E)
#plt.show()


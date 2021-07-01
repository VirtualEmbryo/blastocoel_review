#!/usr/bin/env python3
# blastula.py

import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arc

from ipywidgets import interactive

import ipywidgets as widgets

import blastula_class as blc


def calc_l_ext(Re, R, angle, d) :
    return np.cos(angle)*d + np.sqrt(R**2 - d**2 * np.sin(angle)**2)

def calc_l_in(Re, R, angle, d) :
    return np.cos(angle)*d - np.sqrt(R**2 - d**2 * np.sin(angle)**2)

def calc_psi_ext(Re, R, angle, d) :
    l = calc_l_ext(Re, R, angle, d)
    return np.arcsin( l *np.sin(angle)/R )

def calc_psi_in(Re, R, angle, d) :
    l = calc_l_in(Re, R, angle, d)
    return np.arcsin( l *np.sin(angle)/R )

c_list = ['#780096','#FF36B8','#FF0700','#62FFFF','#009BFF','#0000FF','#D9E34C','#63DB00','#00B400']


def plot_interactive_RaRbNbd(Ra, Rb, Nb, d, Re): 
    e = blc.Embryo(R_embryo=Re, Nb=Nb)
    e.plot_embryo(d=d, Ra=Ra, Rb=Rb, centers=False, circles=False, print_areas=True)
    #print(e.Nb, e.cell_area, Ra, Rb, d)
    #s = str(e.Nb) + '\t' + str(e.cell_area) + '\t' + str(Ra) + '\t' + str(Rb) + '\t' + str(d) + '\n'
    psi_a = calc_psi_ext(Re, Ra, angle=pi/Nb, d=Re-Ra)
    H = Ra+Rb+d
    L = 2*Ra*np.sin(psi_a)
    r = L/H
    print('=====================')
    print('Aspect ratio = ', "{:2.4f}".format(r))
    #save = True
    #if save :
    #    filename = 'nocavity_areas.dat'
    #f.open(filename, 'a').write(s).close()
    #    f.write()
    #    f.close()
    #    print('Saved !')
    #return s

def main() :
    interactive_plot = interactive(plot_interactive_RaRbNbd, Nb=(0, 30), d=(-3.00, 3.), Ra=(0., 1.), Rb=(0., 1.), Re=(0., 5.))
    output = interactive_plot.children[-1]
    output.layout.height = '450px'
    #f.write(interactive_plot.children[-1].outputs[-1]['text'], 'a')

    #interactive_plot
    plt.show(block=True)

    #interactive_plot.children[-1].outputs[-1]['text']

main()

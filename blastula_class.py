import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arc

from ipywidgets import interactive

import ipywidgets as widgets

def plot_vector(pos, magnitude, angle, ax=None, color='k', ls='-', lw=1) :
    X, Y = ((pos[0], pos[0]+magnitude*np.cos(angle)), (pos[1], pos[1]+magnitude*np.sin(angle))) 
    if ax is None :
        plt.plot(X, Y, color=color, linestyle=ls, linewidth=lw)
    else :
        ax.plot(X, Y, color=color, linestyle=ls, linewidth=lw)
    return;

class Embryo:
    def __init__(self, R_embryo, Nb):
        self.R_embryo = R_embryo
        self.Nb = Nb
        self.theta_embryo = pi/Nb
        
    def calc_area(self, Ra, Rb, d) :
        self.total_area = pi*self.R_embryo**2
        
        blastomere = Blastomere(Ra, Rb, d, self.theta_embryo, self.R_embryo, Theta=0.)
        db = blastomere.calc_db(self.R_embryo, Ra, d)
        psi_b = blastomere.calc_psi_in(self.R_embryo, Rb, self.theta_embryo, db)
        
        self.cell_area = self.Nb * blastomere.calc_surf_blastomere(self.R_embryo, self.theta_embryo, Ra, Rb, d, print_areas=False)
        self.cavity_area = self.Nb * blastomere.calc_surf_cavity(Rb, psi_b, theta_e=self.theta_embryo)
        
        return {'total':self.total_area, 'cell':self.cell_area, 'cavity':self.cavity_area}
        
    def plot_embryo(self, d, Ra, Rb, lines=False, centers=False, circles=False, print_areas=False, figsize=(5, 5), savefig=False, savename='pic.eps', cavity_circle=False) :
        fig, ax = plt.subplots(figsize=figsize)
        if lines :
            ax.plot((0, self.R_embryo*np.cos(self.theta_embryo)), (0, self.R_embryo*np.sin(self.theta_embryo)), color='k', linestyle='--', alpha=0.1)
            ax.plot((0, self.R_embryo*np.cos(self.theta_embryo)), (0, -self.R_embryo*np.sin(self.theta_embryo)), color='k', linestyle='--', alpha=0.1)
    
        draw_circle_embryo = plt.Circle((0, 0), radius=self.R_embryo, fill=0, linestyle='--', linewidth=2)
        ax.add_artist(draw_circle_embryo)
    
        self.blastomere_list = []
        
        for n in range(self.Nb) :
            Theta_n = n*2*pi/self.Nb
            self.blastomere = Blastomere(Ra, Rb, d, self.theta_embryo, self.R_embryo, Theta=Theta_n)
            self.blastomere_list += [self.blastomere]
            res = self.blastomere.plot_blastomere(ax=ax, circles=circles, center=centers, print_areas=print_areas)
            if res is None :
                break ;
            
        try :
        #if 1:
            
            if print_areas :
                self.calc_area(Ra, Rb, d)
                print('Total embryo area = ', self.total_area)
                print('Total blastomere area = ', self.cell_area)
                print('Single blastomere area = ', self.cell_area/self.Nb)
                print('Cavity volume = ', self.cavity_area)
                #print('Left space = ', self.total_area - (self.cell_area + self.cavity_area))
        except :
            print('Problem in calculating the areas') 
            pass
            
        # Optional circle (for cavity)
        if cavity_circle :
            R_cavity = self.R_embryo-(Ra+Rb+d)
            draw_circle_cavity = plt.Circle((0, 0), radius=R_cavity, fill=0, linestyle='--', color='green', linewidth=2, alpha=0.5)
            ax.add_artist(draw_circle_cavity)
            print('R_cavity = ', R_cavity)

        ax.axis('equal')
        ax.set_xlim(-self.R_embryo*1.1, self.R_embryo*1.1)
        ax.set_ylim(-self.R_embryo*1.1, self.R_embryo*1.1)
    
        ax.axis('off')
        if savefig :
            plt.savefig(savename, format='eps')
        
        plt.show()

class Blastomere:
    def __init__(self, Ra, Rb, d, theta_embryo, R_embryo, Theta=0.) :
        self.Ra = Ra
        self.Rb = Rb
        self.d = d
        self.theta_embryo = theta_embryo
        self.R_embryo = R_embryo
        self.Theta = Theta
        
    def check_variables_blastomere(self):
        # Check conditions
        ## Minimal radii
        Ra_min = self.R_embryo*np.sin(self.theta_embryo)/(1.+np.sin(self.theta_embryo))
        Rb_min = (self.R_embryo-(self.Ra+self.d))*np.sin(self.theta_embryo)
        if self.Ra <= Ra_min or self.Rb <= Rb_min :
            print('Radii too small !')
            if self.Ra <= Ra_min :
                print('Ra too small (Ra > '+ "{:2.8f}".format(Ra_min)+ ')')
            elif self.Rb <= Rb_min :
                print('Rb too small (Rb > '+ "{:2.8f}".format(Rb_min)+ ')')
            return ;
        ## Minimal distance
        R_cavity = self.R_embryo - (self.Ra+self.Rb+self.d)
        if R_cavity < 0 :
            print('Cavity disappears ! Check for the condition R_embryo (= '+"{:2.2f}".format(self.R_embryo)+') >= Ra+Rb+d (='+"{:2.2f}".format(R_cavity)+')')
            return ;

        ## Negative distance
        self.da = self.R_embryo - self.Ra
        self.db = self.R_embryo - (self.Ra+self.d)
        self.la = self.calc_l_ext(self.R_embryo, self.Ra, self.theta_embryo, self.da)
        self.lb = self.calc_l_in(self.R_embryo, self.Rb, self.theta_embryo, self.db)
        if self.la < self.lb :
            print('Negative lateral length L<0 ! Check for the condition l_a (= '+"{:2.2f}".format(self.la)+') >= lb (='+"{:2.2f}".format(self.lb)+')')
            return ;

        return self.da, self.db, self.la, self.lb
    
    def calc_blastomere(self) :
        try :
            da, db, la, lb = self.check_variables_blastomere()
            print('Test done')
        except :
            pass
            return ;

        self.psi_a = self.calc_psi_ext(self.R_embryo, self.Ra, self.theta_embryo, self.da)*180/pi
        self.psi_b = self.calc_psi_in(self.R_embryo, self.Rb, self.theta_embryo, self.db)*180/pi

        A = self.calc_surf_blastomere(self.R_embryo, self.theta_embryo, self.Ra, self.Rb, self.d)

        return A
        
    def calc_psi_ext(self, Re, R, angle, d) :
        l = self.calc_l_ext(Re, R, angle, d)
        return np.arcsin( l *np.sin(angle)/R )

    def calc_l_ext(self, Re, R, angle, d) :
        return np.cos(angle)*d + np.sqrt(R**2 - d**2 * np.sin(angle)**2)

    def calc_psi_in(self, Re, R, angle, d) :
        l = self.calc_l_in(Re, R, angle, d)
        return np.arcsin( l *np.sin(angle)/R )

    def calc_l_in(self, Re, R, angle, d) :
        return np.cos(angle)*d - np.sqrt(R**2 - d**2 * np.sin(angle)**2)
    
    def calc_da(self, Re, Ra) :
        return Re-Ra

    def calc_db(self, Re, Ra, d) :
        return Re-(Ra+d)

    def calc_lap(self, Re, Ra, d, theta_e) :
        da = self.calc_da(Re, Ra)
        return da*np.cos(theta_e) + np.sqrt(Ra**2 - da**2 * np.sin(theta_e)**2)

    def calc_lam(self, Re, Ra, d, theta_e) :
        da = self.calc_da(Re, Ra)
        return da*np.cos(theta_e) - np.sqrt(Ra**2 - da**2 * np.sin(theta_e)**2)

    def calc_lbp(self, Re, Ra, d, theta_e) :
        db = self.calc_db(Re, Ra, d)
        return db*np.cos(theta_e) + np.sqrt(Rb**2 - db**2 * np.sin(theta_e)**2)

    def calc_lbm(self, Re, Ra, d, theta_e) :
        db = self.calc_db(Re, Ra, d)
        return db*np.cos(theta_e) - np.sqrt(Rb**2 - db**2 * np.sin(theta_e)**2)

    def calc_ratio_gamma_a(self, Re, Ra, d, theta_e) :
        """gamma_L/gamma_a"""
        lap = self.calc_lap(Re, Ra, d, theta_e)
        return np.sin(np.arcsin(lap*np.sin(theta_e)/Ra) - theta_e)

    def calc_ratio_gamma_b(self, Re, Rb, d, theta_e) :
        """gamma_L/gamma_b"""
        lbm = self.calc_lbm(Re, Rb, d, theta_e)
        return np.sin(np.arcsin(lbm*np.sin(theta_e)/Rb) - theta_e)

    def calc_surf_cap(self, radius, angle) :
        return 0.5*(angle-np.sin(angle))*radius**2

    def calc_Hab(self, d, Ra, Rb, psi_a, psi_b) :
        ha = Ra*(1. - np.cos(psi_a))
        hb = Rb*(1. - np.cos(psi_b))
        Hab = d + Ra + Rb - (ha + hb)
        return Hab

    def calc_surf_blastomere(self, R_embryo, theta_embryo, Ra, Rb, d, print_areas=False) :
        da = R_embryo - Ra
        db = R_embryo - (Ra+d)

        psi_a = self.calc_psi_ext(R_embryo, Ra, theta_embryo, da)
        psi_b = self.calc_psi_in(R_embryo, Rb, theta_embryo, db)
        beta = pi-psi_b

        a, b = Ra*np.sin(psi_a), Rb*np.sin(beta)

        Hab = self.calc_Hab(d, Ra, Rb, psi_a, psi_b)

        A1 = Hab*b
        A2 = 0.5*Hab*(a-b)
        A3 = self.calc_surf_cap(Rb, 2*psi_b)
        A4 = self.calc_surf_cap(Ra, 2*psi_a)

        self.Alist = [A1, A2, A3, A4]
        self.area = 2*A1 + 2*A2 + A3 + A4
        return self.area

    def calc_surf_cavity(self, Rb, psi_b, theta_e) :
        A_cap = self.calc_surf_cap(Rb, 2*psi_b)
        A_tr = Rb**2 * np.sin(psi_b)**2 / np.tan(theta_e)
        return A_tr - A_cap

    def plot_blastomere(self, ax=None, circles=False, lines=True, center=False, print_areas=False, savefig=False, savename='pic.eps') :
        if ax is None : 
            fig, ax = plt.subplots(1, 1, figsize=(5, 5))
            ax.axis('equal')
            ax.set_xlim(0., self.R_embryo*1.1)
            ax.set_ylim(-self.R_embryo*0.55, self.R_embryo*0.55)
            
        
        try :
            self.da, self.db, self.la, self.lb = self.check_variables_blastomere()
        except : 
            return None
        

        # Apical Side
        self.psi_a = self.calc_psi_ext(self.R_embryo, self.Ra, self.theta_embryo, self.da)*180/pi

        self.x_apical_up = self.la*np.cos(self.theta_embryo+self.Theta)
        self.x_apical_down = self.la*np.cos(-self.theta_embryo+self.Theta)
        self.y_apical_up = self.la*np.sin(self.theta_embryo+self.Theta)
        self.y_apical_down = self.la*np.sin(-self.theta_embryo+self.Theta)

        draw_circle_apical = plt.Circle((self.da*np.cos(self.Theta), self.da*np.sin(self.Theta)), radius=self.Ra, fill=0, linestyle='-', color='#CB1423', linewidth=2, alpha=0.1)
        draw_arc_apical = Arc((self.da*np.cos(self.Theta), self.da*np.sin(self.Theta)), width=2*self.Ra, height=2*self.Ra, angle=self.Theta*180/pi, theta1=(-self.psi_a)%360, theta2=+self.psi_a%360, linestyle='-', color='#CB1423', linewidth=2)

        if center :
            ax.scatter(self.da*np.cos(self.Theta), self.da*np.sin(self.Theta), color='#CB1423')

        # Basal side
        self.psi_b = self.calc_psi_in(self.R_embryo, self.Rb, self.theta_embryo, self.db)*180/pi
        #psi_b = calc_psi_ext(R_embryo, Rb, theta_embryo, db)*180/pi

        self.x_basal_up = self.lb*np.cos(self.theta_embryo+self.Theta)
        self.x_basal_down = self.lb*np.cos(-self.theta_embryo+self.Theta)
        self.y_basal_up = self.lb*np.sin(self.theta_embryo+self.Theta)
        self.y_basal_down = self.lb*np.sin(-self.theta_embryo+self.Theta)

        draw_circle_basal = plt.Circle((self.db*np.cos(self.Theta), self.db*np.sin(self.Theta)), radius=self.Rb, fill=0, linestyle='-', color='#36509B', linewidth=2, alpha=0.1)
        draw_arc_basal = Arc((self.db*np.cos(self.Theta), self.db*np.sin(self.Theta)), width=2*self.Rb, height=2*self.Rb, angle=self.Theta*180/pi, theta1=(-self.psi_b+180)%360, theta2=(+self.psi_b+180)%360, linestyle='-', color='#36509B', linewidth=2)
        if center :
            ax.scatter(self.db*np.cos(self.Theta), self.db*np.sin(self.Theta), color='#36509B')


        # Drawings
        ax.add_artist(draw_arc_apical)
        ax.add_artist(draw_arc_basal)

        if lines :
            ax.plot((self.x_basal_up, self.x_apical_up), (self.y_basal_up, self.y_apical_up), color='k', linewidth=2)
            ax.plot((self.x_basal_down, self.x_apical_down), (self.y_basal_down, self.y_apical_down), color='k', linewidth=2)

        if circles :
            ax.add_artist(draw_circle_apical)
            ax.add_artist(draw_circle_basal)

        # Calculate the blastomere area :
        self.A_blastomere = self.calc_surf_blastomere(self.R_embryo, self.theta_embryo, self.Ra, self.Rb, self.d)
        
        if savefig :
            plt.savefig(savename, format='eps')
        
        return [self.x_apical_up, self.y_apical_up, self.psi_a*pi/180], [self.x_basal_up, self.y_basal_up, self.psi_b*pi/180], self.A_blastomere

    def optimize_volume(self, A, val_name, val_list, print_val=False) :
        a = np.zeros(len(val_list))
        
        def give_imin(list_values, target) :
            i_min = np.argmin(np.abs(list_values-target))
            return i_min
        
        if val_name == 'd' :
            for i in range(len(val_list)) :
                d = val_list[i]
                a[i] = self.calc_surf_blastomere(self.R_embryo, self.theta_embryo, self.Ra, self.Rb, d)
            i_min = give_imin(a, A)
            self.d = val_list[i_min]
                
        elif val_name == 'Ra' :
            for i in range(len(val_list)) :
                Ra = val_list[i]
                a[i] = self.calc_surf_blastomere(self.R_embryo, self.theta_embryo, Ra, self.Rb, self.d)
            i_min = give_imin(a, A)
            self.Ra = val_list[i_min] 
                
        elif val_name == 'Rb' :
            for i in range(len(val_list)) :
                Rb = val_list[i]
                a[i] = self.calc_surf_blastomere(self.R_embryo, self.theta_embryo, self.Ra, Rb, self.d)
            i_min = give_imin(a, A)
            self.Rb = val_list[i_min]
                
        elif val_name == 'Nb' :
            for i in range(len(val_list)) :
                theta = pi/val_list[i]
                a[i] = self.calc_surf_blastomere(self.R_embryo, theta, self.Ra, self.Rb, self.d)
            i_min = give_imin(a, A)
            
        return val_list[i_min], a[i_min]
            
#
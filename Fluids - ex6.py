import scipy as sp
import matplotlib
import matplotlib.pyplot as py

#Some global settings for matplotlib plots
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['font.weight'] = 'bold'

#Homogeneous liquid phase diagram
def evaluateMiniumum(rho,beta,mu):
    '''Evaluates the equation for the homogeneous density profile
    which minimizes the grand potential'''
    if sp.all((rho>=0)==(rho<=1))==False:
        raise Exception("Density is unphysical")
    elif beta<0.0:
        raise Exception("Negative temperature")
    return rho-(1.0-rho)*sp.exp(beta*(mu+5*rho))

def main():
    rho = sp.linspace(0,1,100)  #Density 

    #Evaluate the density profiles for the various parameter regimes
    muVals = [-3.0,-2.5,-2.0]   #These are the values of chemical potential we want
    energies_lowT = [evaluateMiniumum(rho,1.0,k) for k in muVals]
    energies_highT = [evaluateMiniumum(rho,float(2/3),k) for k in muVals]

    leg = [muVals]+[muVals]
    #Plotting
    figs,ax = py.subplots(2,1)
    titles = ['Low temperature','High temperature']
    for axis in ax:
        axis.set_xlabel('Density')
    for curve in energies_lowT:
        ax[0].plot(rho,curve)
    for curve in energies_highT:
        ax[1].plot(rho,curve)
    for plot in zip(ax,titles,leg):
        plot[0].set_title(plot[1])
        plot[0].legend(plot[2])
        plot[0].plot(rho,sp.zeros(len(rho)),'k--')
    py.show()

if __name__== "__main__":
    main()

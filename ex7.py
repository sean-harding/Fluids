import scipy as sp
import matplotlib
import matplotlib.pyplot as py
import time as t
import math

#Some global settings for matplotlib plots
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['font.weight'] = 'bold'

'''
> Here I will implement an iterative solver for the bulk phase diagram.
We wish to identify the regions which correspond to the bulk liquid and gas phases
as well as the critical points and coexistence curves. 

> The bulk phases should correspond to a single unique solution to the self-consistent 
density equation.

> By definition, the coexistence region of the phase diagram will have multiple solutions,
such that the solution is unstable with respect to perturbations in the initial guess for the
density. i.e. at coexisence, we should find multiple solutions whereby if we start with a high
density guess, we converge to the high density solution and vice versa for a low density guess

> I will calculate the phase diagram numerically and then plot the analytic results on top.
'''
def getF(rho,beta,mu):
    '''Calculates RHS of equation 48'''
    if sp.all((rho>=0)==(rho<=1))==False:
        raise Exception("Density is unphysical")
    elif beta<0.0:
        raise Exception("Negative temperature")
    return (1.0-rho)*sp.exp(beta*(mu+5*rho))
    #return (1.0-rho)*sp.exp(beta*(4*rho+rho+mu))
def iterate(rho_0,t,mu,mixing=0.01,cTol=10**-6,iterMax = 10000):
    '''Iteratively solves x = f(x)'''
    convergence = 100
    iteration = 0
    dens = rho_0
    while (convergence>cTol) and (iteration<iterMax):
        dens_new = (1-mixing)*dens + mixing*getF(dens,1/t,mu)   #t is temperature, 
        convergence = abs(dens_new-dens)
        #print(convergence)
        iteration = iteration + 1
        dens = dens_new
        #Under iteration, the density can become arbitrarily large
        #when we approach the liquid phase. Hence if the density
        #diverges we break the loop and set the density equal to 1.
        #This might give us sharp unphysical features in the liquid phase region, so we should be careful
        #With i.e. mixing parameters etc. to smooth out any unphysical jumps
        if dens>1.0:
            dens = 1.0
            break
        elif dens<0.0:
            dens = 0.0
            break
        if iteration==iterMax-1:
            print("Not converged after {} iterations for parameters (T={},mu={})".format(iteration,t,mu))
    return dens

sampleRate = 100                             #We will have sampleRate**2 points to evaluate
dx = 5/sampleRate
dy = (2-1e-04)/sampleRate 
muVals = sp.linspace(-5,0,sampleRate)        #Values of chemical potential, in units of the interaction energy
tVals = sp.linspace(1e-03,2.0,sampleRate)  #Values of temperature 1/beta in units of epsilon. We use grown-up units so Kb=1 

density = sp.zeros(shape=[sampleRate,sampleRate],dtype='float32')

time_0 = t.perf_counter()
for n in range(sampleRate):
    for m in range(sampleRate): 
        densities = sp.array([iterate(k,tVals[n],muVals[m]) for k in [0.05,0.95]])  #The code works best if we start from the extreme liquid or gas limits
        if sp.any(abs(densities-densities[0])>0.1)==True:     #If starting from different initial points leads to different
            density[n,m] = sp.nan                             #convergence, we are in the coexistence region
            pass
        else:
            density[n,m] = sp.average(densities)                            
time_1 = t.perf_counter()

print('Elapsed: {}s'.format(time_1-time_0))

stable = sp.ma.masked_where(sp.isnan(density),density)  #Here we mask the NaNs I have stored to designate the coexistence region so matplotlib doesn't plot them
x,y = sp.meshgrid(muVals,tVals)                         #x,y, define the corners of a rectangular grid where we have evaluated z at each corner. (Maybe not meshgrid)
contourlevels = py.MaxNLocator(nbins=100).tick_values(stable.min(), stable.max())   #This bins the values of the density into "nbins" bins 
cmap = py.get_cmap('plasma')
norm = matplotlib.colors.BoundaryNorm(contourlevels, ncolors=cmap.N, clip=True)

#Plot raw data
figs,ax = py.subplots(2,1)
titles = ['Raw','Contour']

for axis in ax:
    axis.set_xlabel('Chemical Potential')
    axis.set_ylabel('Temperature')

#The way this works is that it draws squares with corners at the gridpoints specified by x and y.
#It then gives the value of the function evaluated at the upper left corner to the square. Hence we should remove
#the last row and column of our data because this would be asigned to a square outside of our grid.
raw = ax[0].pcolormesh(x,y,stable[:-1,:-1],cmap=cmap,norm=norm)
figs.colorbar(raw,ax=ax[0])

#Plot contour plot

#Contour plots work in the opposite way to the colormesh plots. We need to shift the co-ordinate axes so that the points
#we are plotting are in the centre of each square and not the corners.
contour = ax[1].contourf(x[:-1,:-1]+dx/2,y[:-1,:-1]+dy/2,stable[:-1,:-1],levels=contourlevels,cmap=cmap)
figs.colorbar(contour,ax=ax[1])

#Calculate the two branches of the spinodal line, and the chemical potential
#Plotting the binodal is slightly more effort; in principle it could be done with a function that can find
#density as a function of T by solving the binodal equation self-consistently, and then this goes
#into the equation for mu(T) like as for the spinodal line.
splus = lambda T:(1+sp.sqrt(1-4*T/5))/2
sminus = lambda T:(1-sp.sqrt(1-4*T/5))/2
chem = lambda rho,T:T*math.log(rho/(1-rho))-5*rho 

spin_plus = map(splus,tVals)
mu_plus = map(chem,spin_plus,tVals)
spin_minus = map(sminus,tVals)
mu_minus = map(chem,spin_minus,tVals)

mu_plus = sp.array(list(mu_plus))
mu_minus = sp.array(list(mu_minus))
f1 = sp.where(sp.isreal(mu_plus))
f2 = sp.where(sp.isreal(mu_minus))

for plot in zip(ax,titles):
    plot[0].set_title(plot[1])
    plot[0].plot(mu_plus[f1],tVals[f1],'r--',linewidth=5)
    plot[0].plot(mu_minus[f2],tVals[f2],'r--',linewidth=5)
py.show()

print('Done')
import scipy as sp
import matplotlib.pyplot as py
from functools import partial

'Part 1: Solve the simple two-body problem with one-body fixed'

#For motion in a plane we have four coupled differential equations
#for the (x,y) co-ordinate and (vx,vx) velocities/momenta
#To solve this problem we implement the following algorithm:
#1.) Initialize with (x0,y0), (vx0,vy0)
#2.) We have set up the RK4 equations which evalutate kn = fn(x0,y0,vx0,vy0)
#    need to now map the f's across various inputs to get the kn at each step

class body:
    ''' Defines a massive body object with various properties as defined in __init__'''
    def __init__(self,initState,M):
        self.state = initState      #Current state-point
        self.M = M                  #Mass
        #self.F = [lambda var:0]*2   #External acting on body, i.e. not due to the other bodies in motion. Initialized to zer0
    def getR(self):
        '''Return position co-ordinates'''
        return (self.state[1],self.state[2])
    def getV(self):
        '''Return velocity co-ordinates'''
        return (self.state[3],self.state[4])
    def Ek(self):
        '''Returns kinetic energy'''
        return 0.5*self.M*sum(v**2 for v in self.getV)
    def propagate(self,F,dx):
        '''Perform one RK4 update for a single body in motion around another'''
        while True:
            var = self.state[1:]  
            K1 = [f(var) for f in F]
            var = [s+dx*k/2 for (s,k) in zip(var,K1)]
            K2 = [f(var) for f in F]
            var = [s+dx*k/2 for (s,k) in zip(var,K2)]
            K3 = [f(var) for f in F]
            var = [s+dx*k for (s,k) in zip(var,K3)]
            K4 = [f(var) for f in F]
            #Update the state. First list is the timestep, second list is the RK update. Does not yet update internally
            yield [self.state[0]+dx] + [x+dx*(k1+2*(k2+k3)+k4)/6 for (x,k1,k2,k3,k4) in zip(self.state[1:],K1,K2,K3,K4)],dx

    def propagateMulti(self,bodies,dx):
        '''Perform one RK4 update for a body in motion relative to another'''
        #'var' is x,y,vx,vy
        def force(k,bodies,var):
            '''Calculates component k of the force on self due to the bodies in bodies'''
            #First construct the distance squared between the body and all other bodies
            sepVectors = ([a-b for a,b in zip((var[0],var[1]),b.getR())] for b in bodies)
            #Now I construct the components of force acting on self
            return -1*sum(b.M*s[k]/sum(a**2 for a in s)**1.5 for b,s in zip(bodies,sepVectors))
        F = [lambda var:var[2],lambda var: var[3]]
        F.extend([partial(force,k,bodies) for k in (0,1)])
        while True:
            var = self.state[1:]                        #t,x,y,vx,vy
            K1 = [f(var) for f in F]
            var = [s+dx*k/2 for (s,k) in zip(var,K1)]
            K2 = [f(var) for f in F]
            var = [s+dx*k/2 for (s,k) in zip(var,K2)]
            K3 = [f(var) for f in F]
            var = [s+dx*k for (s,k) in zip(var,K3)]
            K4 = [f(var) for f in F]
            #Update the state. First list is the timestep, second list is the RK update.
            yield [self.state[0]+dx] + [x+dx*(k1+2*(k2+k3)+k4)/6 for (x,k1,k2,k3,k4) in zip(self.state[1:],K1,K2,K3,K4)],dx
'Part 1: Simple orbit around a massive object'
mSun = 1                                #Use units of mSun*G = 1

earth = body([0,0,7,0.36,0],mSun)         #Set up initial state of the FF. These parameters give a decently stable orbit
sun = body([0,0,0,0,0],mSun)

trajectory = []                         #Store all position co-ordinates. Could also choose to store only a few
dev = []                                #Honestly, should probably use scipy arrays here
k=0
newState = earth.propagateMulti([sun],10**-2)
while len(dev)<5:
    newPoint,dx = next(newState)           #Calculate next state
    if (earth.getR()[0]<0 and earth.getR()[1]>0) and newPoint[1]>0:
        dev.append(newPoint[1:3])
    earth.state = newPoint               #Propagate one timestep and update the state
    trajectory.append(earth.getR())
    if k%10==0:                         #Store each 10th co-ordinate point
        trajectory.append(earth.getR())
    k+=1
trajectory = sp.array(trajectory)
#dev = [(d-7)*100/7 for d in dev]
#py.plot(dev)
#py.ylabel('Percent deviation ')
py.plot(trajectory[:,0],trajectory[:,1])
py.plot(0,0,'o')
py.show()

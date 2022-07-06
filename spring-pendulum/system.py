from vpython import *
import numpy as np
import matplotlib.pyplot as plt
import random
import sympy as sym
import os


class System():
    def __init__(self,numPend: int, len: float, g:float, 
                k: float, mass: float, 
                b:float, T: int) -> None:
        """Construct the system and initialize position values

        Args:
            numPend (int): Number of pendulums to be used in the system
            len (float): Length of the pendulums
            g (float): Gravitational constant being used
            k (float): Coefficient of restitution of the springs
            mass (float): Mass of the pendulums
            b (float): Air resistance coefficient 

        Returns:
            None
        """
        self.numPend = numPend
        self.len = len
        self.g = g
        self.k = k
        self.m = mass
        self.b = b
        self.gamma = self.b/self.m
        self.sep = len*2  # separation between adjacent pendulums
        self.thetas = np.zeros(self.numPend)  # initial values for each angle
        self.xs = np.zeros(self.numPend)
        self.ys = np.zeros(self.numPend)
        self.zs = np.zeros(self.numPend)
        self.rs = np.zeros((self.numPend,3))  # array of positions as vectors
        self.nails = np.zeros((self.numPend,3))  # top of the string of pendulums
        self.t = 0
        self.fps = 100
        self.numIter = T*self.fps 
        self.dt_theta = np.zeros(self.numPend)
        self.dtdt_theta = np.zeros(self.numPend)
        self.kineticEs = np.zeros(self.numIter)
        self.gravitEs = self.numIter
        self.elasticEs = self.numIter
        self.totEs = self.numIter
        """initialization of values"""
        for i in range(self.numPend):
            self.thetas[i] = random.randint(-45,45) # in degrees
            self.thetas[i] = self.thetas[i]*np.pi/180
            self.xs[i] = self.len*np.sin(self.thetas[i])+self.sep*(i+1)
            self.ys[i] = self.len*np.cos(self.thetas[i]) 
            self.rs[i] = [self.xs[i], self.ys[i], self.zs[i]]
            self.nails[i] = [self.sep*(i+1), 0, 0]
        return None
    
    
    def create_visuals(self) -> None:
        """Create animated objects

        Returns:
            None: Only makes visuals vpython
        """
        self.roof = box(color=vector(.39,.3,.3), 
            pos=vector((self.numPend+1)*self.sep/2, 0,0),
            length=(self.numPend+1)*self.sep, height=0.03, width=0.5)
        
        self.string = curve(vector(0,0,0))
        for i in range(self.numPend):
            self.string.append(vector(self.nails[i][0], self.nails[i][1], self.nails[i][2]),
                               vector(self.xs[i], -self.ys[i], self.zs[[i]]),
                               vector(self.nails[i][0], self.nails[i][1], self.nails[i][2])) 
            self.string.append(vector((self.numPend+1)*self.sep, 0, 0))
        
        self.ball = [None]*self.numPend
        for i in range(self.numPend):
            self.ball[i] = sphere(pos=vector(self.xs[i], -self.ys[i], self.zs[i]),
                             radius=0.04)
        """for the springs connecting adjacent pendulums"""
        self.spring = [None]*self.numPend
        for i in range(1,self.numPend):
            self.spring[i] = helix(color=vector(.8,.8,.8),thickness=np.sqrt(
            self.k)*.005,coils=11,radius=.01)
            self.spring[i].pos = self.ball[i-1].pos  # direction of spring
            self.spring[i].axis = self.ball[i].pos - self.ball[i-1].pos 
            
        self.gd = graph(width=600, height=350,
            title='<b>Energies</b>',
            xtitle='<i>t</i>', ytitle='<i>E(t)</i>',
            foreground=color.black, background=color.white,
            xmin=0, xmax=self.numIter/self.fps)
        self.kineticC = gcurve(color = color.red, graph = self.gd,
        label = "Kinetic Energy")
        self.potC = gcurve(color = color.blue, graph = self.gd,
        label = "Potential Energy")
        self.totC = gcurve(color = color.black, graph = self.gd, 
        label = "Total Energy")
        return None

   
    def set_camera(self) -> None:
        """Set the camera position
        """
        scene.camera.pos = vector((self.numPend+1)*self.sep/2,-0.5*self.len,(
            (self.numPend+1)*self.sep/16)**0.5)


    def update_values(self) -> None:
        """Compute the angular acceleration and update the values of theta using
        Euler integration

        Returns:
            None: Updates values internally in the System class
        """    
        
        
        len_vec = [np.zeros(3)]*self.numPend
        d = [np.zeros(3)]*(self.numPend-1)
        dnorm = [np.zeros(3)]*(self.numPend-1)  # Normalized distance vector
        I = self.m*self.len**2  # moment of inertia
        grav_force = np.array([0, self.m*self.g, 0])
                
        for i in range(0, self.numPend): 
            len_vec[i] = self.rs[i] - self.nails[i]
            if i != 0:
                d[i-1] = self.rs[i] - self.rs[i-1]
                dnorm[i-1] = d[i-1]/np.linalg.norm(d[i-1])
            """torque by gravitational force"""
            torque = np.cross(len_vec[i], grav_force)
            """torque from left spring"""
            if i != 0:
                torque = torque-np.cross(
                    len_vec[i], self.k*(np.linalg.norm(d[i-1])-self.sep)*dnorm[i-1])
            """torque from right spring"""
            if i != self.numPend-1:
                torque = torque+np.cross(
                    len_vec[i], self.k*(np.linalg.norm(d[i])-self.sep)*dnorm[i])
            self.dtdt_theta[i] = -torque[2]/I
            

            """torque by air resistance (force inversely proportional to velocity)"""
            self.dtdt_theta[i] = self.dtdt_theta[i] - self.gamma*self.dt_theta[i]*self.len/I
            
        for i in range(self.numPend):
            self.dt_theta[i] = self.dt_theta[i] + 1/self.fps * self.dtdt_theta[i]
            self.thetas[i] = self.thetas[i] + 1/self.fps * self.dt_theta[i]

        for i, theta in enumerate(self.thetas): 
            self.xs[i] = self.len*np.sin(theta)+self.sep*(i+1)
            self.ys[i] = self.len*np.cos(theta) 
            self.rs[i] = [self.xs[i], self.ys[i], self.zs[i]]
            
        return None


    def calc_energies(self) -> None:
        """Calculate all types of energy in the system

        Returns:
            None 
        """
        self.kinetic = 0
        self.gravitE = 0
        self.elastic = 0
        self.totE = 0
        
        for i, omega in enumerate(self.dt_theta):
            self.kinetic += self.m*(omega*self.len)**2 / 2
            self.gravitE += (-self.ys[i])*self.m*self.g
            if i != 0:
                dis = self.rs[i] - self.rs[i-1]
                self.elastic += (np.linalg.norm(dis) - self.sep)**2 *self.k/2
        self.totE += self.kinetic+self.gravitE+self.elastic
        return None
    
    
    def plot_energies(self) -> None:
        """Plot the energies in animated graph
        """
        self.kineticC.plot(self.t/self.fps, self.kinetic)
        self.potC.plot(self.t/self.fps, self.gravitE+self.elastic)
        self.totC.plot(self.t/self.fps, self.totE)
        return None


    def update_visuals(self) -> None:
        """Modify the visuals according to the modified values

        Returns:
            None
        """
        for i in range(self.numPend):
            self.ball[i].pos = vector(self.xs[i], -self.ys[i], self.zs[i])            
            self.string.modify(4*i+2, self.ball[i].pos)
            if (i != 0):
                self.spring[i].pos = self.ball[i-1].pos
                self.spring[i].axis = self.ball[i].pos-self.ball[i-1].pos
        
        return None   


    def animate(self) -> None:
        """Implement the entire animation process (updating steps repeatedly)
        """
        self.set_camera()
        self.create_visuals()
        sleep(0.5)
        while(self.t <= self.numIter):
            rate(self.fps)
            self.update_values()
            self.update_visuals()
            self.calc_energies()
            self.plot_energies()
            self.t += 1
        # os.popen('import -window 0x3a00003 frames/vp'+str(self.t).zfill(4)+'.gif')
        return None

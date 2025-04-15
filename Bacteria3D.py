'''
This script contains the Bacterium3D class. This class is used to decribe a 
bacterium that is present within a clinostat. 
'''

# Imports #####################################################################

# modules
import numpy as np

# external files
import Functions as f

#np.random.seed(66)

###############################################################################


class Bacteria3D(object):
    '''
    Class used to describe bacteria modelled as point particles in fluid moving
    as a solid body within a clinostat.
    '''
    
    def __init__(self, mass, position, radius, swimming_vel, organism_density, medium_density):
        """
        Initialises a point particle in 3D space

        :param mass: float, mass of the bacteria
        :param position: [3] float array w/ position
        :param radius: float, radius of bacteria assumed to be shape sphere
        """
        
        self.mass = float(mass) # bouyant bacterial mass in kg
        self.pos = np.array(position, float) # bacterial position in [x,y,z], each component in m
        self.rad = float(radius) # bacterial radius in m
        self.organism_density = float(organism_density)
        self.bouyant_mass(medium_density) # calculate the bouyant density of the bacterium in the clinostat media
        
        # initialising direction of swimming
        self.swim = float(swimming_vel) # swimming speed in m/s, float
        self.swim_direction = f.initialise_swimming_direction() # random initial unit direction for swimming velocity, 3D array of mag 1
        #self.swim_direction = np.array([0, 1, 0])
        #print(self.swim_direction)
        #self.swim_direction = np.array([1, 2, 3])
        self.swim_vel = self.swim*self.swim_direction # swimming velocity in m/s, 3D array [vx, vy, vz]
        
    
    def __str__(self):
        """
        XYZ-compliant string. The format is
        <label>    <x>  <y>  <z>
        """
        xyz_string = f"{self.pos[0]}  {self.pos[1]} {self.pos[2]}" # xyz coords in m from origin
            
        return xyz_string
    
    
    def bouyant_mass(self, medium_density):
        '''
        Calculate the bouyant mass of a bacterium inside clinostat medium
        
        :param medium_density: float, density in kg/m^3 for the same liquid
        '''
        
        # bouyant mass assuming a spherical bacterium of constant density 
        self.bm = (4/3)*np.pi*medium_density*(self.rad**3)*((self.organism_density/medium_density) - 1)
        print(self.bm)
        
    
    # storing old terminal velocity just in case needed
    def terminal_vel(self, viscosity_coeff, density, g):
        '''
        Calculates velocity at t = 0 for each bacteria instance.
        
        It is assumed that organism is at terminal velocity at t = 0 seconds, i.e sedimentation
        speed through the culture medium is equal to terminal velocity.
        
        :param viscosity_coeff: float, viscosity coefficient in PaS for liquid in sim.
        :param density: float, density in kg/m^3 for the same liquid
        :param g: float, gavitational constant in kg/m^2 for desired environment
        '''   
        
        VTy = self.bm*g/(6*np.pi*viscosity_coeff*self.rad) # magnitude of terminal velocity
        #print(VTy)
        
        self.term_vel = np.array([0, -VTy, 0]) # negative comes from coordinate definition, positive y goes up 
        print(self.term_vel)
        
        
    def centripetal_force(self, viscosity_coeff, fluid_density, omega, planar_position, dt, status):
        '''
        Calculate velocity due to centripetal (fugal) force
        
        :param viscosity_coeff: float, viscosity coefficient in PaS for liquid in sim.
        :param fluid_density: float, density in kg/m^3 for the same liquid
        :param omega: float, density of object in clinostat, kg/m^3
        :param planar_position: numpy array [3], xyz position of bacterium with z set to 0.
        :param status: String, True means centripetal force is on, False means its off
        '''   
        
        # centripetal force is on
        if status == 'True':
            
            # viscously overdamped dynamics 
            self.centripetal_vel = (self.bm*(omega**2)*planar_position)/(6*np.pi*viscosity_coeff*self.rad)
         
        # centripetal force is off thus set to zero
        elif status == 'False':
            self.centripetal_vel = [0, 0, 0]
        
        
    def rotational_vel(self, omega):
        '''
        This function calculates the rotational velocity of the bacteria
        at its current position.
        
        :param omega: float, rotational speed of clinostat in rad/s
        '''
        x = self.pos[0] # current x position of bacteria
        y = self.pos[1] # current y position of bacteria
        
        self.rot_vel = np.array([-y*omega, x*omega, 0]) # updating the current rotational velocity of the bacteria 
    
    
    @staticmethod
    def tumble_probability(dt, tumbling_rate):
        '''
        Calculates if a bacterium is going to tumble in the current timestep.
        This allows implementation of run and tumble motion into bacteriums motion.
        
        Assuming that the rate of tumbling is once per second
        
        :param dt: float, timestep of simulation in s
        :param tumbling_rate: integer, number of tumbles per second in bacteriums motion
        '''
    
        # counting number of decimal places in timestep float
        #dps = str(0.00001)[::-1].find('.')
        dps = str(dt)[::-1].find('.')
        #print(dps)
        
        # generating random number between 0 & 1 with same timestep as dt        
        random_number = round(np.random.uniform(0, 1-dt), dps) # 1-dps ensures that when tumbling_rate is 0, you dont get any tumbles at all
        #print(random_number)
        #tumbling_rate = 1 # tumble per second
        
        # threshold that allows bacterium to tumble
        tumble_prob = 1 - (tumbling_rate*dt) # unitless
        #print(tumble_prob)
        # if randomly generated number greater than or equal to theshold
        # bacterium tumbles
        if random_number >= tumble_prob:
            #print('Tumbling........................')
            does_bacterium_tumble = 1
            #print(random_number)
            #print(tumble_prob)
        
            
        else: # bacterium doesnt tumble
            #print('No Tumbles')
            does_bacterium_tumble = 0  
        
        return does_bacterium_tumble # 1 = does tumble, 0 = does not tumble
        
    
    def update_swimming_vel(self, omega, rotational_diffusion_coefficient, dt, does_bacterium_tumble):
        '''
        This function updates the swimming velocity of the bacterium.
        
        :param omega: float, rotational speed of clinostat in rad/s
        :param rotational_diffusion_coefficient: float, rotational diffusion coefficient in m^2/s
        :param dt: float, timestep of simulation in s
        '''
        # If tumble_probabilty is true then the bacterium tumbles in a random direction
        if does_bacterium_tumble == 1:
            #print('tumbling.....................................')
            
            # random new tumbling direction
            new_direction = f.initialise_swimming_direction()           
            magnitude = np.linalg.norm(new_direction)
        
            # updating the swimming direction and then the swimming velocity with that vector
            self.swim_direction = new_direction/magnitude # normalising for unit vector  
            self.swim_vel = self.swim*self.swim_direction
       
        # if tumble probability is false it swims in its new direction
        else:
            #print('swimming')
            # xyz components of current swimming direction unit vector
            ex = self.swim_direction[0]
            ey = self.swim_direction[1]
            
            #print(self.swim_direction)
            #ez = self.swim_direction[2]
        
            # DEFINING RATE OF CHANGE OF SWIMMING UNIT VECTOR #####
        
            # rotation term in rate of change of swimming direction
            rotation = omega*np.array([-ey, ex, 0])
        
            # diffusion term in rate of change of swimming direction
            coeff = np.sqrt((2*rotational_diffusion_coefficient)/dt) # coefficient on second term of rate of change vector
            noise = np.random.normal(0, 1, size=3) # different from noise vector in diffusion velocity (avoids coupling)
            delta = np.identity(3) # rank 2 tensor
            outer_product = np.outer(self.swim_direction, self.swim_direction) # outer product of swimming unit vectors
        
            diffusion_term = coeff*((delta - outer_product)@noise)
        
            # rate of change of the swimming unit vector
            dedt = rotation + diffusion_term
            #print(dedt)
        
            # calculating new direction and 
            new_direction = (dedt*dt) + self.swim_direction
            magnitude = np.linalg.norm(new_direction)
        
            # updating the swimming direction and then the swimming velocity with that vecot
            self.swim_direction = new_direction/magnitude # normalising for unit vector
            #print(self.swim_direction)
            self.swim_vel = self.swim*self.swim_direction
            #print(np.linalg.norm(self.swim))
        


    def update_vel(self, dt, diffusion_coefficient):
        '''
        Calculates velocity at t = 0 for each bacteria instance.
        
        :param dt: float, timestep of simulation
        :param diffusion_coefficient: float, diffusion coefficent of bacteria in medium, in m^2/s
        :param noise: [3] float array, random 3D noise vector taken from a normal distribution (cartesian coords)
        :param rot_vel: [3] float array, 3D vector representing the rotational velocity of a bacteria (cartesian coords)
        
        '''   
        # assuming to be at terminal velocity once in system (can add factor from 23/9/23 notes in not assuming this)
        # currently only using Vt, Vd and Vr but Vs will be implemented
  
        # generating a noise vector from a normal distribution
        noise = np.random.normal(0, 1, size=3)
        #print(np.linalg.norm(noise))
        
        # diffusion
        diffusion = noise*np.sqrt(2*diffusion_coefficient/dt)
        
        self.vel = self.term_vel + diffusion + self.rot_vel + self.swim*self.swim_direction + self.centripetal_vel
        #print(self.term_vel)
       # print(self.vel, self.term_vel, diffusion, self.rot_vel, self.swim*self.swim_direction, self.centripetal_vel)
        
        #print(np.linalg.norm(0.5*self.mass*(self.vel**2)))
        
        #print((self.mass*np.abs(self.pos[1])*9.8) + np.linalg.norm(0.5*self.mass*(self.vel**2)))
 
        
    def update_pos(self, dt):
        """
        Updates the position of a Particle 3D instance
        
        :param dt: float, timestep
        """
        
        self.pos += self.vel*dt # adding new position vector to previous position vector
        
        
    @staticmethod # meaning it could be written as an independant function 
    def new_b3d(file_handle, medium_density):
        """
        Initialises a Particle3D instance given an input file handle.
        
        The input file should contain one line per particle in the following
        format:
        <mass>  <radius>  <x> <y> <z>  <vs> 
        
        :param file_handle: Readable file handle in the above format
        :param medium_density: float, density in kg/m^3 for the same liquid
        :return: Particle3D instance
        """
        data = file_handle.readline() # file_handle will contain information about bacteria
        p = data.split() # getting info for each individual bacteria
        
        mass = float(p[0]) # bacterial mass in kg
        radius = float(p[1]) # radius of bacteria in m
        position = np.array([p[2] ,p[3] ,p[4]], float) # x,y,z component in m from origin
        swimming_vel = float(p[5]) # swimming velocity of bacteria in m/s
        organism_density = float(p[6])

        return Bacteria3D(mass, position, radius, swimming_vel, organism_density, medium_density) # instance of Bacteria3D class
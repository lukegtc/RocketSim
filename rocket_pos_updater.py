from class_file import *
import numpy as np
import matplotlib.pyplot as plt
from Atmosphere import *
from math import *
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

#rocket created
rocket1 = Rocket(prop_mass=418700.,mass_empty=27200.)
#prop created
propellant = Propellant(mass=rocket1.prop_mass,radius=rocket1.radius,height=rocket1.height)

#engine Created
engine1 = Engine(mdot=845)
engine2 = Engine(mdot=845)
#rocket thrust force
rocket1.thrust = np.array([0,0,7605e3])
actual_latitude = 0
latitude = (actual_latitude-90)*pi/180
longitude = 0*pi/180
#exit velocity
ve_tot = 700
ve = np.array([0,0,ve_tot])
dt = 0.1
t = 0
engine1.thrust = ve*engine1.mdot
#atmosphere values
layers = [0.,11.,20.,32.,47.,51.,71.,86.]
avals = [-6.5,0,1,2.8,0,-2.8,-2.]
base_temp = 288.15
base_press = 101325.
r = 287.
#Planet created
g0=9.80665
GM = 3.986004418e14
  # 28 degrees latitude
earth = Planet(6371e3,g0,layers,avals,base_temp,base_press,r,GM)
engine1.pos = np.array([0.,0.,0. + earth.radius]) #Changed




running =True
zpos = []
ypos = []
#zpos = [earth.radius*sin(latitude),earth.radius*sin(latitude)]
#ypos = [-earth.radius*cos(latitude),-earth.radius*cos(latitude)]
xpos = []
time1 = []
rads = []
time = []
ttlforce = []
dragforce = []
densitys = []
temps = []
skin_temps = []
presss = []
vels = []
cg_vals = []
moments= []
rocket_vel = np.array([0.,0.,0.]) #Dont forget to change some stuff about relative velocity and position ya know
off_pad = False
#values for circular orbit
orbital_alt = 500e3
orbit_vel = sqrt(earth.GM/(earth.r+orbital_alt))
radius = earth.radius
skin_temp = 0
angular_vel = np.array([0.,0.,0.])
total_gimbal = np.array([0.,0.,0.])
while running:
    tot_prop_cg = propellant.prop_cg(rocket1.prop_mass)
    propellant.tot_prop_cg = tot_prop_cg
    rocket1.tot_cgpos = rocket1.tot_cgpos1(tot_prop_cg)




    #estimate from NASA estimated vals
    radius = sqrt(engine1.pos[0]**2 + engine1.pos[1]**2 + (engine1.pos[2])**2)

    #calculates pressure, temperature and density using the atmossolver equation from the atmosphere file
    press,temp,density=atmossolver(earth.g0,earth.layers,earth.avals,radius-earth.radius,earth.r,earth.base_temp,earth.base_press)


    engine1.isp = ve/earth.g0
    #engine1 gimballing
    #Relation of gimballed engine to pitch, roll and yaw angle of rocket
    #This is the gimbal of the engine
    # if rocket1.prop_mass<= 0.25*418700:
    #
    #     engine1.gimbal = np.array([pi/2/1800,0,0])
    #
    # else:
    #     engine1.gimbal[0] = 0
    # if engine1.gimbal[2]>= -np.pi/2:
    #     engine1.gimbal-=np.array([0,0,pi/2/1800])
    #
    # else:
    #     engine1.gimbal[2] = -pi/2
    matrices = Matrices()
    roll_mat = matrices.roll_matrix(engine1.gimbal[0]) #around x axis
    pitch_mat = matrices.pitch_matrix(engine1.gimbal[1]) #around the y axis
    yaw_mat = matrices.yaw_matrix(engine1.gimbal[2]) #around the z axis

    #this is the gimbal of the rocket


    #moment of inertia
    moi_mat = np.array([[(rocket1.prop_mass+rocket1.mass_empty+rocket1.payload_mass)*(3*(rocket1.diameter/2)**2 +
                                rocket1.height**2)/12 + (rocket1.mass_empty+rocket1.payload_mass + rocket1.prop_mass)*
                                  (rocket1.tot_cgpos[1]**2 + rocket1.tot_cgpos[2]**2),0.,0.],
                  [0.,(rocket1.prop_mass+rocket1.mass_empty+rocket1.payload_mass)*(3*(rocket1.diameter/2)**2 + rocket1.height**2)/12 +
                   (rocket1.mass_empty+rocket1.payload_mass + rocket1.prop_mass)*(rocket1.tot_cgpos[0]**2 + rocket1.tot_cgpos[2]**2),0.],
                  [0.,0.,0.5*(rocket1.prop_mass+rocket1.mass_empty+rocket1.payload_mass)*(rocket1.diameter/2)**2 +
                   (rocket1.mass_empty+rocket1.payload_mass + rocket1.prop_mass)*(rocket1.tot_cgpos[1]**2 + rocket1.tot_cgpos[0]**2)]])

    moment_from_gimbal_engine = engine1.engine2bodygimbal(rocket1.tot_cgpos,pitch_mat,roll_mat,yaw_mat)
    #print(moment_from_gimbal_engine)

    angular_accel = rocket1.angular_accel_mat(moi_mat,moment_from_gimbal_engine) #CHECK THIS

    angular_vel += angular_accel*dt
    total_gimbal +=angular_vel*dt

    tot_roll_mat = matrices.roll_matrix(total_gimbal[0])
    tot_pitch_mat = matrices.pitch_matrix(total_gimbal[1])
    tot_yaw_mat = matrices.yaw_matrix(total_gimbal[2])

    #fuel mass change


    thrust1 = engine1.mdot * ve
    thrust1 = np.matmul(yaw_mat,np.matmul(pitch_mat,np.matmul(roll_mat,thrust1)))
    thrust2 = engine2.mdot*ve
    rocket1.prop_mass=rocket1.prop_mass - dt*engine1.mdot-dt*engine2.mdot*8
    if rocket1.prop_mass <=0:
        engine1.mdot = 0
        engine2.mdot = 0
        rocket1.prop_mass = 0
    #Matrix stuff that was in the classes but now i gotta take it out and into this loop fml

















    #Drag force change
    if radius<earth.radius:
        drag_force = np.array([0.,0.,0.])
        rocket_vel = np.array([0.,0.,0.])
    else:
        if np.linalg.norm(rocket_vel) != 0.:
            drag_force = -0.5*density*(np.linalg.norm(rocket_vel-np.array([0.,0.,0.]))**2)*rocket1.area*rocket1.cD*rocket_vel/(np.linalg.norm(rocket_vel))
            #fix something with the very last bit of the drag equation because its giving the wrong drag vectors i think invert the lat long stuff

        #    print(rocket_vel,-np.array([0.,465.1*cos(latitude),465.1*sin(latitude)]))
        else:
            drag_force = np.array([0,0,0])

    # if np.linalg.norm(rocket_vel) == 0:
    #     norm_vel = 1
    # else:
    #     norm_vel = np.linalg.norm(rocket_vel)
    #summation of all forces acting on the rocket

    rotated_thrust = np.matmul(tot_yaw_mat,np.matmul(tot_pitch_mat,np.matmul(tot_roll_mat,thrust2*8+thrust1)))
    long_yaw_mat = matrices.yaw_matrix(longitude)
    lat_roll_mat = matrices.roll_matrix(latitude)
    lat_pitch_mat = matrices.pitch_matrix(-pi/2)
    total_force =(rocket1.mass_empty+rocket1.prop_mass)*(engine1.pos/radius)*-9.80665+rotated_thrust+drag_force #Fix coordinate stuff



    rocket_accel = total_force/(rocket1.mass_empty+rocket1.prop_mass)

    rocket_vel +=rocket_accel*dt
    engine1.pos +=(np.array([465.1*cos(latitude),465.1*sin(latitude),0.])+rocket_vel)*dt
    temps.append(temp)

    t+=dt
    #change in temp due to atmosphere
    nose_radius = 0.7 #don't know if this is true but eh close enough for F9
    if radius<= 6371e3 + earth.layers[-1]*10**3:
        #added temp per dt
        fric_temp = 1.83*10**(-4)*25*np.linalg.norm(rocket_vel)**3*sqrt(density/nose_radius)/rocket1.specific_heat/(rocket1.mass_empty+rocket1.prop_mass) #the 25 sqm is the SA of the nose cone exposed to the air friction
        skin_temp+=fric_temp
    else:
        if skin_temp >=2:
            skin_temp-=  0.1 #honestly have no clue if this is how fast heat dissipates

    tot_skin_temp = temp + skin_temp


    if skin_temp >= 1923.15:  #This reentry temp is totally made up based on space shuttle values
        running = False
        print('Burnt up on reentry or ascent')

    if t > 7000 or radius<earth.radius:
        running = False
        break


    plotter_pos = np.matmul(long_yaw_mat,np.matmul(lat_pitch_mat,np.matmul(lat_roll_mat, (engine1.pos)))) #

    #plotting section
    zpos.append(plotter_pos[2])
    ypos.append(plotter_pos[1])
    xpos.append((plotter_pos[0]))
    ttlforce.append((np.linalg.norm(total_force)))
    skin_temps.append(tot_skin_temp)
    dragforce.append(np.linalg.norm(drag_force))
    time1.append(t)
    time.append(t)
    moments.append(moment_from_gimbal_engine[0])
    densitys.append(density)
    temps.append(temp)
    presss.append(press)
    rads.append(radius-6371e3)
    vels.append(np.linalg.norm(rocket_vel))
    cg_vals.append(rocket1.tot_cgpos[2])
mpl.use('Qt5Agg')
fig = plt.figure()
# Make data
ax = fig.gca(projection='3d')
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 6371e3 * np.outer(np.cos(u), np.sin(v))
y = 6371e3 * np.outer(np.sin(u), np.sin(v))
z = 6371e3 * np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface
ax.plot_surface(x, y, z, color='b')

ax.plot(xpos,ypos,zpos,label = 'position',color = 'r')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()
plt.title('ttlforce')
plt.plot(time,rads)
plt.show()
fig, ax = plt.subplots()
plt.title('Y-Z Plot')
circle1 = plt.Circle((0,0),6371e3,color='b')
ax.add_artist(circle1)
plt.plot(ypos,zpos,color = 'r')
plt.xlabel('Y')
plt.ylabel('Z')
plt.show()

fig1, ax1 = plt.subplots()
plt.title('X-Z Plot')
circle2 = plt.Circle((0,0),6371e3, color='b')

ax1.add_artist(circle2)
plt.plot(xpos,zpos,color = 'r')
plt.xlabel('X')
plt.ylabel('Z')
plt.show()

fig2, ax2 = plt.subplots()
plt.title('X-Y Plot')
circle3 = plt.Circle((0,0),6371e3, color='b')

ax2.add_artist(circle3)
plt.plot(xpos,ypos,color = 'r')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

plt.title('moments')
plt.plot(time,moments)
#print(max(rads))
plt.show()
plt.title('Drag Force')
plt.plot(time,dragforce,color = 'r')
plt.show()
plt.title('Skin Temp')
plt.plot(time,skin_temps)
plt.show()

# plt.title('Density')
# plt.plot(time,densitys)
# plt.show()
# plt.title('temp')
# plt.plot(time,temps)
# plt.show()




#
# plt.show()




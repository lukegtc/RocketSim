import vpython
from rocket_pos_updater import *
from Stage_2_Calc import xpos_s2,ypos_s2,zpos_s2#,check
earth_location = vpython.vector(0,0,0)
rocket_location = vpython.vector(xpos[0],ypos[0],zpos[0])
vpython.sphere(pos = earth_location, radius= earth.radius, color = vpython.color.blue)
rocket = vpython.cylinder(pos = rocket_location,axis = vpython.vector(0,0,40.9000),radius = stage1.radius,color = vpython.color.red)
vpython.cylinder(pos=vpython.vector(xpos1[0],ypos1[0],zpos1[0]), axis = vpython.vector(1,0,0),size = vpython.vector(1000,50,50),color = vpython.color.orange)
vpython.cylinder(pos=vpython.vector(xpos1[0],ypos1[0],zpos1[0]), axis = vpython.vector(0,1,0),size = vpython.vector(1000,50,50),color = vpython.color.purple)
vpython.cylinder(pos=vpython.vector(xpos1[0],ypos1[0],zpos1[0]), axis = vpython.vector(0,0,1),size = vpython.vector(1000,50,50),color = vpython.color.white)
xpos1 = xpos1+xpos_s2
ypos1 = ypos1+ypos_s2
zpos1 =zpos1+ zpos_s2
i = 0
running1 = True
vpython.scene.camera.follow(rocket)

while running1:
    vpython.rate(100)
    #print(xpos[i],ypos[i],zpos[i])
    rocket.pos.x = xpos1[i]
    rocket.pos.y = ypos1[i]
    rocket.pos.z = zpos1[i]
    #vpython.cylinder(pos = vpython.vector(xpos[i],ypos[i],zpos[i]),axis = vpython.vector(vels2[i][0],vels2[i][1],vels2[i][2]),size = vpython.vector(40.9,rocket1.radius*10,rocket1.radius*10),color = vpython.color.white)
    #vpython.cylinder(pos = vpython.vector(xpos[i],ypos[i],zpos[i]),size = vpython.vector(40.9,rocket1.radius*10,rocket1.radius*100),color = vpython.color.white)


    color1 = vpython.color.green
    if i <len(xpos):
        vpython.sphere(pos=vpython.vector(xpos[i], ypos[i], zpos[i]), radius=stage1.radius * 100,
                       color=color1)
    if stage1.prop_mass<=0. and i <= len(xpos1):
        vpython.sphere(pos=vpython.vector(xpos1[i], ypos1[i], zpos1[i]), radius=stage1.radius * 100,color=vpython.color.orange)
    #print(vpython.scene.camera.pos)
    i+=1
    if i == len(xpos) and check == 0:
        running1 = False
    if i == len(xpos1) and check ==1:
        running1 = False

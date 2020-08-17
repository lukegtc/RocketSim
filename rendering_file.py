import vpython
from rocket_pos_updater import *
earth_location = vpython.vector(0,0,0)
rocket_location = vpython.vector(xpos[0],ypos[0],zpos[0])
vpython.sphere(pos = earth_location, radius= earth.radius, color = vpython.color.blue)
rocket = vpython.cylinder(pos = rocket_location,axis = vpython.vector(0,0,40.9000),radius = rocket1.radius,color = vpython.color.red)
i = 0
running1 = True
vpython.scene.camera.follow(rocket)
while running1:
    vpython.rate(10)
    #print(xpos[i],ypos[i],zpos[i])
    rocket.pos.x = xpos[i]
    rocket.pos.y = ypos[i]
    rocket.pos.z = zpos[i]
    vpython.cylinder(pos = vpython.vector(xpos[i],ypos[i],zpos[i]),axis = vpython.vector(vels2[i][0],vels2[i][1],vels2[i][2]),size = vpython.vector(40.9,rocket1.radius*10,rocket1.radius*10),color = vpython.color.white)


    print(vpython.scene.camera.pos)
    i+=1
    if i == len(points):
        running1 = False
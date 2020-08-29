# RocketSim
This code allows for the 3-D plotting of a rocket's trajectory if it is given certain parameters such as the exit velocity of the engines used, the mass flow of the engines used,
the amount of fuel in the 1st and 2nd stages and the final payload that must be put into orbit.
  This code also allows one to set a gimballing rate that could be executed at a specific moment in order to perform a gravity turn, for example. This code also makes use of vpython
to three dimensionally plot the trajectory of the 1st ande 2nd stages. Overall, there are still some changes required to make the code more user friendly, such as consolidating all of the
inputs one could make into one file as opposed to being spread out throughout the entire program.
The class_file.py contains all of the relevant classes.
The rocket_pos_updater.py file contains the initialization of all of the classes and where the gimballing can be input. This folder also creates all of the plots.
The Stage_2_calc.py file calculates the position of the 2nd stage once it is detached from the 1st stage.
The atmosphere.py file contains the simulated atmosphere.
The rendering_file.py file renders the trajectory of the 1st stage and the 2nd stage and plots them in vpython (you have to zoom in in order to see the trajectories as the blue sphere rendered is the Earth.)


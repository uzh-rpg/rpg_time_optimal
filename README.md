# Time-Optimal Planning for Quadrotor Waypoint Flight
This is an **example implementation** of the paper "Time-Optimal Planning for Quadrotor Waypoint Flight" published in Science Robotics the "Robotics and Perception Group" at the Dep. of Informatics (University of Zurich), and Dep. of Neuroinformatics (ETH and University of Zurich).

Make sure you check out the [video](https://www.youtube.com/watch?v=ZPI8U1uSJUs) below, read the [paper](http://rpg.ifi.uzh.ch/docs/ScienceRobotics21_Foehn.pdf), the [instructions](#instructions), and cite us accordingly.


# Video
[![Time-Optimal-Planning](https://img.youtube.com/vi/ZPI8U1uSJUs/0.jpg)](https://www.youtube.com/watch?v=ZPI8U1uSJUs)

# Paper
Find the Science Robotics article [here](https://robotics.sciencemag.org/content/6/56/eabh1221) or read [our full version](http://rpg.ifi.uzh.ch/docs/ScienceRobotics21_Foehn.pdf) of the paper and cite us accordingly:

**Time-Optimal Planning for Quadrotor Waypoint Flight**  
Philipp Foehn, Angel Romero, Davide Scaramuzza  
2021, Science Robotics, Volume 6, Issue 56  
DOI: 10.1126/scirobotics.abh1221  
  
  
  
Bibtex:
```
@article {foehn2021CPC,
	author = {Foehn, Philipp and Romero, Angel and Scaramuzza, Davide},
	title = {Time-Optimal Planning for Quadrotor Waypoint Flight},
	volume = {6},
	number = {56},
	elocation-id = {eabh1221},
	year = {2021},
	doi = {10.1126/scirobotics.abh1221},
	publisher = {Science Robotics},
	URL = {https://robotics.sciencemag.org/content/6/56/eabh1221},
	eprint = {https://robotics.sciencemag.org/content/6/56/eabh1221.full.pdf},
	journal = {Science Robotics}
}
```

# Instructions
1. Make sure you have Python3 running.
4. Clone this repository with `git clone git@github.com:uzh-rpg/rpg_time_optimal.git`.
5. Navigate into the root folder of the clone repository `cd rpg_time_optimal`.
2. Install the requirements `pip3 install -r requirements.txt`
3. Download CasADi from [the official website](https://web.casadi.org) or with `pip install casadi`.
6. Run the example with `python3 example/optimization.py`.
7. Show some sparkly plots with `python3 example/plotting.py`.

This will create output `.csv` files with the trajectory and some `.pdf` visualizing the result.

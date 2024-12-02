# Distributed Control of Multi-Agent Quadrupedal Robots with Real-Time Nonlinear MPC and Collision Avoidance
- Distributed and Nonlinear Model Predictive Control using Reference Models for MultiAgent A1 Robots (Lab Website: https://www.kavehakbarihamed.com/)
- Uses SNOPT as a nonlinear solver for NMPC

## Code Organization
The bulk of the code and computations is run through the class "LocoWrapper.cpp". However, the main function for initiating RaiSim and the simulations is in src/robots/A1_Sim.cpp. Note that other files in the robots folder may be outdated and are largely used for troublshooting.

## Where to git clone
Clone this repository under raisim_workspace

## Prerequisite
- In bash, the following must be defined
```sh
export WORKSPACE=${HOME}/RAISIM_WORKSPACE
export LOCAL_INSTALL=${HOME}/RAISIM_WORKSPACE/build
```
- Be sure to add your license to the /rsc folder, and change the license name in src/robots/A1_Sim.cpp


## Building and Running Code
### Build Code From Scratch
```sh
cd $WORKSPACE
cd A1_Robot && mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$LOCAL_INSTALL
make -j8
```

### Build and Run Code Using Script
```sh
cd $WORKSPACE/Distributed_NMPC_A1
./run_Raisim.sh
```

### Run Example from Command Line
```sh
cd $WORKSPACE
cd A1_Robot/build
./A1_Sim 
```

Finally, there are several different checkerboard colors provided in the material file (rsc/material/myMaterials.material) and can be changed in the A1_Sim.cpp file. To generate a different cherckerboard color pattern, use the MATLAB script (rsc/material/ColoredChecker.m) and create a new material in myMaterials.material following the same convention as the other checkerboard materials.

# Citation
This work has been featured as  in ASME Journal of Dynamic Systems, Measurements, and Control in November 2024: 
[Imran, B. M., Fawcett, R. T., Kim, J., Leonessa, A., and Hamed, K. A. (October 23, 2024). "A Distributed Layered Planning and Control Algorithm for Teams of Quadrupedal Robots: An Obstacle-Aware Nonlinear Model Predictive Control Approach." ASME. J. Dyn. Sys., Meas., Control. May 2025; 147(3): 031006. https://doi.org/10.1115/1.4066632](https://asmedigitalcollection.asme.org/dynamicsystems/article-abstract/147/3/031006/1206862/A-Distributed-Layered-Planning-and-Control?redirectedFrom=fulltext)

If you benefit from the code or work, please remember to cite us:
```
@article{10.1115/1.4066632,
    author = {Imran, Basit Muhammad and Fawcett, Randall T. and Kim, Jeeseop and Leonessa, Alexander and Hamed, Kaveh Akbari},
    title = {A Distributed Layered Planning and Control Algorithm for Teams of Quadrupedal Robots: An Obstacle-Aware Nonlinear Model Predictive Control Approach},
    journal = {Journal of Dynamic Systems, Measurement, and Control},
    volume = {147},
    number = {3},
    pages = {031006},
    year = {2024},
    month = {10},
    issn = {0022-0434},
    doi = {10.1115/1.4066632},
    url = {https://doi.org/10.1115/1.4066632},
    eprint = {https://asmedigitalcollection.asme.org/dynamicsystems/article-pdf/147/3/031006/7390960/ds\_147\_03\_031006.pdf},
}
```

# YouTube video
The YouTube link for the experiments and simulations is [here](https://youtu.be/hwEhA7JCXAU?feature=shared).


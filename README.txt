Sat-Tracker v1.0

Python-based Interactive Program to track a satellite, and 
check if you will be able to see the satellite in sky.

By Team Celestial Observers
Team Members (all from 1st Year):
    Bhuvan Narula (author)
    Ayush
    Suhana

https://github.com/bhuvannarula/Sat-Tracker

*PLEASE READ ALL NOTES BELOW TILL END BEFORE RUNNING PROGRAM*

Requirements:
- Anaconda
- Anaconda Library
    basemap
- Python Libraries
    astropy
    sgp4
    matplotlib==3.2
    itertools (build-in)
- tkinter (comes preinstalled)
Note: it is recommended to make a private environment to install packages and run Program
Program has been tested SUCCESSFULLY on Ubuntu 18.04. Windows and MacOS should be supported as well if anaconda is used.


Virtual Environment Setup (optional):
- Install Anaconda from Anaconda website
- Create a virtual environment
    conda create --name {env_name}
    - type 'y' and press enter
- Use the Environment
    conda activate {env_name}

Installing Required Libraries
- Install basemap
    conda install basemap
- Install Python Libraries
    python3 -m pip install sgp4 astropy
    python3 -m pip install -U matplotlib==3.2


Pre-Run Setup:
'satellite.txt' contains some sample satellite orbit elements. 
The satellite orbit data in first 3 lines will be used. (so, STARLINK-2409 will be tracked if file is not modified)
To change the satellite to be tracked, add the satellite orbit data to the first 3 lines.
Format of Satellite Orbit Data (TLE format) is as below:
{Satellite Name}
{Line 1}
{Line 2}
Any data after the first 3 lines will be ignored.
You can refer to celestrak.com to find orbit data


How to Run:
0. (optional) Activate Virtual Environment by 'conda activate {env_name}' (if created).
1. Execute sat-tracker.py file through terminal. GUI should open
2. Then operate program GUI

IMPORTANT NOTES
- Program is slow and unoptimised. It will take a minute or two for calculations.
- Effect of any Pollution on Visibility is not considerd.

Notes on GUI:
- When 'Track Satellite' button is clicked, current window of program will close and a new window will open.
DO NOT PANIC OR DEDUCT MARKS thinking that program crashed. Give it a minute or two to respond. Please
- Program may seem to freeze. This is NORMAL. Wait for 2 min. It will respond.
- No 'Back' button has been implemented due to shortage of time (tkinter is complicated). If more time is provided, it can be implemented.
So, to go back to main menu, please QUIT and RELAUNCH the program. main Menu should open.
- 4th Point of the Rulebook Requirements was unable to make it to final code due to lack of time.

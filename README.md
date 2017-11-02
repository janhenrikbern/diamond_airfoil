# diamond_airfoil
Calculates and plots pressure distribution across an diamond shaped airfoil at specified free stream pressure and mach, wedge angle, and angle of attack. Furthermore, it saves the property changes in a .csv-file. 

Example output:
![Airfoil Plot Example](https://raw.githubusercontent.com/janhenrikbern/diamond_airfoil/master/example.png)

## System Requirements
- Python 3.6
- Anaconda (for Python 3.6)

**Or**
- numpy
- matplotlib
- seaborn
- scipy

## Installation
In order for the program to run, `aerodynamics.py` and `diamond_airfoil.py` have to be in the same folder. The easiest way to do that is to clone the repositiory. 

## How to use it
The only script that needs to be called is `diamond_airfoil.py`. 

`$ python3 diamond_airfoil.py 3 5 -aoa 6 -ipres 1.5`

It's command line input arguments are: 
* Mach number of the incoming flow. 
* Half-wedge angle of the incoming flow in degrees. 
* (Optional) Angle of attach of flow to airfoil in degrees. The default is 0 degrees.
* (Optional) The freestream pressure of the incoming flow in atm. The default is 1 atm.  

For additional documentation use `-h`.
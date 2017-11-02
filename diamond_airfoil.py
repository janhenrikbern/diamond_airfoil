#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 15 12:06:03 2017

@author: Jan Bernhard
"""
from aerodynamics import Compression, Expansion
import math, csv
import numpy as np
from argparse import ArgumentParser
from matplotlib.transforms import Affine2D
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import seaborn; seaborn.set()

def input_argument():
	# Read commandline input. 
	parser = ArgumentParser()
	parser.add_argument("mach_input", 
		help = 'Mach_Input: Mach number of the freestream flow.', 
		type = float)
	parser.add_argument("wedge_angle", 
		help = 'Wedge_Angle [Degrees]: Angle from horizontal to airfoil boundary.',
		type = float)
	parser.add_argument("-aoa","--angle_of_attack", 
		help = 'Angle_Of_Attack [Degrees]: Angle from horizontal to airfoil cordline.',
		type = float)
	parser.add_argument("-ipres","--initial_pressure", 
		help = 'Initial_Pressure [atm]: Pressure of freestream flow.',
		type = float)

	args = parser.parse_args()
	base = [args.mach_input,args.wedge_angle]
	base.append(args.angle_of_attack) if args.angle_of_attack is not None else base.append(0.)
	base.append(args.initial_pressure) if args.initial_pressure is not None else base.append(1.)
	
	return base


def main(Mach_input, wedge_angle, angle_of_attack = 0, initial_pressure = 1):

	if angle_of_attack != 0:
		theta = [wedge_angle - angle_of_attack, wedge_angle + angle_of_attack]
	else:
		theta = [wedge_angle]

	regions = []
	for i,th in enumerate(theta):
	# Array format: region_ 1 = [Mach_1,rho_ratio_1,pressure_ratio_1,temp_ratio_1]
		region_1 = Compression(Mach_input,th).getValues()
		regions.append(region_1)
		region_2 = Expansion(region_1[0],2*wedge_angle).getValues()
		regions.append(region_2)
		if i < 1:
			region_3 = Compression(region_2[0],th+2*wedge_angle).getValues()
		else:
			region_3 = Compression(region_2[0],th-2*wedge_angle).getValues()
		regions.append(region_3)

	# Plot Airfoil cross-section at angle of attack
	fig = plt.figure()
	ax = fig.add_subplot(212)

	y_coor = .5*math.tan(math.radians(wedge_angle))
	verts = np.array([[0,0],[.5,y_coor],[1.0,0],[.5,-y_coor],[0,0]])
	diamond_airfoil = Polygon(verts,True)
	transform = Affine2D().rotate_deg(-angle_of_attack) + ax.transData
	diamond_airfoil.set_transform(transform)
	ax.add_patch(diamond_airfoil)

	plt.xlim([-.3,1.3])
	plt.ylim([-1,1])
	plt.xlabel('Normalized positon along cord [x/c]')
	plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
	
	# Plot pressure distribution along airfoil surfaces
	ax2 = fig.add_subplot(211)

	# Calculate the exact pressure
	pressure_ratios = []
	if len(regions) < 4:
		[pressure_ratios.append(regions[i][2]) for i in range(2)]
		[pressure_ratios.append(regions[i][2]) for i in range(2)]
	else:
		[pressure_ratios.append(regions[i][2]) for i in range(2)]
		[pressure_ratios.append(regions[i][2]) for i in (3,4)]

	p_distribution = []
	pressure = []
	for i,p in enumerate(pressure_ratios):
		if i in (0,2):
			pressure.append(initial_pressure * p)
		else:
			pressure.append(pressure[i-1] * p)

	# Create plot line
	p_real = [(p-initial_pressure)/initial_pressure for p in pressure]

	# add free-stream region to start and end of plot
	p_plot1 = [0 for i in range (30)]
	p_plot2 = p_plot1.copy()
	for i,p in enumerate(p_real):
		if i < 2:
			p_plot1.extend([p for i in range(50)])
		else:
			p_plot2.extend([p for i in range(50)])

	tail = [0 for i in range (30)]
	p_plot1.extend(tail)
	p_plot2.extend(tail)
	x_axis = [0.01*i - .3 for i in range(160)]
	plt1, plt2 = plt.plot(x_axis,p_plot1,'g--',x_axis,p_plot2,'b:')
	plt.ylabel('(P - P_i) / P_i')
	plt.xlim([-.3,1.3])
	plt.legend([plt1,plt2],['Along Top of Airfoil','Along Bottom of Airfoil'])

	plt.show()

	# Save properties to .csv
	with open("diamond_airfoil_at_%s_AoA.csv"%(angle_of_attack),'w') as file:
		writer = csv.writer(file)
		writer.writerow(['Mach','density_ratio','pressure_ratio','temperature_ratio','actual pressure [atm]'])
		counter = 0
		for i,properties in enumerate(regions):
			if i not in (0,1,3,4):
				continue
			temp = [regions[i][0],regions[i][1],regions[i][2],regions[i][3],pressure[counter]]
			temp = map(lambda x: round(x,4),temp)
			counter = counter + 1
			writer.writerow(temp)
		


if __name__ == '__main__':
	Mach,Wedge,AoA,Pressure = input_argument()
	main(Mach,Wedge,AoA,Pressure)

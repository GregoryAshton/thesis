import pylab as py 
import numpy as np 
from mpl_toolkits.mplot3d import Axes3D
import sys
from math import pi

# Set the defaults
py.rcParams['axes.color_cycle'] = ['k', 'r', 'cyan']
py.rcParams['font.size'] = 16
#py.rcParams['lines.linewidth'] = 2
py.rcParams['axes.grid']=True
py.subplots_adjust(left=0.13, right=0.9, top=0.9, bottom=0.12,wspace=0.0)

sys.path.append("/home/ga7g08/nsmod/")
sys.path.append("/home/greg/nsmod/")
from nsmod.Useful_Tools import ThreeD_Sphere



if "bi" in sys.argv:
	I1 = 1.0
	I2 = 1.0 
	I3 = 1.8

	#fig = py.figure()
	#ax = fig.add_subplot(111, projection='3d')

	n=3000

	E=2.0
	#p=7
	#M2_list=py.linspace(2.0*E*I1,2.0*E*I3,p)
	#M2_list = [0.5*(2*E*I2+2.0*E*I1)] #,2*E*I2,0.5*(2.0*E*I3+2*E*I2)]
	#M2_list = [0.45*(2*E*I2+2.0*E*I1),0.5*(2*E*I2+2.0*E*I1),0.55*(2*E*I2+2.0*E*I1)] # low
	#M2_list = [2*E*I2]
	M2_list = [0.39*(2.0*E*I3+2*E*I2),0.5*(2.0*E*I3+2*E*I2),0.62*(2.0*E*I3+2*E*I2)]

	colors=["r","g","b","r","y","purple"]

	ax=py.subplot(111, projection='3d')
	elevation = 15.0
	azimuth = -40.0
	ax.view_init(elevation,azimuth) 
	# Draw an ellipse
	for i in range(len(M2_list)):

		M=py.sqrt(M2_list[i])
		print i,"m=", M
		m1_list = py.linspace(-M,M,n)

		x=[] ; y_pos=[] ; z_pos=[] ;y_neg=[] ; z_neg=[]
		for m1 in m1_list:
			m2_pos = py.sqrt(
					(2*E-pow(m1,2)/I1 - (pow(M,2) - pow(m1,2))/I3)/(pow(I2,-1)-pow(I3,-1)))

			m3_pos = py.sqrt(
					(2*E-pow(m1,2)/I1 - (pow(M,2) - pow(m1,2))/I2)/(pow(I3,-1)-pow(I2,-1)))
			m2_neg = -py.sqrt(
					(2*E-pow(m1,2)/I1 - (pow(M,2) - pow(m1,2))/I3)/(pow(I2,-1)-pow(I3,-1)))

			m3_neg = -py.sqrt(
					(2*E-pow(m1,2)/I1 - (pow(M,2) - pow(m1,2))/I2)/(pow(I3,-1)-pow(I2,-1)))

			x.append(m1) 
			y_pos.append(m2_pos)
			z_pos.append(m3_pos)
			y_neg.append(m2_neg)
			z_neg.append(m3_neg)
		
		#Plot the data
		ThreeD_Sphere(ax,elevation,azimuth,x,y_pos,z_pos,ls="-",lw=2)
		ThreeD_Sphere(ax,elevation,azimuth,x,y_pos,z_neg,ls="-",lw=2)
		ThreeD_Sphere(ax,elevation,azimuth,x,y_neg,z_pos,ls="-",lw=2)
		ThreeD_Sphere(ax,elevation,azimuth,x,y_neg,z_neg,ls="-",lw=2)		
	
	if "sphere" in sys.argv:
		# Draw a sphere
		u_offset = 0.0*pi
		v_offset = 0.1*pi
		u = np.linspace(0+u_offset, 2 * np.pi + u_offset, 100)
		v = np.linspace(0+v_offset, np.pi + v_offset, 100)

		M = py.sqrt(M2_list[2])
		x = M * np.outer(np.cos(u), np.sin(v))
		y = M * np.outer(np.sin(u), np.sin(v))
		z = M * np.outer(np.ones(np.size(u)), np.cos(v))
		ax.plot_surface(x, y, z,  rstride=7, cstride=9, color='green',alpha=0.9,lw=0.2,shade=True)
		ax.set_axis_off()


	ax.set_aspect('equal',"box")

	mx=3.8
	ax.plot(py.linspace(-mx,mx,n),py.zeros(n),py.zeros(n),color="k",ls="-",lw=1.2)
	ax.text(mx*1.2,0.05,0.05,r"$\hat{y}$")
	ax.plot(py.zeros(n),py.linspace(-mx,mx,n),py.zeros(n),color="k",ls="-",lw=1.2)
	ax.text(0.05,-mx*1.2,0.05,r"$\hat{x}$")
	ax.plot(py.zeros(n),py.zeros(n),py.linspace(-mx,mx,n),color="k",ls="-",lw=1.2)
	ax.text(0.0,0.0,mx*1.1,r"$\hat{z}$")




	# Draw an ellipse
	u = np.linspace(0, 2 * np.pi, 100)
	v = np.linspace(0, np.pi, 100)

	x = py.sqrt(2*E*I1) * np.outer(np.cos(u), np.sin(v))
	y = py.sqrt(2*E*I2) * np.outer(np.sin(u), np.sin(v))
	z = py.sqrt(2*E*I3) * np.outer(np.ones(np.size(u)), np.cos(v))
	ax.plot_surface(x, y, z,  rstride=8, cstride=9, color='white',alpha=0.7,lw=0.1,shade=True)
	ax.set_axis_off()
	#ax.set_xlim(0,1.5)
	#ax.set_ylim(-1.5,1.5)
	#ax.set_zlim(-1.5,1.5)


	py.show()

if "tri" in sys.argv:
	# Viewing plotform and init. axis
	elevation = 15.0
	azimuth = 50.0
	ax=py.subplot(111, projection='3d')
	ax.set_aspect('equal',"box")


	# Define some arbitary numbers for the moment of inertia
	I1 = 0.8
	I2 = 1.0 
	I3 = 1.5

	n=5000 # Number of points used

	E=2.0 # Energy 
	#M2_list=py.linspace(2.0*E*I1,2.0*E*I3,p)
	#M2_list = [0.5*(2*E*I2+2.0*E*I1)] #,2*E*I2,0.5*(2.0*E*I3+2*E*I2)]
	M2_list = [0.45*(2*E*I2+2.0*E*I1),0.5*(2*E*I2+2.0*E*I1),0.55*(2*E*I2+2.0*E*I1)] # low
	#M2_list = [2*E*I2]
	M2_list = [0.42*(2.0*E*I3+2*E*I2),0.5*(2.0*E*I3+2*E*I2),0.55*(2.0*E*I3+2*E*I2)]


	# Draw the intersections
	for i in range(len(M2_list)):

		M=py.sqrt(M2_list[i])
		m1_list = py.linspace(-M,M,n)
		x=[] ; y_pos=[] ; z_pos=[] ;y_neg=[] ; z_neg=[]
		for m1 in m1_list:
			m2_pos = py.sqrt(
					(2*E-pow(m1,2)/I1 - (pow(M,2) - pow(m1,2))/I3)/(pow(I2,-1)-pow(I3,-1)))

			m3_pos = py.sqrt(
					(2*E-pow(m1,2)/I1 - (pow(M,2) - pow(m1,2))/I2)/(pow(I3,-1)-pow(I2,-1)))
			m2_neg = -py.sqrt(
					(2*E-pow(m1,2)/I1 - (pow(M,2) - pow(m1,2))/I3)/(pow(I2,-1)-pow(I3,-1)))

			m3_neg = -py.sqrt(
					(2*E-pow(m1,2)/I1 - (pow(M,2) - pow(m1,2))/I2)/(pow(I3,-1)-pow(I2,-1)))

			x.append(m1) 
			y_pos.append(m2_pos)
			z_pos.append(m3_pos)
			y_neg.append(m2_neg)
			z_neg.append(m3_neg)
		
		# Plot the data
		ThreeD_Sphere(ax,elevation,azimuth,x,y_pos,z_pos,ls="-",lw=2)
		ThreeD_Sphere(ax,elevation,azimuth,x,y_pos,z_neg,ls="-",lw=2)
		ThreeD_Sphere(ax,elevation,azimuth,x,y_neg,z_pos,ls="-",lw=2)
		ThreeD_Sphere(ax,elevation,azimuth,x,y_neg,z_neg,ls="-",lw=2)	


	# Plot and label the axis
	mx=3.4
	C="k" ; LW=0.8
	ax.plot(py.linspace(-mx,mx,n),py.zeros(n),py.zeros(n),color=C,alpha=0.8,lw=LW)
	ax.text(mx*1.2,0.05,0.05,r"$e_{1}$")
	ax.plot(py.zeros(n),py.linspace(-mx,mx,n),py.zeros(n),color=C,alpha=0.8,lw=LW)
	ax.text(0.05,mx*1.2,0.05,r"$e_{2}$")
	ax.plot(py.zeros(n),py.zeros(n),py.linspace(-mx,mx,n),color=C,alpha=0.8,lw=LW)
	ax.text(0.05,0.05,mx*1.2,r"$e_{3}$")



	# Draw an ellipsoid
	u = np.linspace(0, 2 * np.pi, 100)
	v = np.linspace(0, np.pi, 100)

	x = py.sqrt(2*E*I1) * np.outer(np.cos(u), np.sin(v))
	y = py.sqrt(2*E*I2) * np.outer(np.sin(u), np.sin(v))
	z = py.sqrt(2*E*I3) * np.outer(np.ones(np.size(u)), np.cos(v))
	ax.plot_surface(x, y, z, rstride=7, cstride=9, color='white',alpha=0.7,lw=0.1,shade=True)
	ax.set_axis_off()
	py.show()

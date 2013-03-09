#!/usr/bin/env python
#
# ========== Py_Astro_Utils.PyAstroCoords ===========
#
# Author: wujinnnnn@qq.com 
# Purpose: Provide basic functions to manipulate astronomical coordinate systems.
#
#
# PyAstroCoords is free software; you can redistribute it and/or modify it under the terms 
# of the GNU General Public License as published by the Free Software Foundation; 
# either version 2 of the License, or (at your option) any later version.
# PyAstroCoords is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
#
#



import numpy as np
from copy import deepcopy
import time
import pdb

def rect2Sphere(*args):
	'''
	rectangular coordinates to spherical coordinates
	x = r*cos(theta)*cos(fai)
	y = r*sin(theta)*cos(fai)
	z = r*sin(fai)
	output angle in degree
	
	examples:
		for single coordinates 
			input: np.array([x,y,z])
			
			rect2Sphere(np.array([1,1,1])) 
			output = (1.7320508075688772, 45.0, 35.264389682754668)
			it is (r, theta, fai)

		for array coordinates
			input: np.array([x1,x2,x3.....]), np.array([y1,y2,y3.....]), np.array([z1,z2,z3.......])
			
			rect2Sphere(np.array([1,-2,1]),np.array([1,-2,0]),np.array([1,-2,0]))
			output = 
			(array([ 1.73205081,  3.46410162,  1.        ]),
			 array([  45.,  225.,    0.]),
			 array([ 35.26438968, -35.26438968,   0.        ]))
			they are 
			(np.array([r1,r2,r3....]),
			 np.array([theta1,theta2,theta3....]),
			 np.array([fai1,fai2,fai3....]))
	'''
	if len(args)==1:
		x=args[0][0]
		y=args[0][1]
		z=args[0][2]
	elif len(args)==3:
		x=args[0]
		y=args[1]
		z=args[2]

	r=np.sqrt(np.power(x,2)+np.power(y,2)+np.power(z,2))
	xx = x/r
	yy = y/r
	zz = z/r
	alpha = np.arcsin(zz)
	delta = np.arctan(yy/xx)
	delta = delta + (xx<0)*np.pi
	delta = delta + (np.logical_and(xx>0,yy<0))*np.pi*2
	return r,rad2Dec(delta),rad2Dec(alpha)

def sphere2Rect(*args):
	'''
	spherical coordinates to rectangular coordinates
	x = r*cos(theta)*cos(fai)
	y = r*sin(theta)*cos(fai)
	z = r*sin(fai)
	input angle in degree

	examples:
		for single coordinates
			input: np.array([r, theta, fai])
			   or  np.mat([r, theta, fai])
		
			sphere2Rect(np.array([1,45,0]))
			output = (0.70710678118654757, 0.70710678118654746, 0.0)
			it is (x, y, z)

		for array coordinates
			input: np.array([r1,r2,r3....]), np.array([theta1, theta2, theta3....]), np.array([fai1, fai2, fai3.....]) or the same np.mat
		
			sphere2Rect(np.array([1,2,3]), np.array([0,45,90]), np.array([30,0,-30]))
			output = 
			(array([  8.66025404e-01,   1.41421356e+00,   1.59081032e-16]),
			 array([ 0.        ,  1.41421356,  2.59807621]),
			 array([ 0.5,  0. , -1.5]))
			they are
			(np.array([x1,x2,x3]),
			 np.array([y1,y2,y3]),
			 np.array([z1,z2,z3]))
	'''
	if len(args)==1:
		r=args[0][0]
		theta=dec2Rad(args[0][1])
		fai=dec2Rad(args[0][2])
	elif len(args)==3:
		r=args[0]
		theta=dec2Rad(args[1])
		fai=dec2Rad(args[2])

	x = np.multiply(np.multiply(r, np.cos(theta)), np.cos(fai))
	y = np.multiply(np.multiply(r, np.sin(theta)), np.cos(fai))
	z = np.multiply(r, np.sin(fai))
	return x,y,z

def dec2Rad(Dec):
	'''
	degree to radian
	input:	Dec(can be np.array)
	output:	Rad
	'''
	return Dec/180.*np.pi

def rad2Dec_postive(Rad):
	'''
	radian to degree
	input:	Rad(can be np.array)
	output:	Dec(between 0 and 360)
	'''
	temp = (Rad<0)
	return (Rad + temp*2*np.pi)/np.pi*180.

def rad2Dec(Rad):
	'''
	radian to degree
	input:	Rad(can be np.array)
	output:	Dec(be negative if Rad is negative)
	'''
	return Rad/np.pi*180.

def rotz(theta):
	'''
	input: theta (in radian);
	output: Rotation matrix M.
	If (x,y,z) is a set of coordinates in the original coordinate system
	R = M * np.matrix([[x], 
	                   [y],
	                   [z])
	Where R is the coordinates in a new system obtained by rotating the original system about z-axis by an angle theta counterclockwisely.
	'''

	Rot = np.matrix([[np.cos(theta),np.sin(theta),0],[-np.sin(theta),np.cos(theta),0],[0,0,1]])
	return Rot

def rotx(theta):
	'''
	input: theta (in radian);
	output: Rotation matrix M.
	If (x,y,z) is a set of coordinates in the original coordinate system
	R = M * np.matrix([[x], 
	                   [y],
	                   [z])
	Where R is the coordinates in a new system obtained by rotating the original system about x-axis by an angle theta counterclockwisely.
	'''
	Rot = np.matrix([[1,0,0],[0,np.cos(theta),np.sin(theta)],[0,-np.sin(theta),np.cos(theta)]])
	return Rot

def roty(theta):
	'''
	input: theta (in radian);
	output: Rotation matrix M.
	If (x,y,z) is a set of coordinates in the original coordinate system
	R = M * np.matrix([[x], 
	                   [y],
	                   [z])
	Where R is the coordinates in a new system obtained by rotating the original system about y-axis by an angle theta counterclockwisely.
	'''
	Rot = np.matrix([[np.cos(theta),0,-np.sin(theta)],[0,1,0],[np.sin(theta),0,np.cos(theta)]])
	return Rot


def fx():
	'''
	return the matrix
	np.mat('-1, 0, 0;
	         0, 1, 0;
			 0, 0, 1;')
	'''

	return np.mat('-1,0,0;0,1,0;0,0,1')

def fy():
	'''
	return the matrix
	np.mat(' 1, 0, 0;
	         0,-1, 0;
			 0, 0, 1;')
	'''

	return np.mat('1,0,0;0,-1,0;0,0,1')

def fz():
	'''
	return the matrix
	np.mat(' 1, 0, 0;
	         0, 1, 0;
			 0, 0,-1;')
	'''

	return np.mat('1,0,0;0,1,0;0,0,1')

def dec22Dec(Dec2):
	'''
	change degree to dec2
	dec2 is the degree: np.array([ XXdegree, XXminute, XXsecond])
	if dec2 is negative, it should be like np.array([-XXdegree, -XXminute, -XXsecond])

	examples:
		for single dec2
			input: np.array([degree, minute, second])

			dec22Dec(np.array([45,30,00])) # 45 degree 30 minute 0 second
			output = 
			45.0
			it is the degree

		for array dec2s
			input: np.array([[degree1, minute1, second1],
			                 [degree2, minute2, second2],
							 [degree3, minute3, second3],
							 ...)
np.array([45.51]))= 
			array([  45.34166667, -128.50833333,  249.675   ])
			they are the degrees 
	'''
	if len(Dec2.shape)==2:
		Dec = Dec2[:,0]+Dec2[:,1]/60.+Dec2[:,2]/3600.
		return Dec
	elif len(Dec2.shape)==1:
		Dec = Dec2[0]+Dec2[1]/60.+Dec2[2]/3600.
		return Dec

def dec2Dec2(Dec):
	'''
	change degree to dec2
	dec2 is the degree: np.array([ XXdegree, XXminute, XXsecond])
	if dec2 is negative, it should be like np.array([-XXdegree, -XXminute, -XXsecond])

	examples:
		for single degree
			input: np.array([degree])

			dec2Dec2(np.array([45.51]))
			output = 
			array([[ 45.,  30.,  36.]])
			it is the dec2 representing 45 degree 30 minute 36 second

		for array degrees
			input: np.array([degree1, degree2, degree3....])

			dec2Dec2(np.array([-256.32, -135.5, 10.3234, 256.3333333333]))
			output = 
			array([[-256.        ,  -19.        ,  -12.        ],
			       [-135.        ,  -30.        ,    0.        ],
			       [  10.        ,   19.        ,   24.24      ],
			       [ 256.        ,   19.        ,   59.99999988]])
			they are 4 dec2s
	'''
	Dec = np.float64(Dec)
	degree = np.fix(Dec)
	#print 'degree=',degree
	_point = (Dec-degree)*60.0
	minute = np.fix(_point)
	#print 'minute=',minute
	arc = (_point-minute)*60.0
	#print 'arc=',arc
	return np.array([degree,minute,arc]).T   

def time2Dec(Time):
	'''
	to show the angle, sometimes we use 0~24h(we say it hourAngle) to represent 0~360(we say it degree) degree, so 1h means 15 degrees
	hourAngle is like np.array([XXhour, XXminute, XXsecond])
	if negative, it is like np.array([-XXhour, -XXminute, -XXsecond]
	(caution: the "minute" and "second" in hourAngle have different meaning compared with that in dec2:dec2 is the degree: np.array([ XXdegree, XXminute, XXsecond]))
	
	this function change hourAngle to degree
	examples:
		for single hourAngle
			input: np.array([XXhour, XXminute, XXsecond])

			time2Dec(np.array([10, 30, 1]))
			output = 
			157.50416666666666
			it is the degree

		for array hourAngles
			input: np.array([[hour1, minute1, second1],
			                 [hour2, minute2, second2],
							 [hour3, minute3, second3]
							 .....)

			time2Dec(np.array([[-10,-30,-2],
			                   [-1,-40,-2],
							   [3,20,12]]))	
			output = 
			array([-157.50833333,  -25.00833333,   50.05      ])
			they are 3 degrees
	'''	
	return dec22Dec(Time)/24.0*360.0

def dec2Time(dec):
	'''
	to show the angle, sometimes we use 0~24h(we say it hourAngle) to represent 0~360(we say it degree) degree, so 1h means 15 degrees
	hourAngle is like np.array([XXhour, XXminute, XXsecond])
	if negative, it is like np.array([-XXhour, -XXminute, -XXsecond]
	(caution: the "minute" and "second" in hourAngle have different meaning compared with that in dec2:dec2 is the degree: np.array([ XXdegree, XXminute, XXsecond]))
	
	this function chage degree to hourAngle
	examples:
		for single degree
			input: np.array([degree])

			dec2Time(np.array([15.51]))
			output = 
			array([[ 1. ,  2. ,  2.4]])

		for array degrees
			input: np.array([degree1, degree2, degree3....])

			dec2Time(np.array([-30.1, -15.5, 300.5]))
			array([[ -2.00000000e+00,  -1.00000000e+00,  -6.00000000e+01],
			       [ -1.00000000e+00,  -1.00000000e+00,  -6.00000000e+01],  # it should be -1, -2, 0, the error comes from the truncation error
				   [  2.00000000e+01,   2.00000000e+00,   5.96855898e-12]]) # it should be 20, 2, 0,  the error comes from the truncation error
				                                                              this error will not affect calculation....
	'''
	return dec2Dec2(dec/360.0*24)
	

def HourAngle2Horizontal(inputHourAngle, inputDE,  latitude=39.93):
	'''
	HourAngle coordinate system to Horizontal coordinate system
	HourAngle2Horizontal(inputHourAngle, inputDE,  latitude=39.93)
		inputHourAngle is the hourAngle of the object (should be degrees)
		inputDE is the latitude of the object (should be degrees)
		latitude is the geographic latitude of the observatory(default is 39.93 in Beijing)

	examples:
	for single input:
	input: np.array([XXhour, XXminute, XXsecond]), np.array([XXdegree, XXminute, XXsecond]), latitude
		HourAngle2Horizontal(time2Dec(np.array([22,32,19])),dec22Dec(np.array([-23, -43, -03])))
		output = (the answer of)
		(dec22Dec(array([ 158.        ,   10.        ,   40.69303984])), # Az
		 dec22Dec(array([ 23.        ,   8.        ,  58.66820603])))    # Alt

		HourAngle2Horizontal(time2Dec(np.array([22,32,19])),dec22Dec(np.array([-23, -43, -03])),120)
		output = (the answer of)
		(dec22Dec(array([ 147.        ,   23.        ,   58.90998699])), # Az
		 dec22Dec(array([-50.        , -37.        , -30.68584123])))    # Alt

	 for array input:
		HourAngle2Horizontal(time2Dec(np.array([[22,32,19],[23,42,19]])),dec22Dec(np.array([[-23, -43, -03],[-23,0,-32]])))
		output = (the answer of)
		(dec22Dec(array([[ 158.        ,   10.        ,   40.69303984],
		       [ 175.        ,   26.        ,    8.90077392]])), # Azs
		 dec22dec(array([[ 23.        ,   8.        ,  58.66820603],
		        [ 26.        ,  55.        ,  33.89142343]])))   # Alts
	'''
	hourAngle = inputHourAngle
	delta = inputDE
	x,y,z = sphere2Rect((hourAngle<999)*1,hourAngle,delta)
	rotmat = fy()*rotz(np.pi/2)*rotx(dec2Rad(90-latitude))*rotz(np.pi/2)*fy()
	if len(x.shape)==0:
		x=np.array([x])
		y=np.array([y])
		z=np.array([z])
	
	vectors = np.matrix([x,y,z])
	rotVectors = rotmat * vectors
	_x = rotVectors[0,:].A1
	_y = rotVectors[1,:].A1
	_z = rotVectors[2,:].A1

	_r, _Az, _Alt = rect2Sphere(_x,_y,_z)
	
	return _Az, _Alt 

def Horizontal2HourAngle(inputAz, inputAlt, latitude=39.93):	
	'''
	Horizontal coordinate system to HourAngle coordinate system
	Horizontal2HourAngle(inputAz, inputAlt, latitude=39.93)
		inputAz is the azimuth of the object (should be degrees)
		inputAlt is the altitude of the object (should be degrees)
		latitude is the geographic latitude of the observatory(default is 39.93 in Beijing)

	examples:
	for single input:
		input: np.array([XXdegree, XXminute, XXsecond]), np.array([XXdegree, XXminute, XXsecond]), latitude
		Horizontal2HourAngle(dec22Dec(np.array([158,10,40.7])), dec22Dec(np.array([23, 8, 58.7])))
		output = (the answer of)
		(time2Dec(array([ 22.        ,  32.        ,  19.00116365])), # hourAngle
		dec22dec(array([-23.        , -43.        ,  -2.97177883])))  # latitude

	for array input:
		Horizontal2HourAngle(dec22Dec(np.array([[158,10,40.7],[175,26,8.9]])),dec22Dec(np.array([[23,8,58.7],[26,55,33.9]])))
		output = (the answer of)
		(time2Dec(array([[ 22.        ,  32.        ,  19.00116365],
		       [ 23.        ,  42.        ,  18.99999131]])),
		 time2Dec(array([[-23.        , -43.        ,  -2.97177883],
			   [-23.        ,  -0.        , -31.99139656]])))
	'''
	Az = inputAz
	Alt = inputAlt

	x,y,z = sphere2Rect((Az<999)*1,Az,Alt)
	rotmat = fy()*rotz(-np.pi/2)*rotx(dec2Rad(latitude-90))*rotz(-np.pi/2)*fy()

	if len(x.shape)==0:
		x=np.array([x])
		y=np.array([y])
		z=np.array([z])
	
	vectors = np.matrix([x,y,z])
	rotVectors = rotmat * vectors
	_x = rotVectors[0,:].A1
	_y = rotVectors[1,:].A1
	_z = rotVectors[2,:].A1

	_r, _t, _DE = rect2Sphere(_x,_y,_z)
	return _t, _DE

def genTimeSeries(timeChar):
	'''
	genTimeSeries(timeChar,times=1)
	input: timeChar = '%Y-%m-%d-%H-%M-%S'
	output: np.array([[Y, m, d, H, M, S])

	'''
	return np.tile(np.array([int(timeChar.split('-')[i]) for i in range(6)]),(1,1))

def Equatorial2Horizontal(inputAlpha, inputDelta, longitude=116.4, latitude=39.93, Time=[]):
	'''
	Equatorial coordinate system to Horizontal coordinate system
		Equatorial2Horizontal(inputAlpha, inputDelta, longitude=116.4, latitude=39.93, *args)
		inputAlpha is the longitude of the object (should be degrees)
		inputDelta is the latitude of the object (should be degrees)
		longitude is the geographic longitude of the observatory (default is 116.4 in Beijing)
		latitude is the geographic latitude of the observatory(default is 39.93 in Beijing)
		Time is the timeAndData list generated by function genTimeSeries, which represent the time, default is NOW.

	It use Equatorial2HourAngle and HourAngle2Horizontal
	for the input format, please see Equatorial2HourAngle.__doc__
	for the output result, please see HourAngle2Horizontal.__doc__
	'''
	hourAngle, delta = Equatorial2HourAngle(inputAlpha, inputDelta, longitude,Time)
	Az, Alt = HourAngle2Horizontal(hourAngle, delta, latitude)
	return Az, Alt

def Horizontal2Equatorial(inputAz, inputAlt, longitude=116.4, latitude=39.93, Time=[]):
	'''
	Horizontal coordinate system to Equatorial coordinate system
		Horizontal2Equatorial(inputAz, inputAlt, longitude=116.4, latitude=39.93, *args)
		inputAz is the azimuth of the object (should be degrees)
		inputAlt is the altitude of the object (should be degrees)
		longitude is the geographic longitude of the observatory (default is 116.4 in Beijing)
		latitude is the geographic latitude of the observatory(default is 39.93 in Beijing)
		Time is the timeAndData list generated by function genTimeSeries, which represent the time, default is NOW.

	It use Horizontal2HourAngle and HourAngle2Equatorial
	for the input format, please see Horizontal2HourAngle.__doc__
	for the output result, please see HourAngle2Equatorial.__doc__
	'''
	hourAngle, delta = Horizontal2HourAngle(inputAz, inputAlt, latitude)
	alpha, delta = HourAngle2Equatorial(hourAngle,delta,longitude,Time)

	return alpha, delta

def Equatorial2HourAngle(inputAlpha,inputDelta,longitude=116.4,Time=[]):
	'''
	Equatorial coordinate system to HourAngle coordinate system
		Equatorial2HourAngle(inputAlpha,inputDelta,longitude=116.4,Time=[])
		inputAlpha is the longitude of the object (should be degrees)
		inputDelta is the latitude of the object (should be degrees)
		longitude is the geographic longitude of the observatory (default is 116.4 in Beijing)
		Time is the timeAndData list generated by function genTimeSeries, which represent the time, default is NOW.

	output: hourAngle,latitude

	for single input:
		input: np.array([XXhour, XXminute, XXsecond]), np.array([XXdegree, XXminute, XXsecond]), longitude, genTimeSeries(timeChar)
			Equatorial2HourAngle(time2Dec(np.array([22,58,21])),dec22Dec(np.array([-29,-33,-14]))) # get the current HourAngle coordinates in Beijing 
			output = (the answer of)
			(time2Dec(array([ -1.        , -34.        , -47.54860318])), dec22Dec(array([-29, -33, -14])))

			Equatorial2HourAngle(time2Dec(np.array([22,58,21])),dec22Dec(np.array([-29,-33,-14])),-30) # get the current HourAngle coordinates at w 30 
			output = (the answer of )
			(time2Dec(array([-11.        , -17.        , -44.11327565])), dec22Dec(array([-29, -33, -14])))

			Equatorial2HourAngle(time2Dec(np.array([22,58,21])),dec22Dec(np.array([-29,-33,-14])),116.4, 
								 genTimeSeries('2012-10-29-19-10-25')) # get the HourAngle coordinates in Beijing at 2012-10-29, 19:10:25
			output = (the answer of )
			(time2Dec(array([ -1.        , -29.        , -35.69711348])), dec22Dec(array([-29, -33, -14])))

	for array input:
		Equatorial2HourAngle(time2Dec(np.array([[22,58,21],[1,38,11]])),dec22Dec(np.array([[-29,-33,-14],[-57,-10,-19]])))
		output = (the answer of)
		(time2Dec(array([[ -1.        , -22.        , -39.56088124],
		       [ 19.        ,  57.        ,  30.43911876]])),
		 dec22Dec(array([[-29, -33, -14],
			   [-57, -10, -19]])))
		
		Equatorial2HourAngle(time2Dec(np.array([[22,58,21],[1,38,11]])),dec22Dec(np.array([[-29,-33,-14],[-57,-10,-19]])),116.4,genTimeSeries('2012-10-29-19-10-25'))
		output = (the answer of)
		(time2Dec(array([[ -1.        , -29.        , -35.69711348],
		       [ 19.        ,  50.        ,  34.30288652]])),
		 dec22Dec(array([[-29, -33, -14],
		       [-57, -10, -19]])))

	'''

	if len(Time)==0:
		dateAndTime = np.tile(np.array([int(time.strftime('%Y-%m-%d-%H-%M-%S').split('-')[i]) for i in range(6)]),(1,1))
	elif len(Time)==1:
		dateAndTime = Time

	alpha = inputAlpha
	delta = inputDelta
	
	lst = LST(dateAndTime[:,0], dateAndTime[:,1], dateAndTime[:,2], dateAndTime[:,3] -8 + dateAndTime[:,4]/60. + dateAndTime[:,5]/3600., longitude)
	t = lst - alpha
	#t = dec2Time(t)

	return t, delta 

def HourAngle2Equatorial(inputHourAngle,inputDelta,longitude=116.4,Time=[]):
	'''
	Equatorial coordinate system to HourAngle coordinate system
		HourAngle2Equatorial(inputHourAngle,inputDelta,longitude=116.4,Time=[])
		inputHourAngle is the hourAngle of the object (should be degrees)
		inputDelta is the latitude of the object (should be degrees)
		longitude is the geographic longitude of the observatory (default is 116.4 in Beijing)
		Time is the timeAndData list generated by function genTimeSeries, which represent the time, default is NOW.

	output: longitude,latitude

	inputHourAngle is like np.array([XXhour, XXminute, XXsecond])
	inputDelta is link np.array([XXdegree, XXminute, XXsecond])
	the input examples are similiar to that of Equatorial2HourAngle, see Equatorial2HourAngle.__doc__
	'''

	if len(Time)==0:
		dateAndTime = np.tile(np.array([int(time.strftime('%Y-%m-%d-%H-%M-%S').split('-')[i]) for i in range(6)]),(1,1))
		print 'dateAndTime=', dateAndTime
	elif len(Time)==1:
		dateAndTime = Time

	hourAngle = inputHourAngle
	delta = inputDelta
	
	lst = LST(dateAndTime[:,0], dateAndTime[:,1], dateAndTime[:,2], dateAndTime[:,3] -8 + dateAndTime[:,4]/60. + dateAndTime[:,5]/3600., longitude)
	theta = (lst - hourAngle) % 360
	#theta = dec2Time(theta)
	return theta, delta



def LST(y, m, d, ut1, EL): 
	'''
	calculate the Local sidereal time
	input: year,month,day(do not have to be integer),ut1,Longitude of the site (all input can be array)
	output: Local sidereal time in degree
	examples:
		LST(2012,10,29,19+53./60 -8 ,116.4) # the Local sidereal time in Beijing(+8) at 2012-10-29 19:53:00
		output = 
		332.86374268340001 # (degree)
	'''
	def J0(year, month, day):
		'''temp function of LST'''
		j0 = 367.0*year - np.fix(7.0*(year + np.fix((month + 9)/12.0))/4.0)+ np.fix(275.0*month/9) + day + 1721013.5
		return j0

	y = np.float64(y)
	m = np.float64(m)
	d = np.float64(d)
	ut1 = np.float64(ut1)
	EL = np.float64(EL)

	dd = np.fix(d)
	time = (d-dd)*24.0
	ut = (ut1 + time)%24
	j0 = J0(y, m, dd)
	j = (j0 - 2451545.0)/36525.0
	g0 = 100.4606184 + (36000.77004*j) + (0.000387933*j**2)- 2.583e-8*j**3
	g0 = g0%360
	gst = g0 + (360.98564724*ut/24.0)
	lst = gst + EL
	lst = lst - 360.0*np.fix(lst/360.0)
	return lst
 
def juliandate(year,month,day,hour,min,sec):
	'''
	JULIANDATE Calculate Julian date.
	Limitations: 
	This function is valid for all common era (CE) dates in the Gregorian
	calendar.
	The calculation of Julian date does not take into account leap seconds.
	Examples:
		Calculate Julian date for October 10, 2004 at 12:21:00 pm:
		juliandate(2004,10,10,12,21,0)
		output = 
		2453289.0145833334
	'''

	year = np.float64(year)
	month = np.float64(month)
	day = np.float64(day)
	hour = np.float64(hour)
	min = np.float64(min)
	sec = np.float64(sec)

	month_LT_2 = month<=2;
	year = year -1.0 * month_LT_2
	month = month + 12.0 * month_LT_2

	jd =  (np.floor( 365.25*(year + 4716.0)) + np.floor( 30.6001*( month + 1.0)) + 2.0 -
		np.floor( year/100.0 ) + np.floor( np.floor( year/100.0 )/4.0 ) + day - 1524.5 + 
		(hour + min/60.0 + sec/3600.0)/24.0)
	return jd
 
def StereographicProjection(inputLong,inputLat,scale,Dtype="sphere",*Data):
	"""
		inputLong, inputLat, Should be in degree!
	"""
	long = dec2Rad(inputLong)
	lat = dec2Rad(-inputLat)
	
	if Dtype=="sphere":	
		alpha = Data[0]
		delta = Data[1]
		N = len(alpha)
		initX, initY, initZ = sphere2Rect(np.ones((N)), alpha, delta)
	elif Dtype=="rect":
		if len(Data)==1:
			initX, initY, initZ = Data[0][0],Data[0][1],Data[0][2]
		else:
			initX, initY, initZ = Data[0],Data[1],Data[2]
	else:
		print "Error Format!\n"

	initCords = np.mat([initX, initY, initZ])
	newCords = rotx(np.pi/2-lat) *rotz(long-np.pi/2) * initCords
	temp, newAlpha, newDelta = rect2Sphere(newCords[0].A1, newCords[1].A1, newCords[2].A1)
	#pdb.set_trace()	
	newAlpha = dec2Rad(newAlpha)
	newDelta = dec2Rad(newDelta)

	Rs = scale/np.tan( (np.pi/2-newDelta)/2)
	x = Rs * np.cos(newAlpha)
	y = Rs * np.sin(newAlpha)
	return x,y

def xyzAngularDistance(Apoints, Bpoints):
	ARpos = rect2Sphere(Apoints)
	BRpos = rect2Sphere(Bpoints)
	DeltaAlpha = dec2Rad(np.abs(ARpos[1]-BRpos[1]))
	DeltaDelta = dec2Rad(np.abs(ARpos[2]-BRpos[2]))
	return rad2Dec_postive(np.arccos(np.cos(DeltaAlpha) * np.cos(DeltaDelta)))
	
def sphereAngularDistance(pointA,pointB):
	'''
	points = np.array([alpha,
					   delta]) # all in degree
	'''
	DeltaAlpha = dec2Rad(np.abs(pointA[0]-pointB[0]))
	DeltaDelta = dec2Rad(np.abs(pointA[1]-pointB[1]))
	return rad2Dec_postive(np.arccos(np.cos(DeltaAlpha) * np.cos(DeltaDelta)))


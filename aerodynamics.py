import math
from scipy.optimize import fsolve

class Compression():

	def __init__(self, M,theta,gamma=1.4):
		self.M = M
		self.theta = math.radians(theta)
		self.gamma = gamma

	def getWaveAngle(self, M):
		'''
		Returns wave angle (beta).
		Solves equation 9.23 from the book for beta.
		'''
		gamma = self.gamma
		theta = self.theta
		func = lambda beta: ((2*(1/math.tan(beta)) * (M**2 * math.sin(beta)**2 -1))/\
							(M**2*(gamma+math.cos(2*beta)) + 2)) - math.tan(theta)
		beta = fsolve(func, .707)
		beta = min(beta) # weak solution
		return beta

	def getMachNormal(self, M):
		beta = self.getWaveAngle(M)
		M_normal = M*math.sin(beta)
		return M_normal

	def getMach2(self,M,theta,gamma):
		'''
		Returns Mach number after Oblique shock. 
		M 		= the initial Mach number before the wedge.
		theta	= the wedge angle from horizontal.
		gamma 	= the ratio of specific heats. (1.4 @ air at sea level conditions)
		'''
		M_normal = self.getMachNormal(M)
		M_normal2 = math.sqrt(\
			(1+((gamma-1)/2)*M_normal**2)/\
			(gamma*M_normal**2 - (gamma-1)/2))

		beta = self.getWaveAngle(M)
		M_2 = M_normal2/(math.sin(beta - theta))
		return M_2

	def getPropertyRatios(self,M,gamma):
		'''
		Returns ratios for density (rho), pressure, and temperature
		as tuple in the respective order.

		Ratios are equivalent to:
		
		 (final property value)
		-------------------------
		 (initial property value)
		'''
		M_normal = self.getMachNormal(M)

		rho_ratio = (gamma+1)*M_normal**2/\
					(2+(gamma-1)*M_normal**2)
		pressure_ratio = 1 + (2*gamma*(M_normal**2 - 1)/\
						(gamma+1))
		temp_ratio = pressure_ratio/rho_ratio

		return(rho_ratio,pressure_ratio,temp_ratio)

	def getValues(self):
		M_2 = self.getMach2(self.M,self.theta,self.gamma)
		rho_ratio,pressure_ratio,temp_ratio = self.getPropertyRatios(self.M,self.gamma)
		return (M_2,rho_ratio,pressure_ratio,temp_ratio)

class Expansion():

	def __init__(self,M,theta,gamma = 1.4):
		self.M = M
		self.theta = math.radians(theta)
		self.gamma = gamma


	def getPrantlMeyer(self,M,gamma):
		vM = math.sqrt((gamma+1)/(gamma-1))*\
			math.atan(math.sqrt((gamma-1)*(M**2 - 1)/(gamma+1))) \
			- math.atan(math.sqrt(M**2 - 1))
		return vM

	def getInversePrantlMeyer(self,vM):
		'''
		Returns Mach number from from Prantl-Meyer function output v(M).
		Uses the Hall approximation for the inverse of the Prantl-Meyer function:
		http://www.pdas.com/pm.pdf

		The maximum error of Hall approximation for gamma = 1.4 is less than 1%.
		'''
		# Approximation coefficients for gamma  = 1.4:
		A,B,C,D,E = (1.3604,0.0962,-0.5127,-0.6722,-0.3278)
		# Input angle (vM) over the maximum turning angle.
		x = (vM / (math.pi * 0.5 * (math.sqrt(6) - 1)))**(2/3)
		# Hall approximation:
		M = (1 + A*x + B*x**2 + C*x**3) /\
			(1 + D*x + E*x**2)
		return M

	def getMach2(self,M,theta,gamma):
		vM1 = self.getPrantlMeyer(M,gamma)
		vM2 = theta + vM1
		M2 = self.getInversePrantlMeyer(vM2)
		return M2


	def getPropertyRatios(self,M2,M1,gamma):
		'''
		Returns ratios for density (rho), pressure, and temperature
		as tuple in the respective order.

		Ratios are equivalent to:
		
		 (final property value)
		-------------------------
		 (initial property value)
		'''
		temp_ratio = (1+((gamma-1)/2)*M1**2)/\
					(1+((gamma-1)/2)*M2**2)
		pressure_ratio = ((1+((gamma-1)/2)*M1**2)/\
					(1+((gamma-1)/2)*M2**2))**(gamma/(gamma-1))
		rho_ratio = pressure_ratio/temp_ratio
		return(rho_ratio,pressure_ratio,temp_ratio)

	def getValues(self):
		M_2 = self.getMach2(self.M,self.theta,self.gamma)
		rho_ratio,pressure_ratio,temp_ratio \
            = self.getPropertyRatios(M_2,self.M,self.gamma)
		return (M_2,rho_ratio,pressure_ratio,temp_ratio)

def main():

	# Test compression, should be M2 = 1.21,p2 = 2.82, T2 = 399.7
	M1,p1,T1 = (2.0,1.0,288.0) #Mach number, atm, K
	theta = 20 # degrees
	M2,rho_ratio,pressure_ratio,temp_ratio = Compression(M1,theta).getValues()
	p2 = p1*pressure_ratio
	T2 = T1*temp_ratio
    
	print('Compression:')
	print(M2,p2,T2)

	# Test expansion, should be M2 = 2.0, p2 = 0.469, T2 = 232
	M1,p1,T1 = (1.5,1.0,288.0) #Mach number, atm, K
	theta = 15 # degrees
	M2,rho_ratio,pressure_ratio,temp_ratio = Expansion(M1,theta).getValues()
	p2 = p1*pressure_ratio
	T2 = T1*temp_ratio
    
	print('Expansion:')
	print(M2,p2,T2)
    
if __name__ == '__main__':
	main()

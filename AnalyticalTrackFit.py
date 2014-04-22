from sympy import *
from ROOT import *

###############################################################################################################################
#
#        perform the analytical track fit taking into account the multiple scattering
#
###############################################################################################################################

R = Rational
F = Float
I = Integer
S = symbols

def Theta_scatering(dx,X0,impulsion,i) :

    return (F(13.6)/impulsion)*sqrt(dx[i]/X0[i])*(I(1)+F(0.038)*log(dx[i]/X0[i]))

def A(x,i) :

    return I(1)/(x[i+1]-x[i])


def Chi2(x,y,p,eps,N,dx,X0,impulsion) :

    res = I(0)
    for index in range(N) :
        res += eps[index]*(y[index]-p[index])**2

    for index in range(1,N-1) :

        numerator = (A(x,index)+A(x,index-1))*p[index]-A(x,index-1)*p[index-1]-A(x,index)*p[index+1]
        denominator = Theta_scatering(dx,X0,impulsion,index)
        res += (numerator/denominator)**2

    return res


def BuildMatrix(aChi2,N)  :

    aMatrix = [[0 for x in range(N)] for x in range(N)]
    for i in range(N) :
        for j in range(N) :
            aMatrix[i][j]=R(1,2)*diff(diff(aChi2,p[i]),p[j])

    return Matrix(aMatrix)

def SolveForPlane(invMatrix,eps,y,i,N) :
    res = 0
    print "Solution for plane %i"%i
    for j in range(N):
        res += invMatrix[i,j]*eps[j]*y[j]
        diff
    return res


def PrintResolution(invMatrix,N) :
    for i in range(N) :
        print "Total Resolution at plane %i = %f"%(i,sqrt(invMatrix[i,i]))

x = [0,59.,119.5,194.,404.,534.,664.]
y = []
p = []
sigma_sp = [3.5e-3,3.5e-3,3.5e-3,oo,3.5e-3,3.5e-3,3.5e-3]
eps =[1/s**2 for s in sigma_sp]
dx = [0.05,0.05,0.05,1.6,0.05,0.05,0.05]
X0 = [96.,96.,96.,142.,96.,96.,96.]
N = 7
impulsion = 4000

for i in range(N) :
    #x.append(S("x_%i"%i))
    y.append(S("y_%i"%i))
    p.append(S("p_%i"%i))
    #eps.append(S("eps_%i"%i))
    #dx.append(S("dx_%i"%i))
    #X0.append(S("X0_%i"%i))


aChi2 = Chi2(x,y,p,eps,N,dx,X0,impulsion)
print aChi2

Mat = BuildMatrix(aChi2,N)
invMat = Mat.inv()

for i in range(N) :
    print SolveForPlane(invMat,eps,y,i,N)

PrintResolution(invMat,N)

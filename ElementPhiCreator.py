import sympy as sym
import math

#Defining the phi functions.
x,y = sym.symbols("x y")
def phiCreatorForElements(x1,y1,x2,y2,x3,y3):#Function to create phi with given x1,x2,x3,y1,y2,y3 for CST.
    #print("The phi values:",x1,y1,x2,y2,x3,y3,A)
    A = (1/2)*(math.sqrt((y2-y1)**2+(x2-x1)**2))*(math.sqrt((y3-y2)**2+(x3-x2)**2))
    phi1 = (1/(2*A))*((x2*y3-x3*y2)+(y2-y3)*x+(x3-x2)*y)
    phi2 = (1/(2*A))*((x3*y1-x1*y3)+(y3-y1)*x+(x1-x3)*y)
    phi3 = (1/(2*A))*((x1*y2-x2*y1)+(y1-y2)*x+(x2-x1)*y)
    phi = [phi1,phi2,phi3]
    coordinates=(x1,y1,x2,y2,x3,y3)
    return phi,coordinates

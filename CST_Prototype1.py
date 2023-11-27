#Meshing program for assignment Question 1
#FEM project.
import numpy as np
from ElementPhiCreator import phiCreatorForElements
import sympy as sym
import pandas as pd

#Defining the question.
E = 200*10**9 #Youngs Modulus
mu = 0.3 #Poissons Ratio
h = 0.01 #Thickness of element

Quest_coordinates = np.array([[0,0],[4,0],[4,3]])
q = Quest_coordinates

#Meshing
#Creating the Node List; NL
Q_base = Quest_coordinates[1][0]-Quest_coordinates[0][0]
Q_height = Quest_coordinates[2][1]-Quest_coordinates[1][1]
n = int(input("Enter number of divisions along horizontal: "))
#n=2
m=n
baseDivLen = Q_base/n
heightDivLen = Q_height/m

NoNodes = ((n+1)*(m+2))/2
NL_matrix = []
#NL = np.zeros((int(NoNodes),2))#2 since 2D problem.

NodeCount = 1
st = 0 #Start index for the horizontal's x coordinate.
for i in range(0,m+1):
    temp = []
    for j in range(st,n+1):
        
        tempNodeDetails = (NodeCount,(q[0,0]+j*baseDivLen),(q[0,1]+i*heightDivLen))
        temp.append(tempNodeDetails)
        NodeCount+=1
    NL_matrix.append(temp)
    st += 1#st=st+1

noOfHorDiv = n
EL_matrix = []
elementNo = 1

for i in range(0,m): #For row below the ith division use i itself and for the upper row use i+1
    
    #print("Vertical Div: "+str(i))
    if i==m-1:
        #print("Last vertical div.\n")
        temp = [NL_matrix[i][0],NL_matrix[i][1],NL_matrix[i+1][0]]
        EL_matrix.append(temp)
        elementNo+=1
        #print(elementNo)
        break
    #elif i==0:
        #temp = NL_matrix[i][0],NL_matrix[i][1],NL_matrix[i+1][0]
        #EL_matrix.append(temp)
    for j in range(0,noOfHorDiv):
        if j==0:
            #print("Horizontal div "+str(j))
            temp = [NL_matrix[i][0],NL_matrix[i][1],NL_matrix[i+1][j]]
            EL_matrix.append(temp)
            #print(elementNo)
            if i>0:
                elementNo+=1
            continue
            pass
        #For the next horizontal division we need to define 2 elements a rectangular division at one time
        #Next element in that vertical division
        #print("Horizontal div "+str(j))
        temp1 = [EL_matrix[elementNo-1][1],EL_matrix[elementNo-1][2],NL_matrix[i+1][j]]
        EL_matrix.append(temp1)
        elementNo += 1
        #print(elementNo)
        #Second one
        temp2 = [EL_matrix[elementNo-1][0],NL_matrix[i][j+1],EL_matrix[elementNo-1][2]]
        EL_matrix.append(temp2)
        elementNo += 1
        #print(elementNo)
        
    noOfHorDiv -= 1
  

    
#Forming the golbal stiffness matrix. Assembling.
#Creating the stress-strain constants.
#Plane stress
C11 = (E/(1-mu**2))
C12 = ((E*mu)/(1-mu**2))
C21 = ((E*mu)/(1-mu**2))
C22 = (E/(1-mu**2))
C33 = (E*(1-mu)/(1-mu**2))

SJ = np.zeros((2*int(NoNodes),2*int(NoNodes)))
LoadVector = np.zeros((2*int(NoNodes),1))
x,y,t = sym.symbols("x y t")

elementStiffnessOfAll = []
elementLoadVectorOfAll = []

for i in range(0,elementNo):
    #The local stiffness matrix will be created first for 1 element here.
    print(i)
    node1 = EL_matrix[i][0]
    node2 = EL_matrix[i][1]
    node3 = EL_matrix[i][2]
    MemDofNo = (node1[0]*2-1),(node1[0]*2),(node2[0]*2-1),(node2[0]*2),(node3[0]*2-1),(node3[0]*2)
    #Adding tx ty values.
    if node2[1]==4.0:
        EL_matrix[i].append((50*10**3,0))
    else:
        EL_matrix[i].append((0,0))
    
    phis,coordinates = phiCreatorForElements(node1[1], node1[2], node2[1], node2[2], node3[1], node3[2])
    x1,y1,x2,y2,x3,y3=coordinates #For this element
    #Creating the member stiffness matrix.
    elementStiffness = np.zeros((6,6))
    elementLoadVector = np.zeros((6,1))
    
    #First Set of filling up elemental stiffness matrix
    for r in range(1,7):#Will run from 1 to 6, specified this way cause CST have 6 DOF, and it prolly wont change.
        if r in (1,2):
            w=phis[0]
            #print("For "+str(r),w)
        elif r in (3,4):
            #print("For "+str(r),w)
            w=phis[1]
        elif r in (5,6):
            #print("For "+str(r),w)
            w=phis[2]
        if r%2 != 0:
            a=C11
            b=C33
            countForCol = 0
            for c in range(0,6,2):
                #print("Odd Col ",c,"\n")
                f1 = h*(a*sym.diff(w,x)*sym.diff(phis[countForCol],x))
                f2 = h*(b*sym.diff(w,y)*sym.diff(phis[countForCol],y))
                #print(f1,f2)
                result = sym.integrate(f1,(y,y1,((y3-y1)/(x3-x1))*x),(x,x1,x3))+sym.integrate(f2,(y,y1,((y3-y1)/(x3-x1))*x),(x,x1,x3))
                elementStiffness[r-1][c]=elementStiffness[r-1][c]+result
                countForCol+=1
        else:
            a=C33
            b=C22
            countForCol = 0
            for c in range(1,6,2):
                #print("Even Col ",c,"\n")
                f1 = h*(a*sym.diff(w,x)*sym.diff(phis[countForCol],x))
                f2 = h*(b*sym.diff(w,y)*sym.diff(phis[countForCol],y))
                #print(f1,f2)
                result = sym.integrate(f1,(y,y1,((y3-y1)/(x3-x1))*x),(x,x1,x3))+sym.integrate(f2,(y,y1,((y3-y1)/(x3-x1))*x),(x,x1,x3))
                elementStiffness[r-1][c]=elementStiffness[r-1][c]+result
                countForCol+=1
                pass
            pass
        #Calculating and assembling the load vector
        if EL_matrix[i][3][0]==50*10**3: # or EL_matrix[i][3][0] != 50 
            #For n number of elements it seems that the edge triangles are the only ones which have forces acting on them. i.e. the seg 2 of all elements where load is acting is used like class work question. Bit of an bummer, needs future improvement.
            #Basically I am keeping the format of line integral same here like the single element. But this may be not be the case everywhere.
            #Calculating the line integral for CST.
            tx = EL_matrix[i][3][0]
            ty = EL_matrix[i][3][1]
            #Seg 1
            if r%2!=0:
                T=0#tx=0 for segment 1
            else:
                T=ty
            a=t
            b=node1[2] #Always y constant.
            dx_dt=sym.diff(a,t)
            dy_dt=sym.diff(b,t)
            seg1 = sym.integrate(h*w.subs({x:a,y:b})*T*(dx_dt**2 + dy_dt**2)**0.5, (t, node1[1], node2[1]))
            #Seg2
            if r%2!=0:
                T=tx
            else:
                T=ty
            a=node2[1]
            b=t
            dx_dt=sym.diff(a,t)
            dy_dt=sym.diff(b,t)
            seg2 = sym.integrate(h*w.subs({x:a,y:b})*T*(dx_dt**2 + dy_dt**2)**0.5, (t, node2[2], node3[2]))
            #Seg3
            if r%2!=0:
                T=0#tx=0 for segment 2
            else:
                T=ty
            a=t
            b=((y3-y1)/(x3-x1))*a
            dx_dt=sym.diff(a,t)
            dy_dt=sym.diff(b,t)
            seg3 = sym.integrate(h*w.subs({x:a,y:b})*T*(dx_dt**2 + dy_dt**2)**0.5, (t, node1[1], node3[1]))
            
            elementLoadVector[r-1] = elementLoadVector[r-1]+seg1+seg2+seg3
        pass
    #Second Set for filling up elemental matrix
    for r in range(1,7):#Will run from 1 to 6, specified this way cause CST have 6 DOF, and it prolly wont change.
        if r in (1,2):
            w=phis[0]
        elif r in (3,4):
            w=phis[1]
        elif r in (5,6):
            w=phis[2]
        
        if r%2!=0:
            a=C12
            b=C33
            countForCol = 0
            for c in range(1,6,2):
                f1 = h*(a*sym.diff(w,x)*sym.diff(phis[countForCol],y))
                f2 = h*(b*sym.diff(w,y)*sym.diff(phis[countForCol],x))
                result = sym.integrate(f1,(y,y1,((y3-y1)/(x3-x1))*x),(x,x1,x3))+sym.integrate(f2,(y,y1,((y3-y1)/(x3-x1))*x),(x,x1,x3))
                elementStiffness[r-1][c]=elementStiffness[r-1][c]+result
                countForCol+=1
        else:
            a=C33
            b=C21
            countForCol = 0
            for c in range(0,6,2):
                f1 = h*(a*sym.diff(w,x)*sym.diff(phis[countForCol],y))
                f2 = h*(b*sym.diff(w,y)*sym.diff(phis[countForCol],x))
                result = sym.integrate(f1,(y,y1,((y3-y1)/(x3-x1))*x),(x,x1,x3))+sym.integrate(f2,(y,y1,((y3-y1)/(x3-x1))*x),(x,x1,x3))
                elementStiffness[r-1][c]=elementStiffness[r-1][c]+result
                countForCol+=1
                pass
            pass
        pass
    elementStiffnessOfAll.append(elementStiffness)
    elementLoadVectorOfAll.append(elementLoadVector)
    
    for r in range(0,elementLoadVector.shape[0]):
        LoadVector[MemDofNo[r]-1,0] = LoadVector[MemDofNo[r]-1,0] + elementLoadVector[r,0] 
    
    for c in range(0,elementStiffness.shape[0]): #I am giving the size of member stiffness matrix as the no of rows and columns to traverse, but the format may change for frame and truss. Look for it.
        for r in range(0,elementStiffness.shape[0]):
            SJ[MemDofNo[r]-1,MemDofNo[c]-1] = SJ[MemDofNo[r]-1,MemDofNo[c]-1] + elementStiffness[r,c]#This using of row and column numbers based on the size of k_temp may back fire, or may not, as the size of member stiffness matrix will always be same as that of the no of dofs for a memeber available
            pass
        pass
    
df = pd.DataFrame(SJ)
df.to_excel("9_Element\AssembledCheck9.xlsx", index=False)
df2 = pd.DataFrame(LoadVector)
df2.to_excel("9_Element\AssembledLoadVector9.xlsx", index=False)

#Prototype of applying boundary conditions using compensation method.
dispVector = []
countU = 1
countV = 1
for i in range(1,2*int(NoNodes)+1):
    if i%2 != 0:
        dispVector.append(sym.symbols("u"+str(countU)))
        countU+=1
    elif i%2==0:
        dispVector.append(sym.symbols("v"+str(countV)))
        countV+=1

#Since this program is for our own project I am not generalising it and applying as per the known values.
modifiedSJ = SJ
equationsToSolve = []
for i in range(0,2*(n+1)):
    modifiedSJ[i][i] += 10**20 #Check if it add or multiply

df3 = pd.DataFrame(modifiedSJ)
df3.to_excel("4_Element\ModifiedSJ_Test_4eles.xlsx", index=False)

for r in range(0,modifiedSJ.shape[0]):
    temp=0
    for c in range(0,modifiedSJ.shape[1]):
        temp += modifiedSJ[r][c]*dispVector[c]
    equationsToSolve.append(sym.Eq(temp,float(LoadVector[r])))

FinalDisp = sym.solve(equationsToSolve,dispVector)

#Displacement
#df_disp = pd.DataFrame(FinalDisp)
#df_disp.to_excel("4_Element\Displacement_4elements.xlsx", index = False)

#Reactions
Reactions_SJ = SJ@np.array(list(FinalDisp.values()))

df4 = pd.DataFrame(Reactions_SJ)

df4.to_excel("4_Element\ReactionsCheck_4eles.xlsx", index=False)

Reactions_modSJ = modifiedSJ@np.array(list(FinalDisp.values()))

df5 = pd.DataFrame(Reactions_modSJ)

df5.to_excel("4_Element\ModdedR_4eles.xlsx",index=False)

#Add the tx and ty values in the EL_matrix. Fix the application of boundary conditions and calculations of the line integral.
#Application of boundary conditions is through compensation method, increasing stiffness of the boudary conditioned diagonal values.
 
  
'''
Daniel Zheng, 17/10/26 meow!!
'''
################[SET THESE BEFORE STARTING AND KEEP THE SAME]##################
filepath = 'C:\\Uni work\\Electric Propulsions\\IterativeBase - Copy.fem'
isolatorratio = (1.5) #bottom thickness / side thickness (NOT IN CODE RN), if changed baseline file may need to be manually editted
###############################################################################

import femm
import numpy as np
from scipy.optimize import differential_evolution #not used yet will be used for iterative section later

def open_femm():
    femm.openfemm()
    femm.opendocument(filepath)
    
def run_femm():
    femm.mi_createmesh()
    femm.mi_analyse(0)
    femm.mi_loadsolution()
    
def check_nodes():
    '''prints each of the shield nodes that are moved in console'''
    for i in shieldnodes:
        print(i.__repr__())
    
class ShieldingNodes():
    def __init__(self, name, x, y, offset, swh, swt, swg): 
        self.name = name
        self.x = x
        self.y = y
        self.offset = offset #x displacement from channel
        self.swh = swh #shield wall height (outer length)
        self.swt = swt #shield wall thickness
        self.swg = swg #shiedl wall gap (x displacement from channel)
        
    def __str__(self): #what to show to people
        return f"{self.name}, x = {self.x}, y = {self.y}"
    
    def __repr__(self): #what to show for debugging, prob anything else thats not x or y
        return f"ShieldingNode, {self.name}, x = {self.x}, y = {self.y}, height = {self.swh}, thickness = {self.swt}, gap = {self.swg}"
    
    def move(self, dx, dy):
        self.x += dx
        self.y += dy
        
#Shielding points (sp = shielding points) INITIAL DIMENSIONS FROM BASE FEMM MODEL, draw rectangle tool not used incase of future geometry iteration thing
sp1 = ShieldingNodes("sp1", -94.25, 5, 2, 72, 2, 2)
sp2 = ShieldingNodes("sp2", -94.25, 77, 2, 72, 2, 2)
sp3 = ShieldingNodes("sp3", -92.25, 77, 2, 72, 2, 2)
sp4 = ShieldingNodes("sp4", -92.25, 49, 2, 72, 2, 2)
sp5 = ShieldingNodes("sp5", -92.25, 26, 2, 72, 2, 2)
sp6 = ShieldingNodes("sp6", -71.75 ,26, 2, 72, 2, 2)
sp7 = ShieldingNodes("sp7", -71.75, 49, 2, 72, 2, 2)
sp8 = ShieldingNodes("sp8", -71.75, 77, 2, 72, 2, 2)
sp9 = ShieldingNodes("sp9", -69.75, 77, 2, 72, 2, 2)
sp10 = ShieldingNodes("sp10", -69.75, 5, 2, 72, 2, 2)

shieldnodes = [sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, sp9, sp10]

def redraw_shielding(shieldthickness, shieldheight, shieldgap):
    '''Redrawst the mag shielding walls given: a thickness, height (from point 1 to 2), and gap (displacement in x from the channel)'''
    #checks inputs are valid
    if shieldthickness<0 or shieldheight<0 or shieldgap<0: #dimensions cant be negative (maybe change for height? small did but i think that would be bad for mag field)
        raise ValueError("Values must be >= 0")
    if shieldthickness+shieldgap >= 6.75 and shieldheight >= 73: #intersects with the top plate
        raise ValueError("Physically impossible scenario")
    if shieldthickness+shieldgap > 7.4: #intersects with coil
        raise ValueError("Physically impossible scenario")
    if shieldheight < 49: #i just dont wanna move the isolator lol
        raise ValueError("Shield height must be more than 49! (lazy programming?)")
    
    oldthickness = sp1.swt
    oldgap = sp1.swg
    gapchange = shieldgap - oldgap
    thicknesschange = shieldthickness - oldthickness
    femm.mi_clearselected() #incase any selected before
    #selects and deletes all nodes
    for i in shieldnodes:
        femm.mi_selectnode(i.x, i.y)
    femm.mi_deleteselectednodes()
    #moves nodes to correvt x coordniate
    for i in (sp1, sp2, sp3, sp4, sp5):
        i.move(-gapchange, 0)
    for i in (sp6, sp7, sp8, sp9, sp10):
        i.move(gapchange, 0)
    #makes thickness correct
    for i in (sp1, sp2):
        i.move(-thicknesschange, 0)
    for i in (sp9, sp10):
        i.move(thicknesschange, 0)
    #makes y coordinate correct    
    for i in (sp2, sp3, sp8, sp9):
        i.y = shieldheight+5 
    #update all additional info 
    for i in shieldnodes:
        i.swh = shieldheight
        i.swt = shieldthickness
        i.swg = shieldgap
    #channel doesnt move so 29 can be used as a reference point
    sp5.y = (29-(shieldgap*isolatorratio))
    sp6.y = (29-(shieldgap*isolatorratio))
    #draw new lines+nodes
    for i, j in zip(shieldnodes[:-1], shieldnodes[1:]): #combines combinations 1 and 2, 2 and 3, 3 and 4...
        femm.mi_drawline(i.x, i.y, j.x, j.y)  
    #finish drawing theuster (isolator and connect to coils)
    femm.mi_addsegment(-99, 5, sp1.x, sp1.y)
    femm.mi_addsegment(-66.3493, 5, sp10.x, sp10.y)
    femm.mi_addsegment(sp4.x, sp4.y, -90.25, 49)
    femm.mi_addsegment(sp7.x, sp7.y, -73.75, 49)
    femm.mi_addsegment(sp1.x, sp1.y, sp10.x, sp10.y)
    
def change_current(newcurrent):
    if newcurrent <0:
        raise ValueError("Value must be more than 0")
    femm.mi_modifycircprop("Outer Coil", 1, newcurrent)
    femm.mi_modifycircprop("Inner Coil", 1, newcurrent)
    
def change_turns(innerturns, outerturns):
    if innerturns <0 or outerturns <0:
        raise ValueError("Value must be more than 0")
    #select and adjust inner coil
    femm.mi_selectlabel(-62.411, 41.968)
    femm.mi_setblockprop("24 AWG", 1, 0, 'Inner Coil', 0, 2, innerturns)
    #select and adjust outer coil 
    femm.mi_selectlabel(-100.15, 18.60)
    femm.mi_setblockprop("24 AWG", 1, 0, 'Outer Coil', 0, 2, outerturns)
    
########################################ITERATIVE PART########################################

#DO LATER: THESE ARE INPUTTED (maybe idk if its necessary to input)
middlerp = (-82, 32) #middle reference point
leftrp = (-87.2, 32) #left reference point
rightrp = (-76.8, 32) #right reference point

#these retun a tuple (B_x, B_y)
def get_fdM(y): 
    '''finds the mag flux density at the centre of the channel at a given Y coordinate'''
    return(femm.mo_getb(middlerp[0], y))
def get_fdL(y): 
    '''finds the mag flux density at the inner left of the channel at a given Y coordinate'''
    return(femm.mo_getb(leftrp[0], y))
def get_fdR(y):
    '''finds the mag flux density at the inner right of the channel at a given Y coordinate'''
    return(femm.mo_getb(rightrp[0], y))

#NEED TO WORK ON X AND Y COMPONENTS DIFFERENCE
def get_total_middlefd():
    ycoords = np.linspace(32, 79.5, 1000) #CHECK 79.5 IS CORRECT FOR TOP OF THE CHANNEL
    totalx = 0 
    totaly = 0
    for i in ycoords:
        totalx += get_fdM(i)[0]
        totaly += get_fdM(i)[1]
    return(totalx, totaly)

def get_total_sidefd():
    ycoords = np.linspace(32, 79.5, 1000) #CHECK 79.5 IS CORRECT FOR TOP OF THE CHANNEL
    totalx = 0 
    totaly = 0
    for i in ycoords:
        totalx += (get_fdL(i)[0] + get_fdR(i)[0])
        totaly += (get_fdL(i)[1] + get_fdR(i)[1])
    return(totalx, totaly)

def get_xuniformity():
    '''finds root mean square error for the X component of flux density at 1000 Y coordinate'''
    ycoords = np.linspace(32, 79.5, 1000) #CHECK 79.5 IS CORRECT FOR TOP OF THE CHANNEL
    uniformity = 0
    for i in ycoords:
        lp = get_fdL(i)[0]
        rp = get_fdR(i)[0]
        mp = get_fdM(i)[0]
        uniformity += np.sqrt(((lp-mp)**2 + (rp-mp)**2)/3)
    return(uniformity)

def get_yuniformity():
    '''finds root mean square error for the Y component of flux density at 1000 y coordinates'''
    ycoords = np.linspace(32, 79.5, 1000) #CHECK 79.5 IS CORRECT FOR TOP OF THE CHANNEL
    uniformity = 0
    for i in ycoords:
        lp = get_fdL(i)[1]
        rp = get_fdR(i)[1]
        mp = get_fdM(i)[1]
        uniformity += np.sqrt(((lp-mp)**2 + (rp-mp)**2)/3)
    return(uniformity)

leakrpL = (-90.26, 49.5)
leakrpR = (-73.7, 49.5)

def get_leak():
    '''finds the total mag fluxs at 1000 y levels for the gaps between the mag screen and channel, returns a tupple for (X,Y) components'''
    total_xleak = 0
    total_yleak = 0
    ycoords = np.linspace(leakrpL[1], sp1.swh-0.5, 1000) #sp1.swh so 1000 points only between the reference point and height of the shielding
    for i in ycoords:
        Lleak = femm.mo_getb(leakrpL[0], i)
        Rleak = femm.mo_getb(leakrpR[0], i)
        total_xleak += (Rleak[0] + Lleak[0])
        total_yleak += (Rleak[1] + Lleak[1])
        return((total_xleak, total_yleak))

##################################################################################################################

runcount = 0
def get_score(x):
    '''scores a given design using an objective function, lower score better
    input: Shield thickness, shield height, shield gap, inner turns, outer turns, current'''
    try: #normally run this, CHANGES VARIABLES
        redraw_shielding(x[0], x[1], x[2])
        change_turns(x[3], x[4])
        change_current(x[5])
    except ValueError: #when a valueerror is raised from a function like redraw_shielding() itll do this instead
        return 1e6
    
    print(check_nodes())
    run_femm()
    
    uniformity = get_yuniformity() + get_xuniformity()
    leakx = get_leak()[0]
    leaky = get_leak()[1]
    side_Xfd = get_total_sidefd()[0]
    side_Yfd = get_total_sidefd()[1]
    middle_Xfd = get_total_middlefd()[0]
    middle_Yfd = get_total_middlefd()[1]

    score = uniformity + side_Xfd + side_Yfd + middle_Xfd + middle_Yfd + (1e6*leaky) + (1e6*leakx)
    return(score)

########################################OPTIMISATION PART########################################
#Order: Shield thickness, shield height, shield gap, inner turns, outer turns, current
bounds = [(0.1, 6.75), (0.1,80), (0.1, 6.75), (1,1000000), (1, 1000000), (1,10)]    

def optimise():
    results = differential_evolution(get_score, bounds, integrality=(False, False, False, True, True, False))
    print(results)
    
    
    
    
    
    
    
    

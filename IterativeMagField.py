'''
Daniel Zheng, 17/10/26 meow!!
'''
################[SET THESE BEFORE STARTING AND KEEP THE SAME]##################
filepath = 'C:\\Uni work\\Electric Propulsions\\IterativeBase - Copy.fem'
isolatorratio = (1.5)
###############################################################################

import femm
import numpy as np
from scipy.optimize import differential_evolution

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
'''PLAN: 
    start from reference points and get sum of values of points in periodic spacings upwards until top of channel and we can use the values at left, middle, and right to see how good the mag field is
    penalties and bonus points can then be given for bottom 5 middle points etc'''
middlerp = (-82, 32) #middle reference point
leftrp = (-87.2, 32) #left reference point
rightrp = (-76.8, 32) #right reference point

#x1,x2,y1,y2
region1 = (-82.5, -81.5, 78, 86)
region2L = (-87.2, -86.2, 40, 82)
region2R = (-77.8, -76.8, 40, 82)

class lowfluxdensity(Exception):
    pass

def get_region2_score(x1,x2,y1,y2):   #imagine a rectangle: start left to right, bottom to top
    '''Gets a regional score,not just a point like before, minimising Br at channel inner walls'''
    xsegments = np.linspace(x1,x2,5)
    ysegments = np.linspace(y1,y2,100)
    totalr = 0
    for i in xsegments:
        for j in ysegments:
            Br = femm.mo_getb(i, j)[0]
            totalr += (Br/0.03)**2       #divide by 0.03 so design score later is better, dimensionless!
    return(totalr/500)  #500 from 5*100, returns a tuple (Br, Bz)

def get_region1_score(x1,x2,y1,y2):   #imagine a rectangle: start left to right, bottom to top
    '''Gets a regional score for how close to the target BR is it (0.03), best score is Bz is minimised'''
    xsegments = np.linspace(x1,x2,5)
    ysegments = np.linspace(y1,y2,100)
    totalr = 0
    totalz = 0
    for i in xsegments:
        for j in ysegments:
            Br, Bz = femm.mo_getb(i, j)
            totalz += (Bz/0.03)**2
            totalr += ((Br-0.03)/0.03)**2       #divide by 0.03 so design score later is better, dimensionless!
    return((totalr/500)+(totalz/500)) 

def check_max_fd():
    '''makes sure the flux is atleast 0.02T (documentation said 0.03T max so 0.02 for now) and rejects design if it isnt'''
    ycoords = np.linspace(78, 86, 50)
    values = []
    for i in ycoords:
        Br, Bz = femm.mo_getb(-82, i)
        values.append(abs(Br))
    if max(values) < 0.02:
        raise lowfluxdensity()

def get_parallel_score(): #if field lines stay parallel to the channel wall
    '''Measures how parallel flux is to the channel and gives a penalty if it isnt perfectly parallel'''
    total = 0
    ycoords = np.linspace(40, 79.5, 100)
    for i in ycoords:
        Lflux = femm.mo_getb(leftrp[0], i)
        Rflux = femm.mo_getb(rightrp[0], i)
        bonus = ((abs(Lflux[1])/np.hypot(*Lflux)) + (abs(Rflux[1])/np.hypot(*Rflux))) / 2 #find cos from trig but we leave in tersmf of cos so its between 0 and 1
        total += bonus
    return(1 - (total/100))  #so 0 best score, max penalty = 1
        
    
########################################OPTIMISATION PART########################################
runcount = 0 
def get_score(x):
    '''where objective function thingy is, lower score = more desirable
    Order for input list: Shield thickness, shield height, shield gap, inner turns, outer turns, current'''
    global runcount
    runcount += 1 
    try: #normally run this, CHANGES DIMENSIONS
        redraw_shielding(x[0], x[1], x[2])
        change_turns(x[3], x[4])
        change_current(x[5])
        run_femm()
        check_max_fd() #reject if too low (<0.02T)
        
        region1_score = get_region1_score(*region1)
        region2_score = (get_region2_score(*region2L) + get_region2_score(*region2R))/2
        parallel = get_parallel_score()
        
        penalty = parallel #add others later?
        
        score = region1_score + region2_score + penalty

        #score = uniformity + (5*normalval) - (10*depthscore) + (5*fd_error) # - middlefdR + middlefdZ
        
        print("\n########################################")
        print("Run: ", runcount)
        print("########################################")
        print("Shield thickness: ", sp1.swt)
        print("Shield height: ", sp1.swh)
        print("Shield gap: ", sp1.swg)
        print("Inner coil turns: ", x[3])
        print("Outer coil turns: ", x[4])
        print("Current: ", x[5])
        
        print("region1_score: ", region1_score)
        print("region2_score: ", region2_score)
        print("Stickiness penalty: ", stickiness)
        print("########################################")
        print("Score: ", score)
        
        return(score)    
    
    except (ValueError, lowfluxdensity): #when a valueerror is raised from a function like redraw_shielding() itll do this instead
        print("\n Run", runcount, "rejected")
        return 1e10 #large score so programme will not use this design 

#Order: Shield thickness, shield height, shield gap, inner turns, outer turns, current
bounds = [(0.5, 6.75), (49.5,72.95), (0.5, 6.75), (1,100000), (1, 100000), (1,10)]    
def optimise():
    results = differential_evolution(get_score, bounds, integrality=(False, False, False, True, True, False),popsize=10, maxiter=100, polish=False, disp=True)
    print("\n\n")
    print("################################### Optimisation Finished ###################################")
    print("Nodes:")
    print(check_nodes())
    print("#############################################################################################")
    print("Shield thickness =", results.x[0])
    print("Shield height    =", results.x[1])
    print("Shield gap       =", results.x[2])
    print("Inner turns      =", (results.x[3]))
    print("Outer turns      =", (results.x[4]))
    print("Current          =", results.x[5]) 
    print("#############################################################################################")
    print("Best score: ", results.fun)


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
    def __init__(self, name, x, y, swh, swt, swg): 
        self.name = name
        self.x = x
        self.y = y
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
        
class PlateNodes():
    def __init__(self, name, x, y, in_out, thickness, length):
        self.name = name
        self.x = x
        self.y = y
        self.in_out = in_out
        self.thickness = thickness
        self.length = length
    def __str__(self):
        return f"{self.name}, x = {self.x}, y = {self.y}, {self.in_out}"
    def __repr__(self):
        return f"{self.name}, x = {self.x}, y = {self.y}, {self.in_out}, thickness = {self.thickness}, length = {self.length}"
    def move(self, dx, dy):
        self.x += dx
        self.y += dy
    
#PLAN: add new nodes between isolator points and top shield height and between 2 top nodes. delete 2 lines and connect 2 new nodes to form chamfered geometry 
class MagCircuitChamfer():
    def __init__(self, name, x, y, chamferwidth, chamferlength, which): 
        self.name = name
        self.x = x
        self.y = y
        self.chamferwidth = chamferwidth 
        self.chamferlength = chamferlength
    def move(self, dx, dy):
        self.x += dx
        self.y += dy
    
#Shielding points (sp = shielding points) INITIAL DIMENSIONS FROM BASE FEMM MODEL, draw rectangle tool not used incase of future geometry iteration thing
sp1 = ShieldingNodes("sp1", -37.995, -21, 42, 3, 1)
sp2 = ShieldingNodes("sp2", -37.995, 21, 42, 3, 1) 
sp3 = ShieldingNodes("sp3", -35, 21, 42, 3, 1)
sp4 = ShieldingNodes("sp4", -35, -13, 42, 3, 1)
sp5 = ShieldingNodes("sp5", -35, -18.5, 42, 3, 1)
sp6 = ShieldingNodes("sp6", -18 , -18.5, 42, 3, 1)
sp7 = ShieldingNodes("sp7", -18, -13, 42, 3, 1)
sp8 = ShieldingNodes("sp8", -18, 21, 42, 3, 1)
sp9 = ShieldingNodes("sp9", -15.005, 21, 42, 3, 1)
sp10 = ShieldingNodes("sp10", -15.005, -21, 42, 3, 1)
shieldnodes = [sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, sp9, sp10]

#Top plate points
tpp1 = PlateNodes("tpp1", -45, 24, "out", 2, 10)
tpp2 = PlateNodes("tpp2", -45, 26, "out", 2, 10)
tpp3 = PlateNodes("tpp3", -35, 26, "out", 2, 10)
tpp4 = PlateNodes("tpp4", -35, 24, "out", 2, 10)
tpp5 = PlateNodes("tpp5", -38, 24, "out", 2, 10)
tpp6 = PlateNodes("tpp1", -15, 24, "in", 2, 18)
tpp7 = PlateNodes("tpp2", -18, 24, "in", 2, 18)
tpp8 = PlateNodes("tpp3", -18, 26, "in", 2, 18)
tpp9 = PlateNodes("tpp4", -0, 26, "in", 2, 18)
tpp10 = PlateNodes("tpp5", -0, 24, "in", 2, 18)
movablenodes = [tpp2,tpp3,tpp4,tpp7,tpp8,tpp9]
outerplatenodes = [tpp1,tpp2,tpp3,tpp4,tpp5]
innerplatenodes = [tpp6,tpp7,tpp8,tpp9,tpp10]

#Chamfer points (points that make up the mag screening as well!)
#S for screen P for plate
Pleftwidth = MagCircuitChamfer("leftwidth", -36, 21, 1, 10, "Screen")
Pleftlength = MagCircuitChamfer("leftlength", -38, 11, 1,10, "Screen")
Prightwidth = MagCircuitChamfer("rightwidth", -17, 21, 1, 10, "Screen")
Prightlength = MagCircuitChamfer("rightlength", -15, 11, 1, 10, "Screen")
platechamfernodes = [Pleftwidth, Pleftlength, Prightwidth, Prightlength]

Sleftwidth = MagCircuitChamfer("leftwidth", -36, 21, 1, 10, "Screen")
Sleftlength = MagCircuitChamfer("leftlength", -38, 11, 1,10, "Screen")
Srightwidth = MagCircuitChamfer("rightwidth", -17, 21, 1, 10, "Screen")
Srightlength = MagCircuitChamfer("rightlength", -15, 11, 1, 10, "Screen")
screenchamfernodes = [Sleftwidth, Sleftlength, Srightwidth, Srightlength]

def chamfer_top_plate(width, length):
    if width < 0 or length < 0:
        raise ValueError("Must be more than 0")
    if width > abs(tpp3.y-tpp4.y) or length > abs(tpp4.x-tpp5.x):
        raise ValueError("cant do that twin")
        
    femm.mi_clearselected()
    Pleftwidth.x = tpp3.x
    Pleftwidth.y = tpp3.y - width
    Pleftlength.x = tpp4.x - length
    Pleftlength.y = tpp4.y

    Prightwidth.x = tpp8.x
    Prightwidth.y = tpp8.y - width
    Prightlength.x = tpp7.x + length
    Prightlength.y = tpp7.y
    
    for i in platechamfernodes:
        femm.mi_addnode(i.x,i.y)
        
    femm.mi_selectnode(tpp4.x, tpp4.y)
    femm.mi_selectnode(tpp7.x, tpp7.y)
    femm.mi_deleteselectednodes()
    
    femm.mi_addsegment(Pleftwidth.x, Pleftwidth.y, Pleftlength.x, Pleftlength.y)
    femm.mi_addsegment(Prightwidth.x, Prightwidth.y, Prightlength.x, Prightlength.y)

def chamfer_screening(width, length):
    femm.mi_clearselected()
    Sleftwidth.x = sp3.x - width
    Sleftwidth.y = sp3.y 
    Sleftlength.x = sp2.x
    Sleftlength.y = sp2.y - length
    
    Srightwidth.x = sp8.x + width
    Srightwidth.y = sp8.y 
    Srightlength.x = sp9.x
    Srightlength.y = sp9.y - length
    
    for i in screenchamfernodes:
        femm.mi_addnode(i.x,i.y)
    
    femm.mi_selectnode(sp2.x, sp2.y)
    femm.mi_selectnode(sp9.x, sp9.y)
    femm.mi_deleteselectednodes()
    
    femm.mi_addsegment(Sleftwidth.x, Sleftwidth.y, Sleftlength.x, Sleftlength.y)
    femm.mi_addsegment(Srightwidth.x, Srightwidth.y, Srightlength.x, Srightlength.y)
    
def undo_chamfers():
    femm.mi_clearselected()
    #Screen
    femm.mi_selectnode(Sleftwidth.x, Sleftwidth.y)
    femm.mi_selectnode(Sleftlength.x, Sleftlength.y)
    femm.mi_selectnode(Srightwidth.x, Srightwidth.y)
    femm.mi_selectnode(Srightlength.x, Srightlength.y)
    femm.mi_deleteselectednodes()
    
    femm.mi_addnode(sp2.x, sp2.y)
    femm.mi_addnode(sp9.x, sp9.y)
    
    femm.mi_addsegment(sp1.x, sp1.y, sp2.x, sp2.y)
    femm.mi_addsegment(sp2.x, sp2.y, sp3.x, sp3.y)
    femm.mi_addsegment(sp8.x, sp8.y, sp9.x, sp9.y)
    femm.mi_addsegment(sp9.x, sp9.y, sp10.x, sp10.y)
    
    #Plate
    femm.mi_selectnode(Pleftwidth.x, Pleftwidth.y)
    femm.mi_selectnode(Pleftlength.x, Pleftlength.y)
    femm.mi_selectnode(Prightwidth.x, Prightwidth.y)
    femm.mi_selectnode(Prightlength.x, Prightlength.y)
    femm.mi_deleteselectednodes()
    
    femm.mi_addnode(tpp4.x, tpp4.y)
    femm.mi_addnode(tpp7.x, tpp7.y)
    
    femm.mi_addsegment(tpp3.x, tpp3.y, tpp4.x, tpp4.y)
    femm.mi_addsegment(tpp4.x, tpp4.y, tpp5.x, tpp5.y)
    femm.mi_addsegment(tpp7.x, tpp7.y, tpp8.x, tpp8.y)
    femm.mi_addsegment(tpp6.x, tpp6.y, tpp7.x, tpp7.y)
    
    
    

def redraw_top_plate(outerthickness, innerthickness, outerlength, innerlength):
    if outerthickness<0 or innerthickness<0:
        raise ValueError("Values must be more than 0")
    if outerlength<7 or innerlength<15:
        raise ValueError("Physically impossible scenario")   #Check this, idk if the top plate must cover the coils and stuff. It current assumes it does
    
    femm.mi_clearselected()
    oldouterlength = tpp1.length
    oldinnerlength = tpp6.length
    innerlengthchange = innerlength - oldinnerlength
    outerlengthchange = outerlength - oldouterlength
    #selects and deletes old nodes
    for i in movablenodes:
        femm.mi_selectnode(i.x, i.y)
    femm.mi_deleteselectednodes()
    #corrects outer length
    for i in (tpp3,tpp4):
        i.move(outerlengthchange, 0)
    #correct outer thickness
    for i in (tpp2, tpp3):
        i.y = outerthickness+24
    #correct inner length
    for i in (tpp7,tpp8):
        i.move(-innerlengthchange, 0)
    #correct inner thickness
    for i in (tpp8, tpp9):
        i.y = innerthickness+24
    #draw new lines and nodes
    for i, j in zip(innerplatenodes[:-1], innerplatenodes[1:]):
        femm.mi_drawline(i.x, i.y, j.x, j.y) 
    for i, j in zip(outerplatenodes[:-1], outerplatenodes[1:]):
        femm.mi_drawline(i.x, i.y, j.x, j.y)
    #update data  
    for i in outerplatenodes:
        i.length = outerlength
        i.thickness = outerthickness
    for i in innerplatenodes:
        i.length = innerlength
        i.thickness = innerthickness
    
def redraw_shielding(shieldthickness, shieldheight, shieldgap):
    '''Redrawst the mag shielding walls given: a thickness, height (from point 1 to 2), and gap (displacement in x from the channel)'''
    #checks inputs are valid
    if shieldthickness<0 or shieldheight<0 or shieldgap<0: #dimensions cant be negative (maybe change for height? small did but i think that would be bad for mag field)
        raise ValueError("Values must be >= 0")
    if shieldthickness+shieldgap > 3.995: #intersects with coil
        raise ValueError("Physically impossible scenario")
    if shieldheight < 8: #i just dont wanna move the isolator lol
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
        i.y = shieldheight-21
    #update all additional info 
    for i in shieldnodes:
        i.swh = shieldheight
        i.swt = shieldthickness
        i.swg = shieldgap
    #channel doesnt move so 29 can be used as a reference point
    sp5.y = (-17-(shieldgap*isolatorratio))
    sp6.y = (-17-(shieldgap*isolatorratio))
    #draw new lines+nodes
    for i, j in zip(shieldnodes[:-1], shieldnodes[1:]): #combines combinations 1 and 2, 2 and 3, 3 and 4...
        femm.mi_drawline(i.x, i.y, j.x, j.y)  
    #finish drawing theuster (isolator and connect to coils)
    femm.mi_addsegment(-38, -21, sp1.x, sp1.y)
    femm.mi_addsegment(-15, -21, sp10.x, sp10.y)
    femm.mi_addsegment(sp4.x, sp4.y, -34, -13)
    femm.mi_addsegment(sp7.x, sp7.y, -19, -13)
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
    femm.mi_clearselected()
    femm.mi_selectlabel(-10, -8)
    femm.mi_setblockprop("24 AWG", 1, 0, 'Inner Coil', 0, 2, innerturns)
    #select and adjust outer coil 
    femm.mi_clearselected()
    femm.mi_selectlabel(-42, -8)
    femm.mi_setblockprop("24 AWG", 1, 0, 'Outer Coil', 0, 2, outerturns)
    
   
########################################ITERATIVE PART########################################
middlerp = (-26.5, -4) #middle reference point
leftrp = (-31, -4) #left reference point
rightrp = (-22, -4) #right reference point

class lowfluxdensity(Exception):
    pass

class lowpenetration(Exception):
    pass

def nd_target(max_value, y, y_peak, constant): #constant = 2c^2, from general normal distribution form 
    '''just the normal distribution, purpose of this function is to scale the target flux density at differeny y coordinates'''
    return (max_value * np.exp(-(y - y_peak)**2 / constant))

#x1, x2, y1, y2, imagine the rectangle formed by these limits
region1 = (-27.5, -25.5, -4, 28.5)
region2L = (-30.95, -30.5, 1, 26)   #CHECK REGION Y LEVELS: mag field curves -> does it go fully in, where should we stop measuring
region2R = (-22.05, -22.5, 1, 26)
region3 = () #we might actually not need this

def get_region1_score(x1,x2,y1,y2):
    '''maximise Br and minimises Bz, returning a dimensionless score'''
    xsegments = np.linspace(x1,x2,3)
    ysegments = np.linspace(y1,y2,50)
    scoreR = 0 
    scoreZ = 0
    for i in xsegments:
        for j in ysegments:
            Br, Bz = femm.mo_getb(i, j)
            target = nd_target(0.03, j, region1[3], 200) #TEST THIS LAST VARIABLE (2c^2 value). thought process for 162: channel width=9 so 2c^2 = 2(9^2)
            scoreR += ((Br-target)/0.03)**2
            scoreZ += (Bz/0.03)**2
    return((scoreR/150 + scoreZ/150)/2)

def get_region2_score(x1,x2,y1,y2):
    '''minimises Br, returning a dimensionless score'''
    xsegments = np.linspace(x1,x2,3)
    ysegments = np.linspace(y1,y2,50)
    scoreR = 0 #idk if we need to maximise Bz here, prob not right? as region1 maximised
    for i in xsegments:
        for j in ysegments:
            Br, Bz = femm.mo_getb(i, j)
            scoreR += (Br/0.03)**2
    return(scoreR/150)

def get_parallel_score(): #if field lines stay parallel to the channel wall
    '''Measures how parallel flux is to the channel and gives a penalty if it isnt perfectly parallel'''
    total = 0
    ycoords = np.linspace(1, 25, 50)
    for i in ycoords:
        Lflux = femm.mo_getb(leftrp[0], i)
        Rflux = femm.mo_getb(rightrp[0], i)
        bonus = ((abs(Lflux[1])/np.hypot(*Lflux)) + (abs(Rflux[1])/np.hypot(*Rflux))) / 2 #find cos from trig but we leave in tersmf of cos so its between 0 and 1
        total += bonus
    return(1 - (total/50))
########################################OPTIMISATION PART########################################
runcount = 0 

#ISSUE: DIFFERENTIAL EVOLUTION CHANGES THE VALUES SO WHEN IT TRIED TO UNDO CHAMFER TI FUCKS UP}

def get_score(x):
    '''where objective function thingy is, lower score = more desirable
    Order for input list: Shield thickness, shield height, shield gap, inner turns, outer turns, current'''
    global runcount
    runcount += 1 
    try: #normally run this, CHANGES DIMENSIONS
        redraw_shielding(x[0], x[1], x[2])
        chamfer_screening(x[10], x[11])
        chamfer_top_plate(x[12], x[13])
        change_turns(x[3], x[4])
        change_current(x[5])
        run_femm()
        
        parallel = get_parallel_score()
        region1_score = get_region1_score(*region1)
        region2_score = (get_region2_score(*region2L) + get_region2_score(*region2R))/2
        
        score = region1_score + region2_score + parallel

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
        print("########################################")
        print("Score: ", score)
        
        undo_chamfers()
        return(score)    
    
    except (ValueError, lowfluxdensity, lowpenetration): #when a valueerror is raised from a function like redraw_shielding() itll do this instead
        print("\n Run", runcount, "rejected")
        return 1e10 #large score so programme will not use this design 

#Order: Shield thickness, shield height, shield gap, inner turns, outer turns, current, outerthickness, innerthickness, outerlength, innerlength, screenchamferwidth, screenchamferlength, platechamferwidth, platechamferlength
bounds = [(0.5, 2.99), (13,44.95), (0.5, 1), (1,10000), (1, 10000), (1,10), (0.5,10), (0.5,10), (7, 10.9), (15,18.9), (0.1,  (abs(sp3.x - sp2.x))), (0.1, (abs(sp1.swh-21))), (0.1, (abs(tpp3.y-tpp4.y))), (0.1,(abs(tpp4.x-tpp5.x)))]   
def optimise():
    results = differential_evolution(get_score, bounds, integrality=(False, False, False, True, True, False, False, False, False, False, False, False, False, False),popsize=10, maxiter=100, polish=False, disp=True)
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


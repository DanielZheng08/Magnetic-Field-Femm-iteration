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

#Chamfer points (points that make up the mag screening as well!)
#CHANGE TS
leftwidth = MagCircuitChamfer("leftwidth", -93.25, 77, 1, 10, "Screen")
leftlength = MagCircuitChamfer("leftlength", -92.25, 67, 1,10, "Screen")
rightwidth = MagCircuitChamfer("rightwidth", -70.75, 77, 1, 10, "Screen")
rightlength = MagCircuitChamfer("rightlength", -69.75, 67, 1, 10, "Screen")

#FOR TOP PLATE 'which' IS "top plate"!!
    
#Shielding points (sp = shielding points) INITIAL DIMENSIONS FROM BASE FEMM MODEL, draw rectangle tool not used incase of future geometry iteration thing
#CORRECTED
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
#CORRECTED
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

def redraw_top_plate(outerthickness, innerthickness, outerlength, innerlength):
    #CORRECTED
    if outerthickness<0 or innerthickness<0:
        raise ValueError("Values must be more than 0")
    if outerlength<6.5 or innerlength<27.5:
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
    if shieldheight < 44.995: #i just dont wanna move the isolator lol
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
        i.y = shieldheight-20
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
middlerp = (-26.5, 14) #middle reference point
leftrp = (-31, 14) #left reference point
rightrp = (-22, 14) #right reference point

#x1,y1,x2,y2 (Top left and Top right)
region1 = (-31, 23.5, -22, 28)
region2L = (-31, -9, -22, 23.5)
region2R = ()
region3 = (-31, -14, -22, -9)

class lowfluxdensity(Exception):
    pass

class lowpenetration(Exception):
    pass

########################################OPTIMISATION PART########################################
runcount = 0 
def get_score(x):
    '''where objective function thingy is, lower score = more desirable
    Order for input list: Shield thickness, shield height, shield gap, inner turns, outer turns, current'''
    global runcount
    runcount += 1 
    try: #normally run this, CHANGES DIMENSIONS
        redraw_shielding(x[0], x[1], x[2])
        redraw_top_plate(x[6], x[7], x[8], x[9])
        change_turns(x[3], x[4])
        change_current(x[5])
        run_femm()
        
        check_max_fd()
        
        region1_score = get_region1_score(*region1)
        region2_score = get_region2_score(*region2L) + get_region2_score(*region2R))/2
        parallel = get_parallel_score()
        depth = get_penetration_depth_score()
        
        penalty = parallel #add others later?
        
        score = region1_score + region2_score + penalty
        
        if parallel < 0.2:
            score += depth
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
        print("Parallel penalty: ", parallel)
        print("Penetration bonus: ", depth)
        if parallel < 0.2:
            print("parallel bonus given")
        else:
            print("no bonus")
        print("########################################")
        print("Score: ", score)
        
        return(score)    
    
    except (ValueError, lowfluxdensity, lowpenetration): #when a valueerror is raised from a function like redraw_shielding() itll do this instead
        print("\n Run", runcount, "rejected")
        return 1e10 #large score so programme will not use this design 

#Order: Shield thickness, shield height, shield gap, inner turns, outer turns, current, outerthickness, innerthickness, outerlength, innerlength
bounds = [(0.5, 6.75), (49.5,72), (0.5, 6.75), (1,10000), (1, 10000), (1,10), (0.5,10), (0.5,10), (6.5, 15), (27.4,34.7)]    
def optimise():
    results = differential_evolution(get_score, bounds, integrality=(False, False, False, True, True, False, False, False, False, False),popsize=10, maxiter=100, polish=False, disp=True)
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


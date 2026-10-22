'''
Daniel Zheng, 17/10/26
'''
################[SET THESE BEFORE STARTING AND KEEP THE SAME]##################
isolatorratio = (1.5) #bottom thickness / side thickness (NOT IN CODE RN), if changed baseline file may need to be manually editted
###############################################################################

import femm
import numpy as np
from scipy.optimize import differential_evolution #not used yet will be used for iterative section later

def open_femm():
    femm.openfemm()
    femm.opendocument('C:\\Uni work\\Electric Propulsions\\IterativeBase - Copy.fem')

def check_nodes():
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
    #checks inputs are valid
    if shieldthickness<0 or shieldheight<0 or shieldgap<0:
        raise ValueError("Values must be >= 0")
    if shieldthickness+shieldgap > 7.4:
        raise ValueError("Physically impossible scenario")
    if shieldheight < 49:
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
    for i, j in zip(shieldnodes[:-1], shieldnodes[1:]):
        femm.mi_drawline(i.x, i.y, j.x, j.y)  
    #finish drawing theuster (isolator and connect to coils)
    femm.mi_addsegment(-99, 5, sp1.x, sp1.y)
    femm.mi_addsegment(-66.3493, 5, sp10.x, sp10.y)
    femm.mi_addsegment(sp4.x, sp4.y, -90.25, 49)
    femm.mi_addsegment(sp7.x, sp7.y, -73.75, 49)
    femm.mi_addsegment(sp1.x, sp1.y, sp10.x, sp10.y)
    
def change_current(newcurrent):
    femm.mi_modifycircprop("Outer Coil", 1, newcurrent)
    femm.mi_modifycircprop("Inner Coil", 1, newcurrent)
    
def change_turns(innerturns, outerturns):
    #select and adjust inner coil
    femm.mi_selectlabel(-62.411, 41.968)
    femm.mi_setblockprop("24 AWG", 1, 0, 'Inner Coil', 0, 2, innerturns)
    #select and adjust outer coil 
    femm.mi_selectlabel(-100.15, 18.60)
    femm.mi_setblockprop("24 AWG", 1, 0, 'Outer Coil', 0, 2, outerturns)

def test123():
    return "hello world"
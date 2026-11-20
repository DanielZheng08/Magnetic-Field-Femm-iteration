'''
Daniel Zheng, 17/10/26
'''
################[SET THESE BEFORE STARTING]##################
filepath = "C:\\Uni work\\Electric Propulsions\\Multiple Sims\\IterativeBase2.fem"
portval=48017
isolatorratio = 1.5
c = 200
targetflux = 0.025
ypeak = 30.5
ymin = -4
#############################################################

import femm
import femm as femm_module
femm = femm_module.FemmSession(port=portval)
import numpy as np
from scipy.optimize import differential_evolution

def open_femm():
    femm.opendocument(filepath)

def run_femm():
    femm.mi_createmesh()
    femm.mi_analyse(0)
    femm.mi_loadsolution()

def check_nodes():
    for i in shieldnodes:
        print(i.__repr__())

class screenNodes():
    def __init__(self, name, x, y, swh, swt, swg):
        self.name = name
        self.x = x  # radial coordinate r
        self.y = y  # axial coordinate z
        self.swh = swh
        self.swt = swt
        self.swg = swg

    def __str__(self):
        return f"{self.name}, r = {self.x}, z = {self.y}"

    def __repr__(self):
        return (f"screenNode, {self.name}, r = {self.x}, z = {self.y}, "
                f"height = {self.swh}, thickness = {self.swt}, gap = {self.swg}")

    def move(self, dx, dy):
        self.x += dx
        self.y += dy

class PlateNodes():
    def __init__(self, name, x, y, in_out, thickness, length, oldx, oldy):
        self.name = name
        self.x = x
        self.y = y
        self.in_out = in_out
        self.thickness = thickness
        self.length = length
        self.oldx = oldx
        self.oldy = oldy

    def __str__(self):
        return f"{self.name}, r = {self.x}, z = {self.y}, {self.in_out}"

    def __repr__(self):
        return (f"{self.name}, r = {self.x}, z = {self.y}, {self.in_out}, "
                f"thickness = {self.thickness}, length = {self.length}")

    def move(self, dx, dy):
        self.oldx = self.x
        self.oldy = self.y
        self.x += dx
        self.y += dy

class MagCircuitChamfer():
    def __init__(self, name, x, y, chamferwidth, chamferlength, which):
        self.name = name
        self.x = x
        self.y = y
        self.chamferwidth = chamferwidth
        self.chamferlength = chamferlength
        self.which = which

    def move(self, dx, dy):
        self.x += dx
        self.y += dy

#screening points
sp1 = screenNodes("sp1", 37.995, -21, 42, 3, 1)
sp2 = screenNodes("sp2", 37.995,  21, 42, 3, 1)
sp3 = screenNodes("sp3", 35.000,  21, 42, 3, 1)
sp4 = screenNodes("sp4", 35.000, -13, 42, 3, 1)
sp5 = screenNodes("sp5", 35.000, -18.5, 42, 3, 1)
sp6 = screenNodes("sp6", 18.000, -18.5, 42, 3, 1)
sp7 = screenNodes("sp7", 18.000, -13, 42, 3, 1)
sp8 = screenNodes("sp8", 18.000,  21, 42, 3, 1)
sp9 = screenNodes("sp9", 15.005,  21, 42, 3, 1)
sp10 = screenNodes("sp10", 15.005, -21, 42, 3, 1)
shieldnodes = [sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, sp9, sp10]

#top plate points
tpp1 = PlateNodes("tpp1", 45, 24, "out", 2, 10, 45, 24)
tpp2 = PlateNodes("tpp2", 45, 26, "out", 2, 10, 45, 26)
tpp3 = PlateNodes("tpp3", 35, 26, "out", 2, 10, 35, 26)
tpp4 = PlateNodes("tpp4", 35, 24, "out", 2, 10, 35, 24)
tpp5 = PlateNodes("tpp5", 38, 24, "out", 2, 10, 38, 24)
tpp6 = PlateNodes("tpp6", 15, 24, "in", 2, 18, 15, 24)
tpp7 = PlateNodes("tpp7", 18, 24, "in", 2, 18, 18, 24)
tpp8 = PlateNodes("tpp8", 18, 26, "in", 2, 18, 18, 26)
tpp9 = PlateNodes("tpp9", 0, 26, "in", 2, 18, 0, 26)
tpp10 = PlateNodes("tpp10", 0, 24, "in", 2, 18, 0, 24)
movablenodes = [tpp2, tpp3, tpp4, tpp7, tpp8, tpp9]
outerplatenodes = [tpp1, tpp2, tpp3, tpp4, tpp5]
innerplatenodes = [tpp6, tpp7, tpp8, tpp9, tpp10]

#plate chamfer points
Pleftwidth = MagCircuitChamfer("leftwidth", 36, 21, 1, 10, "Plate")
Pleftlength = MagCircuitChamfer("leftlength", 38, 11, 1, 10, "Plate")
Prightwidth = MagCircuitChamfer("rightwidth", 17, 21, 1, 10, "Plate")
Prightlength = MagCircuitChamfer("rightlength", 15, 11, 1, 10, "Plate")
platechamfernodes = [Pleftwidth, Pleftlength, Prightwidth, Prightlength]

#screen chamfer points
Sleftwidth = MagCircuitChamfer("leftwidth", 36, 21, 1, 10, "Screen")
Sleftlength = MagCircuitChamfer("leftlength", 38, 11, 1, 10, "Screen")
Srightwidth = MagCircuitChamfer("rightwidth", 17, 21, 1, 10, "Screen")
Srightlength = MagCircuitChamfer("rightlength", 15, 11, 1, 10, "Screen")
screenchamfernodes = [Sleftwidth, Sleftlength, Srightwidth, Srightlength]

def chamfer_screening(width, length):
    max_width = min(abs(sp3.x - sp2.x), abs(sp9.x - sp8.x)) - 0.01
    max_length = min(abs(sp2.y - sp1.y), abs(sp9.y - sp10.y)) - 0.01

    if width <= 0 or length <= 0:
        raise ValueError("Screen chamfer dimensions must be positive")
    if width > max_width or length > max_length:
        raise ValueError("Screen chamfer too deep")

    femm.mi_clearselected()

    Sleftwidth.x = sp3.x + width
    Sleftwidth.y = sp3.y
    Sleftlength.x = sp2.x
    Sleftlength.y = sp2.y - length

    Srightwidth.x = sp8.x - width
    Srightwidth.y = sp8.y
    Srightlength.x = sp9.x
    Srightlength.y = sp9.y - length

    femm.mi_selectnode(sp2.x, sp2.y)
    femm.mi_selectnode(sp9.x, sp9.y)
    femm.mi_deleteselectednodes()

    for i in screenchamfernodes:
        femm.mi_addnode(i.x, i.y)

    femm.mi_addsegment(sp1.x, sp1.y, Sleftlength.x, Sleftlength.y)
    femm.mi_addsegment(Sleftlength.x, Sleftlength.y, Sleftwidth.x, Sleftwidth.y)
    femm.mi_addsegment(Sleftwidth.x, Sleftwidth.y, sp3.x, sp3.y)

    femm.mi_addsegment(sp8.x, sp8.y, Srightwidth.x, Srightwidth.y)
    femm.mi_addsegment(Srightwidth.x, Srightwidth.y, Srightlength.x, Srightlength.y)
    femm.mi_addsegment(Srightlength.x, Srightlength.y, sp10.x, sp10.y)

    femm.mi_clearselected()

def chamfer_top_plate(width, length):
    if width <= 0 or length <= 0:
        raise ValueError("Plate chamfer dimensions must be positive")

    max_width = min(abs(tpp3.y - tpp4.y), abs(tpp8.y - tpp7.y)) - 0.01
    max_length = min(abs(tpp4.x - tpp5.x), abs(tpp6.x - tpp7.x)) - 0.01

    if width > max_width or length > max_length:
        raise ValueError("Plate chamfer too deep")

    femm.mi_clearselected()

    Pleftwidth.x = tpp3.x
    Pleftwidth.y = tpp3.y - width
    Pleftlength.x = tpp4.x + length
    Pleftlength.y = tpp4.y

    Prightwidth.x = tpp8.x
    Prightwidth.y = tpp8.y - width
    Prightlength.x = tpp7.x - length
    Prightlength.y = tpp7.y

    femm.mi_selectnode(tpp4.x, tpp4.y)
    femm.mi_selectnode(tpp7.x, tpp7.y)
    femm.mi_deleteselectednodes()

    for node in platechamfernodes:
        femm.mi_addnode(node.x, node.y)

    # Outer plate chamfer.
    femm.mi_addsegment(tpp3.x, tpp3.y, Pleftwidth.x, Pleftwidth.y)
    femm.mi_addsegment(Pleftwidth.x, Pleftwidth.y, Pleftlength.x, Pleftlength.y)
    femm.mi_addsegment(Pleftlength.x, Pleftlength.y, tpp5.x, tpp5.y)

    # Inner plate chamfer (opposite radial direction to outer plate).
    femm.mi_addsegment(tpp6.x, tpp6.y, Prightlength.x, Prightlength.y)
    femm.mi_addsegment(Prightlength.x, Prightlength.y, Prightwidth.x, Prightwidth.y)
    femm.mi_addsegment(Prightwidth.x, Prightwidth.y, tpp8.x, tpp8.y)

    femm.mi_clearselected()

def undo_screen_chamfers():
    femm.mi_clearselected()
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

def undo_plate_chamfers():
    femm.mi_clearselected()
    femm.mi_selectnode(Pleftwidth.x, Pleftwidth.y)
    femm.mi_selectnode(Pleftlength.x, Pleftlength.y)
    femm.mi_selectnode(Prightwidth.x, Prightwidth.y)
    femm.mi_selectnode(Prightlength.x, Prightlength.y)
    femm.mi_deleteselectednodes()

    femm.mi_addnode(tpp4.x, tpp4.y)
    femm.mi_addnode(tpp7.x, tpp7.y)

    femm.mi_addsegment(tpp3.x, tpp3.y, tpp4.x, tpp4.y)
    femm.mi_addsegment(tpp4.x, tpp4.y, tpp5.x, tpp5.y)
    femm.mi_addsegment(tpp6.x, tpp6.y, tpp7.x, tpp7.y)
    femm.mi_addsegment(tpp7.x, tpp7.y, tpp8.x, tpp8.y)

    femm.mi_clearselected()

def redraw_top_plate(outerthickness, innerthickness, outerlength, innerlength):
    if outerthickness < 0 or innerthickness < 0:
        raise ValueError("Values must be more than 0")
    if outerlength < 7.01 or innerlength < 15.01:
        raise ValueError("Physically impossible scenario, length too small")
    if outerlength > 10.99 or innerlength > 18.99:
        raise ValueError("Physically impossible scenario, length too large")

    femm.mi_clearselected()
    oldouterlength = tpp1.length
    oldinnerlength = tpp6.length
    innerlengthchange = innerlength - oldinnerlength
    outerlengthchange = outerlength - oldouterlength

    for i in movablenodes:
        femm.mi_selectnode(i.x, i.y)
    femm.mi_deleteselectednodes()

    for i in (tpp3, tpp4):
        i.move(-outerlengthchange, 0)

    for i in (tpp2, tpp3):
        i.y = outerthickness + 24

    for i in (tpp7, tpp8):
        i.move(innerlengthchange, 0)

    for i in (tpp8, tpp9):
        i.y = innerthickness + 24

    for i in movablenodes:
        femm.mi_addnode(i.x, i.y)

    for i, j in zip(innerplatenodes[:-1], innerplatenodes[1:]):
        femm.mi_addsegment(i.x, i.y, j.x, j.y)
    for i, j in zip(outerplatenodes[:-1], outerplatenodes[1:]):
        femm.mi_addsegment(i.x, i.y, j.x, j.y)

    for i in outerplatenodes:
        i.length = outerlength
        i.thickness = outerthickness
    for i in innerplatenodes:
        i.length = innerlength
        i.thickness = innerthickness

    femm.mi_clearselected()

def redraw_screen(shieldthickness, shieldheight, shieldgap):
    if shieldthickness < 0 or shieldheight < 0 or shieldgap < 0:
        raise ValueError("Values must be >= 0")
    if shieldthickness + shieldgap > 3.995:
        raise ValueError("Physically impossible scenario")
    if shieldheight < 8:
        raise ValueError("Shield height must be more than 8")

    oldthickness = sp1.swt
    oldgap = sp1.swg
    gapchange = shieldgap - oldgap
    thicknesschange = shieldthickness - oldthickness

    femm.mi_clearselected()
    for i in shieldnodes:
        femm.mi_selectnode(i.x, i.y)
    femm.mi_deleteselectednodes()

    for i in (sp1, sp2, sp3, sp4, sp5):
        i.move(gapchange, 0)

    for i in (sp6, sp7, sp8, sp9, sp10):
        i.move(-gapchange, 0)

    for i in (sp1, sp2):
        i.move(thicknesschange, 0)

    for i in (sp9, sp10):
        i.move(-thicknesschange, 0)

    for i in (sp2, sp3, sp8, sp9):
        i.y = shieldheight - 21

    for i in shieldnodes:
        i.swh = shieldheight
        i.swt = shieldthickness
        i.swg = shieldgap

    sp5.y = -17 - shieldgap * isolatorratio
    sp6.y = -17 - shieldgap * isolatorratio

    for i in shieldnodes:
        femm.mi_addnode(i.x, i.y)

    for i, j in zip(shieldnodes[:-1], shieldnodes[1:]):
        femm.mi_addsegment(i.x, i.y, j.x, j.y)

    femm.mi_addsegment(38, -21, sp1.x, sp1.y)
    femm.mi_addsegment(15, -21, sp10.x, sp10.y)
    femm.mi_addsegment(sp4.x, sp4.y, 34, -13)
    femm.mi_addsegment(sp7.x, sp7.y, 19, -13)
    femm.mi_addsegment(sp1.x, sp1.y, sp10.x, sp10.y)
    femm.mi_clearselected()

def change_current(newcurrent):
    if newcurrent < 0:
        raise ValueError("Value must be more than 0")
    femm.mi_clearselected()
    femm.mi_modifycircprop("Outer Coil", 1, newcurrent)
    femm.mi_modifycircprop("Inner Coil", 1, newcurrent)
    femm.mi_clearselected()

def change_turns(innerturnsaddition, outerturns):
    if innerturnsaddition < 0 or outerturns < 0:
        raise ValueError("Value must be more than 0")

    femm.mi_clearselected()
    femm.mi_selectlabel(10, -8)
    femm.mi_setblockprop("24 AWG", 1, 0, "Inner Coil", 0, 2,
                         innerturnsaddition + outerturns)

    femm.mi_clearselected()
    femm.mi_selectlabel(42, -8)
    femm.mi_setblockprop("24 AWG", 1, 0, "Outer Coil", 0, 2, outerturns)
    femm.mi_clearselected()

########################################ITERATIVE PART########################################
middlerp = (26.5, -4)
leftrp = (31, -4)   # outer wall
rightrp = (22, -4)  # inner wall

class lowfluxdensity(Exception):
    pass

class lowpenetration(Exception):
    pass

def nd_target(max_value, y, y_peak, constant):
    return max_value * np.exp(-(y - y_peak) ** 2 / constant)

region1 = (25.5, 27.5, ymin, ypeak)
region2L = (22.05, 22.55, 1, 23.5)
region2R = (30.95, 39.9, 1, 23.5)

def get_uniformity_score(y1, y2):
    ysegments = np.linspace(y1, y2, 20)
    difference = 0
    for i in ysegments:
        BrL, BzL = femm.mo_getb(region2L[0], i)
        BrR, BzR = femm.mo_getb(region2R[0], i)
        difference += ((BrL + BrR) / targetflux) ** 2
        difference += ((BzL - BzR) / targetflux) ** 2
    return difference / (2 * len(ysegments))

def get_region1_score(x1, x2, y1, y2):
    xsegments = np.linspace(x1, x2, 3)
    ysegments = np.linspace(y1, y2, 50)
    scoreR = 0
    scoreZ = 0
    for i in xsegments:
        for j in ysegments:
            Br, Bz = femm.mo_getb(i, j)
            target = nd_target(targetflux, j, region1[3], c)
            scoreR += ((Br - target) / targetflux) ** 2
            scoreZ += (Bz / targetflux) ** 2
    return (scoreR / 150 + scoreZ / 150) / 2

def get_region2_score(x1, x2, y1, y2):
    xsegments = np.linspace(x1, x2, 3)
    ysegments = np.linspace(y1, y2, 50)
    scoreR = 0
    for i in xsegments:
        for j in ysegments:
            Br, Bz = femm.mo_getb(i, j)
            scoreR += (Br / targetflux) ** 2
    return scoreR / 150

def get_region3_score(x1, x2, y1, y2):
    xsegments = np.linspace(x1, x2, 3)
    ysegments = np.linspace(y1, y2, 50)
    scoreR = 0
    for i in xsegments:
        for j in ysegments:
            Br, Bz = femm.mo_getb(i, j)
            scoreR += (Br / targetflux) ** 2
    return scoreR / 150
    

def check_max_fd():
    ycoords = np.linspace(22, ypeak, 10)
    values = []
    for i in ycoords:
        Br, Bz = femm.mo_getb(26.5, i)
        values.append(abs(Br))
    if max(values) < 0.02:
        raise lowfluxdensity()

def get_penetration_depth_score():
    y_exit = 26
    y_anode = ymin
    minimum_Br = 0.005  # 5 mT
    ycoords = np.linspace(y_exit, y_anode, 50)
    depth = 0
    for i in ycoords:
        Br, Bz = femm.mo_getb(middlerp[0], i)
        if abs(Br) >= minimum_Br:
            depth = y_exit - i
        else:
            break
    return 1 - depth / (y_exit - y_anode)

########################################OPTIMISATION PART########################################
runcount = 0
region1 = (25.5, 27.5, ymin, ypeak)
region2L = (22.05, 22.55, 1, 23.5)
region2R = (30.95, 39.9, 1, 23.5)
region3L = (tpp8.x, 19, tpp8.y, Prightwidth.y) #LEFT SIDE
region3R = (34, tpp3.x, tpp3.y, Pleftwidth.y) #RIGHT SIDE

def get_score(x):
    plate_chamfered = False
    screen_chamfered = False
    global runcount
    runcount += 1

    try:
        redraw_screen(x[0], x[1], x[2])
        # redraw_top_plate(x[10], x[11], x[12], x[13])  #disabled for now, too mych variables
        chamfer_screening(x[6], x[7])
        screen_chamfered = True
        chamfer_top_plate(x[8], x[9])
        plate_chamfered = True
        change_turns(x[3], x[4])
        change_current(x[5])

        run_femm()
        check_max_fd()

        region1_score = get_region1_score(*region1)
        inner_wall_score = abs(get_region2_score(*region2L))
        outer_wall_score = abs(get_region2_score(*region2R))
        region2_score = max(inner_wall_score, outer_wall_score)
        inner_grazing_score = abs(get_region3_score(*region3L))
        outer_grazing_score = abs(get_region3_score(*region3R))
        region3_score = max(inner_grazing_score, outer_grazing_score)
        #bonus+penaltys
        depth = get_penetration_depth_score()
        uniformity = get_uniformity_score(23.5, 10)

        #score
        score = region1_score + region2_score + region3_score + depth
        
        print("########################################")
        print("Run:", runcount)
        print("########################################")
        print("Shield thickness:", sp1.swt)
        print("Shield height:", sp1.swh)
        print("Shield gap:", sp1.swg)
        print("Inner coil turns:", x[3] + x[4])
        print("Outer coil turns:", x[4])
        print("Current:", x[5])
        print("Screen Chamfer   =", x[6], " x ", x[7])
        print("Plate Chamfer    =", x[8], " x ", x[9])
        print("########################################")
        print("Region 1 score:", region1_score)
        print("Inner wall score:", inner_wall_score)
        print("Outer wall score:", outer_wall_score)
        print("Region 2 score:", region2_score)
        print("Depth score:", depth)
        print("Uniformity score:", uniformity)
        print("########################################")
        print("Total Score:", score)
        return score

    except (ValueError, lowfluxdensity, lowpenetration) as e:
        print("Run", runcount, "rejected:", repr(e))
        return 1e10

    finally:
        if plate_chamfered:
            undo_plate_chamfers()
        if screen_chamfered:
            undo_screen_chamfers()

#order: shield thickness, shield height, shield gap, inner turns addition, outer turns, current, screen chamfer width/length, plate chamfer width/length.
bounds = [(0.5, 2.99), (30, 44.95), (0.1, 1), (1, 500), (100, 500), (1, 10), (0.1, abs(sp3.x - sp2.x) - 0.01), (0.1, abs(sp1.swh - 21)), (0.1, abs(tpp3.y - tpp4.y) - 0.01), (0.1, abs(tpp4.x - tpp5.x) - 0.01)]

def optimise():
    results = differential_evolution( get_score, bounds, strategy='rand1bin', integrality=(False, False, False, True, True, False, False, False, False, False), popsize=15, maxiter=150, polish=True, disp=True)

    print("################################### Optimisation Finished ###################################")
    print("Nodes:")
    check_nodes()
    print("#############################################################################################")
    print("Shield thickness =", results.x[0])
    print("Shield height    =", results.x[1])
    print("Shield gap       =", results.x[2])
    print("Inner turns      =", results.x[3] + results.x[4])
    print("Outer turns      =", results.x[4])
    print("Current          =", results.x[5])
    print("Screen Chamfer   =", results.x[6], " x ", results.x[7])
    print("Plate Chamfer    =", results.x[8], " x ", results.x[9])
    print("#############################################################################################")
    print("Best score: ", results.fun)
    redraw_screen(results.x[0], results.x[1], results.x[2])

    chamfer_screening(results.x[6], results.x[7])
    chamfer_top_plate(results.x[8], results.x[9])

    change_turns(results.x[3], results.x[4])
    change_current(results.x[5])

import io
import os
from tkinter import *
import datetime
import random
from matplotlib.ticker import FuncFormatter
import math
from numba import prange
import decimal
import numpy
import pandas as pd
import json
import matplotlib.pyplot as plt
from PIL import Image
import ghostscript
import moviepy.editor as mp

continuedVer = '1.18_3'  # Number of the set of initial conditions and settings

# Technical model introductions
nol0 = decimal.Decimal('0')  # Zero of type Decimal
continued = False  # Renewable mode for distributed computing
saveGif = True;
showVisualDelay = 100;  # Skipping steps for the next visualization step
unlimetedSteps = True; # Compute until all available model state transition probabilities are exhausted
unlimetedLimits = False; # infinite lattice allowed
averBoardHLimit = True; # surround the workspace with border cells with "H"
startIcellsFromCenter = True; # Fill the initial state of the grid around the geometric center
t0 = decimal.Decimal('5') # Initial time
xlimits = [0, 60 * 24 * 3] # Limits of X-axis
T = 300; # Limit number of steps (when unlimetedSteps = False)


# Input data
X = 30;
Y = 30;
N_I = X * Y * 0.10 # Initial number of I-cells to fill (total number of cells multiplied by I fraction)
N_D = X * Y * 0.000
N_F = X * Y * 0.0
k1 = decimal.Decimal('0') / 90  # H->I
k1minus = decimal.Decimal('0') / 90  # I->H
k2 = decimal.Decimal('0.0033333333333333335')  # I->D
k4 = decimal.Decimal('0.0')  # I->F
k4minus = decimal.Decimal('0.0')  # F->I
k5 = decimal.Decimal('0.0')  # F->H
k7 = decimal.Decimal('0.08888888888888889')  # IH->HH
k8 = decimal.Decimal('0.044444444444444446')  # HI->II
k9 = decimal.Decimal('0.08888888888888889')  # HD->ID
k10 = decimal.Decimal('0.0033333333333333335')  # ID->DD
k11 = decimal.Decimal('0.0033333333333333335')  # II->DI



# global variables of the model implementation
R = nol0 # global variable with R
out_data_array = [] # global variable with output data
speeds_dict = {} # global variable with process speeds
images = [] # global variable - list of frames for writing the visualization
for i in prange(Y):
    for j in prange(X):
        speeds_dict[f'{i}_{j}'] = {'i': i, 'j': j}


# Saving the imaging frame
def addCanvasToArr(canvas):
    global images

    # producing postscript image
    canvasImage = canvas.postscript(colormode='color')

    # use PIL to convert to image
    img = Image.open(io.BytesIO(canvasImage.encode('utf-8')))
    images.append(img)


# Initializing canvas
def create_TK(width, height):
    root = Tk()
    canvas = Canvas(root, width=width, height=height)
    canvas.pack()
    return root, canvas


# Building the visualization grid
def build_board(canvas, st, X, Y, x=False, y=False):
    fill = '#FECD72'
    outline = '#825100'
    if x and y:
        canvas.create_rectangle(x * st, y * st, x * st + st, y * st + st, fill=fill, outline=outline)

    for i in range(0, X):
        for j in range(0, Y):
            canvas.create_rectangle(i * st, j * st, i * st + st, j * st + st, fill=fill, outline=outline)


# Dictionary generation with output data. Visualization on the grid
def prepareOutData_1step(root=False, canvas=False, st=False, X=False, Y=False, board=False, visibable=False):
    all = X * Y
    fetaD = 0
    fetaI = 0
    fetaH = 0
    fetaF = 0
    minor = st * 0.95
    If visibable:
        canvas.delete("all")
        build_board(canvas, st, X, Y)
    # I is red, H is blue, D is black, F is green


    outline = '#000'
    for i in range(Y):
        for j in range(X):
            value = board[i][j]
            if value == "H":
                color = 'blue';fetaH = fetaH + 1
            elif value == "D":
                color = 'black';
                fetaD = fetaD + 1
            elif value == "F":
                color = 'green';fetaF = fetaF + 1
            elif value == "I":
                color = 'red';fetaI = fetaI + 1

            if visibable:
                x1, y1, x2, y2 = j * st + minor, i * st + minor, j * st + st - minor, i * st + st - minor
                canvas.create_oval(x1, y1, x2, y2, fill=color, outline=outline)
    if visibable:
        print(f'H+I: {(fetaD + fetaI) / all}')
        root.update()
    fetas = {"fetaH": fetaH / all, "fetaI": fetaI / all, "fetaD": fetaD / all, "fetaF": fetaF / all}
    return fetas


# Filling cells for the initial state of the model grid
def start_status():
    board = [['H' for i2 in prange(X)] for i in prange(Y)]
    I = 0
    D = 0
    F = 0
    if startIcellsFromCenter:
        xCenter = X / 2
        yCenter = Y / 2
        radius = math.sqrt(N_I / math.pi)
        pointsInRadius = list()
        for x in range(0, X):
            for y in range(0, Y):
                if math.sqrt((x - xCenter) ** 2 + (y - yCenter) ** 2) < radius:
                    pointsInRadius.append({"x": x, "y": y})
    while I < N_I:
        if startIcellsFromCenter:
            if not len(pointsInRadius):
                break
            onePoint = pointsInRadius.pop()
            x = onePoint['x'];
            y = onePoint['y'];
        else:
            x = random.randint(0, X - 1)
            y = random.randint(0, Y - 1)
        if board[y][x] == "H":
            I = I + 1
            board[y][x] = "I"
    while (D < N_D) and ((I + D) < X * Y):

        x = random.randint(0, X - 1)
        y = random.randint(0, Y - 1)
        if board[y][x] == "H":
            D = D + 1
            board[y][x] = "D"
    while (F < N_F) and ((I + D + F) < X * Y):
        x = random.randint(0, X - 1)
        y = random.randint(0, Y - 1)
        if board[y][x] == "H":
            F = F + 1
            board[y][x] = "F"
    return board


# boundary conditions for 2x2
def get_event_2x2_limits(x, y, dx, dy):
    x2 = x + dx
    y2 = y + dy
    if unlimetedLimits:
        # infinite lattice
        if x2 == X:
            x2 = 0;
        elif y2 == Y:
            y2 = 0;
        return x, y, dx, dy, x2, y2
    else:
        # Boundary conditions.
        if x2 == X:
            x2 = False;
        elif y2 == Y:
            y2 = False;
        return x, y, dx, dy, x2, y2


# A study of 1 scenario with 2 nodes in 2 directions
def get_event_2x2(x, y, dx, dy, board):
    ret = {'speed': nol0}
    x, y, dx, dy, x2, y2 = get_event_2x2_limits(x, y, dx, dy)
    ret['y'] = y;
    ret['x'] = x;
    if y2 and x2:
        secondCell = board[y2][x2]
        ret['y2'] = y2;
        ret['x2'] = x2;
    else:
        if averBoardHLimit:
            secondCell = "H"
        else:
            return []
    if board[y][x] == 'I':
        if secondCell == 'I':
            ret['speed'] = k11  #
            ret['yx'] = "D"  #
            if y2 and x2:
                ret['y2x2'] = "I"  #
            return [ret]
    return []


# boundary conditions for 2x4
def get_event_2x4_limits(x, y, dx, dy):
    x2 = x + dx
    y2 = y + dy
    if unlimetedLimits:
        if x2 == X:
            x2 = 0;
        elif y2 == Y:
            y2 = 0;
        elif x2 < 0:
            x2 = X - 1;
        elif y2 < 0:
            y2 = Y - 1;
        return x, y, dx, dy, x2, y2
    else:
        if x2 == X:
            x2 = False;
        elif y2 == Y:
            y2 = False;
        elif x2 < 0:
            x2 = False;
        elif y2 < 0:
            y2 = False;
        return x, y, dx, dy, x2, y2


# Study of 1 scenario with 2 nodes on 4 sides
def get_event_2x4(x, y, dx, dy, board):
    ret = {'speed': nol0}
    x, y, dx, dy, x2, y2 = get_event_2x4_limits(x, y, dx, dy)
    ret['y'] = y;
    ret['x'] = x;
    if y2 and x2:
        secondCell = board[y2][x2]
        ret['y2'] = y2;
        ret['x2'] = x2;
    else:
        if averBoardHLimit:
            secondCell = "H"
        else:
            return []
    if board[y][x] == 'H':
        if secondCell == 'I':
            ret['speed'] = k8  
            ret['yx'] = "I"  
            if y2 and x2:
                ret['y2x2'] = "I"  
            return [ret]
        if secondCell == 'D':
            ret['speed'] = k9  
            ret['yx'] = "I"  
            if y2 and x2:
                ret['y2x2'] = "D"  
            return [ret]
    if board[y][x] == 'I':
        if secondCell == 'H':
            ret['speed'] = k7  
            ret['yx'] = "H"  
            if y2 and x2:
                ret['y2x2'] = "H"  
            return [ret]
        if secondCell == 'D':
            ret['speed'] = k10  
            ret['yx'] = "D"  
            if y2 and x2:
                ret['y2x2'] = "D"  
            return [ret]
    return []


# Calculation of R
def get_r_event(R):
    Ep_minus1 = nol0
    E = decimal.Decimal(str(numpy.random.uniform()))

    for key in speeds_dict:
        if speeds_dict[key]['R_'] == 0: continue
        position_events = speeds_dict[key]['events_']
        for ev in position_events:
            ER = E * R
            if (ER > Ep_minus1) and ((ER) <= (Ep_minus1 + ev['speed'])):
                return ev
            Ep_minus1 = Ep_minus1 + ev['speed']
    return False


# get all possible events for one node of the lattice
def get_event_1(x, y, board):
    ret = []
    if board[y][x] == 'H':
        ret.append({"yx": "I", "speed": k1, "x": x, "y": y});

    elif board[y][x] == 'I':
        ret.append({"yx": "H", "speed": k1minus, "x": x, "y": y})
        ret.append({"yx": "D", "speed": k2, "x": x, "y": y})
        ret.append({"yx": "F", "speed": k4, "x": x, "y": y});

    elif board[y][x] == 'F':
        ret.append({"yx": "H", "speed": k5, "x": x, "y": y});
        ret.append({"yx": "I", "speed": k4minus, "x": x, "y": y});
    return ret


# Calculations of elementary processes and interactions for 1 lattice cell
def inv_1_point(i, j, board):
    R_ = nol0;
    events = []
    # generation of velocities for the 4 directions added to the list of events and calculation of R
    # we check for all 4 directions=> reverse reactions are taken into account when relative to the next point
    left = get_event_2x4(x=j, y=i, dx=-1, dy=0, board=board)
    right = get_event_2x4(x=j, y=i, dx=+1, dy=0, board=board)
    up = get_event_2x4(x=j, y=i, dx=0, dy=-1, board=board)
    down = get_event_2x4(x=j, y=i, dx=0, dy=+1, board=board)

    # nodes are the same - only 1 time should be taken into account
    down2 = get_event_2x2(x=j, y=i, dx=0, dy=+1, board=board)
    right2 = get_event_2x2(x=j, y=i, dx=+1, dy=0, board=board)

    # node out of touch with the neighboring nodes
    this_point = get_event_1(x=j, y=i, board=board)

    for i in left: R_ = R_ + i['speed']
    for i in right: R_ = R_ + i['speed']
    for i in right2: R_ = R_ + i['speed']
    for i in up: R_ = R_ + i['speed']
    for i in down: R_ = R_ + i['speed']
    for i in down2: R_ = R_ + i['speed']
    for i in this_point: R_ = R_ + i['speed']

    # add events
    events = events + left
    events = events + right
    events = events + up
    events = events + down
    events = events + down2
    events = events + right2
    events = events + this_point
    return events, R_


# change the board to a new state
def change_board(board, t, changed_points):
    global speeds_dict
    global R
    # step 2
    # Calculation of speeds of elementary events. At the current time t1 and calculate the total speed R
    events = []

    # Acceleration of the algorithm - working only with changed grid cells and their environment
    if changed_points:
        for point in changed_points:
            i = point['y']
            j = point['x']
            events_, R_ = inv_1_point(i, j, board)
            events = events + events_
            R = R - speeds_dict[f'{i}_{j}']['R_']
            speeds_dict[f'{i}_{j}']['R_'] = R_;
            speeds_dict[f'{i}_{j}']['events_'] = events_;
            R = R + speeds_dict[f'{i}_{j}']['R_']

    if not changed_points:
        for i in prange(Y):
            for j in prange(X):
                events_, R_ = inv_1_point(i, j, board)
                events = events + events_
                R = R + R_
                speeds_dict[f'{i}_{j}']['R_'] = R_;
                speeds_dict[f'{i}_{j}']['events_'] = events_;

    # Step 3 One of the possible elementary events is randomly chosen with a probability proportional to its speed.
    # Changes the state of the grid according to the chosen event
    ev = get_r_event(R=R)
    # Step 4 Calculates the time step. The time moment t2 of the system exit is calculated
    # from the current state: t2=t1-ln(E)/R, where E is a random variable uniformly distributed on the interval (0,1).
    E = decimal.Decimal(str(random.random()))

    # if there are no available grid states to change, then the algorithm terminates
    if R == 0:
        return False
    dt = decimal.Decimal(str(math.log(1 / E)) / R
    t = t + dt
    try:
        board[ev['y']][ev['x']] = ev['yx']]
    except:
        return False


    # Get the closest cell environment
    def get_near_points(x, y):
        x1 = x + 1
        y1 = y + 0
        x2 = x - 1
        y2 = y + 0
        x3 = x + 0
        y3 = y + 1
        x4 = x + 0
        y4 = y - 1
        if x1 == X: x1 = 0
        if x2 < 0: x2 = X - 1
        if y3 == Y: y3 = 0;
        if y4 < 0: y4 = Y - 1
        return [{'x': x1, 'y': y1}, {'x': x2, 'y': y2}, {'x': x3, 'y': y3}, {'x': x4, 'y': y4}]


    # Get only changed points
    def get_changed_points(ev):
        points = []
        points = points + get_near_points(x=ev['x'], y=ev['y'])
        if ev.get('y2x2'):
            points = points + get_near_points(x=ev['x2'], y=ev['y2'])
        else:
            points = points + [{'x': ev['x'], 'y': ev['y']}]
        return points

    if ev.get('y2x2'):
        board[ev['y2']][ev['x2']] = ev['y2x2']]

    changed_points = get_changed_points(ev)
    return board, t, changed_points

# Visualize on the graph at the end of the program and write to the data repository
def showGraph(FetaD_array, filename_suffix=""):
    def prepareConfigsToSave():
        configs = {"saveGif": saveGif, "continued": continued, "showVisualDelay": showVisualDelay,
                   "unlimetedSteps": unlimetedSteps,
                   "unlimetedLimits": unlimetedLimits, "averBoardHLimit": averBoardHLimit,
                   "startIcellsFromCenter": startIcellsFromCenter, "t0": t0, "continuedVer": continuedVer,
                   "xlimits": xlimits, "X": X, "Y": Y, "T": T, "N_I": N_I, "N_D": N_D, "N_F": N_F, "k1": k1,
                   "k1minus": k1minus, "k2": k2, "k4": k4, "k4minus": k4minus, "k5": k5, "k7": k7, "k8": k8, "k9": k9,
                   "k10": k10, "k11": k11}
        for key in configs:
            if isinstance(configs[key], decimal.Decimal):
                configs[key] = float(configs[key])

        return json.dumps(configs)

    def formatOx(x, pos):
        delta = datetime.timedelta(minutes=float(x))
        deltaStr = str(delta).replace(" days, ", 'd')
        return deltaStr

    dF = pd.DataFrame(FetaD_array)
    if not os.path.isdir('out'):
        os.mkdir("out")
    with open(f'out/outData_{continuedVer}{filename_suffix}.csv', 'w') as ff:
        ff.write(prepareConfigsToSave() + '\n')
    dF.to_csv(f'out/outData_{continuedVer}{filename_suffix}.csv', index=False, mode='a')
    fig = plt.figure(figsize=(18, 6), dpi=200)
    yD = dF['fetaD'].tolist()
    yF = dF['fetaF'].tolist()
    yI = dF['fetaI'].tolist()
    yH = dF['fetaH'].tolist()
    yID = (dF['fetaI'] + dF['fetaD']).tolist()
    x = dF['t'].tolist()

    ax = fig.add_subplot(111)
    ax.plot(x, yD, color="black", label='D')
    ax.plot(x, yI, color="red", label='I')
    ax.plot(x, yH, color="blue", label='H')
    ax.plot(x, yF, color="green", label='F')
    ax.plot(x, yID, '.', color="orange", label='I+D')
    ax.grid()
    ax.set_title(f"t(0)={t0}, N={X}x{Y}, I(0)={N_I / (X * Y)} ")
    ax.set_xlabel("t", fontsize=9, color='blue')
    ax.set_ylabel("", fontsize=9, color='orange')
    ax.legend()
    ax.set_xlim(xlimits)

    ax.xaxis.set_major_formatter(FuncFormatter(formatOx))

    fig.show()
    fig.savefig(f'out/outDataGraph_{continuedVer}{filename_suffix}.png')
    print(dF)

    
# Change the format of t for the visualization
def tToint(t):
    delta = datetime.timedelta(minutes=float(t))
    ret = f "t={int(t)}min ({delta})"
    return ret


# TK event processing to pause the algorithm
def pauseUI(root):
    pause_var2 = StringVar()
    root.bind('<Button-1>', lambda e: pause_var2.set(1))
    root.wait_variable(pause_var2)
    pause_var2.set(0)
    root.bind('<Button-1>', lambda e: pauseUI(root))


# main() function
def main():
    global speeds_dict
    t = t0
    step = 0
    changed_points = []
    # Initial drawing of the visualization objects
    st = (700 * 2) / (X + Y)
    root, canvas = create_TK(width=st * X, height=st * Y)

    # Event handling - closing of the rendering window
    def on_closing():
        root.event_generate("<Key>")
        showGraph(out_data_array, filename_suffix="_interact")
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_closing)
    root.title('Monte Karlo method. t=0')

    # отрисовка решетки
    build_board(canvas, st, X, Y)

    # step 1 Fill in the initial state of the grid
    board = start_status()

    # drawing particles in the initial state and saving a rendering frame
    out_data_1step = prepareOutData_1step(root, canvas, st, X, Y, board=board, visibable=True)
    out_data_1step['t'] = t
    out_data_array.append(out_data_1step)
    root.title(f'Monte Carlo method. t={tToint(t)}, step={step}')
    root.bind('<Button-1>', lambda e: pauseUI(root))
    If saveGif:
        addCanvasToArr(canvas=canvas)

    while unlimetedSteps or t <= 300: # Algorithm runs until time T, display every 10, create animation Gif
        step = step + 1
        try:
            board, t, changed_points = change_board(board=board, t=t, changed_points=changed_points)
        except:
            if saveGif:
                if not os.path.isdir('out2'):
                    os.mkdir('out2')
                print('saving gif')
                images[0].save(f'out2/evolution_{continuedVer}_board.gif',
                               save_all=True,
                               append_images=images[1:],
                               duration=1000,
                               loop=0)
                clip = mp.VideoFileClip(filename=f'out2/evolution_{continuedVer}_board.gif', audio=False,
                                        target_resolution=(1000, 1000))
                print('saving mp4')
                clip.write_videofile(f'out2/evolution_{continuedVer}_board.mp4')
            print('no possible events for the current system state.\n Press any key')
            pause_var = StringVar()
            root.bind('<Key>', lambda e: pause_var.set(1))
            showGraph(out_data_array)
            root.wait_variable(pause_var)
            exit()

        # end of basic algorithm ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        # Drawing visualization
        if step % showVisualDelay == 0:
            root.title(f'Monte Carlo method. t={tToint(t)} min, step={step}')
            out_data_1step = prepareOutData_1step(root, canvas, st, X, Y, board=board, visibable=True)
            out_data_1step['t'] = t
            out_data_array.append(out_data_1step)

            if saveGif:
                addCanvasToArr(canvas=canvas)

    root.mainloop()


if __name__ == '__main__':
    main()

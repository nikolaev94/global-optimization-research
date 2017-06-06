
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.mlab import griddata
import re
import sys


def get_trial_points(filename):
    trial_points = {}

    x_values = []
    y_values = []

    x_ref = 0
    y_ref = 0

    x_calc = 0
    y_calc = 0

    file = open(filename, 'r')
    file.readline()
    txt_lines = file.read().split('\n')

    file.close()

    for txt_line in txt_lines[1:]:
        if txt_line == '----------':
            break
        match = re.search(r'.+\((-?\d+\.\d+);(-?\d+\.\d+);', txt_line)
        if match:
            x_values.append(match.group(1))
            y_values.append(match.group(2))
            continue

        match_calc_sln = re.search(r'^Calculated.+\((-?\d+\.\d+);(-?\d+\.\d+);', txt_line)
        match_ref_sln = re.search(r'Reference.+\((-?\d+\.\d+);(-?\d+\.\d+);', txt_line)

        if match_calc_sln:
            x_calc = match_calc_sln.group(1)
            y_calc = match_calc_sln.group(2)
        elif match_ref_sln:
            x_ref = match_ref_sln.group(1)
            y_ref = match_ref_sln.group(2)

    trial_points['x_calc'] = x_calc
    trial_points['y_calc'] = y_calc

    trial_points['x_ref'] = x_ref
    trial_points['y_ref'] = y_ref

    trial_points['x'] = x_values
    trial_points['y'] = y_values

    return trial_points


if len(sys.argv) < 2:
    sys.exit(1)

contour_plot_data_txt = sys.argv[1]

trials_data_txt = ''

if len(sys.argv) == 3:
    trials_data_txt = sys.argv[2]

trials_data = {}

if trials_data_txt:
    trials_data = get_trial_points(trials_data_txt)


contour_plot_data = np.genfromtxt(contour_plot_data_txt)

x = contour_plot_data[:,0]
y = contour_plot_data[:,1]
z = contour_plot_data[:,2]

print x
print y
print z

xi = np.linspace(min(x), max(x))
yi = np.linspace(min(x), max(y))

X, Y = np.meshgrid(xi, yi)
Z = griddata(x, y, z, xi, yi, interp='linear')

plt.figure()

plt.contour(X, Y, Z, 50)

plt.plot(trials_data['x'], trials_data['y'], 'bo', markersize=1.8)

plt.plot(trials_data['x_ref'], trials_data['y_ref'], 'ro', markersize=2.8)

plt.plot(trials_data['x_calc'], trials_data['y_calc'], 'go', markersize=2.8)



plt.show()


import matplotlib.pyplot as plt
import re
import sys


def get_errors_per_trials_values(filename):
    trials = []
    avg_errors = []
    max_errors = []

    file = open(filename, 'r')
    file.readline()
    csv_data = file.read().split('\n')

    file.close()

    for csv_line in csv_data:
        match = re.search(r'(\d+);(0\.\d+|1{1});(0\.\d+|1{1})', csv_line)
        if match:
            trials.append(match.group(1))
            avg_errors.append(match.group(2))
            max_errors.append(match.group(3))

    return {'x' : trials, 'y_avg' : avg_errors, 'y_max' : max_errors}

def get_portions_per_trials_values(filename):
    trials = []
    portions = []

    file = open(filename, 'r')
    file.readline()
    csv_data = file.read().split('\n')

    file.close()

    for csv_line in csv_data:
        match = re.search(r'(\d+);(0\.\d+|0{1}|1{1})', csv_line)
        if match:
            trials.append(match.group(1))
            portions.append(match.group(2))

    return {'x' : trials, 'y' : portions}



if len(sys.argv) < 3:
    sys.exit(1)

first_csv = sys.argv[1]
second_csv = sys.argv[2]

matchobj_first = re.search(r'^(.+)_(dynamic|simult)_(portions|errors)\.csv$', first_csv)
matchobj_second = re.search(r'^(.+)_(dynamic|simult)_(portions|errors)\.csv$', second_csv)

if not matchobj_first or not matchobj_second:
    sys.exit(1)

if matchobj_first.group(1) != matchobj_second.group(1):
    print "can't compare when problem classes or dimesions are different"
    sys.exit(1)

if matchobj_first.group(2) == matchobj_second.group(2):
    print "can't compare when solving strategy is equal"
    sys.exit(1)

if matchobj_first.group(3) != matchobj_second.group(3):
    print "can't compare errors with portions: different CSV formats"
    sys.exit(1)

if matchobj_first.group(3) == 'portions':
    dynamic_portions_values = {}
    simultaneous_portions_values = {}

    if matchobj_first.group(2) == 'dynamic':
        dynamic_portions_values = get_portions_per_trials_values(first_csv)
        simultaneous_portions_values = get_portions_per_trials_values(second_csv)
    else:
        simultaneous_portions_values = get_portions_per_trials_values(first_csv)
        dynamic_portions_values = get_portions_per_trials_values(second_csv)

    plt.plot(dynamic_portions_values['x'], dynamic_portions_values['y'], 'b', label='Dynamic')
    plt.plot(simultaneous_portions_values['x'], simultaneous_portions_values['y'], 'g', linestyle='--', label='Simultaneous')

    plt.title('Solved problem portions per trials')
    fig = plt.gcf()
    fig.canvas.set_window_title(matchobj_first.group(1) + '_portions')

    plt.legend()
    plt.grid(True)
    plt.show()
else: # comparing errors
    dynamic_errors_values = {}
    simultaneous_errors_values = {}

    if matchobj_first.group(2) == 'dynamic':
        dynamic_errors_values = get_errors_per_trials_values(first_csv)
        simultaneous_errors_values = get_errors_per_trials_values(second_csv)
    else:
        simultaneous_errors_values = get_errors_per_trials_values(first_csv)
        dynamic_errors_values = get_errors_per_trials_values(second_csv)

    plt.figure(1)
    plt.plot(dynamic_errors_values['x'], dynamic_errors_values['y_avg'], 'b', label='Dynamic')
    plt.plot(simultaneous_errors_values['x'], simultaneous_errors_values['y_avg'], 'g', linestyle='--', label='Simultaneous')

    plt.title('Average errors per trials')
    fig = plt.gcf()
    fig.canvas.set_window_title(matchobj_first.group(1) + '_avg_errors')

    plt.legend()
    plt.grid(True)

    plt.figure(2)
    plt.plot(dynamic_errors_values['x'], dynamic_errors_values['y_max'], 'b', label='Dynamic')
    plt.plot(simultaneous_errors_values['x'], simultaneous_errors_values['y_max'], 'g', linestyle='--', label='Simultaneous')

    plt.title('Max errors per trials')
    fig = plt.gcf()
    fig.canvas.set_window_title(matchobj_first.group(1) + '_max_errors')

    plt.legend()
    plt.grid(True)

    plt.show()



'''
first = get_axis_values(sys.argv[1])

second = get_axis_values(sys.argv[2])
'''
'''
plt.plot(first['x'], first['y'], 'b', label='Dynamic')
plt.plot(second['x'], second['y'], 'g', label='Simult')

plt.legend()

plt.grid(True)

plt.show()
'''
#plt.plot([1,2,3,4])
#plt.ylabel('some numbers')
#plt.show()

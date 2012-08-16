import csv
from sympy import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)
rc('font', family='serif')

calcium_concentations = ['0.001', '0.01', '0.1', '1.0', '5.0', '10.0', '100.0']


for ca in calcium_concentations:
    data_file = 'ca-' + ca + '.csv'
    data_read = csv.reader(open(data_file, 'rb'), delimiter=',')
    data = []

    for row in data_read:
        data.append(row)

    data = np.array(data)

    print data[-1, 10], data[-1, 48], data[-1, 7], data[-1, 6]

    plt.plot(data[1:, 0], data[1:, 10])
    plt.plot(data[1:, 0], data[1:, 48])
    plt.plot(data[1:, 0], data[1:, 7])
    plt.plot(data[1:, 0], data[1:, 6])
    plt.legend(['N', 'P', 'XB_{\mathrm PreR}', 'XB_{\mathrm PostR}'], loc='best')
    plt.xlabel('Time (mS)')
    plt.title('[Ca$^{2+}$] = ' + ca + ' $\mu$M')
    plot_file =  'ca-' + ca + '.pdf'
    plt.savefig(plot_file)
    plt.close()

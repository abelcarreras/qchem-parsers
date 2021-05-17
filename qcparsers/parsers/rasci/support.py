import numpy as np
import re


def get_occupied_electrons(configuration, structure):

    alpha_e = np.add.reduce([int(c) for c in configuration['alpha']])
    beta_e = np.add.reduce([int(c) for c in configuration['beta']])
    hole = 0 if configuration['hole'] == '' else 1
    part = 0 if configuration['part'] == '' else 1

    return (structure.number_of_electrons + structure.charge - (alpha_e + beta_e + part - hole))//2


def get_rasci_occupations_list(configuration, structure, total_orbitals):
    occupied_orbitals = get_occupied_electrons(configuration, structure)
    n_extra = total_orbitals - occupied_orbitals - len(configuration['alpha'])
    vector_alpha = [1] * occupied_orbitals + [int(c) for c in configuration['alpha']] + [0] * n_extra

    n_extra = total_orbitals - occupied_orbitals - len(configuration['beta'])
    vector_beta = [1] * occupied_orbitals + [int(c) for c in configuration['beta']] + [0] * n_extra

    if configuration['hole'] != '':
        if np.add.reduce(vector_alpha) > np.add.reduce(vector_beta):
            vector_alpha[int(configuration['hole']) - 1] = 0
        else:
            vector_beta[int(configuration['hole']) - 1] = 0

    if configuration['part'] != '':
        if np.add.reduce(vector_alpha) < np.add.reduce(vector_beta):
            vector_alpha[int(configuration['part']) - 1] = 1
        else:
            vector_beta[int(configuration['part']) - 1] = 1

    return {'alpha': vector_alpha, 'beta': vector_beta}

def read_simple_matrix(header, output, maxchar=10000, foot='-------'):
    matrix_list = []
    for m in re.finditer(header, output):
        section_state = output[m.end():m.end() + maxchar]  # 10000: assumed to max of section
        section_state = section_state[:section_state.find(foot)]
        dim = len(section_state.split('\n')[1].split())
        matrix = section_state.split('\n')[1:dim + 1]
        matrix = [[float(n) for n in l.split()] for l in matrix]
        matrix_list.append(matrix)

    return matrix_list


def read_soc_matrix(lines, dimensions):
    # for line in lines:
    #     print(line)

    col_per_line = 5
    matrix = []
    for ib in range(dimensions[0]):
        real = []
        complex = []
        for j in range((dimensions[1] - 1) // col_per_line + 1):
            real += lines[j*dimensions[0] + 1 * (j+1) + ib][11:].split()[0::2]
            complex += lines[j*dimensions[0] + 1 * (j+1) +ib][11:].split()[1::2]

        row = [float(r) + float(c[:-1]) * 1j for r, c in zip(real, complex)]
        matrix.append(row)

    return matrix

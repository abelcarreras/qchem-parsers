import numpy as np


def basis_format(basis_set_name,
                 atomic_numbers,
                 atomic_symbols,
                 shell_type,
                 n_primitives,
                 atom_map,
                 p_exponents,
                 c_coefficients,
                 p_c_coefficients):
    """
    Function to generate standard basis dictionary

    :param basis_set_name: the name of basis set
    :param atomic_numbers: atomic numbers
    :param atomic_symbols: the symbols
    :param shell_type: types of shell (check typeList)
    :param n_primitives: number of primitives
    :param atom_map: map shell - atom
    :param p_exponents: exponents of basis functions
    :param c_coefficients: coefficients of basis functions
    :param p_c_coefficients: coefficients of P functions in SP shells
    :return:
    """

    # print(n_primitives)

    typeList = {'0': ['s', 1],
                '1': ['p', 3],
                '2': ['d', 6],
                '3': ['f', 10],
                '-1': ['sp', 4],
                '-2': ['d_', 5],
                '-3': ['f_', 7]}

    atomic_numbers = [int(an) for an in atomic_numbers]
    atom_map = np.array(atom_map, dtype=int)
    # print(atom_map)
    basis_set = {'name': basis_set_name,
                 'primitive_type': 'gaussian'}

    shell_type_index = [0] + np.cumsum([typeList['{}'.format(s)][1]
                                        for s in shell_type]).tolist()
    prim_from_shell_index = [0] + np.cumsum(np.array(n_primitives, dtype=int)).tolist()

    # print(shell_type_index)
    # print(prim_from_shell_index)

    atoms_data = []
    for iatom, atomic_number in enumerate(atomic_numbers):
        symbol = str(atomic_symbols[iatom])

        shell_from_atom_counts = np.unique(atom_map, return_counts=True)[1]
        shell_from_atom_index = np.unique(atom_map, return_index=True)[1]
        # print(shell_from_atom_counts)
        # print('atom_indexes', shell_from_atom_index)
        # print('atom_number', iatom)
        # print('shells index', shell_from_atom_index[iatom])
        # print('number of shells', shell_from_atom_counts[iatom])

        shells_data = []
        for ishell in range(shell_from_atom_counts[iatom]):
            st = typeList['{}'.format(shell_type[shell_from_atom_index[iatom] + ishell])]
            # print(st, ishell)
            ini_prim = prim_from_shell_index[shell_from_atom_index[iatom] + ishell]
            fin_prim = prim_from_shell_index[shell_from_atom_index[iatom] + ishell+1]
            # print(ini_prim)
            # print(fin_prim)

            shells_data.append({
                'shell_type': st[0],
                'functions': st[1],
                'p_exponents': p_exponents[ini_prim: fin_prim],
                'con_coefficients': c_coefficients[ini_prim: fin_prim],
                'p_con_coefficients': p_c_coefficients[ini_prim: fin_prim],
            })

        atoms_data.append({'shells': shells_data,
                           'symbol': symbol,
                           'atomic_number': atomic_number})

    basis_set['atoms'] = atoms_data

    return basis_set


def reformat_input(array):
    flat_list = []
    for sublist in array:
        for item in sublist:
            if len(item) > 2:
                flat_list.append(item)
            else:
                flat_list.append(item)
    return flat_list


def vect_to_mat(vector):
    n = int(np.sqrt(0.25 + 2 * len(vector)) - 0.5)

    k = 0
    matrix = np.zeros([n, n])
    for i in range(n):
        for j in range(0, i+1):
            matrix[i, j] = vector[k + j]
            matrix[j, i] = vector[k + j]
        k += i+1

    return matrix


def get_all_nato(output):
    import re

    nato_coefficients_list = []
    nato_occupancies_list = []

    for m in re.finditer('Alpha NATO coefficients', output):
        n_elements = int(output[m.end():m.end() + 100].replace('\n', ' ').split()[2])
        nbas = int(np.sqrt(n_elements))
        data = output[m.end(): m.end() + n_elements*50].split()[3:n_elements+3]
        nato_coefficients_list.append({'alpha': np.array(data, dtype=float).reshape(nbas, nbas).tolist()})

    for m in re.finditer('Alpha Natural Orbital occupancies', output):
        n_elements = int(output[m.end():m.end() + 100].replace('\n', ' ').split()[2])
        data = output[m.end(): m.end() + n_elements*50].split()[3:n_elements+3]
        nato_occupancies_list.append({'alpha': np.array(data, dtype=float).tolist()})

    for i, m in enumerate(re.finditer('Beta NATO coefficients', output)):
        n_elements = int(output[m.end():m.end() + 100].replace('\n', ' ').split()[2])
        nbas = int(np.sqrt(n_elements))
        data = output[m.end(): m.end() + n_elements*50].split()[3:n_elements+3]
        nato_coefficients_list[i]['beta'] = np.array(data, dtype=float).reshape(nbas, nbas).tolist()

    for i, m in enumerate(re.finditer('Beta Natural Orbital occupancies', output)):
        n_elements = int(output[m.end():m.end() + 100].replace('\n', ' ').split()[2])
        data = output[m.end(): m.end() + n_elements*50].split()[3:n_elements+3]
        nato_occupancies_list[i]['beta'] = np.array(data, dtype=float).tolist()

    return nato_coefficients_list, nato_occupancies_list


def get_all_nto(output):
    import re

    nto_coefficients_list = []
    nto_occupancies_list = []

    for m in re.finditer('Natural Transition Orbital U coefficients', output):
        n_elements = int(output[m.end():m.end() + 100].replace('\n', ' ').split()[2])
        nbas = int(np.sqrt(n_elements))
        data = output[m.end(): m.end() + n_elements*50].split()[3:n_elements+3]
        nto_coefficients_list.append({'U': np.array(data, dtype=float).reshape(nbas, nbas).tolist()})

    for m in re.finditer('Natural Transition Orbital occupancies', output):
        n_elements = int(output[m.end():m.end() + 100].replace('\n', ' ').split()[2])
        data = output[m.end(): m.end() + n_elements*50].split()[3:n_elements+3]
        nto_occupancies_list.append(np.array(data, dtype=float).tolist())

    for i, m in enumerate(re.finditer('Natural Transition Orbital V coefficients', output)):
        n_elements = int(output[m.end():m.end() + 100].replace('\n', ' ').split()[2])
        nbas = int(np.sqrt(n_elements))
        data = output[m.end(): m.end() + n_elements*50].split()[3:n_elements+3]
        nto_coefficients_list[i]['V'] = np.array(data, dtype=float).reshape(nbas, nbas).tolist()

    return nto_coefficients_list, nto_occupancies_list


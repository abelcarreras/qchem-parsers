import numpy as np
import pickle


class BasisSet():
    """
    This object contains the Basis Set data of the molecule
    """

    def __init__(self,
                 basis_set_name,
                 atomic_numbers,
                 atomic_symbols,
                 shell_type,
                 n_primitives,
                 atom_map,
                 p_exponents,
                 c_coefficients,
                 p_c_coefficients):

        """
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
        self._basis_set = {'name': basis_set_name,
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
                fin_prim = prim_from_shell_index[shell_from_atom_index[iatom] + ishell + 1]
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

        self._basis_set['atoms'] = atoms_data

    def __hash__(self):
        return hash(pickle.dumps(self._basis_set, protocol=2))

    def __eq__(self, other):
        return hash(other) == hash(self)

    def get_dictionary(self):
        return self._basis_set

    def get_qc_input_txt(self):
        """
        Return basis in plane text in the format of Q-chem/Gaussian input

        :return: the basis set
        """

        basis_txt = ''

        for atom in self._basis_set['atoms']:
            basis_txt += atom['symbol'] + '\n'
            for shell in atom['shells']:
                basis_txt += '{} {} {}\n'.format(shell['shell_type'].upper(), len(shell['p_exponents']), 1.00)
                for p, c, pc in zip(shell['p_exponents'], shell['con_coefficients'], shell['p_con_coefficients']):
                    if shell['shell_type'].upper() in ['SP']:
                        basis_txt += '{:15.10e} {:15.10e} {:15.10e} \n'.format(p, c, pc)
                    else:
                        basis_txt += '{:15.10e} {:15.10e} \n'.format(p, c)

            basis_txt += '****\n'
        return basis_txt

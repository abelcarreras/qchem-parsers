from qcparsers.abstractions.molecule import Molecule
from qcparsers.parsers.fchk.support import get_all_nato, get_all_nto, reformat_input, basis_format, vect_to_mat
import numpy as np
from qcparsers.abstractions.basis import BasisSet

def parser_fchk(output):

    def convert_to_type(item_type, item):
        item_types = {'I': int,
                      'R': float}

        if type(item) is list:
            return [item_types[item_type](e) for e in item]
        else:
            return item_types[item_type](item)

    key_list = ['Charge', 'Multiplicity', 'Number of alpha electrons', 'Number of beta electrons',
                'Atomic numbers', 'Current cartesian coordinates', 'Number of basis functions', 'Shell types',
                'Number of primitives per shell', 'Shell to atom map', 'Primitive exponents',
                'Contraction coefficients', 'P(S=P) Contraction coefficients', 'Alpha MO coefficients',
                'Beta MO coefficients', 'Coordinates of each shell', 'Overlap Matrix',
                'Core Hamiltonian Matrix', 'Alpha Orbital Energies', 'Beta Orbital Energies',
                'Total SCF Density', 'Alpha NATO coefficients', 'Beta NATO coefficients',
                'Alpha Natural Orbital occupancies', 'Beta Natural Orbital occupancies',
                'Natural Transition Orbital occupancies', 'Natural Transition Orbital U coefficients',
                'Natural Transition Orbital V coefficients'
                ]

    basis_set = output.split('\n')[1].split()[-1]
    words_output = output.replace('\n', ' ').split()

    data = {}
    nw = len(words_output)
    for key in key_list:
        wc = len(key.split())
        for i in range(nw):
            word = ' '.join(words_output[i:i+wc])
            if word == key:
                item_type = words_output[i+wc]
                if words_output[i + wc + 1] == 'N=':
                    n_elements = int(words_output[i + wc + 2])
                    data[word] = convert_to_type(item_type, words_output[i + wc + 3: i + wc + n_elements + 3])
                else:
                    data[word] = convert_to_type(item_type, words_output[i + wc + 1])
                break

    bohr_to_angstrom = 0.529177249

    coordinates = np.array(data['Current cartesian coordinates']).reshape(-1, 3) * bohr_to_angstrom
    structure = Molecule(coordinates=coordinates.tolist(),
                         atomic_numbers=data['Atomic numbers'],
                         multiplicity=data['Multiplicity'],
                         charge=data['Charge'])

    if not 'P(S=P) Contraction coefficients' in data:
        data['P(S=P) Contraction coefficients'] = np.zeros_like(data['Contraction coefficients']).tolist()

    #basis = basis_format(basis_set_name=basis_set,

    basis = BasisSet(basis_set_name=basis_set,
                         atomic_numbers=structure.get_atomic_numbers(),
                         atomic_symbols=structure.get_symbols(),
                         shell_type=data['Shell types'],
                         n_primitives=data['Number of primitives per shell'],
                         atom_map=data['Shell to atom map'],
                         p_exponents=data['Primitive exponents'],
                         c_coefficients=data['Contraction coefficients'],
                         p_c_coefficients=data['P(S=P) Contraction coefficients'])

    nbas = data['Number of basis functions']

    final_dict = {'structure': structure,
                  'basis': basis,
                  'number_of_electrons': {'alpha': data['Number of alpha electrons'],
                                          'beta': data['Number of beta electrons']}
                  }

    if 'Alpha MO coefficients' in data:
        final_dict['coefficients'] = {'alpha': np.array(data['Alpha MO coefficients']).reshape(nbas, nbas).tolist()}
        final_dict['mo_energies'] = {'alpha': data['Alpha Orbital Energies']}

    if 'Beta MO coefficients' in data:
        final_dict['coefficients']['beta'] = np.array(data['Beta MO coefficients']).reshape(nbas, nbas).tolist()
        final_dict['mo_energies']['beta'] = data['Beta Orbital Energies']

    if 'Total SCF Density' in data:
        final_dict['scf_density'] = vect_to_mat(data['Total SCF Density']).tolist()

    if 'Core Hamiltonian Matrix' in data:
        final_dict['scf_density'] = vect_to_mat(data['Core Hamiltonian Matrix']).tolist()

    if 'Overlap Matrix' in data:
        final_dict['overlap'] = vect_to_mat(data['Overlap Matrix']).tolist()

    if 'Alpha NATO coefficients' in data:
        final_dict['nato_coefficients'] = {'alpha': np.array(data['Alpha NATO coefficients']).reshape(nbas, nbas).tolist()}
        final_dict['nato_occupancies'] = {'alpha': data['Alpha Natural Orbital occupancies']}

    if 'Beta NATO coefficients' in data:
        final_dict['nato_coefficients'].update({
            'beta': np.array(data['Beta NATO coefficients']).reshape(nbas, nbas).tolist()})
        final_dict['nato_occupancies'].update({'beta': data['Beta Natural Orbital occupancies']})

    # check multiple NATO (may be improved)
    if 'Alpha NATO coefficients' in data:
        nato_coefficients_list, nato_occupancies_list = get_all_nato(output)
        if len(nato_occupancies_list) > 1:
            final_dict['nato_coefficients_multi'] = nato_coefficients_list
            final_dict['nato_occupancies_multi'] = nato_occupancies_list

    if 'Natural Transition Orbital occupancies' in data:
        nat_coefficients_list, nat_occupancies_list = get_all_nto(output)
        if len(nat_occupancies_list) > 1:
            final_dict['nto_coefficients_multi'] = nat_coefficients_list
            final_dict['nto_occupancies_multi'] = nat_occupancies_list

    return final_dict


if __name__ == '__main__':
    a = parser_fchk(open('fchk_1.out').read())
    basis = a['basis']

    print(basis.get_qc_format())
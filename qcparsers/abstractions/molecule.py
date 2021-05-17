import numpy as np


class Molecule:
    """
    This object contains the geometrical data of the molecule
    """
    def __init__(self,
                 coordinates,
                 symbols=None,
                 atomic_numbers=None,
                 charge=0,
                 multiplicity=1,
                 name=None):
        """
        :param coordinates: List containing the cartesian coordinates of each atom in Angstrom
        :param symbols: Symbols of the atoms within the molecule
        :param atomic_numbers: Atomic numbers of the atoms within the molecule
        :param charge: charge of the molecule
        :param multiplicity: multiplicity of the molecule
        :param name: name of the molecule
        """

        self._coordinates = np.array(coordinates)
        self._atomic_numbers = atomic_numbers
        self._symbols = symbols
        self._charge = charge
        self._multiplicity = multiplicity
        self._name = name

        self._atomic_masses = None
        self._number_of_atoms = None

        if atomic_numbers is not None:
            self._symbols = [atom_data[i][1] for i in atomic_numbers]

    def __hash__(self):
        return hash((np.array_str(self._coordinates, precision=8),
                     tuple(self._symbols),
                     self._charge,
                     self._name,
                     self._multiplicity))

    def __eq__(self, other):
        return hash(other) == hash(self)

    def __str__(self):
        return self.get_xyz()

    def get_coordinates(self, fragment=None):
        """
        gets the cartesian coordinates

        :param fragment: list of atoms that are part of the fragment

        :return: coordinates list
        """

        if fragment is None:
            return np.array(self._coordinates).tolist()
        else:
            return np.array(self._coordinates)[fragment].tolist()

    def get_positions(self, fragment=None):
        return self.get_coordinates(fragment)

    def export_ase_atoms(self):
        from ase import Atoms
        name_ase = ''.join(self.get_symbols())
        charge_per_atom = self.charge/self.get_number_of_atoms()
        charges = [charge_per_atom] * self.get_number_of_atoms()
        return Atoms(name_ase, positions=self.get_coordinates(), charges=charges)

    @property
    def name(self):
        """
        returns the name
        :return: structure name
        """
        return self._name


    @property
    def charge(self):
        """
        returns the charge
        :return: the charge
        """
        return self._charge

    @charge.setter
    def charge(self, charge):
        self._charge = charge

    @property
    def multiplicity(self):
        """
        returns the multiplicity

        :return: the multiplicity
        """
        return self._multiplicity

    @multiplicity.setter
    def multiplicity(self, multiplicity):
        self._multiplicity = multiplicity

    @property
    def number_of_electrons(self):
        """
        returns the total number of electrons

        :return: number of total electrons
        """
        return int(np.sum(self.get_atomic_numbers()) + self.charge)

    @property
    def alpha_electrons(self):
        """
        returns the alpha electrons

        :return: number of alpha electrons
        """
        alpha_unpaired = self.multiplicity // 2 + 1 if (self.number_of_electrons % 2) else self.multiplicity // 2
        return self.number_of_electrons // 2 + alpha_unpaired

    @property
    def beta_electrons(self):
        """
        returns the number of beta electrons

        :return: number of beta electrons
        """
        return self.number_of_electrons - self.alpha_electrons

    def get_atomic_numbers(self):
        """
        get the atomic numbers of the atoms of the molecule

        :return: list with the atomic numbers
        """
        if self._atomic_numbers is None:
            self._atomic_numbers = [[data[1].upper() for data in atom_data].index(element.upper())
                                    for element in self.get_symbols()]
        return self._atomic_numbers

    def get_symbols(self):
        """
        get the  atomic element symbols of the atoms of the molecule

        :return: list of symbols
        """
        if self._symbols is None:
            self._symbols = np.array(atom_data)[self.get_atomic_numbers()].T[1]
        return np.array([i for i in self._symbols if i != "X"], dtype=str)

    def get_number_of_atoms(self):
        """
        get the number of atoms

        :return: number of atoms
        """

        return np.array(self.get_coordinates()).shape[0]

    def get_atomic_masses(self):
        """
        get the atomic masses of the atoms of the molecule

        :return: list of atomic masses
        """
        if self._atomic_masses is None:

            try:
                masses_string = np.array(atom_data)[:, 3:4][[np.where(np.array(atom_data)==element)[0][0]
                                                             for element in self.get_symbols()]]
                self._atomic_masses = np.array(masses_string, dtype=float).T[0]
            except TypeError:
                print('Error reading element labels')
                exit()
        return self._atomic_masses

    def get_valence_electrons(self):
        """
        gets number of valence electrons

        :return: number of valence electrons
        """
        valence_electrons = 0
        for number in self.get_atomic_numbers():
            if 2 >= number > 0:
                valence_electrons += np.mod(number, 2)
            if 18 >= number > 2:
                valence_electrons += np.mod(number-2, 8)
            if 54 >= number > 18:
                valence_electrons += np.mod(number-18, 18)
            if 118 >= number > 54:
                valence_electrons += np.mod(number-54, 32)
            if number > 118:
                raise Exception('Atomic number size not implemented')

        valence_electrons -= self.charge

        return valence_electrons

    def get_xyz(self, title=''):
        """
        generates a XYZ formatted file

        :param title: title of the molecule
        :return: string with the formatted XYZ file
        """
        txt = '{}\n{}\n'.format(self.get_number_of_atoms(), title)
        for s, c in zip(self.get_symbols(), self.get_coordinates()):
            txt += '{:2} '.format(s) + '{:15.10f} {:15.10f} {:15.10f}\n'.format(*c)

        return txt


atom_data = [
    [  0, "X", "X", 0], # 0
    [  1, "H", "Hydrogen", 1.00794], # 1
    [  2, "He", "Helium", 4.002602], # 2
    [  3, "Li", "Lithium", 6.941], # 3
    [  4, "Be", "Beryllium", 9.012182], # 4
    [  5, "B", "Boron", 10.811], # 5
    [  6, "C", "Carbon", 12.0107], # 6
    [  7, "N", "Nitrogen", 14.0067], # 7
    [  8, "O", "Oxygen", 15.9994], # 8
    [  9, "F", "Fluorine", 18.9984032], # 9
    [ 10, "Ne", "Neon", 20.1797], # 10
    [ 11, "Na", "Sodium", 22.98976928], # 11
    [ 12, "Mg", "Magnesium", 24.3050], # 12
    [ 13, "Al", "Aluminium", 26.9815386], # 13
    [ 14, "Si", "Silicon", 28.0855], # 14
    [ 15, "P", "Phosphorus", 30.973762], # 15
    [ 16, "S", "Sulfur", 32.065], # 16
    [ 17, "Cl", "Chlorine", 35.453], # 17
    [ 18, "Ar", "Argon", 39.948], # 18
    [ 19, "K", "Potassium", 39.0983], # 19
    [ 20, "Ca", "Calcium", 40.078], # 20
    [ 21, "Sc", "Scandium", 44.955912], # 21
    [ 22, "Ti", "Titanium", 47.867], # 22
    [ 23, "V", "Vanadium", 50.9415], # 23
    [ 24, "Cr", "Chromium", 51.9961], # 24
    [ 25, "Mn", "Manganese", 54.938045], # 25
    [ 26, "Fe", "Iron", 55.845], # 26
    [ 27, "Co", "Cobalt", 58.933195], # 27
    [ 28, "Ni", "Nickel", 58.6934], # 28
    [ 29, "Cu", "Copper", 63.546], # 29
    [ 30, "Zn", "Zinc", 65.38], # 30
    [ 31, "Ga", "Gallium", 69.723], # 31
    [ 32, "Ge", "Germanium", 72.64], # 32
    [ 33, "As", "Arsenic", 74.92160], # 33
    [ 34, "Se", "Selenium", 78.96], # 34
    [ 35, "Br", "Bromine", 79.904], # 35
    [ 36, "Kr", "Krypton", 83.798], # 36
    [ 37, "Rb", "Rubidium", 85.4678], # 37
    [ 38, "Sr", "Strontium", 87.62], # 38
    [ 39, "Y", "Yttrium", 88.90585], # 39
    [ 40, "Zr", "Zirconium", 91.224], # 40
    [ 41, "Nb", "Niobium", 92.90638], # 41
    [ 42, "Mo", "Molybdenum", 95.96], # 42
    [ 43, "Tc", "Technetium", 0], # 43
    [ 44, "Ru", "Ruthenium", 101.07], # 44
    [ 45, "Rh", "Rhodium", 102.90550], # 45
    [ 46, "Pd", "Palladium", 106.42], # 46
    [ 47, "Ag", "Silver", 107.8682], # 47
    [ 48, "Cd", "Cadmium", 112.411], # 48
    [ 49, "In", "Indium", 114.818], # 49
    [ 50, "Sn", "Tin", 118.710], # 50
    [ 51, "Sb", "Antimony", 121.760], # 51
    [ 52, "Te", "Tellurium", 127.60], # 52
    [ 53, "I", "Iodine", 126.90447], # 53
    [ 54, "Xe", "Xenon", 131.293], # 54
    [ 55, "Cs", "Caesium", 132.9054519], # 55
    [ 56, "Ba", "Barium", 137.327], # 56
    [ 57, "La", "Lanthanum", 138.90547], # 57
    [ 58, "Ce", "Cerium", 140.116], # 58
    [ 59, "Pr", "Praseodymium", 140.90765], # 59
    [ 60, "Nd", "Neodymium", 144.242], # 60
    [ 61, "Pm", "Promethium", 0], # 61
    [ 62, "Sm", "Samarium", 150.36], # 62
    [ 63, "Eu", "Europium", 151.964], # 63
    [ 64, "Gd", "Gadolinium", 157.25], # 64
    [ 65, "Tb", "Terbium", 158.92535], # 65
    [ 66, "Dy", "Dysprosium", 162.500], # 66
    [ 67, "Ho", "Holmium", 164.93032], # 67
    [ 68, "Er", "Erbium", 167.259], # 68
    [ 69, "Tm", "Thulium", 168.93421], # 69
    [ 70, "Yb", "Ytterbium", 173.054], # 70
    [ 71, "Lu", "Lutetium", 174.9668], # 71
    [ 72, "Hf", "Hafnium", 178.49], # 72
    [ 73, "Ta", "Tantalum", 180.94788], # 73
    [ 74, "W", "Tungsten", 183.84], # 74
    [ 75, "Re", "Rhenium", 186.207], # 75
    [ 76, "Os", "Osmium", 190.23], # 76
    [ 77, "Ir", "Iridium", 192.217], # 77
    [ 78, "Pt", "Platinum", 195.084], # 78
    [ 79, "Au", "Gold", 196.966569], # 79
    [ 80, "Hg", "Mercury", 200.59], # 80
    [ 81, "Tl", "Thallium", 204.3833], # 81
    [ 82, "Pb", "Lead", 207.2], # 82
    [ 83, "Bi", "Bismuth", 208.98040], # 83
    [ 84, "Po", "Polonium", 0], # 84
    [ 85, "At", "Astatine", 0], # 85
    [ 86, "Rn", "Radon", 0], # 86
    [ 87, "Fr", "Francium", 0], # 87
    [ 88, "Ra", "Radium", 0], # 88
    [ 89, "Ac", "Actinium", 0], # 89
    [ 90, "Th", "Thorium", 232.03806], # 90
    [ 91, "Pa", "Protactinium", 231.03588], # 91
    [ 92, "U", "Uranium", 238.02891], # 92
    [ 93, "Np", "Neptunium", 0], # 93
    [ 94, "Pu", "Plutonium", 0], # 94
    [ 95, "Am", "Americium", 0], # 95
    [ 96, "Cm", "Curium", 0], # 96
    [ 97, "Bk", "Berkelium", 0], # 97
    [ 98, "Cf", "Californium", 0], # 98
    [ 99, "Es", "Einsteinium", 0], # 99
    [100, "Fm", "Fermium", 0], # 100
    [101, "Md", "Mendelevium", 0], # 101
    [102, "No", "Nobelium", 0], # 102
    [103, "Lr", "Lawrencium", 0], # 103
    [104, "Rf", "Rutherfordium", 0], # 104
    [105, "Db", "Dubnium", 0], # 105
    [106, "Sg", "Seaborgium", 0], # 106
    [107, "Bh", "Bohrium", 0], # 107
    [108, "Hs", "Hassium", 0], # 108
    [109, "Mt", "Meitnerium", 0], # 109
    [110, "Ds", "Darmstadtium", 0], # 110
    [111, "Rg", "Roentgenium", 0], # 111
    [112, "Cn", "Copernicium", 0], # 112
    [113, "Uut", "Ununtrium", 0], # 113
    [114, "Uuq", "Ununquadium", 0], # 114
    [115, "Uup", "Ununpentium", 0], # 115
    [116, "Uuh", "Ununhexium", 0], # 116
    [117, "Uus", "Ununseptium", 0], # 117
    [118, "Uuo", "Ununoctium", 0], # 118
    ]

if __name__ == '__main__':
    mol = Molecule([[1, 0, 0],
                    [2, 0, 0],
                    [0, 0, 1]],
                   symbols=['H', 'H', 'O'],
                   charge=1)

    print(mol)
    print('--------------------')
    asemol = mol.export_ase_atoms()
    print(asemol)
    print(asemol.get_positions())
    print(asemol.symbols)
    print('ic', asemol.get_initial_charges())
    print(asemol.get_chemical_symbols())
    print(asemol.get_global_number_of_atoms())
    print(asemol)

    print(asemol.get_distances(0, 1))

    import ase.io.gaussian
    with open('test.gj', 'w') as f:
        #ase.io.gromacs.write_gromacs(f, asemol)
        ase.io.gaussian.write_gaussian_in(f, asemol)
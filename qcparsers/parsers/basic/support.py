

def get_orbital_energies(orbitals_section):
    # print(orbitals_section)
    occupied = orbitals_section.find('-- Occupied --')
    virtual = orbitals_section.find('-- Virtual --')

    occupied_section = orbitals_section[occupied:virtual]
    virtual_section = orbitals_section[virtual:]

    occupied_energies = [float(energy) if energy != '********' else None for energy in ' '.join(occupied_section.split('\n')[1::2]).split()]
    virtual_energies = [float(energy) if energy != '********' else None for energy in ' '.join(virtual_section.split('\n')[1::2]).split()]

    return occupied_energies + virtual_energies

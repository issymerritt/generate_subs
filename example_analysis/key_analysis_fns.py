import numpy as np
from copy import copy, deepcopy
from tqdm import tqdm
import periodictable
from typing import List
import scipy.constants as cst
import sys
ThreeStringList = List[str]

def hartree_to_eV(energy_h):
    return energy_h*cst.physical_constants['Hartree energy in eV'][0]

def hartree_to_kJmol(energy_h):
    return energy_h*cst.physical_constants['Hartree energy'][0]*cst.Avogadro/1000 

def eV_to_hartree(energy_eV):
    return energy_eV/cst.physical_constants['Hartree energy in eV'][0]

def read_freqs(filename):
    freqs_list = []
    with open(filename,'r') as t:
        react_lines = t.read().splitlines()
    last = 0
    for idx, line in enumerate(react_lines):
        if "Harmonic frequencies" in line:
            last = idx
    for line in react_lines[last:]:
        if "Frequencies --" in line:
            for a in line.split()[-3:]:
                freqs_list += [float(a)]
    return np.array(freqs_list)

def read_energy(filename):
    with open(filename,'r') as t:
        react_lines = t.read().splitlines()
        for line in react_lines:
            if "SCF Done:" in line:
                energy = line.split()[-5]
    return float(energy)

def check_geom_conv_and_freq(filename, ts=False):
    converged = False
    neg_freqs = False
    freqs_list = []
    with open(filename,'r') as t:
        react_lines = t.read().splitlines()
        for line in react_lines:
            if "Frequencies --" in line:
                for a in line.split()[-3:]:
                    freqs_list += [float(a)]
            elif "Optimization completed" in line:
                converged = True
    nb_allowed_neg = 0 if ts == False else 1
    if len(freqs_list) == 0:
        neg_freqs = None
    elif freqs_list[nb_allowed_neg] < 0:
        neg_freqs = True
    return converged, neg_freqs

def get_polarity_polarization_info(filename):
    dipole_mom, isotrop_pol, anisotrop_pol = False, False, False
    with open(filename,'r') as q:
        a = q.read().splitlines()
    for idx, line in enumerate(a):
        if "Electric dipole moment (input orientation):" in line:
            dm_idx = idx + 3
        elif "Dipole polarizability, Alpha (input orientation)" in line:
            iso_idx = idx + 4
            aniso_idx = idx + 5
    dipole_mom = float(a[dm_idx].split()[2].replace('D','E'))
    isotrop_pol = float(a[iso_idx].split()[1].replace('D','E'))
    anisotrop_pol = float(a[aniso_idx].split()[1].replace('D','E'))
    return dipole_mom, isotrop_pol, anisotrop_pol


class Key_Geometry(object):
    def __init__(self, filename : str, ts : bool=False, verbose : bool=False):
        if verbose:
            print(f'Parsing data for {filename}.')
        with open(filename,'r') as t:
            self.filelines = t.read().splitlines()  
        self.name = filename
        self._set_coords(self.filelines)
        if not check_geom_conv_and_freq(filename, ts)[0]:
            raise Exception(f'{filename} not converged/too many negative frequencies')
        else:
            if verbose:
                print('Optimisation completed successfully, correct number of negative frequencies (Min = 0, TS = 1).')
        self.freqs = read_freqs(filename)
        self.nb_atoms = len(self.coords)
        self.scf_energy = read_energy(filename)
        if len(self.freqs) != self.nb_atoms*3 - 6:
            print(f'Incorrect nb of frequencies for file {filename}')
        self._set_nb_elecs()
        self.dipole_moment, self.iso_polarisability, self.aniso_polarisability = get_polarity_polarization_info(filename)
        self._get_excited_states(self.filelines, verbose)
        self._get_orbital_energies(self.filelines)
        del self.filelines

    def molecule_report(self):
        print(f'---- {self.name} ----')
        print(f'Dipole Moment: {self.dipole_moment}')
        print(f'Isotropic Polarisability: {self.iso_polarisability}')
        print(f'Anisotropic Polarisability: {self.aniso_polarisability}')
        print(f'Excited States: {self.excited_states}')
        print(f'Oscillator Strengths: {self.es_oscillator_strengths}')
        try:
            coords = self.coords
            try:
                freqs = self.freqs
            except AttributeError:
                print(f'Frequencies not found for {self.name}')
            else:
                print(f'Frequencies and Coordinates available')
        except AttributeError:
            print(f'Coordinates not found for {self.name}')
        try:
            homo = self.occ_orbs_energies[-1]
            lumo = self.virt_orbs_energies[0]
            print(f'HOMO energy: {homo}')
            print(f'LUMO energy: {lumo}')
        except:
            print('Orbital energies not extracted')

    def _set_states_hartree(self, excited_states : List[float]):
        self.all_states = [self.scf_energy]
        for state in excited_states:
            self.all_states.append(self.scf_energy + eV_to_hartree(state))
    
    def _get_orbital_energies(self, filelines : List[str]):
        occ_orbital_en = []
        virt_orbital_en = []
        for idx,line in enumerate(filelines):
            if "Population analysis using the SCF Density" in line:
                final_popline = idx
        for line in filelines[final_popline:]:
            if "Alpha  occ. eigenvalues" in line:
                eigenvals = [float(a) for a in line.split()[4:]]
                occ_orbital_en.extend(eigenvals)
            elif "Alpha virt. eigenvalues" in line:
                eigenvals = [float(a) for a in line.split()[4:]]
                virt_orbital_en.extend(eigenvals) 
        self.occ_orbs_energies = occ_orbital_en
        self.virt_orbs_energies = virt_orbital_en

    def _set_coords(self, filelines : List[str]):
        for idx, line in enumerate(filelines):
            if "Input orientation" in line:
                geom_idx = idx + 5
        xyz_lines_raw = []
        for line in filelines[geom_idx:]:
            if '---------------------------------------' in line:
                break
            xyz_lines_raw += [line]
        xyz_lines = []
        for atom_nb_line in xyz_lines_raw:
            xyz_lines += [f'{periodictable.elements[float(atom_nb_line.split()[1])]}   {atom_nb_line.split()[3]}  {atom_nb_line.split()[4]}  {atom_nb_line.split()[5]}']
        self.coords = {}
        for idx, line in enumerate(xyz_lines):
            self.coords[idx+1] = Coord(line)

    def _set_nb_elecs(self):
        self.electrons = 0
        for i in self.coords:
            self.electrons += periodictable.elements.symbol(self.coords[i].at_type).number

    def _get_excited_states(self, filelines : List[str], verbose : bool=False):
        self.excited_states = []
        self.es_oscillator_strengths = []
        for line in filelines:
            if "Excited State   1" in line:
                self.excited_states = [float(line.split()[4])]
                self.es_oscillator_strengths = [float(line.split()[8][2:])]
            elif "Excited State" in line:
                self.excited_states += [float(line.split()[4])]
                self.es_oscillator_strengths += [float(line.split()[8][2:])]
        if not self.excited_states or not self.es_oscillator_strengths:
            if verbose == True :
                print(f'Excited State Calculation not located for {self.name}')
        else:
            self._set_states_hartree(self.excited_states)
    
    def write_mol(self,outfile='out.xyz',comment_string='Molecule written by Molecule.py'):
        write_string = str(self.nb_atoms) + '\n' + f'{comment_string} \n'
        with open(outfile,'w') as t:
            t.write(write_string)
            for atom in self.coords:
                coordstring = self.coords[atom].at_type
                coordstring += '      ' + f'{self.coords[atom].xyz_coords[0]:.10f}' + '       ' + f'{self.coords[atom].xyz_coords[1]:.10f}' + '       ' + f'{self.coords[atom].xyz_coords[2]:.10f}' + '\n'
                t.write(coordstring)
                
    def gen_xyz_string(self):
        xyzstr = ''
        for coord in self.coords:
            xyzstr += self.coords[coord].at_type 
            for axis in range(0,3):
                xyzstr += ' '
                xyzstr += f'{self.coords[coord].xyz_coords[axis]:.8f}'
            xyzstr += '; '
        return xyzstr

    def bond_length(self,at1,at2):
        atom1 = self.coords[at1].xyz_coords
        atom2 = self.coords[at2].xyz_coords
        bond = atom1 - atom2
        return np.linalg.norm(bond)
    
    def angle(self,at1,at2,at3):
        v1 = self.coords[at1].xyz_coords - self.coords[at2].xyz_coords
        v2 = self.coords[at3].xyz_coords - self.coords[at2].xyz_coords
        norm1 = np.linalg.norm(v1)
        norm2 = np.linalg.norm(v2)
        return np.degrees(np.arccos(np.dot(v1,v2)/norm1/norm2))

    def dihedral(self,at1,at2,at3,at4):
        """Praxeolitic formula
        1 sqrt, 1 cross product"""

        b0 = -1.0*(self.coords[at2].xyz_coords - self.coords[at1].xyz_coords)
        b1 = self.coords[at3].xyz_coords - self.coords[at2].xyz_coords
        b2 = self.coords[at4].xyz_coords - self.coords[at3].xyz_coords
        b1 /= np.linalg.norm(b1)
        v = b0 - np.dot(b0, b1)*b1
        w = b2 - np.dot(b2, b1)*b1
        x = np.dot(v, w)
        y = np.dot(np.cross(b1, v), w)
        return np.degrees(np.arctan2(y, x))

class Coord(object):
    def __init__(self,xyz_line):
        self.at_type = xyz_line.split()[0]
        self._xyz_list = [float(a) for a in xyz_line.split()[1:4]]
        self.xyz_coords = np.array(self._xyz_list)
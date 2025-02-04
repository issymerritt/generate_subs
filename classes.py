from typing import Optional, List

import numpy as np
import periodictable


class Coord:
    def __init__(self, xyz_line: str):
        self.at_type = xyz_line.split()[0]
        self._xyz_list = [float(a) for a in xyz_line.split()[1:4]]
        self.xyz_coords = np.array(self._xyz_list)
        self.label: Optional[str] = None

    def set_label(self, lab_name: str):
        self.label = lab_name


class Molecule:
    def __init__(self, coord_list: Optional[List[Coord]] = None, filename: Optional[str] = None):
        self.coords = {}

        if not coord_list and not filename:
            raise Exception('Either list of Coord objects (coord_list = ) or xyz input file (filename = ) must be given')
        elif coord_list and filename:
            raise Exception('Only one of atom coordinates (coord_list = ) or xyz input file (filename = ) should be provided')
        if coord_list:
            self._set_coords_from_list(coord_list)
        elif filename:
            self._set_coords_from_file(filename)

        self.nb_atoms = len(self.coords)
        self._set_nb_electrons()

    def _set_coords_from_list(self, coord_list: List[Coord]):
        for idx, coord in enumerate(coord_list):
            self.coords[idx + 1] = coord
        self.label = 'Molecule created from coordinate list'

    def _set_coords_from_file(self, filename: str):
        with open(f'{filename}.xyz') as fp:
            self._xyz_lines = fp.read().splitlines()
        for idx, line in enumerate(self._xyz_lines[2:]):
            if len(line.split()) > 1:
                self.coords[idx + 1] = Coord(line)
        self.label = self._xyz_lines[1]
        del self._xyz_lines

    def _set_nb_electrons(self):
        self.electrons = 0
        for i in self.coords:
            self.electrons += periodictable.elements.symbol(self.coords[i].at_type).number

    def write_mol(self, outfile: str = 'out.xyz', comment_string: str = 'Molecule written by Molecule.py'):
        write_string = f'{self.nb_atoms} \n {comment_string} \n'
        with open(outfile, 'w') as t:
            t.write(write_string)
            for atom in self.coords:
                coordstring = self.coords[atom].at_type
                coordstring += '      ' + f'{self.coords[atom].xyz_coords[0]:.10f}' + '       ' + f'{self.coords[atom].xyz_coords[1]:.10f}' + '       ' + f'{self.coords[atom].xyz_coords[2]:.10f}' + '\n'
                t.write(coordstring)

    def gen_xyz_string(self) -> str:
        xyzstr = ''
        for coord in self.coords:
            xyzstr += self.coords[coord].at_type
            for axis in range(0, 3):
                xyzstr += ' '
                xyzstr += f'{self.coords[coord].xyz_coords[axis]:.8f}'
            xyzstr += '; '
        return xyzstr

    def bond_length(self, at1: int, at2: int) -> float:
        atom1 = self.coords[at1].xyz_coords
        atom2 = self.coords[at2].xyz_coords
        bond = atom1 - atom2
        return float(np.linalg.norm(bond))

    def angle(self, at1: int, at2: int, at3: int) -> float:
        v1 = self.coords[at1].xyz_coords - self.coords[at2].xyz_coords
        v2 = self.coords[at3].xyz_coords - self.coords[at2].xyz_coords
        norm1 = np.linalg.norm(v1)
        norm2 = np.linalg.norm(v2)
        return np.degrees(np.arccos(np.dot(v1, v2) / norm1 / norm2))

    def dihedral(self, at1: int, at2: int, at3: int, at4: int) -> float:
        """ Praxeolitic formula, 1 sqrt, 1 cross product """
        b0 = -1.0 * (self.coords[at2].xyz_coords - self.coords[at1].xyz_coords)
        b1 = self.coords[at3].xyz_coords - self.coords[at2].xyz_coords
        b2 = self.coords[at4].xyz_coords - self.coords[at3].xyz_coords
        b1 /= np.linalg.norm(b1)
        v = b0 - np.dot(b0, b1) * b1
        w = b2 - np.dot(b2, b1) * b1
        x = np.dot(v, w)
        y = np.dot(np.cross(b1, v), w)
        return np.degrees(np.arctan2(y, x))

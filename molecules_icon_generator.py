#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Luca Monari"
__credits__ = ["Luca Monari"]
__version__ = "1.0"
__email__ = "luca.monari@etu.unistra.fr"

import numpy as np
import cv2
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import math
import itertools
import argparse
import os

# change the Path below with your path tho the icons
atom_icon_dir = "base-icons/"

def load_icons(folder, resize_dim=(300, 300)):
    icon_map = dict()
    for file in os.listdir(folder):
        file_img = cv2.imread(folder + os.sep + file, cv2.IMREAD_UNCHANGED)
        file_img = cv2.resize(file_img, resize_dim, interpolation=cv2.INTER_AREA)
        icon_map[file.split('.')[0]] = file_img
    return icon_map


# LOAD ICONS IMAGES
icon_map = load_icons(atom_icon_dir)


def rotate_image(image, angle):
    image_center = tuple(np.array(image.shape[1::-1]) / 2) # the '//' division gives error in getRotationMatrix2D
    rot_mat = cv2.getRotationMatrix2D(image_center, angle, 1.0)
    result = cv2.warpAffine(image, rot_mat, image.shape[1::-1], flags=cv2.INTER_LINEAR)
    return result


def add_image(src, new, position, overwrite=True):
    # position (x,y)
    new_x = new.shape[1]
    new_y = new.shape[0]
    y_offset = position[1] - new_y // 2
    x_offset = position[0] - new_x // 2
    # overwrite the source array with the new image only if the alpha channel is equal to 0
    for y_index, y in enumerate(range(y_offset, y_offset + new.shape[0])):
        for x_index, x in enumerate(range(x_offset, x_offset + new.shape[1])):
            if new[y_index, x_index][3] != 0:
                # don't overwrite the image if the pixel is coloured (not black nor white)
                if not overwrite and src[y, x][:3] != (255, 255, 255) and src[y, x][:3] != (0, 0, 0):
                    continue
                src[y, x] = new[y_index, x_index]
    return src


def add_bond(src, bond_type, degree, position, length):
    # resize the bond length to make it reach the atoms center
    # resized_bond = cv2.resize(bond_type.copy(), (length, bond_type.shape[0]), interpolation=cv2.INTER_AREA)
    # the resize method fail, thus I extend the image array manually to match the lenght of the bond
    missing_length = length - bond_type.shape[0]
    one_column = bond_type[:, 0]
    missing_array = np.repeat(one_column, missing_length)
    resized_bond = np.concatenate((bond_type, missing_array), axis=1)
    rotated_bond = rotate_image(resized_bond, degree)
    add_image(src, rotated_bond, position)


def icon_print(SMILES, name='molecule_icon', directory=os.getcwd(), rdkit_img=False,
               single_bonds=False, remove_H=False, verbose=False, save=True,
               symbol_img_dict=icon_map, position_multiplier=150):
    mol = Chem.MolFromSmiles(SMILES)
    if not remove_H:
        mol = Chem.AddHs(mol)
    AllChem.Compute2DCoords(mol)
    mol.GetConformer()

    atom_map = dict()
    atom_type_map = dict()
    atom_bond_map = dict()
    for i, atom in enumerate(mol.GetAtoms()):
        positions = mol.GetConformer().GetAtomPosition(i)
        atom_symbol = atom.GetSymbol()
        # 150 is a nice value to obtain big icons
        x = int(positions.x * position_multiplier)
        y = int(positions.y * position_multiplier)
        atom_map[i] = (x, -y)  # y axis of an image is the inverse
        atom_bonds = rdkit.Chem.rdchem.Atom.GetBonds(atom)
        atom_bond_map[i] = len(atom_bonds)
        atom_type_map[i] = atom_symbol
        if verbose:
            print(atom_symbol, x, y)

    # get the maximum position from the atom (positive or negative
    max_pos = max([abs(pos) for pos in list(itertools.chain(*list(atom_map.values())))])
    # to the max position we have to add width of the atom icons
    max_tot_pos = max_pos + list(symbol_img_dict.values())[0].shape[1]
    # multiply for 2 because the max position is considered from the center of the image
    dimension = int(max_tot_pos * 2)
    # base black image, four channel to include alpha channel
    img_rgb = np.zeros((dimension, dimension, 3), np.uint8)
    # make the image white
    img_rgb.fill(255)
    # add alpha channel
    img = cv2.cvtColor(img_rgb, cv2.COLOR_RGB2RGBA)
    # make alpha channel equal to 0
    img[:, :, 3] = np.zeros((dimension, dimension), np.uint8)
    
    # add bond
    bonds_list = list(itertools.combinations(list(atom_map.keys()), 2))
    aromatic_index = set()
    for couple in bonds_list:
        atom1 = couple[0]
        atom2 = couple[1]
        BOND = mol.GetBondBetweenAtoms(atom1, atom2)
        if BOND:
            x1 = atom_map[atom1][0] + dimension // 2
            y1 = atom_map[atom1][1] + dimension // 2
            x2 = atom_map[atom2][0] + dimension // 2
            y2 = atom_map[atom2][1] + dimension // 2
            mid_x = (x1 + x2) // 2
            mid_y = (y1 + y2) // 2
            position = (mid_x, mid_y)
            myradians = math.atan2(y1 - y2, x2 - x1)
            mydegrees = math.degrees(myradians)
            length = int(math.dist([x1, y1], [x2, y2]))
            print(length)
            b_type = BOND.GetBondType()
            bond_img = symbol_img_dict['single']
            if rdkit.Chem.rdchem.BondType.DOUBLE == b_type and not single_bonds:
                bond_img = symbol_img_dict['double']
            elif rdkit.Chem.rdchem.BondType.AROMATIC == b_type and not single_bonds:
                conditions = [atom1 not in aromatic_index, atom2 not in aromatic_index,
                              atom_type_map[atom1] != 'O', atom_type_map[atom2] != 'O',
                              atom_type_map[atom1] != 'S', atom_type_map[atom2] != 'S',
                              atom_type_map[atom1] != 'N' or atom_bond_map[atom1] < 3,
                              atom_type_map[atom2] != 'N' or atom_bond_map[atom2] < 3]
                if all(conditions):
                    bond_img = symbol_img_dict['double']
                    aromatic_index.add(atom1)
                    aromatic_index.add(atom2)
            elif rdkit.Chem.rdchem.BondType.TRIPLE == b_type and not single_bonds:
                bond_img = symbol_img_dict['triple']
            add_bond(img, bond_img, mydegrees, position, length)

    # add atoms (to start from the Hydrogens, the atom index must be reversed)
    for i in reversed(range(len(mol.GetAtoms()))):
        atom = mol.GetAtoms()[i]
        atom = atom.GetSymbol()
        # add dimension to center with respect to the center of the blank image
        atom_x = atom_map[i][0] + dimension // 2
        atom_y = atom_map[i][1] + dimension // 2
        if atom not in symbol_img_dict:
            atom = 'other'
        add_image(img, symbol_img_dict[atom], (atom_x, atom_y))

    if rdkit_img:
        rdkit.Chem.Draw.MolToImageFile(mol, directory + os.sep + name + "_rdkit.png")
    if save:
        cv2.imwrite(directory + os.sep + name + ".png", img)
    print('\033[0;32m' + name + ' completed' + '\033[0;0;m')
    return img


def parse():
    # create a parser for command line
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Produce 2D icons of molecules from smiles.')
    main = parser.add_argument_group('[ Required ]')
    main.add_argument('SMILE',
                      metavar='SMILE string',
                      help='Smile of the molecule to produce the icon')

    optional = parser.add_argument_group('[ Optional ]')
    optional.add_argument("--name",
                          metavar='STR',
                          default='molecule_icon',
                          help='Name of the png output file')
    optional.add_argument("-d", '--directory',
                          metavar='FOLDER',
                          default=os.getcwd(),
                          help='Path to the folder to save the icon file ')
    optional.add_argument("--rdkit_draw",
                          action='store_true',
                          help='Use this flag to save also the rdkit 2D image of the molecule')
    optional.add_argument("-s", "--single_bond",
                          action='store_true',
                          help='Use this flag to draw single bonds only')
    optional.add_argument("--remove_H",
                          action='store_true',
                          help='Use this flag to remove the hydrogens from the structure')
    optional.add_argument("-v", "--verbose",
                          action='store_true',
                          help='Print the 2D coordinates of each atom')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    parsed = parse()
    icon_map = load_icons(atom_icon_dir)
    icon_print(parsed.SMILE, parsed.name, parsed.directory, parsed.rdkit_draw,
               parsed.single_bond, parsed.remove_H, parsed.verbose)

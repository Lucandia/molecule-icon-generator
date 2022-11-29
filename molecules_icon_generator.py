#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Luca Monari"
__credits__ = ["Luca Monari"]
__version__ = "1.0"
__email__ = "luca.monari@etu.unistra.fr"

import numpy as np
import svgwrite
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF
from pdf2image import convert_from_path # require poppler
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import math
import itertools
import argparse
import os
import colorsys
import warnings
warnings.filterwarnings("ignore")  # brute force approach to avoid decompression bomb warning by pdf2image and PIL


color_map = {"H": "#FFFFFF", "D": "#FFFFC0", "T": "#FFFFA0", "He": "#D9FFFF", "Li": "#CC80FF", "Be": "#C2FF00",
             "B": "#FFB5B5", "C": "#909090", "C-13": "#505050", "C-14": "#404040", "N": "#3050F8", "N-15": "#105050",
             "O": "#FF0D0D", "F": "#90E050", "Ne": "#B3E3F5", "Na": "#AB5CF2", "Mg": "#8AFF00", "Al": "#BFA6A6",
             "Si": "#F0C8A0", "P": "#FF8000", "S": "#FFFF30", "Cl": "#1FF01F", "Ar": "#80D1E3", "K": "#8F40D4",
             "Ca": "#3DFF00", "Sc": "#E6E6E6", "Ti": "#BFC2C7", "V": "#A6A6AB", "Cr": "#8A99C7", "Mn": "#9C7AC7",
             "Fe": "#E06633", "Co": "#F090A0", "Ni": "#50D050", "Cu": "#C88033", "Zn": "#7D80B0", "Ga": "#C28F8F",
             "Ge": "#668F8F", "As": "#BD80E3", "Se": "#FFA100", "Br": "#A62929", "Kr": "#5CB8D1", "Rb": "#702EB0",
             "Sr": "#00FF00", "Y": "#94FFFF", "Zr": "#94E0E0", "Nb": "#73C2C9", "Mo": "#54B5B5", "Tc": "#3B9E9E",
             "Ru": "#248F8F", "Rh": "#0A7D8C", "Pd": "#006985", "Ag": "#C0C0C0", "Cd": "#FFD98F", "In": "#A67573",
             "Sn": "#668080", "Sb": "#9E63B5", "Te": "#D47A00", "I": "#940094", "Xe": "#429EB0", "Cs": "#57178F",
             "Ba": "#00C900", "La": "#70D4FF", "Ce": "#FFFFC7", "Pr": "#D9FFC7", "Nd": "#C7FFC7", "Pm": "#A3FFC7",
             "Sm": "#8FFFC7", "Eu": "#61FFC7", "Gd": "#45FFC7", "Tb": "#30FFC7", "Dy": "#1FFFC7", "Ho": "#00FF9C",
             "Er": "#00E675", "Tm": "#00D452", "Yb": "#00BF38", "Lu": "#00AB24", "Hf": "#4DC2FF", "Ta": "#4DA6FF",
             "W": "#2194D6", "Re": "#267DAB", "Os": "#266696", "Ir": "#175487", "Pt": "#D0D0E0", "Au": "#FFD123",
             "Hg": "#B8B8D0", "Tl": "#A6544D", "Pb": "#575961", "Bi": "#9E4FB5", "Po": "#AB5C00", "At": "#754F45",
             "Rn": "#428296", "Fr": "#420066", "Ra": "#007D00", "Ac": "#70ABFA", "Th": "#00BAFF", "Pa": "#00A1FF",
             "U": "#008FFF", "Np": "#0080FF", "Pu": "#006BFF", "Am": "#545CF2", "Cm": "#785CE3", "Bk": "#8A4FE3",
             "Cf": "#A136D4", "Es": "#B31FD4", "Fm": "#B31FBA", "Md": "#B30DA6", "No": "#BD0D87", "Lr": "#C70066",
             "Rf": "#CC0059", "Db": "#D1004F", "Sg": "#D90045", "Bh": "#E00038", "Hs": "#E6002E", "Mt": "#EB0026",
             'other': '#f5c2cb', 'bond': '#575757'}


def hex_to_rgb(color):
    r = int(color[1:3], 16)
    g = int(color[3:5], 16)
    b = int(color[5:], 16)
    return r, g, b


def rgb_to_hex(color):
    r = hex(int(color[0]))[2:]
    g = hex(int(color[1]))[2:]
    b = hex(int(color[2]))[2:]
    check_len = []
    for i in [r, g, b]:
        if len(i) < 1:
            check_len.append('00')
        elif len(i) < 2:
            check_len.append('0'+i)
        else:
            check_len.append(i)
    return '#' + ''.join(check_len)


def shadow_color_correction(color, shadow_light):
    r, g, b = hex_to_rgb(color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    rgb = colorsys.hls_to_rgb(h, l*shadow_light, s)
    return rgb_to_hex(rgb)

def circ_post(degree, size, center):
    angle_rad = math.radians(degree)
    x = int(size*math.cos(angle_rad) + center[0])
    y = int(size*math.sin(angle_rad) + center[1])
    return x, y


def add_atom_svg(src, center, radius, color, shadow=True, shadow_curve=1.125, shadow_deg=45, shadow_light=0.35,
                 thickness=1/3):
    shadow_color = shadow_color_correction(color, shadow_light)
    src.add(src.circle(center=(center[0], center[1]), r=radius, fill=color, stroke=shadow_color, stroke_width=radius*thickness/2))
    if shadow:
        start_shade = circ_post(-shadow_deg, radius, center)
        end_shade = circ_post(-shadow_deg+180, radius, center)
        src.add(src.path(d=f'M{start_shade[0]},{start_shade[1]} A{radius},{radius} 0, 1,1 {end_shade[0]},{end_shade[1]} M{end_shade[0]},{end_shade[1]} A{radius*shadow_curve},{radius*shadow_curve} 0, 0,0 {start_shade[0]},{start_shade[1]} Z',
                         fill=shadow_color, stroke_width=0))
        # patch to cover line that appears in jpeg and png images with pdf2image
        start_patch = circ_post(-shadow_deg, radius-radius*thickness/4, center)
        end_patch = circ_post(-shadow_deg+180, radius-radius*thickness/4, center)
        src.add(src.path(d=f'M{start_patch[0]},{start_patch[1]} L {end_patch[0]},{end_patch[1]}',
                         stroke=color))
    return


def add_bond_svg(src, bond_type, x1, y1, x2, y2, line_thickness, bondcolor='#575757', shadow_light=0.35, shadow=True):
    start = np.array((x1, y1))
    end = np.array((x2, y2))
    d_space = int(line_thickness * 1)
    t_space = int(line_thickness * 2)
    # calculate the degree of the bond line, y axis is reversed in images
    radians = math.atan2(-y1 + y2, x1 - x2) + math.pi/2  # add 90 degree to make the angle perpendicular
    contour_color = shadow_color_correction(bondcolor, shadow_light)

    def dist_point(point, spacer):
        dist_x = int(math.cos(radians) * spacer)
        dist_y = int(math.sin(radians) * spacer)
        # y axis is reversed in images
        pt1 = point + np.array((dist_x, -dist_y))
        pt2 = point - np.array((dist_x, -dist_y))
        return pt1, pt2

    def add_bond(p, q):
        src.add(src.line(p.tolist(), q.tolist(), stroke=bondcolor, stroke_width=line_thickness,
                         stroke_linecap="round"))
    def add_bond_contour(p, q):
        src.add(src.line(p.tolist(), q.tolist(), stroke=contour_color, stroke_width=line_thickness+line_thickness/5,
                         stroke_linecap="round"))

    if bond_type == 2:
        start_1, start_2 = dist_point(start, d_space)
        end_1, end_2 = dist_point(end, d_space)
        if shadow:
            add_bond_contour(start_1, end_1)
            add_bond_contour(start_2, end_2)
        add_bond(start_1, end_1)
        add_bond(start_2, end_2)
    else:
        if shadow:
            add_bond_contour(start, end)
        add_bond(start, end)
    if bond_type == 3:
        start_1, start_2 = dist_point(start, t_space)
        end_1, end_2 = dist_point(end, t_space)
        if shadow:
            add_bond_contour(start_1, end_1)
            add_bond_contour(start_2, end_2)
        add_bond(start_1, end_1)
        add_bond(start_2, end_2)


def icon_print(SMILES, name='molecule_icon', directory=os.getcwd(), rdkit_png=False, rdkit_svg=False, single_bonds=False,
                remove_H=False, verbose=False, save_svg=True, save_png=True, save_jpeg=True, save_pdf=True,
               atom_color=color_map, position_multiplier=160, atom_radius=100, bw=False, shadow=True,
               black=False, thickness=1/4, shadow_light=0.35):
    if black:
        atom_color = {key: '#000000' for key in atom_color}
    elif bw:
        atom_color = {key: '#ffffff' for key in atom_color}
        shadow_light = 0
        shadow = False
        atom_color['bond'] = '#000000'

    if '.svg' in name:
        fullname = directory + os.sep + name
    else:
        fullname = directory + os.sep + name + '.svg'

    mol = Chem.MolFromSmiles(SMILES)
    if not remove_H:
        mol = Chem.AddHs(mol)
    AllChem.Compute2DCoords(mol)
    mol.GetConformer()
    mol = rdkit.Chem.Draw.rdMolDraw2D.PrepareMolForDrawing(mol)

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

    # get the maximum position from the atom (positive or negative)
    max_pos = max([abs(pos) for pos in list(itertools.chain(*list(atom_map.values())))])
    # to the max position we have to add width of the atom icons (*2 because the
    max_tot_pos = max_pos + atom_radius * 2
    # multiply for 2 because the max position is considered from the center of the image
    dimension = int(max_tot_pos * 2)
    # create the svg Drawing
    svg = svgwrite.Drawing(fullname, (str(dimension), str(dimension)))


    # add bond
    bonds_list = list(itertools.combinations(list(atom_map.keys()), 2))
    aromatic_index = set()
    bond_thickness = int(atom_radius * thickness)
    for couple in bonds_list:
        atom1 = couple[0]
        atom2 = couple[1]
        BOND = mol.GetBondBetweenAtoms(atom1, atom2)
        if BOND:
            x1 = atom_map[atom1][0] + dimension // 2
            y1 = atom_map[atom1][1] + dimension // 2
            x2 = atom_map[atom2][0] + dimension // 2
            y2 = atom_map[atom2][1] + dimension // 2
            b_type = BOND.GetBondType()
            bond_type = 1
            if rdkit.Chem.rdchem.BondType.DOUBLE == b_type and not single_bonds:
                bond_type = 2
            elif rdkit.Chem.rdchem.BondType.AROMATIC == b_type and not single_bonds:
                conditions = [atom1 not in aromatic_index, atom2 not in aromatic_index,
                              atom_type_map[atom1] != 'O', atom_type_map[atom2] != 'O',
                              atom_type_map[atom1] != 'S', atom_type_map[atom2] != 'S',
                              atom_type_map[atom1] != 'N' or atom_bond_map[atom1] < 3,
                              atom_type_map[atom2] != 'N' or atom_bond_map[atom2] < 3]
                if all(conditions):
                    bond_type = 2
                    aromatic_index.add(atom1)
                    aromatic_index.add(atom2)
            elif rdkit.Chem.rdchem.BondType.TRIPLE == b_type and not single_bonds:
                bond_type = 3
                aromatic_index.add(atom1)
                aromatic_index.add(atom2)
            add_bond_svg(svg, bond_type, x1, y1, x2, y2, bond_thickness, bondcolor=atom_color['bond'],
                         shadow_light=0.35,shadow=shadow)

    # add atoms (to start from the Hydrogens, the atom index must be reversed)
    for i in reversed(range(len(mol.GetAtoms()))):
        atom = mol.GetAtoms()[i]
        atom = atom.GetSymbol()
        # add dimension to center with respect to the center of the blank image
        atom_x = atom_map[i][0] + dimension // 2
        atom_y = atom_map[i][1] + dimension // 2
        if atom not in atom_color:
            atom = 'other'
        add_atom_svg(svg, (atom_x, atom_y), atom_radius, atom_color[atom], shadow=shadow, shadow_light=shadow_light,
                     thickness=thickness)

    if rdkit_png:
        rdkit.Chem.Draw.MolToFile(mol, directory + os.sep + name + "_rdkit.png")
    if rdkit_svg:
        drawer = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DSVG(300, 300)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        rdkit_svg_text = drawer.GetDrawingText()
        with open(directory + os.sep + name + "_rdkit.svg", 'w') as f:
            f.write(rdkit_svg_text)

    pdf_name = directory + os.sep + name + ".pdf"
    if save_pdf or save_png or save_jpeg:
        save_svg = True
        if save_png or save_jpeg:
            save_pdf = True
    if save_svg:
        svg.save()
    if save_pdf:
        svg.save()
        drawing = svg2rlg(fullname)
        renderPDF.drawToFile(drawing, pdf_name)
    if save_png or save_jpeg:
        pages = convert_from_path(pdf_name)
    if save_png:
        pages[0].save(directory + os.sep + name + '.png', 'PNG')
    if save_jpeg:
        pages[0].save(directory + os.sep + name + '.jpeg', 'JPEG')

    print('\033[0;32m' + fullname + ' completed' + '\033[0;0;m')
    return svg


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
    optional.add_argument("-a", '--atom_multiplier',
                          metavar='FLOAT',
                          type=float,
                          default=1,
                          help='Increase or decrease the atom size in the image')
    optional.add_argument("-p", '--position_multiplier',
                          metavar='FLOAT',
                          type=float,
                          default=1,
                          help='Increase or decrease the image size in the image')
    optional.add_argument("-d", '--directory',
                          metavar='FOLDER',
                          default=os.getcwd(),
                          help='Path to the folder to save the icon file ')
    optional.add_argument("--rdkit_svg",
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
    optional.add_argument("-bw", "--black_and_white",
                          action='store_true',
                          help='Print the image in black and white')
    optional.add_argument("--hide_shadows",
                          action='store_true',
                          help='Hide the shadows of the atoms')
    optional.add_argument('-b', "--black",
                          action='store_true',
                          help='Draw a black icon')
    optional.add_argument('-t', "--thickness",
                          type=float,
                          default=0.35,
                          help='Line thickness compared to atom size, range [0:1]')
    optional.add_argument("--shadow_light",
                          type=float,
                          default=0.35,
                          help='Select how dark the shadow should be in the range [0:1]')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    parsed = parse()
    icon_print(parsed.SMILE, name=parsed.name, directory=parsed.directory, rdkit_svg=parsed.rdkit_svg,
               single_bonds=parsed.single_bond, remove_H=parsed.remove_H, verbose=parsed.verbose, save_png=True,
               bw=parsed.black_and_white, position_multiplier=int(160*parsed.position_multiplier),
               atom_radius=int(100*parsed.atom_multiplier), shadow=not parsed.hide_shadows, black=parsed.black,
               shadow_light=parsed.shadow_light, thickness=parsed.thickness)

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
from pdf2image import convert_from_path  # require poppler
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdCoordGen
from rdkit.Chem import rdDepictor
import math
import itertools
import argparse
import os
import colorsys
import warnings

warnings.filterwarnings("ignore")  # brute force approach to avoid decompression bomb warning by pdf2image and PIL

# dictionary containing the default colour for each atom, according to CPK colour convention
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
             'other': '#f5c2cb', 'Bond': '#575757'}

# dictionary containing the value to multiply to the final radius of each atom
atom_resize = {'All atoms': 1.0, 'H': 1.0, 'D': 1.0, 'T': 1.0, 'He': 1.0, 'Li': 1.0, 'Be': 1.0, 'B': 1.0, 'C': 1.0,
               'C-13': 1.0,
               'C-14': 1.0, 'N': 1.0, 'N-15': 1.0, 'O': 1.0, 'F': 1.0, 'Ne': 1.0, 'Na': 1.0, 'Mg': 1.0, 'Al': 1.0,
               'Si': 1.0, 'P': 1.0, 'S': 1.0, 'Cl': 1.0, 'Ar': 1.0, 'K': 1.0, 'Ca': 1.0, 'Sc': 1.0, 'Ti': 1.0, 'V': 1.0,
               'Cr': 1.0, 'Mn': 1.0, 'Fe': 1.0, 'Co': 1.0, 'Ni': 1.0, 'Cu': 1.0, 'Zn': 1.0, 'Ga': 1.0, 'Ge': 1.0,
               'As': 1.0, 'Se': 1.0, 'Br': 1.0, 'Kr': 1.0, 'Rb': 1.0, 'Sr': 1.0, 'Y': 1.0, 'Zr': 1.0, 'Nb': 1.0,
               'Mo': 1.0, 'Tc': 1.0, 'Ru': 1.0, 'Rh': 1.0, 'Pd': 1.0, 'Ag': 1.0, 'Cd': 1.0, 'In': 1.0, 'Sn': 1.0,
               'Sb': 1.0, 'Te': 1.0, 'I': 1.0, 'Xe': 1.0, 'Cs': 1.0, 'Ba': 1.0, 'La': 1.0, 'Ce': 1.0, 'Pr': 1.0,
               'Nd': 1.0, 'Pm': 1.0, 'Sm': 1.0, 'Eu': 1.0, 'Gd': 1.0, 'Tb': 1.0, 'Dy': 1.0, 'Ho': 1.0, 'Er': 1.0,
               'Tm': 1.0, 'Yb': 1.0, 'Lu': 1.0, 'Hf': 1.0, 'Ta': 1.0, 'W': 1.0, 'Re': 1.0, 'Os': 1.0, 'Ir': 1.0,
               'Pt': 1.0, 'Au': 1.0, 'Hg': 1.0, 'Tl': 1.0, 'Pb': 1.0, 'Bi': 1.0, 'Po': 1.0, 'At': 1.0, 'Rn': 1.0,
               'Fr': 1.0, 'Ra': 1.0, 'Ac': 1.0, 'Th': 1.0, 'Pa': 1.0, 'U': 1.0, 'Np': 1.0, 'Pu': 1.0, 'Am': 1.0,
               'Cm': 1.0, 'Bk': 1.0, 'Cf': 1.0, 'Es': 1.0, 'Fm': 1.0, 'Md': 1.0, 'No': 1.0, 'Lr': 1.0, 'Rf': 1.0,
               'Db': 1.0, 'Sg': 1.0, 'Bh': 1.0, 'Hs': 1.0, 'Mt': 1.0, 'other': 1.0}


def hex_to_rgb(color):
    """It takes a hexadecimal color string and returns a rgb tuple.

    Parameters
    ----------
    color : string
        The hexadecimal color code, it starts with a '#'.

    Returns
    -------
    tuple
        the RGB values of the hexadecimal color code.

    """
    r = int(color[1:3], 16)
    g = int(color[3:5], 16)
    b = int(color[5:], 16)
    return r, g, b


def rgb_to_hex(color):  # based on https://stackoverflow.com/questions/3380726/converting-an-rgb-color-tuple-to-a-hexidecimal-string
    """It takes the rgb tuple and returns the hexadecimal string of the color.

    Parameters
    ----------
    color : tuple
        a tuple of three integers, each between 0 and 255, representing the rgb values.

    Returns
    -------
    string
        The hexadecimal color code, it starts with a '#'.

    """
    return '#%02x%02x%02x' % color


def shadow_color_correction(color, light_multiplier):
    """Given a color and a light multiplier, return the hexadecimal color code with the corrected light.

    Parameters
    ----------
    color : tuple
        The rgb color tuple.
    light_multiplier : float
        The value to multiply to decrease the light.

    Returns
    -------
    string
        The hexadecimal color code, it starts with a '#'.

    """
    r, g, b = hex_to_rgb(color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    rgb = colorsys.hls_to_rgb(h, l * light_multiplier, s)
    return rgb_to_hex(rgb)


def circ_post(degree, size, center):
    """This function takes in an angle, size, and center and returns the x-y coordinates of the point
    on the circumference.

    Parameters
    ----------
    degree : float
        The degree angle where you want to calculate the coordinate on the circumference.
    size : float
        The radius of the circumference.
    center : tuple
        Tuple containing the x-y coordinates of the center of the circumference.

    Returns
    -------
    tuple
        A tuple with the x-y coordinates of the point on the circumference at the defined angle.

    """
    angle_rad = math.radians(degree)
    x = int(size * math.cos(angle_rad) + center[0])
    y = int(size * math.sin(angle_rad) + center[1])
    return x, y


def add_atom_svg(src, center, radius, color, shadow=True, shadow_curve=1.125, shadow_deg=45, shadow_light=0.35,
                 thickness=1 / 3):
    """It draws a circle, filled with the color, in the svg object. A shadow can be drawn on the circle.

    Parameters
    ----------
    src : svgwrite.Drawing
        The svg object.
    center : tuple
        The center of the atom in x-y coordinates.
    radius : float
        The radius of the circle.
    color : string
        The hex color to fill the circle.
    shadow : bool, default: True
        Whether to add a shadow or not on the circle.
    shadow_curve : float, default: 1.125
        The curve of the shadow. 1.125 is a good value.
    shadow_deg : float, default: 45
        The angle in the degree of the shadow start.
    shadow_light : float, default: 0.35
        The lightness of the shadow. 0 is black, 1 is white.
    thickness : float, default: 1/3
        The thickness of the border of the circle.

    """
    shadow_color = shadow_color_correction(color, shadow_light)
    src.add(src.circle(center=(center[0], center[1]), r=radius, fill=color, stroke=shadow_color,
                       stroke_width=radius * thickness / 2))
    if shadow:
        start_shade = circ_post(-shadow_deg, radius, center)
        end_shade = circ_post(-shadow_deg + 180, radius, center)
        src.add(src.path(
            d=f'M{start_shade[0]},{start_shade[1]} A{radius},{radius} 0, 1,1 {end_shade[0]},{end_shade[1]} M{end_shade[0]},{end_shade[1]} A{radius * shadow_curve},{radius * shadow_curve} 0, 0,0 {start_shade[0]},{start_shade[1]} Z',
            fill=shadow_color, stroke_width=0))
        # patch to cover line that appears in jpeg and png images with pdf2image
        start_patch = circ_post(-shadow_deg, radius - radius * thickness / 4, center)
        end_patch = circ_post(-shadow_deg + 180, radius - radius * thickness / 4, center)
        src.add(src.path(d=f'M{start_patch[0]},{start_patch[1]} L {end_patch[0]},{end_patch[1]}',
                         stroke=color))


def add_bond_svg(src, bond_type, x1, y1, x2, y2, line_thickness, bondcolor='#575757', shadow_light=0.35):
    """It adds a line as a bond to an SVG image

    Parameters
    ----------
    src : svg.Drawing
        The svg object.
    bond_type : int
        The type of the bond. 1 stands for single bond, 2 stands for double bond, 3 stands for triple bond.
    x1 : float
        x-coordinate of the first atom.
    y1 : float
        y-coordinate of the first atom.
    x2 : float
        x-coordinate of the second atom.
    y2 : float
        y-coordinate of the second atom.
    line_thickness: float
        The thickness of the bond line.
    bondcolor : string, default: '#575757'
        The hex color of the bond.
    shadow_light : float, default: 0.35
        The lightness of the shadow. 0 is black, 1 is white.

    """
    start = np.array((x1, y1))
    end = np.array((x2, y2))
    d_space = line_thickness * 1.5
    t_space = line_thickness * 2.5
    # calculate the degree of the bond line, y axis is reversed in images
    radians = math.atan2(-y1 + y2, x1 - x2) + math.pi / 2  # add 90 degree to make the angle perpendicular
    contour_color = shadow_color_correction(bondcolor, shadow_light)

    def dist_point(point, spacer):
        """This function takes a point and a spacer distance. It returns two points that have a distance
        equal to the spacer on the direction defined by the radians parent variable

        Parameters
        ----------
        point : tuple
            A tuple of the x and y coordinates of the point
        spacer : float
            The distance between the points

        Returns
        -------
        tuple
            Two points with a spacer distance from point on the radian angle direction.

        """
        dist_x = int(math.cos(radians) * spacer)
        dist_y = int(math.sin(radians) * spacer)
        # y-axis is reversed in images
        pt1 = point + np.array((dist_x, -dist_y))
        pt2 = point - np.array((dist_x, -dist_y))
        return pt1, pt2

    def add_bond(p, q):
        """It adds a bond between two points. It uses many parent variables.

        Parameters
        ----------
        p : tuple
            The first point.
        q : tuple
            The second point.

        """
        src.add(src.line(p.tolist(), q.tolist(), stroke=bondcolor, stroke_width=line_thickness,
                         stroke_linecap="round"))

    def add_bond_contour(p, q):
        """It takes two points, p and q, and adds a bigger and darker bond as a contour

        Parameters
        ----------
        p : tuple
            The first point.
        q : tuple
            The second point.

        """
        src.add(src.line(p.tolist(), q.tolist(), stroke=contour_color, stroke_width=line_thickness + line_thickness,
                         stroke_linecap="round"))

    if bond_type == 2:
        start_1, start_2 = dist_point(start, d_space)
        end_1, end_2 = dist_point(end, d_space)
        add_bond_contour(start_1, end_1)
        add_bond_contour(start_2, end_2)
        add_bond(start_1, end_1)
        add_bond(start_2, end_2)
    else:
        add_bond_contour(start, end)
        add_bond(start, end)
    if bond_type == 3:
        start_1, start_2 = dist_point(start, t_space)
        end_1, end_2 = dist_point(end, t_space)
        add_bond_contour(start_1, end_1)
        add_bond_contour(start_2, end_2)
        add_bond(start_1, end_1)
        add_bond(start_2, end_2)


def parse_structure(mol, position_multiplier, verbose=False):
    """It takes a molecule object and a position multiplier, and returns a list of dictionary with the atoms and bonds
    maps.

    Parameters
    ----------
    mol : mol object
        A molecule object from the rdkit library.
    position_multiplier : int, default: 160
        This is the distance between atoms.
    verbose : bool, optional
        Prints out the atoms coordinates.

    Returns
    -------
    tuple
        A list of dictionaries: the atom_map, atom_type_map, atom_bond_map

    """
    atom_map = dict()
    atom_type_map = dict()
    atom_bond_map = dict()
    for i, atom in enumerate(mol.GetAtoms()):
        positions = mol.GetConformer().GetAtomPosition(i)
        atom_symbol = atom.GetSymbol()
        x = int(positions.x * position_multiplier)
        y = int(positions.y * position_multiplier)
        z = int(positions.z * position_multiplier)
        atom_map[i] = (x, -y)  # y-axis of an image is the inverse
        atom_bonds = rdkit.Chem.rdchem.Atom.GetBonds(atom)
        atom_bond_map[i] = len(atom_bonds)
        atom_type_map[i] = atom_symbol
        if verbose:
            print(atom_symbol, x, y, z)
    return atom_map, atom_type_map, atom_bond_map


def icon_print(SMILES, name='molecule_icon', directory=os.getcwd(), rdkit_png=False, rdkit_svg=False,
               single_bonds=False, remove_H=False, verbose=False, save_svg=True, save_png=False, save_jpeg=False,
               save_pdf=False, atom_color=color_map, position_multiplier=160, atom_radius=100, bw=False, shadow=True,
               black=False, thickness=1/4, shadow_light=0.35, nice_conformation=True, radius_multi=atom_resize):
    """This function takes a SMILES string and returns an icon of the molecule, in format PNG, SVG, JPEG, and PDF.

    Parameters
    ----------
    SMILES : string
        The SMILES string of the molecule you want to draw.
    name : string, default: 'molecule_icon'
        The name of the file to be saved.
    directory : string, default: os.getcwd()
        The directory to save the image in.
    rdkit_png : bool, optional
        If True, will use RDKit to generate a PNG image of the default structure.
    rdkit_svg : bool, optional
        If True, will use RDKit to generate an SVG image of the default structure.
    single_bonds : bool, optional
        If True, all bonds will be single bonds.
    remove_H : bool, optional
        Remove all non-chiral hydrogen from the molecule.
    verbose : bool, optional
        Prints out the atoms coordinates.
    save_svg : bool, default: True
        Save the SVG icon format.
    save_png : bool, default: False
        Save the SVG, PDF and PNG icon formats.
    save_jpeg : bool, default: False
        Save the SVG, PDF and JPEG icon formats.
    save_pdf : bool, default: False
        Save the SVG and PDF icon formats.
    atom_color : dictionary, default: color_map
        a dictionary of atom colors. The keys are the atom symbols, and the values are the hex colors.
    position_multiplier : int, default: 160
        This is the distance between atoms.
    atom_radius : int, default: 100
        The radius of the atoms in the icon.
    bw : bool, optional
        Draw a black image with white atoms.
    shadow : bool, optional
        Whether to add a shadow to the image or not.
    black : bool, optional
        If True, the molecule will be drawn in full black.
    thickness : float, default: 1/4
        The thickness of the bonds compared to the atom radius.
    shadow_light : float, default: 0.35
        How light the shadow is. 0.35 is a good value.
    nice_conformation : bool, default: True
        If True, the molecule will be put into a nice conformation.
    radius_multi : dictionary, default: atom_resize
        A dictionary containing the multiplier for each atom. It multiplies the atom radius.

    Returns
    -------
    svg.Drawing
        The svg object.

    """

    if black:
        atom_color = {key: '#000000' for key in atom_color}  # set all colors black
    elif bw:  # set all atoms black
        atom_color = {key: '#ffffff' for key in atom_color}
        shadow_light = 0
        shadow = False
        atom_color['bond'] = '#000000'

    if '.svg' in name:
        fullname = directory + os.sep + name
    else:
        fullname = directory + os.sep + name + '.svg'

    mol = Chem.MolFromSmiles(SMILES)  # read the molecule
    if not remove_H:
        mol = Chem.AddHs(mol)  # add Hydrogens
    else:
        mol = Chem.RemoveHs(mol)  # remove not chiral Hydrogens
    if nice_conformation:
        rdDepictor.SetPreferCoordGen(True)  # rdCoordGen conformation as default
        rdCoordGen.AddCoords(mol)  # better conformation for macrocycle
    else:
        rdDepictor.SetPreferCoordGen(False)  # rdkit conformation default
        AllChem.Compute2DCoords(mol)  # canonical rdkit conformation

    mol.GetConformer()  # load the conformer
    mol = rdkit.Chem.Draw.rdMolDraw2D.PrepareMolForDrawing(mol)  # clean the conformer

    atom_map, atom_type_map, atom_bond_map = parse_structure(mol, position_multiplier, verbose=False)

    # get the maximum position from the atom (positive or negative)
    max_pos = max([abs(pos) for pos in list(itertools.chain(*list(atom_map.values())))])
    # to the max position we have to add width of the atom icons (*2 for the diameter)
    max_tot_pos = max_pos + atom_radius * max(radius_multi.values()) * 2
    # multiply for 2 because the max position is considered from the center of the image
    dimension = int(max_tot_pos * 2)
    # create the svg Drawing
    svg = svgwrite.Drawing(fullname, (str(dimension), str(dimension)))

    # add bond
    bonds_list = list(itertools.combinations(list(atom_map.keys()), 2))  # generate a list of bonds couple
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
            add_bond_svg(svg, bond_type, x1, y1, x2, y2, bond_thickness, bondcolor=atom_color['Bond'],
                         shadow_light=shadow_light)

    # add atoms (to start from the Hydrogens, the atom index must be reversed)
    for i in reversed(range(len(mol.GetAtoms()))):
        atom = mol.GetAtoms()[i]
        atom = atom.GetSymbol()
        # add dimension to center with respect to the center of the blank image
        atom_x = atom_map[i][0] + dimension // 2
        atom_y = atom_map[i][1] + dimension // 2
        if atom not in atom_color:
            atom = 'other'
        corrected_radius = atom_radius * radius_multi[atom]  # resize the atom dimension
        add_atom_svg(svg, (atom_x, atom_y), corrected_radius, atom_color[atom], shadow=shadow,
                     shadow_light=shadow_light,
                     thickness=thickness)

    if rdkit_png:
        rdkit.Chem.Draw.MolToFile(mol, directory + os.sep + name + "_rdkit.png")
    if rdkit_svg:
        drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
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
               bw=parsed.black_and_white, position_multiplier=int(160 * parsed.position_multiplier),
               atom_radius=int(100 * parsed.atom_multiplier), shadow=not parsed.hide_shadows, black=parsed.black,
               shadow_light=parsed.shadow_light, thickness=parsed.thickness)

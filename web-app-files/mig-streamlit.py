# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 02:56:08 2022

@author: monar
"""

import streamlit as st
import cirpy
from cirpy import Molecule
import os
import cv2
st.write(os.getcwd())
st.write(os.listdir())
import importlib
mig = importlib.import_module("molecules_icon_generator")
import molecules_icon_generator as mig


if __name__ == "__main__":
    
    # select the folder with the atom icons:
    atom_icon_dir = "../base-icons"
    icon_map = mig.load_icons(atom_icon_dir)
    
    st.write('''
    # Molecule-icons generator!
    ''')
    
    input_type = st.selectbox("Create your icon by", 
                 ['smiles', 'name', 'cas_number', 'stdinchi'], 
                 help= 'Chose the input info of your moleculs')
    
    input_string = st.text_input('Input informations', "CC(=O)Nc1ccc(cc1)O" )
    # name = st.text_input('create by NAME:', "paracetamol" )
    
    if input_type == 'name':
        input_string = cirpy.resolve(input_string, 'smiles')
    mol = Molecule(input_string)   
    iupac = mol.iupac_name
    smiles = mol.smiles

    single_bonds = st.checkbox('Draw just single_bonds')
    remove_H = st.checkbox('remove all Hydrogens') 
    rdkit_draw = st.checkbox('show rdkit structure') 
    
    filename = 'molecular-icon' + '.png'
    
    image = mig.icon_print(smiles, name = 'molecular-icon', rdkit_img = rdkit_draw, 
                            single_bonds = single_bonds, remove_H = remove_H, save=True,
                            symbol_img_dict = icon_map)
    
    im_rgba = cv2.cvtColor(image, cv2.COLOR_BGRA2RGBA)
    img_list = [im_rgba]
    caption_list = 'Iupac name: ' + iupac
    
    if rdkit_draw:
        rdkit_img = cv2.imread(os.getcwd() + os.sep + "molecular-icon_rdkit.png", cv2.IMREAD_UNCHANGED)
        img_list.append(rdkit_img)
        caption_list.append('Rdkit 2D conformation')
        
    st.image(img_list, caption = caption_list, channels = 'RGBA', use_column_width=True,)
    
    with open(os.getcwd() + os.sep + filename, "rb") as file:
        btn = st.download_button( label="Download icon",
                                 data=file,
                                 file_name=filename,
                                 mime="image/png" )

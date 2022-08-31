# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 02:56:08 2022

@author: monar
"""

import streamlit as st
from molecules_icon_generator import icon_print 
import cirpy
from cirpy import Molecule
import os
import cv2


if __name__ == "__main__":
    
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
    
    options = st.multiselect(
     'Additional options',
     ['just single bonds', 'remove H'] )
    
    single_bonds  = "just single bonds" in options
    reomve_H  = "remove H" in options
    
    filename = 'molecular-icon' + '.png'
    
    image = icon_print(smiles, name = 'molecular-icon', rdkit_img = False, 
                   single_bonds = single_bonds, remove_H = reomve_H, save=True)
    
    im_rgba = cv2.cvtColor(image, cv2.COLOR_BGRA2RGBA)
    st.image(im_rgba, caption = 'Iupac name: ' + iupac, channels = 'RGBA')
    
    st.write(f'''
    Image saved in:
    {os.getcwd() + os.sep + filename}
    '''    )
    
    with open(os.getcwd() + os.sep + filename, "rb") as file:
        btn = st.download_button( label="Download icon",
                                 data=file,
                                 file_name=filename,
                                 mime="image/png" )



    

# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 02:56:08 2022

@author: monar
"""

import streamlit as st
import cirpy
import rdkit.Chem as Chem
from cirpy import Molecule
import base64
import molecule_icon_generator as mig
import warnings
import json
import os
import shutil
import time

# brute force approach to avoid decompression bomb warning by pdf2image and PIL
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
warnings.filterwarnings("ignore")
warnings.simplefilter('ignore', Image.DecompressionBombWarning)

emoji_license = """The emoji are published under the Creative Commons Share Alike 
                    License 4.0 (CC BY-SA 4.0). Icons containing emojis are distributed under the
                    CC BY-SA 4.0 license (https://creativecommons.org/licenses/by-sa/4.0/#)"""

def render_svg(svg):
    """Renders the given svg string."""
    b64 = base64.b64encode(svg.encode('utf-8')).decode("utf-8")
    html = r'<img width="300px" height="300px" src="data:image/svg+xml;base64,%s"/>' % b64
    st.write(html, unsafe_allow_html=True)


def upload_setting_button():
    """Allow to upload setting"""
    st.session_state['upload_setting'] = True
    return


if __name__ == "__main__":
    # initialize session state
    if 'color_dict' not in st.session_state:
        st.session_state['color_dict'] = mig.color_map.copy()
    if 'resize_dict' not in st.session_state:
        st.session_state['resize_dict'] = mig.atom_resize.copy()
    if 'reset_color' not in st.session_state:
        st.session_state['reset_color'] = False
    if 'reset_size' not in st.session_state:
        st.session_state['reset_size'] = False
    if 'last_atom_size' not in st.session_state:
        st.session_state['last_atom_size'] = None
    if 'last_atom_color' not in st.session_state:
        st.session_state['last_atom_color'] = None
    if 'upload_setting' not in st.session_state:
        st.session_state['upload_setting'] = False
    if 'emoji_dict' not in st.session_state:
        st.session_state['emoji_dict'] = dict()
    # if 'sizes_percentage' not in st.session_state:
    #     st.session_state['sizes_percentage'] = 1

    loading_err = KeyError("""The app encountered a problem in initializing the data. 
Try to reload the page. If the problem persists, contact monarluca@gmail.com""")
    # loading the color, resize and emoji dictionary
    if 'color_dict' in st.session_state:
        new_color = st.session_state['color_dict']
    else:
        st.exception(loading_err)
        print([i for i in st.session_state])
        st.session_state['color_dict'] = mig.color_map.copy()
        new_color = st.session_state['color_dict']
    if 'resize_dict' in st.session_state:
        resize = st.session_state['resize_dict']
    else:
        st.exception(loading_err)
        print([i for i in st.session_state])
        st.session_state['resize_dict'] = mig.atom_resize.copy()
        resize = st.session_state['resize_dict']
    if 'emoji_dict' in st.session_state:
        emoji = st.session_state['emoji_dict']
    else:
        st.exception(loading_err)
        print([i for i in st.session_state])
        st.session_state['emoji_dict'] = dict()
        emoji = st.session_state['emoji_dict']

    # check if the color/resize dictionary have been reset
    if 'atom_color_select' in st.session_state and 'color_picker' in st.session_state and st.session_state[
        'reset_color']:
        st.session_state.color_picker = new_color[st.session_state.atom_color_select]
        st.session_state['last_atom_color'] = None
        st.session_state['reset_color'] = False
    last_atom_color = st.session_state['last_atom_color']
    if 'atom_size_select' in st.session_state and 'sizes_percentage' in st.session_state and st.session_state[
        'reset_size']:
        st.session_state['last_atom_size'] = None
        st.session_state['reset_size'] = False
    last_atom_size = st.session_state['last_atom_size']

    # setting header, description and citation
    st.set_page_config(page_title="Molecule icons")
    st.header('''
    Molecule Icon Generator!
    ''')
    st.write('''
    Generate icons of molecules from SMILES, Name, Cas-number, Inchi, InChIKey, load your molecule file or convert a
    list of SMILES.
    ''')
    st.warning('A new release is about to be published. If you have saved settings, they could not be compatible with the new release.', icon="⚠️")
    st.markdown('''
For more options and information, check out the 
[GitHub repository](https://github.com/lmonari5/molecule-icon-generator.git).\\
[DOI](https://doi.org/10.5281/zenodo.7388429): 10.5281/ZENODO.7388429.
       ''')

    # check the time of the parsing
    start_time = time.time()
    # select the input type
    input_type = st.selectbox("Create your icon by",
                              ['name', 'smiles', 'load file', 'cas_number', 'stdinchi', 'stdinchikey', 'smiles list'],
                              help='Choose the input info of your molecule. SMILES inputs are faster')
    if input_type != 'smiles' and input_type != 'smiles list' and input_type != 'load file':
        # end of parsing time
        end_time = time.time()
        if end_time - start_time > 1.5:
            # The SMILES of the molecule is parsed from the cirpy library.
            st.warning('If the app is slow, use SMILES input.', icon="⚠️")
    # default input for each input_type except 'load file'
    def_dict = {'name': 'paracetamol',
                'smiles': "CC(=O)Nc1ccc(cc1)O",
                'cas_number': '103-90-2',
                'stdinchi': 'InChI=1S/C8H9NO2/c1-6(10)9-7-2-4-8(11)5-3-7/h2-5,11H,1H3,(H,9,10)',
                'stdinchikey': 'RZVAJINKPMORJF-UHFFFAOYSA-N'}

    # load the molecule input
    if input_type == 'load file':
        load_sdf = True
        smiles_list = False
        mol_file = st.file_uploader('Load mol file', type=['sdf', 'mol2', 'pdb'],
                                    help='Load a file with a molecule structure (supported sdf, mol2 and pdb files)')
        st.write('The icon results are highly dependent on the molecule coordinates.')
    # load the smile list input
    elif input_type == 'smiles list':
        smiles_list = True
        load_sdf = False
        smiles_list_file = st.file_uploader('Load the list of SMILES as txt file', type='txt',
                                            help='Load a text file containing a SMILES in each line (encoding utf-8)')
    else:
        load_sdf = False
        smiles_list = False
        input_string = st.text_input(input_type + ' :', def_dict[input_type],
                                     help=f'Insert the corresponding {input_type} of your molecule')

    # load settings
    load_settings = st.checkbox('Upload previous settings',
                                help='''If you have saved a "molecule_icon_settings.json" file, you can upload it 
                                     and use the same settings. You can save the settings with the button at the 
                                     end of the page'''
                                )  # using checkbox to save space, in case doesn't want to save settings
    if load_settings:
        saved_setting = st.file_uploader("Upload previous settings (optional):", on_change=upload_setting_button,
                                         type='json',
                                         help='''If you have saved a "molecule_icon_settings.json" file, you can upload it 
                                                 and use the same settings. You can save the settings with the button at the 
                                                 end of the page''')
        if saved_setting is not None and st.session_state['upload_setting']:
            # To read file as bytes:
            session_old = json.load(saved_setting)
            for key, value in session_old.items():
                st.session_state[key] = value
            st.session_state['upload_setting'] = False
            st.experimental_rerun()

    # select conformation and output format
    col1, col2 = st.columns(2, gap='medium')
    # select whether to se a 2D or 3D conformation
    with col1:
        dimension = st.selectbox("Build a structure:", ['2D', '3D', '3D interactive'], key='dimension_type',
                                 help='Use the classical 2D visualization, or try to visualize a 3D structure')
    # select the download format
    with col2:
        if dimension == '3D interactive':  # plotly doesn't support saving plot in pdf from the camera
            formats = ('svg', 'png', 'jpeg')
        else:
            formats = ('svg', 'png', 'jpeg', 'pdf')
        forms = [False, False, False, False]
        img_format = st.selectbox('Download file format:', formats, key='img_format',
                                  help="""The native file format is svg. Using png and jpeg formats could slow down 
                                       the app""")
        if dimension != '3D interactive':
            for ind, img_form in enumerate(('svg', 'png', 'jpeg', 'pdf')):
                if img_form == img_format:
                    forms[ind] = True

    # set the parameters for a 2D/3D structure
    conf = False
    dimension_3 = False
    rand_seed = -1
    f_field = None
    activate_emoji = False
    if input_type != 'load file':
        col1, col2, col3 = st.columns(3, gap='medium')
        if '3D' in dimension:
            dimension_3 = True
            with col1:
                f_field = st.selectbox("Select force field", ['UFF', 'MMFF'], key='force_filed',
                                       help="Force fields currently supported: UFF (faster) and MMFF (more accurate)")
            with col2:
                rand_seed = st.number_input('Random seed', min_value=-1, value=1000, key='3D_random_seed',
                                            help='''Choose the random seed to generate the molecule. A value of -1 will 
                                                    generate a random structure every time the app is running''')
            if dimension == '3D':
                with col3:
                    st.write('\n')
                    st.write('\n')
                    activate_emoji = st.checkbox('Use emoji', key='use_emoji', help=emoji_license)
        else:
            with col1:
                conf = not st.checkbox('Switch conformation', key='switch_conf', value=False)
            with col3:
                if dimension == '2D':
                    activate_emoji = st.checkbox('Use emoji', key='use_emoji', help=emoji_license)
    else:
        col1, col2, col3 = st.columns(3, gap='medium')
        with col3:
            activate_emoji = st.checkbox('Use emoji', key='use_emoji', help=emoji_license)

    # add common checkbox
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        remove_H = st.checkbox('Remove Hydrogen', key='removeH')
    with col2:
        rdkit_label = 'Show RDKIT'
        if activate_emoji:
            rdkit_label = 'Show RDKIT index'
        rdkit_draw = st.checkbox(rdkit_label, key='show_rdkit')
    with col3:
        if dimension != '3D interactive':
            h_shadow = st.checkbox('Hide shadows', key='remove_shadow')
    with col4:
        if dimension != '3D interactive':
            single_bonds = st.checkbox('Single bonds', key='single_bonds')

    # catch error when using the cirpy library
    try:
        if input_type == 'smiles':  # if the input is a smile, use it directly ignoring the cirpy smiles
            smiles = input_string
        elif load_sdf:
            pass
        elif smiles_list:
            pass
        else:
            if input_type == 'name':
                input_string = cirpy.resolve(input_string, 'smiles')
            mol = Molecule(input_string)
            smiles = mol.smiles
    except Exception as e:
        st.write(f'''
        The cirpy python library is not able to resolve your input {input_type}.
        You can use SMILES to skip the cirpy library.
        ''')
        st.stop()

    # try to build the mol structure
    try:
        molecules = list()
        if load_sdf and mol_file:
            if '.sdf' in mol_file.name:
                molecule = Chem.MolFromMolBlock(mol_file.read(), sanitize=False)
            elif '.mol2' in mol_file.name:
                molecule = Chem.MolFromMol2Block(mol_file.read(), sanitize=False)
            elif '.pdb' in mol_file.name:
                molecule = Chem.MolFromPDBBlock(mol_file.read(), sanitize=False)
            mig.partial_sanitize(molecule)
            molecules.append(molecule)
            if remove_H:
                molecule = Chem.RemoveHs(molecule)
        elif load_sdf and not mol_file:
            st.stop()
        elif smiles_list and smiles_list_file:
            for line in smiles_list_file.readlines():
                smiles = str(line, 'utf-8').strip()
                if smiles == '': # skip empty lines
                    continue
                try:
                    molecule = mig.parse_structure(smiles, remove_H=remove_H, nice_conformation=conf,
                                                   dimension_3=dimension_3,
                                                   force_field=f_field, randomseed=rand_seed)
                    molecules.append(molecule)
                except:
                    st.write(f'Failed building SMILES: {str(smiles)}')
        elif smiles_list and not smiles_list_file:
            st.stop()
        else:
            molecule = mig.parse_structure(smiles, remove_H=remove_H, nice_conformation=conf, dimension_3=dimension_3,
                                           force_field=f_field, randomseed=rand_seed)
            molecules.append(molecule)
    except Exception as e:
        error_txt = f'''
            Rdkit failed in building the structure of the molecule. Full error:
            {e}.'''
        if input_type != 'smiles':
            error_txt += f' \n  Try to use the SMILES instead of {input_type} as input.'
        elif input_type == 'smiles' and 'H' in smiles:
            error_txt += '  \n  If you have written Hydrogen atoms in the SMILES, try to remove them.'
        error_txt += '  \n  If the problem persists, contact monarluca@gmail.com'
        st.error(error_txt)
        st.stop()

    # add emojis
    if activate_emoji:
        atom_and_index = list(range(molecules[0].GetNumAtoms())) + list(mig.atom_resize.keys())
        col1, col2, col3 = st.columns(3, gap='medium')
        with col1:
            st.write('\n')
            atom_emoji = st.selectbox(
                'Select atom index or element:', atom_and_index,
                help="""The atom index depends on rdkit parsing. You can see the atom indexes using 'Show RDKIT index'.
                     To reset all the emojis, choose 'All atoms' without indicating the unicode""")
        with col2:
            if atom_emoji in emoji:
                def_value = emoji[atom_emoji][0]
            else:
                def_value = ''
            emoji_code = st.text_input(f'Emoji unicode from https://openmoji.org/:', value=def_value,
                                       help='''Insert unicode character according to the open-emoji project
                                                https://openmoji.org/''')
        if atom_emoji == 'All atoms':
            for key in atom_and_index:
                emoji[key] = [emoji_code, 1]  # set coloured because black emoji have transparency.
        else:
            emoji[atom_emoji] = [emoji_code, 1]  # set coloured because black emoji have transparency.
        with col3:
            st.write('\n')
            st.write('\n')
            st.write('\n')
            periodic_emoji = st.button('Emoji periodic table', key='periodic_emoji_but',
                                       help=""""In the emoji periodic table, the atoms are replaced by emojis which
                                        represent the element. It was created by Nicola Ga-stan and Andrew White.
                                        To reset all the emojis, select 'All atoms' in the selection window without
                                        unicode""")
            if periodic_emoji:
                emoji = {k: [v, 1] for k, v in mig.emoji_periodic_table.items()}
                st.session_state['emoji_dict'] = emoji
    else:
        emoji = None

    # change the color of the icon (single atoms, all atoms or bonds)
    col1, col2, col3 = st.columns(3, gap='medium')
    with col1:
        atom_color = st.selectbox(
            'Change the color:',
            sorted(list(mig.color_map.keys())), key='atom_color_select')
    with col2:
        if last_atom_color != atom_color:
            def_value = new_color[atom_color]
        else:
            if 'sizes_percentage' in st.session_state:
                def_value = st.session_state.color_picker
            else:
                def_value = new_color[atom_color]
        new_color[atom_color] = st.color_picker(f' Pick {atom_color} color', def_value,
                                                key="color_picker")
        if atom_color == "All icon":  # set all icon same color
            unicolor = new_color[atom_color]
            for key in new_color:  # have to modify directly new_color, which is saved in session state
                new_color[key] = unicolor
        if atom_color == "All atoms":  # set all atoms same color
            unicolor = new_color[atom_color]
            bond_color = new_color['Bond']
            for key in new_color:  # have to modify directly new_color, which is saved in session state
                if key == 'Bond':
                    continue
                new_color[key] = unicolor
    st.session_state['last_atom_color'] = atom_color
    with col3:
        st.write('\n')
        st.write('\n')
        if st.button('Reset colours', help='Reset colours as default CPK', key='reset_color_but'):
            st.session_state['color_dict'] = mig.color_map.copy()
            new_color = st.session_state['color_dict']
            st.session_state['reset_color'] = True
            st.experimental_rerun()

    # change the size of the icon (single atoms, all atoms or bonds)
    col1, col2, col3 = st.columns(3, gap='medium')
    with col1:
        atom_size = st.selectbox(
            'Change the size:',
            ['All atoms'] + sorted(list(mig.atom_resize.keys())), key='atom_size_select')
    with col2:
        if last_atom_size != atom_size:
            def_value = int(resize[atom_size] * 100)
        else:
            if 'sizes_percentage' in st.session_state:
                def_value = st.session_state.sizes_percentage
            else:
                def_value = 100
        resize[atom_size] = st.number_input(f'{atom_size} radius percentage (%)', min_value=0, value=def_value,
                                            key='sizes_percentage', format="%d",
                                            help='''Increase or decrease the size of one specific atom''') / 100
        st.session_state['last_atom_size'] = atom_size
    with col3:
        st.write('\n')
        st.write('\n')
        if st.button('Reset atoms size', help='Reset size to 100% for all atoms', key='reset_size_but'):
            st.session_state['resize_dict'] = mig.atom_resize.copy()
            resize = st.session_state['resize_dict']
            st.session_state['reset_size'] = True
            st.experimental_rerun()
    icon_size = resize['All atoms'] * 100

    # change multiplier, thickness and shadow darkness
    col1, col2, col3 = st.columns(3)
    with col1:
        pos_multi = st.slider('Image size multiplier', 0, 900, 300, key='size_multi_slider',
                              help='''Multiply the position of the atoms with respect to the 2D structure.
                              A higher multiplier leads to higher resolution.''')
    with col2:
        thickness = st.slider('Thickness', 0.0, 0.6, 0.2, key='thickness_slider',
                              help='''Bond and stroke thicknesses over the atom radius.''')
    with col3:
        if dimension != '3D interactive':
            shadow_light = st.slider('Shadow/outline light', 0.0, 1.0, 1 / 3, key='outline_slider',
                                     help='''Regulate the brightness of the shadow''')
        else:
            resolution = st.slider('3D Graph resolution', 0, 100, 30, key='resolution_slider',
                                   help='''Resolution of the bond and atoms 3D mesh''')

    # correct the size of the image according to rdkit default conformation or coordGen
    if conf:
        img_multi = pos_multi
    else:
        img_multi = pos_multi * 2 / 3

    # Add rotation axis to observe the molecule from different angles
    if dimension == '3D':
        col1, col2, col3 = st.columns(3)
        with col1:
            x_rot = st.slider('x-axis rotation', 0, 360, 0, key='x_rot_slider',
                              help='''Multiply the position of the atoms with respect to the 2D structure.
                                  A higher multiplier leads to higher resolution.''')
        with col2:
            y_rot = st.slider('y-axis rotation', 0, 360, 0, key='y_rot_slider',
                              help='''Bond and stroke thicknesses over the atom radius.''')
        with col3:
            z_rot = st.slider('z-axis rotation', 0, 360, 0, key='z_rot_slider',
                              help='''Regulate the brightness of the shadow''')
        rot_angles = (x_rot, y_rot, z_rot)
    else:
        rot_angles = (0, 0, 0)

    # try to produce the image
    if smiles_list:
        direct ='icons_list'
        if direct in os.listdir():
            shutil.rmtree(direct)
        os.mkdir(direct)
    else:
        direct = os.getcwd()
    for index, mol in enumerate(molecules):
        try:
            if dimension == '3D interactive':
                graph = mig.graph_3d(mol, name=str(index), rdkit_svg=rdkit_draw, radius_multi=resize, directory=direct,
                                     atom_color=new_color, pos_multi=img_multi, atom_radius=icon_size,
                                     thickness=thickness, resolution=resolution)
                filename = direct + os.sep + str(index) + ".html"
                # set camera to download the image format selected
                config = {'toImageButtonOptions': {
                    'label': f'Download {img_format}',
                    'format': img_format,  # one of png, svg, jpeg, webp
                    'filename': 'molecule-icon',
                    'scale': 1  # Multiply title/legend/axis/canvas sizes by this factor
                }
                }
                graph.write_html(filename, config=config)
            else:
                mig.icon_print(mol, name=str(index), rdkit_svg=rdkit_draw, pos_multi=img_multi, directory=direct,
                               single_bonds=single_bonds, atom_radius=icon_size, radius_multi=resize,
                               atom_color=new_color, shadow=not h_shadow,
                               save_svg=forms[0], save_png=forms[1], save_jpeg=forms[2], save_pdf=forms[3],
                               thickness=thickness, shadow_light=shadow_light, rotation=rot_angles, emoji=emoji)
        except Exception as e:
            error_txt = f'''
                The program failed at producing the Image/Graph. Full error:
                {e}.'''
            if dimension != '3D interactive':
                if img_format != 'svg':
                    error_txt += '  \n  Try to use the svg format.'
            error_txt += '  \n  If the problem persists, contact monarluca@gmail.com'
            st.error(error_txt)
            if not smiles_list:
                st.stop()

    # show the download button for single molecule
    if not smiles_list:
        # download the html-graph or the image
        if dimension == '3D interactive':
            with open(filename, "rb") as file:
                btn = st.download_button(label="Download 3D plot",
                                         data=file,
                                         file_name="molecule-icon-graph.html",
                                         help=f'''Download the html graph and open it in your browser to take 
                                              {img_format} snapshots with the camera button''')
        else:
            filename = '0.' + img_format
            with open(filename, "rb") as file:
                btn = st.download_button(label="Download icon",
                                         data=file,
                                         file_name='molecule_icon.' + img_format,
                                         mime=f"image/{img_format}")
    else:
        shutil.make_archive('molecules-icons', 'zip', direct)
        filename = 'molecules-icons.zip'
        with open(filename, "rb") as file:
            btn = st.download_button(label="Download icons zip",
                                     data=file,
                                     file_name='molecules-icons.zip',
                                     mime=f"image/{img_format}")

        # cite
        st.markdown('''
        Thanks for using the Molecules icons generators! 
        Please [cite this work](https://doi.org/10.5281/zenodo.7388429): 
        [![DOI](https://zenodo.org/badge/530035520.svg)](https://zenodo.org/badge/latestdoi/530035520)
        ''')

    # add preview for single image
    if not smiles_list:
        if dimension == '3D interactive':
            show_graph = st.checkbox('Show molecule 3D plot (the app will be slower)', )
            if show_graph:
                st.write(f'''
                    Plotly graph preview (download the graph and open it in your browser to take {img_format} 
                    snapshots with the camera button):
                    ''')
        else:
            st.write('''
                Image SVG preview:
                ''')
        col1, col2 = st.columns(2)
        with col1:
            if dimension == '3D interactive':
                if show_graph:
                    graph.update_layout(height=300)
                    st.plotly_chart(graph, config=config, use_container_width=True)
            else:
                f = open("0.svg", "r")
                svg_text = f.read()
                render_svg(svg_text)
            with col2:
                if rdkit_draw:
                    f = open("0_rdkit.svg", "r")
                    svg_text = f.read()
                    render_svg(svg_text)
                    with open("0_rdkit.svg", "rb") as file:
                        btn = st.download_button(label="Download RDKIT icon",
                                                 data=file,
                                                 file_name='molecule_icon_rdkit.svg',
                                                 mime=f"image/{img_format}")


    # save settings and allow download
    with open('molecule_icon_settings.json', 'w') as settings:
        # cannot save session state as it is, I have to convert it to a dictionary
        session_dict = {key: st.session_state[key] for key in st.session_state if 'but' not in key}
        json.dump(session_dict, settings)
    with open("molecule_icon_settings.json", "rb") as file:
        btn = st.download_button(
            label="Download settings",
            data=file,
            file_name="molecule_icon_settings.json",
            mime="application / json",
            help='''Save the current settings (e.g. atoms color, atoms radius), 
                    so you can easily reload them in you refresh the page!'''
        )

    # donate
    st.markdown('I enjoy working on this project in my free time, especially at night. '
                'If you want to support me with a coffee, just ['
                'click on the caffeine:](https://www.paypal.com/donate/?hosted_button_id=V4LJ3Z3B3KXRY)')

    html_string = '''<form action="https://www.paypal.com/donate" method="post" target="_top">
<input type="hidden" name="hosted_button_id" value="V4LJ3Z3B3KXRY" />
<input type="image" src="https://pics.paypal.com/00/s/MTRiMDhhZmMtMTE0YS00MWNjLWJmYjItM2RhZjAwNmFjZGVh/file.PNG" 
border="0" name="submit" title="PayPal - The safer, easier way to pay online!" alt="Donate with PayPal button" />
<img alt="" border="0" src="https://www.paypal.com/en_IT/i/scr/pixel.gif" width="1" height="1" />
</form>
'''
    st.markdown(html_string, unsafe_allow_html=True)

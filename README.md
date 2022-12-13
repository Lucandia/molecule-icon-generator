# Molecular-icons Generator!
<img src="example/paracetamol.jpeg" width=300 height=300>
Generate nice icons of molecules from SMILES.

This program follows the topology of SMILES chemical structures and creates an icon.
The atoms' colours are inspired by the [CPK colouring convention](https://sciencenotes.org/molecule-atom-colors-cpk-colors/)

## Try the web app:

[Molecules-icons generator web app](https://molecule-icon-generator.streamlit.app/) powered by streamlit

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://molecule-icon-generator.streamlit.app/)

## Cite this work:

```
Lmonari5. Lmonari5/Molecule-Icon-Generator; Zenodo, 2022. https://doi.org/10.5281/ZENODO.7388429.
```
[![DOI](https://zenodo.org/badge/530035520.svg)](https://zenodo.org/badge/latestdoi/530035520)

## Run it on local:

### Step 1: clone the repository

```
git clone https://github.com/lmonari5/molecule-icon-generator.git
```

### Step 2: install packages and requirements

For Linux:

```
xargs -a packages.txt sudo apt-get install 
pip install -r requirements.txt
```

### Step 3: Have Fun

- Run the code from the command line:

 ```
  python molecule_icon_generator.py "CC(=O)Nc1ccc(cc1)O" --name paracetamol --rdkit_svg
 ```

- Or from the python interpreter:

 ```
 import molecule_icon_generator as mig 
 molecule = mig.parse_structure("CC(=O)Nc1ccc(cc1)O", remove_H=False)
 mig.icon_print(molecule, name = 'paracetamol', rdkit_svg = True, single_bonds = False, remove_H = False, verbose=False)
 ```

- Or use the Streamlit functionalities! From the terminal run:

 ```
python -m streamlit run streamlit_app.py
 ```
 
## Donate

I enjoy working on this project in my free time, especially at night. If you want to support me with a coffee, just [click here!](https://www.paypal.com/donate/?hosted_button_id=V4LJ3Z3B3KXRY)


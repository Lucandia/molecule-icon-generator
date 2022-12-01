# Molecular-icons Generator!
<img src="example/paracetamol.jpeg" width=300 height=300>
Generate nice icons of molecules from SMILES.

This program follows the topology of SMILES chemical structures and creates an icon.
The atoms' colours are inspired by the [CPK colouring convention](https://sciencenotes.org/molecule-atom-colors-cpk-colors/)

## Try the web app:

[Molecules-icons generator web app](https://molecule-icon-generator.streamlit.app/) powered by streamlit

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://molecule-icon-generator.streamlit.app/)


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
  python molecules_icon_generator.py "CC(=O)Nc1ccc(cc1)O" --name paracetamol --rdkit_svg
 ```

- Or from the python interpreter:

 ```
 from molecules_icon_generator import icon_print 
 icon_print("CC(=O)Nc1ccc(cc1)O", name = 'paracetamol', rdkit_svg = True, single_bonds = False, remove_H = False, verbose=False)
 ```

## Donate

I enjoy working on this project in my free time, especially at night. If you want to support me with a coffee, just [click on the caffeine:](https://www.paypal.com/donate/?hosted_button_id=V4LJ3Z3B3KXRY)

<form action="https://www.paypal.com/donate" method="post" target="_top">
<input type="hidden" name="hosted_button_id" value="V4LJ3Z3B3KXRY" />
<input type="image" src="https://pics.paypal.com/00/s/MTRiMDhhZmMtMTE0YS00MWNjLWJmYjItM2RhZjAwNmFjZGVh/file.PNG" border="0" name="submit" title="PayPal - The safer, easier way to pay online!" alt="Donate with PayPal button" />
<img alt="" border="0" src="https://www.paypal.com/en_IT/i/scr/pixel.gif" width="1" height="1" />
</form>


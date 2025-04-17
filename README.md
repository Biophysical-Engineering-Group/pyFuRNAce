![Alt Text](https://github.com/Biophysical-Engineering-Group/pyFuRNAce/blob/main/pyfurnace/app/static/logo_text.png?raw=true)

**pyFuRNAce** is an open-source Python package and web-based design engine for creating complex RNA nanostructures using the co-transcriptional RNA origami approach. It streamlines the entire design pipeline — from structural motif assembly to sequence generation and primer design — into an intuitive, user-friendly platform.

**Website:** [pyfurnace.de](http://pyfurnace.de)\
**Documentation & Source:** _Coming Soon_\
**PyPI:** _Coming Soon_

---

## 🚀 Features

- 🧩 **Motif-based assembly:** Build RNA structures using a rich, expandable library of motifs including stems, dovetails, kissing loops, aptamers, and ribozymes.
- 🎨 **GUI & Real-time 3D Visualization:** Interactive blueprint editor and real-time 3D rendering via Streamlit and oxView.
- 🔄 **Integrated Workflow:** Design, generate, convert, and prepare your RNA origami in one unified interface.
- 🧬 **Sequence Generation & Optimization:** Built-in support for sequence folding (Revolvr + ViennaRNA).
- 🧪 **Primer & Template Design:** Includes tools for DNA conversion, promoter addition, and primer calculations.
- 💻 **Python Scripting API:** Automate complex designs or build at scale using a programmable interface.

---

## 📦 Installation

Note: pyFuRNAce requires Python 3.10 or later.

### Note: Local Installation always includes the GUI

### Instal from PyPI 
Note: The PyPi version doesn't include the oxDNA analysis tools (OAT). The OAT package is used to convert 3D structures in oxDNA format (natively supported in pyFuRNAce) from/to PDB files and write oxDNA force files. That's why they have to be installed separately.

```bash
pip install pyfurnace
# install oxDNA analysis tools for pdb conversion
pip install pip install git+https://github.com/lorenzo-rovigatti/oxDNA.git#subdirectory=analysis
```

### Install from GitHub
To install the latest version of pyFuRNAce, clone the repository and install all the required dependencies:

```bash
pip install "git+https://github.com/Biophysical-Engineering-Group/pyFuRNAce.git#egg=pyfurnace[oat]"
```

## 🖥️ Running the Web Application

To run the web application locally, clone the repository and install the required dependencies:

```bash
python -m pyfurnace
```

This will lunch the GUI in your default web browser. 
You can also use the hosted version at [pyfurnace.de](http://pyfurnace.de).
The WebApp is built using Streamlit and can be run locally or on a server. You can access the webapp directly at 
[pyfurnace.streamlit.app](https://pyfurnace.streamlit.app).

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://pyfurnace.streamlit.app/)


## 🎛 Modules

1. Design: Create and edit RNA structure blueprints. Visualize assembled structures in 3D.
Define custom motifs via GUI or scripting.
2. Generate: use inverse folding (Revolvr) to produce RNA sequences matching the target structure. Evaluate folding energies and structural ensemble diversity.
3. Convert: Translate RNA sequences to DNA templates. Add transcriptional promoters (e.g., T7). Analyze sequence properties (e.g., GC content, dimers).
4. Prepare: Design PCR primers with melting temperature calculations. Generate input files for molecular dynamics simulations with oxRNA.

## 🧑‍💻 Using the Python API

```python
import pyfurnace as pf

line1 = [pf.TetraLoop(),
        pf.Stem(7),
        pf.Dovetail(-2, up_cross=False),
        pf.Stem(6),
        pf.KissingDimer(),
        pf.Stem(6),
        pf.Dovetail(-2, up_cross=False),
        pf.Stem(7),
        pf.TetraLoop(True),
        ]

line2 = [pf.TetraLoop(),
        pf.Stem(7),
        pf.Dovetail(-2, down_cross=False),
        pf.Stem(10),
        pf.start_end_stem(),
        pf.Stem(10),
        pf.Dovetail(-2, down_cross=False),
        pf.Stem(7),
        pf.TetraLoop(True),
        ]

origami = pf.Origami(line1, line2, aling='center')

print(origami)
print(origami.structure)
print(origami.sequence)
```

#### -> Output:
```
╭CGNKNKNNN──SS──NNKNNNAA┼─NNNNNN╯╭─ANNKNNN──SS──NKNKNNNUU╮
│  ┊┊┊┊┊┊┊  ┊┊  ┊┊┊┊┊┊  │ ┊┊┊┊┊┊ │  ┊┊┊┊┊┊  ┊┊  ┊┊┊┊┊┊┊  │
╰UUNKNKNNN──SS╮╭NNKNNNA─╯╭NNNNNN─┼AANNKNNN──SS╮╭NKNKNNNGC╯
              ││         ╰───────╯            ││
          ╭───╯│                         ╭────╯│
          │╭───╯                         │╭────╯
          ↓↓                             ↓↓
╭CGNNKNNKN╯╰SS──NKNNNKNNNN─3 5─NNKNNNKNNN╯╰SS──NNKNNNNUU╮
│  ┊┊┊┊┊┊┊  ┊┊  ┊┊┊┊┊┊┊┊┊┊     ┊┊┊┊┊┊┊┊┊┊  ┊┊  ┊┊┊┊┊┊┊  │
╰UUNNKNNKN──SS──NKNNNKNNNN─────NNKNNNKNNN──SS──NNKNNNNGC╯

((((((((((((((((((..[[[[[[.))))))))(((((((....)))))))(((((((((....)))))))))))))))))))(((((((((((((((((((....)))))))(((((((((....)))))))))((((((..]]]]]].))))))))))))))))))
NNKNNNKNNNSSNNNKNNAANNNNNNANNKNNNSSNKNKNNNUUCGNNNKNKNSSNNKNNNNUUCGNNNNKNNSSNNNKNNNKNNNNNNKNNNKNSSNKNNKNNUUCGNNKNNKNSSNNNKNKNUUCGNKNKNNNSSNNKNNNAANNNNNNANNNKNNSSNKNNNKNNNN
```


### Example: generate a simple DAE RNA origami with 3 helices (120° angle)
```python
ori = pf.simple_origami(helices=[120], use_angles=True)
```

## 📚 Examples

Explore tutorials and example notebooks in the /examples directory (or website when available).

## 📜 License

Code is licensed under the GNU General Public License v3.0 ([GPL-3.0](https://www.gnu.org/licenses/gpl-3.0.en.html))

[![License: GPL-3.0](https://img.shields.io/badge/License-GPL%20v3-lightgrey.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

## 🧠 Citation

If you use pyFuRNAce in your research, please cite:
Monari L., Braun I., Poppleton E., Göpfrich K. (2025). PyFuRNAce: An integrated design engine for RNA origami. (In preparation)

## 🙏 Acknowledgements

Supported by the ERC Starting Grant “ENSYNC”, DFG, HFSP, and Max Planck Society. Developed by Luca Monari, Ina Braun, Erik Poppleton, and Kerstin Göpfrich.

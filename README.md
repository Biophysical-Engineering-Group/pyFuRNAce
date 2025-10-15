![Alt Text](https://github.com/Biophysical-Engineering-Group/pyFuRNAce/blob/main/pyfurnace/app/static/logo_text.png?raw=true)

**pyFuRNAce** is an open-source Python package and web-based design engine for creating complex RNA nanostructures using the co-transcriptional RNA origami approach. It streamlines the entire design pipeline — from structural motif assembly to sequence generation and primer design — into an intuitive, user-friendly platform.

[![Downloads/MONTH](https://static.pepy.tech/personalized-badge/pyfurnace?period=month&units=INTERNATIONAL_SYSTEM&left_color=GREY&right_color=GREEN&left_text=Downloads/Month)](https://pepy.tech/projects/pyfurnace) [![PyPI Downloads](https://static.pepy.tech/personalized-badge/pyfurnace?period=total&units=INTERNATIONAL_SYSTEM&left_color=GREY&right_color=RED&left_text=Downloads)](https://pepy.tech/projects/pyfurnace)

#### pyFuRNAce design API coverage:
[![codecov](https://codecov.io/gh/Biophysical-Engineering-Group/pyFuRNAce/graph/badge.svg?token=EQX6C990WP)](https://codecov.io/gh/Biophysical-Engineering-Group/pyFuRNAce)

<!-- **Documentation & Source:** _Coming Soon_\ -->
**WebApp:** [pyfurnace.de](http://pyfurnace.de)\
**GitHub:** [Biophysical-Engineering-Group/pyFuRNAce](https://github.com/Biophysical-Engineering-Group/pyFuRNAce)\
**PyPI:** [pyfurnace](https://pypi.org/project/pyfurnace/)\
**Script API examples:** [Code Examples](https://github.com/Biophysical-Engineering-Group/pyFuRNAce/tree/main/examples)\
**Documentation:** [Read the Docs](https://pyfurnace.readthedocs.io/en/latest/)

![Alt Text](https://github.com/Biophysical-Engineering-Group/pyFuRNAce/blob/main/vid/demo_1min.gif?raw=true)

---

## 🚀 Features

- 🧩 **Motif-based assembly:** Build RNA structures using a rich, expandable library of motifs including stems, dovetails, kissing loops, aptamers, and ribozymes.
- 🎨 **GUI & Real-time 3D Visualization:** Interactive blueprint editor and real-time 3D rendering via Streamlit and oxView.
- 🔄 **Integrated Workflow:** Design, generate, convert, and prepare your RNA origami in one unified interface.
- 🧬 **Sequence Generation & Optimization:** Built-in support for sequence folding (Revolvr + ViennaRNA).
- 🧪 **Primer & Template Design:** Includes tools for DNA conversion, promoter addition, and primer calculations.
- 💻 **Python Scripting API:** Automate complex designs or build at scale using a programmable interface.

### System requirements

- Python 3.10 or later
- Operating System: Linux, macOS, or Windows
- dependencies are listed in `requirements.txt`; with stable versions recommended in `stable_requirements.txt`

---

## 📦 Installation

Typical installation time: 5 minutes

### Pyfurnace

Note: pyFuRNAce requires Python 3.10 or later, and the local installation always includes the GUI via Streamlit.

#### Install from PyPI
Install or upgrade the latest stable version of pyFuRNAce from PyPI using pip (and OAT from GitHub):

```bash
pip install --upgrade pyfurnace
```

#### Install from GitHub
Install the latest development version of pyFuRNAce (and OAT) directly from the GitHub repository:

```bash
pip install "git+https://github.com/Biophysical-Engineering-Group/pyFuRNAce.git"
```

### Extra Dependencies

**Note**: The **pyFuRNAce installation does not include the oxDNA analysis tools (OAT)**.

The OAT package is used to convert 3D structures from/to PDB files and write oxDNA force files.
To install the OAT package, you can use the following command (the **git** command is required, you can install it via anaconda with `conda install git`):

```bash
pip install "git+https://github.com/lorenzo-rovigatti/oxDNA.git#subdirectory=analysis"
```

To run sequence generation, a `Perl` interpreter is required by the Revolvr script from ROAD. It is usually installed by default on most systems. If not, you can install it via Anaconda with `conda install bioconda-legacy::perl`.

## 🖥️ Running the Web Application

To run the web application locally after installation, use the following command:

```bash
pyfurnace
```

Or alternatively `python -m pyfurnace`.

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
                        ╭───────╮
╭CGNNNKNNN──SS──NKNNNNAA┼─NNNNNN╯╭─ANNKNNN──SS──NNNKNNNUU╮
│  ┊┊┊┊┊┊┊  ┊┊  ┊┊┊┊┊┊  │ ┊┊┊┊┊┊ │  ┊┊┊┊┊┊  ┊┊  ┊┊┊┊┊┊┊  │
╰UUNNNKNNN──SS╮╭NKNNNNA─╯╭NNNNNN─┼AANNKNNN──SS╮╭NNNKNNNGC╯
              ││         ╰───────╯            ││
          ╭───╯│                         ╭────╯│
          │╭───╯                         │╭────╯
          ↑↓                             ↑↓
╭CGNKNKNNN╯╰SS──NKNNNKNNNN─3 5─NNNKNNNKNN╯╰SS──NKNNKNNUU╮
│  ┊┊┊┊┊┊┊  ┊┊  ┊┊┊┊┊┊┊┊┊┊     ┊┊┊┊┊┊┊┊┊┊  ┊┊  ┊┊┊┊┊┊┊  │
╰UUNKNKNNN──SS──NKNNNKNNNN─────NNNKNNNKNN──SS──NKNNKNNGC╯
((((((((((((((((((..[[[[[[.))))))))(((((((....)))))))(((((((((....)))))))))))))))))))(((((((((((((((((((....)))))))(((((((((....)))))))))((((((..]]]]]].))))))))))))))))))
NNNKNNNKNNSSNNNKNNAANNNNNNANNKNNNSSNNNKNNNUUCGNNNKNNNSSNKNNKNNUUCGNNKNNKNSSNNKNNNKNNNNNNNKNNNKNSSNNNKNKNUUCGNKNKNNNSSNNNKNNNUUCGNNNKNNNSSNKNNNNAANNNNNNANNNNKNSSNKNNNKNNNN
```

### 📚 Examples

Explore tutorials and example notebooks in the examples directory.

## 📜 License

Code is licensed under the GNU General Public License v3.0 ([GPL-3.0](https://www.gnu.org/licenses/gpl-3.0.en.html))

[![License: GPL-3.0](https://img.shields.io/badge/License-GPL%20v3-lightgrey.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)


## 📜 Logo Copyright & Permission
The pyFuRNAce logo is original artwork created by Luca Monari.

Luca Monari retains full copyright ownership of the logo.

Permission is granted to reproduce the logo in academic, scientific, or educational materials that reference or use pyFuRNAce, including publications, presentations, and documentation (whether commercial or non-commercial).

Any other use of the logo requires prior written permission from the copyright holder.

## 🧠 Citation

If you use pyFuRNAce in your research, please cite:
Monari, L., Braun, I., Poppleton, E. & Göpfrich, K. PyFuRNAce: An integrated design engine for RNA origami (2025) [doi:10.1101/2025.04.17.647389](https://doi.org/10.1101/2025.04.17.647389).

## 🙏 Acknowledgements

This work was supported by the ERC Starting Grant ENSYNC (No. 101076997) and Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under CRC 392 and CRC 1638. This work was supported by a Research Grant from HFSP (Ref.-No: RGP003/2023, DOI: https://doi.org/10.52044/HFSP.RGP0032023.pc.gr.168589). The authors thank the Max Planck Society for access to computational resources and the Alfried Krupp von Bohlen und Halbach Foundation. E.P. was supported through state funds approved by the State Parliament of Baden-Württemberg for the Innovation Campus Health + Life Science Alliance Heidelberg Mannheim. We thank Cody Geary for his feedback on the user interface design. We thank Dominic Kempf for his feedback on software development and testing.

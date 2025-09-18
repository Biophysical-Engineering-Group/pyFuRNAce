.. _features:

Features
========

- 🧩 **Motif-based assembly:** Build RNA structures using a rich, expandable library of motifs including stems, dovetails, kissing loops, aptamers, and ribozymes.
- 🎨 **GUI & Real-time 3D Visualization:** Interactive blueprint editor and real-time 3D rendering via Streamlit and oxView.
- 🔄 **Integrated Workflow:** Design, generate, convert, and prepare your RNA origami in one unified interface.
- 🧬 **Sequence Generation & Optimization:** Built-in support for sequence folding (Revolvr + ViennaRNA).
- 🧪 **Primer & Template Design:** Includes tools for DNA conversion, promoter addition, and primer calculations.
- 💻 **Python Scripting API:** Automate complex designs or build at scale using a programmable interface.


.. _installation:

Installation
============

**Note**: pyFuRNAce requires Python 3.10 or later. Local installation includes the GUI via Streamlit.

Install and/or upgrade from PyPI:

.. code-block:: bash

   pip install --upgrade pyfurnace

Install the latest development version from GitHub:

.. code-block:: bash

   pip install "git+https://github.com/Biophysical-Engineering-Group/pyFuRNAce.git"

Perl (**required by Revolvr for sequence optimisation**). Assuming you have an anaconda environment, you can install Perl via conda:

.. code-block:: bash

   conda install bioconda-legacy::perl


Install OAT for oxDNA structure conversion:

.. code-block:: bash

   pip install "git+https://github.com/lorenzo-rovigatti/oxDNA.git#subdirectory=analysis"


.. _webapp:

Web Application
===============

To run locally, use the command in the terminal:

.. code-block:: bash

   pyfurnace

This launches the GUI in your browser. Or alternatively, run the Streamlit app with:

.. code-block:: bash

   python -m pyfurnace

Or use the hosted version: `pyfurnace.de <http://pyfurnace.de>`_ or `pyfurnace.streamlit.app <https://pyfurnace.streamlit.app>`_

.. image:: https://static.streamlit.io/badges/streamlit_badge_black_white.svg
   :target: https://pyfurnace.streamlit.app


.. _modules:

Modules
^^^^^^^

- **Design:** Create and edit RNA blueprints. Visualize in 3D.
- **Generate:** Use inverse folding to produce target-matching sequences.
- **Convert:** Convert to DNA templates, add promoters, analyze properties.
- **Prepare:** Design primers and generate oxDNA simulation inputs.


.. _api_usage:

Using the Python API
====================

.. code-block:: python

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


.. _citation:

Citation
========

If you use pyFuRNAce in your research, please cite:

Monari, L., Braun, I., Poppleton, E. & Göpfrich, K. PyFuRNAce: An integrated design engine for RNA origami. (2025) `doi:10.1101/2025.04.17.647389 <https://doi.org/10.1101/2025.04.17.647389>`_


.. _license:

License
=======

GPL-3.0 License — `GNU General Public License v3.0 <https://www.gnu.org/licenses/gpl-3.0.en.html>`_

.. image:: https://img.shields.io/badge/License-GPL%20v3-lightgrey.svg
   :target: https://www.gnu.org/licenses/gpl-3.0.en.html


.. _acknowledgements:

Acknowledgements
================

This work was supported by the ERC Starting Grant ENSYNC (No. 101076997) and Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under CRC 392 and CRC 1638. This work was supported by a Research Grant from HFSP (Ref.-No: RGP003/2023, DOI: https://doi.org/10.52044/HFSP.RGP0032023.pc.gr.168589). The authors thank the Max Planck Society for access to computational resources and the Alfried Krupp von Bohlen und Halbach Foundation. E.P. was supported through state funds approved by the State Parliament of Baden-Württemberg for the Innovation Campus Health + Life Science Alliance Heidelberg Mannheim. We thank Dr.~Cody Geary for his feedback on the user interface design. We thank Dominic Kempf for his feedback on software development and testing.

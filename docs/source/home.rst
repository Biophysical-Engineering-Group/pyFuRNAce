.. _home:

Home
====

.. image:: https://github.com/Biophysical-Engineering-Group/pyFuRNAce/blob/main/pyfurnace/app/static/logo_text.png?raw=true
   :alt: pyFuRNAce logo

**pyFuRNAce** is an open-source Python package and web-based design engine for creating complex RNA nanostructures using the co-transcriptional RNA origami approach. It streamlines the entire design pipeline — from structural motif assembly to sequence generation and primer design — into an intuitive, user-friendly platform.

**WebApp:** `pyfurnace.de <http://pyfurnace.de>`_  
**GitHub:** `Biophysical-Engineering-Group/pyFuRNAce <https://github.com/Biophysical-Engineering-Group/pyFuRNAce>`_  
**PyPI:** `pyfurnace <https://pypi.org/project/pyfurnace/>`_  
**Script API examples:** `Code Examples <https://github.com/Biophysical-Engineering-Group/pyFuRNAce/tree/main/examples>`_
**Documentation:** `Read the Docs <https://pyfurnace.readthedocs.io/en/latest/>`_

.. image:: https://github.com/Biophysical-Engineering-Group/pyFuRNAce/blob/main/pyfurnace/app/static/walkthrough_1min.gif?raw=true
   :alt: pyFuRNAce demo


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

Install from PyPI:

.. code-block:: bash

   pip install pyfurnace

Install the latest development version from GitHub:

.. code-block:: bash

   pip install "git+https://github.com/Biophysical-Engineering-Group/pyFuRNAce.git"

Install OAT for oxDNA structure conversion:

.. code-block:: bash

   pip install "git+https://github.com/lorenzo-rovigatti/oxDNA.git#subdirectory=analysis"

Perl (required by Revolvr):

.. code-block:: bash

   conda install bioconda-legacy::perl


.. _webapp:

Web Application
===============

To run locally:

.. code-block:: bash

   python -m pyfurnace

This launches the GUI in your browser.

Or use the hosted version: `pyfurnace.de <http://pyfurnace.de>`_ or `pyfurnace.streamlit.app <https://pyfurnace.streamlit.app>`_

.. image:: https://static.streamlit.io/badges/streamlit_badge_black_white.svg
   :target: https://pyfurnace.streamlit.app


.. _modules:

Modules
^^^^^^^

1. **Design:** Create and edit RNA blueprints. Visualize in 3D.
2. **Generate:** Use inverse folding to produce target-matching sequences.
3. **Convert:** Convert to DNA templates, add promoters, analyze properties.
4. **Prepare:** Design primers and generate oxDNA simulation inputs.


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

Supported by the ERC Starting Grant “ENSYNC”, DFG, HFSP, and Max Planck Society.
Developed by Luca Monari, Ina Braun, Erik Poppleton, and Kerstin Göpfrich.


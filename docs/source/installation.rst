.. _installation:

Install and run
===============

.. _anaconda: https://www.anaconda.com/
.. _miniconda: https://www.anaconda.com/docs/getting-started/miniconda/main

**Note**: pyFuRNAce requires Python 3.10 or later. Local installation includes the graphical interface via Streamlit.
If you want to run sequence optimization locally, you also need to install Perl.

**Follow this instructions in the terminal/command line.**

It's recommended to use `Anaconda`_ or `Miniconda`_ to manage create a Python environment:

.. code-block:: bash

   conda create -n pyfurnace python
   conda activate pyfurnace

Install and/or upgrade pyFuRNAce from PyPI:

.. code-block:: bash

   pip install --upgrade pyfurnace

Install the latest development version from GitHub:

Install Perl via Conda (**required by Revolvr for sequence optimisation**):

.. code-block:: bash

   conda install -c bioconda perl


Install OAT for oxDNA structure conversion:

.. code-block:: bash

   pip install "git+https://github.com/lorenzo-rovigatti/oxDNA.git#subdirectory=analysis"


To launch the graphical interface, run:

.. code-block:: bash

   pyfurnace

Alternatively, you can also run `python -m pyfurnace`.

If you want to use pyFuRNAce as a Python package, you can import it in your scripts:

.. code-block:: python

   import pyfurnace as pf

See `API Reference <api.html>`_ for more information.

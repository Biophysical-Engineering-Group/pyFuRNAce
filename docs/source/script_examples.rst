.. _script_examples:

Python Scripting Examples
-------------------------

This page shows practical usage examples of **pyFuRNAce** scripting interfaces.

Two helices with FLAPs (and optimization)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../examples/basic_usage.py
   :language: python
   :linenos:


RNA filament
^^^^^^^^^^^^

.. literalinclude:: ../../examples/rna_filament.py
   :language: python
   :linenos:

.. _barrier-opti:

Folding barrier optimization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A simple folding barrier optimisation is provided in pyFuRNAce, and involves changing the 5′ start position of the structure (see :ref:`barrier-limit`  for details).

It can be accessed via the ``improve_folding_pathway`` method of the ``Origami`` class:

.. code-block:: python

   import pyfurnace as pf
   origami = pf.Origami()
   RENDER_TARGET = origami
   origami.append([]) # Add empty line
   origami = pf.simple_origami(dt_list=[120, 120, 120],
                               kl_columns=3,
                               main_stem=33,
                               use_angles=True) # Create a simple origami

   origami = origami.improve_folding_pathway(kl_delay=150)

Or it can be accessed via the ``Folding Barrier`` view of the graphical user interface.

.. image:: /_static/folding_barrier_view.png
   :alt: Folding Barrier View
   :align: center

For a full folding barrier optimization example, see:

.. _barrier-opti-script:

Full optimisation
+++++++++++++++++

.. literalinclude:: ../../examples/scaffold_barriers_opti_full.py
   :language: python
   :linenos:

Jupyter Notebooks
-----------------

.. toctree::
   :maxdepth: 2

   notebooks/droplets

.. toctree::
   :maxdepth: 2

   notebooks/rna_filament

.. toctree::
   :maxdepth: 2

   notebooks/ROAD_origami

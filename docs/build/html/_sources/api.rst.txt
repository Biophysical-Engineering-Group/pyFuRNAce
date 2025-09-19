..  _api:

Python API Reference
====================

What's the idea?
----------------
PyFuRNAce is inspired by the ROAD RNA origami software (Geary et al., Nat. Chem., 2021), which uses text-based RNA diagrams for RNA nanostructures design.

To facilitate the design of RNA nanostructures, pyFuRNAce aims to provide a simple and flexible Python package for creating, manipulating, and visualizing RNA with text-based diagrams.

The approach of pyFuRNAce is to hierarchically build RNA nanostructures from simple building blocks: strands and motifs.

Therefore, there are three main classes:

1. **Strand**
   A strand is simply a string of ASCII characters, which contains nucleotides (e.g., AUCG) and/or characters to route the strands on the page (e.g., ``-/\|``). If we imagine drawing a strand on a page, each symbol has an ``(x, y)`` position. In pyFuRNAce, each character position is calculated automatically, as the strand is defined by:

   * strand characters (e.g., ``A-UCG|CC—``)
   * start position (e.g., ``(0, 0)``)
   * start direction (e.g., ``pf.RIGHT``)
   * directionality of the sequence (e.g., ``'53'``; optional)

   A ``Strand`` can also have 3D coordinates associated with the sequence, allowing the assembly of 3D models in real time.

2. **Motif**
   The motif is the building block of the RNA origami, exploiting the vast library of RNA aptamers. In pyFuRNAce, a ``Motif`` is simply a collection of strands which can have base pairs. The base pairs define the secondary structure of the motifs. Internally, base pairs are represented as strand–position pairs in a dictionary. In the Python API, the secondary structure can be set and obtained in the more common `dot-bracket notation <https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/io/rna_structures.html>`_. In summary, a motif is made of:

   * strands (e.g., ``[Strand('CCAGG'), Strand('CCGGG', start=(4, 0), direction=pf.LEFT)]``)
   * structure (e.g., ``"((.((&)).))"``)

   A large library of pre-defined motifs and aptamers is provided in pyFuRNAce.

3. **Origami**
   The origami is the final RNA nanostructure. In pyFuRNAce, an ``Origami`` is simply a collection of motifs in a grid. The ``Origami`` concatenates motifs horizontally (in rows) and vertically, automatically joining the vertical junctions. In principle, any RNA nanostructure can be represented as an origami:

   * motifs (e.g., ``[[pf.TetraLoop(), pf.Stem(10), pf.Dovetail(-1)], [pf.KissingLoop180(), pf.Stem(8), pf.Dovetail(-2)]]``)

When your nanostructure blueprint is ready, pyFuRNAce provides functions to optimize the RNA sequence (in the ``generate`` subpackage), design PCR primers, or set up oxDNA simulations (``prepare`` subpackage).

Examples of how to use the pyFuRNAce API are available in the :ref:`tutorials` section.


Design Core Objects
-------------------

Strand
^^^^^^

.. automodule:: pyfurnace.design.core.strand
   :members:
   :undoc-members:
   :show-inheritance:

Motif
^^^^^
.. automodule:: pyfurnace.design.core.motif
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:

Origami
^^^^^^^
.. automodule:: pyfurnace.design.core.origami
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:

Sequence
^^^^^^^^
.. automodule:: pyfurnace.design.core.sequence
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:

Basepair
^^^^^^^^
.. automodule:: pyfurnace.design.core.basepair
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:

Common Variables
^^^^^^^^^^^^^^^^

Useful collections
""""""""""""""""""
.. autodata:: pyfurnace.design.core.symbols.nucleotides
.. autodata:: pyfurnace.design.core.symbols.iupac_code
.. autodata:: pyfurnace.design.core.symbols.base_pairing
.. autodata:: pyfurnace.design.core.symbols.db_pairs
.. autodata:: pyfurnace.design.core.symbols.all_pk_symbols
.. autodata:: pyfurnace.design.core.symbols.accept_symbol
.. autodata:: pyfurnace.design.core.symbols.bp_symbols
.. autodata:: pyfurnace.design.core.symbols.T7_PROMOTER

Translators for strings characters
""""""""""""""""""""""""""""""""""
.. autodata:: pyfurnace.design.core.symbols.pseudo_to_dot
   :no-value:
.. autodata:: pyfurnace.design.core.symbols.pair_db_sym
   :no-value:
.. autodata:: pyfurnace.design.core.symbols.nucl_to_none
   :no-value:
.. autodata:: pyfurnace.design.core.symbols.symb_to_none
   :no-value:
.. autodata:: pyfurnace.design.core.symbols.nucl_to_pair
   :no-value:
.. autodata:: pyfurnace.design.core.symbols.only_nucl
   :no-value:
.. autodata:: pyfurnace.design.core.symbols.horiz_flip
   :no-value:
.. autodata:: pyfurnace.design.core.symbols.verti_flip
   :no-value:
.. autodata:: pyfurnace.design.core.symbols.rotate_90
   :no-value:
.. autodata:: pyfurnace.design.core.symbols.symb_to_road
   :no-value:

Common functions
^^^^^^^^^^^^^^^^
.. automodule:: pyfurnace.design.core.symbols
   :members:
   :undoc-members:
   :exclude-members: nucleotides,iupac_code,base_pairing,db_pairs,all_pk_symbols,accept_symbol,bp_symbols,pseudo_to_dot,pair_db_sym,nucl_to_none,symb_to_none,nucl_to_pair,only_nucl,horiz_flip,verti_flip,rotate_90,symb_to_road,T7_PROMOTER,Node,BasePair

Motifs
------
Stem
^^^^
.. automodule:: pyfurnace.design.motifs.stem
   :members:
   :undoc-members:

Dovetail
^^^^^^^^
.. automodule:: pyfurnace.design.motifs.dovetail
   :members:
   :undoc-members:

Loops
^^^^^^^^^^^^^^^^
.. automodule:: pyfurnace.design.motifs.loops
   :members:
   :undoc-members:

Kissing loops
^^^^^^^^^^^^^^
.. automodule:: pyfurnace.design.motifs.kissing_loops
   :members:
   :undoc-members:

Aptamers
^^^^^^^^
.. automodule:: pyfurnace.design.motifs.aptamers
   :members:
   :undoc-members:

Utilities
---------
Motif Utilities
^^^^^^^^^^^^^^^
.. automodule:: pyfurnace.design.utils.motif_lib
   :members:
   :undoc-members:

Origami Utilities
^^^^^^^^^^^^^^^^^
.. automodule:: pyfurnace.design.utils.origami_lib
   :members:
   :undoc-members:

Sequence Generation
-------------------
.. automodule:: pyfurnace.generate.road
   :members:
   :undoc-members:

Prepare
-------

Primers
^^^^^^^^^^^^^^
.. automodule:: pyfurnace.prepare.utils
   :members:
   :undoc-members:

Oxdna simulation setup
^^^^^^^^^^^^^^^^^^^^^^
.. automodule:: pyfurnace.prepare.oxdna_sim
   :members:
   :undoc-members:

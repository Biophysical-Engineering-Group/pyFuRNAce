.. _gui_examples:

Graphical User Interfaces Examples
----------------------------------

Custom motif: Single Stranded
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Purely single-stranded motifs are discouraged in pyFuRNAce, since they are less predictable. In particular, they could lead to ambiguous scenarios in pyFuRNAce. For example, if you want two stems connected with a single-stranded ``3AAAAAAA5`` strand:

.. code-block:: bash

  5->NNNKNN       NKNNN->3
     ┊┊┊┊┊┊       ┊┊┊┊┊
  3<-NNNKNNAAAAAAANKNNN<-5

pyFuRNAce will try to use the first connection available, which would lead to:

.. code-block:: bash

  5->NNNKNNAAAAAAANKNNN->3
     ┊┊┊┊┊┊       ┊┊┊┊┊
  3<-NNNKNN       NKNNN<-5

Which will throw an error, since the directionality of the ``3AAAAAAA5`` strand is not compatible with the top strands of the stems.

1) Create an Origami with two stems
+++++++++++++++++++++++++++++++++++

.. image:: /_static/two_stems_motif.png
    :alt: Motif with two stems
    :align: center
    :width: 400px

You will have an origami made of two stems:

.. code-block:: bash

   NNNKNNNNNNNKNN
   ┊┊┊┊┊┊┊┊┊┊┊┊┊┊
   NNNKNNNNNNNKNN

You can select the second motif (motif index 1) and click of Custom in the motif menu'.

There are three different ways to create a custom motif:
- Structure converter
- Manual input
- Drawing tool

You can switch between the three during the creation of the custom motif, to use the tools that suit you best.

2. A) Custom structure converter
++++++++++++++++++++++++++++++++

In dot-bracket notation, an unpaired nucleotide is represented with a dot ``.``. So, a single-stranded motif of 7 nucleotides can be represented as ``.......``. In the structure converter, you can add either your structure, your sequence or both.
If a structure is added, the sequence will be composed by any nucleotide (N). If a sequence is added, the structure will be guessed with ViennaRNA. If both are added, the structure will be used to create the motif.

.. code-block:: bash
    Structure: .......
    Sequence: AAAAAAA

.. image:: /_static/custom_motif_structure_converter1.png
    :alt: Custom motif with structure converter
    :align: center
    :width: 400px

If we add a single stranded structure, the custom motif will automatically create a loop:

.. code-block:: bash

  5->AAAAAAA╮
            │
  3<-───────╯

But in this case we want a single-stranded connection. To break the loop, we can use the ViennaRNA dot-bracket cleavage symbol ``&``. Adding it to the end of the dot-bracket structure will break the loop:

.. code-block:: bash

    Structure: .......&
    Sequence: AAAAAAA

**Click convert to apply the changes.**
Now the custom motif will create a purely single-stranded motif:

.. code-block:: bash

  5->AAAAAAA->3

     ───────

Using the flip button at the top will produce the motif that we want:

.. code-block:: bash
     ───────

  3<-AAAAAAA<-5

.. image:: /_static/custom_motif_structure_converter2.png
    :alt: Custom motif with structure converter breaking the loop
    :align: center
    :width: 400px

2. B) Custom manual input
+++++++++++++++++++++++++

The second custom motif creation method is by manual text input. It involves writing the motif as text in the text area.
**Important Note**: to set the directionality of a strand, you need to add **only the 5** symbol at the beginning of the strand.

.. image:: /_static/custom_motif_manual_input1.png
    :alt: Custom motif with manual input
    :align: center
    :width: 400px

You can try to copy-paste the single-stranded motif:

.. code-block::

  5───────

  AAAAAAA5

**Click convert to apply the changes.**

If you wanna draw curves in the strand, you can use the slash symbols ``/`` and ``\``; while the minus symbol ``-`` and the pipe symbol ``|`` can be used to draw straight lines. If you wanna use the ASCII character of ROAD/pyFuRNAce, they can be copy-pasted in the popover at the top right (``Common symbols to copy``).

An example of a curved strand is:

.. code-block::

      ╭╮
  5───╯╰──

   AAAAAAA5

2. C) Custom drawing tool
+++++++++++++++++++++++++

This tools displays the canvas the draw the motif, where all the dots are the available positions.
When you first click on a dot, you are select the starting point of the strand. You can click on another dot in the same line/row to create a straight line. By consecutively clicking on dots in different lines/rows, you can draw you strand. Below the canvas, you can find the strands information such as: starting point, starting direction, characters and directionality. You can edit them to modify the strand in the canvas.

.. image:: /_static/custom_motif_drawing_tool.png
    :alt: Custom motif with drawing tool
    :align: center
    :width: 400px

You can add, select or remove strands with the buttons at the left of the canvas.

3. Complete the structure
+++++++++++++++++++++++++
Once you are satisfied with your custom motif, you can click the green ``Finish editing`` button at the bottom left. You can always go back to edit mode clicking the ``Edit the motif`` button above the motif preview.

Now your origami should look like:

.. image:: /_static/custom_finish_1.png
    :alt: Origami with custom single-stranded motif
    :align: center
    :width: 400px

.. code-block:: bash

   5->NNNNKNN───────NNNNKNN->3
      ┊┊┊┊┊┊┊       ┊┊┊┊┊┊┊
   3<-NNNNKNNAAAAAAANNNNKNN<-5

To add a separtion in the strand at the top, you can select ``Connections`` then ``start_end_stem`` in the motif menu.

Here is the final result:
.. image:: /_static/custom_finish_2.png
    :alt: Origami with custom single-stranded motif and separation
    :align: center
    :width: 400px

.. code-block:: bash

  5NNNNKNN─3 5────────NNNNKNN3
   ┊┊┊┊┊┊┊            ┊┊┊┊┊┊┊
  3NNNNKNN─────AAAAAAANNNNKNN5


3. Equivalent Code
++++++++++++++++++

The equivalent code to create the same origami with the scripting interface is:

.. code-block:: python

  import pyfurnace as pf

  aa_strand = pf.Motif.from_structure(".......&", "AAAAAAA&").flip()

  origami = pf.Origami([pf.Stem(7), pf.start_end_stem(), aa_strand, pf.Stem(7)])

In the motif menu, you can select `Code`, paste the code and click `Run` at the bottom right create the origami from python code.

.. image:: /_static/custom_code.png
    :alt: Origami with custom single-stranded motif and separation code
    :align: center
    :width: 400px

Make an aptamer
^^^^^^^^^^^^^^^

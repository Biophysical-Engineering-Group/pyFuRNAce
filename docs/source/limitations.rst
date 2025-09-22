.. _limitations:

Limitations & Workarounds
=========================

.. _barrier-limit:

PyFuRNAce folding barrier estimation
------------------------------------

To support cotranscriptional folding, the design of RNA origami should carefully avoid kinetic traps and steric hindrance. In the ROAD paper :cite:`geary2021rna`, a folding barrier penalty was introduced to account for a specific type of cotranscriptional steric hindrance. In particular, the RNA origami structures rely on duplexes with adjacent helices connected through 180° kissing loops (which form pseudoknots). For successful cotranscriptional folding, the transcription complex (DNA template, T7 RNA polymerase, and nascent RNA strand) must have sufficient spatial freedom for the RNA to assemble into duplexes. If the neighboring kissing loops are fully paired too early, then one of the strands which makes up the duplex is constrained at both ends, leading to steric hindrance which may prevent proper helix formation. This type of folding barrier has been further studied by Orponen et al. :cite:`orponen2025secondary`.

The folding barrier calculation in ROAD and pyFuRNAce highlights such cases: if an RNA duplex is expected to form only after adjacent kissing loops have already paired, the folding process may be hindered. To approximate this, the algorithm assumes that kissing loops do not pair instantaneously, introducing a default delay of 150 nucleotides during which no barrier is applied. Beyond this, two types of barriers are identified based on the length of the helix:

- **Weak barrier**: ≤ 5 bases (approximately half a helix turn), penalty +1.
- **Strong barrier**: > 5 bases, penalty +2.

To reduce folding barriers, three strategies are useful:

- **Use short dovetails rather than long dovetails,** since short dovetails require less than half a helix turn and thus lead to weaker barriers.
- **Change the 5′ start position of the structure,** which alters strand routing and can reduce the number of kissing loops between complementary sequences. This optimization is available in both the pyFuRNAce GUI and API  (see :ref:`barrier-opti`).
- **Modify the placement of continuous stems and kissing loops in the RNA origami structure,** which also changes routing. Since this requires altering the origami structure itself, it is not implemented as a direct function, but we provide a Python script in the documentation to assist with this (see :ref:`barrier-opti-script`).

Ultimately, cotranscriptional folding depends on many factors beyond this simplified model. While the folding barrier calculation is a useful heuristic for estimating cotranscriptional feasibility, our experimental validation did not reveal a clear correlation between predicted penalties and folding outcomes. At present, the final judgment still rests with the designer and further experimental characterization is needed to fully understand the kinetic landscape of RNA origami folding.

Revolvr
-------

For sequence generation and optimization, we use the Revolvr program from ROAD :cite:`geary2021rna`. Briefly, Revolvr implements a series of optimization steps, first assigning a random sequence which has the complementarity of the target dot-bracket string, followed by mutation steps interleaved with calls to ViennaRNA's RNAfold, until the predicted secondary structure matches the target structure. It then continues to mutate to reduce sequence symmetry (minimizing non-target complementarity and repeated subsequences) and to remove common restriction enzyme sites. Finally, kissing loop (KL) sequences are orthogonalized to prevent unintended multimerization and loop closure. Revolvr introduces a number of limitations which are inherited by pyFuRNAce:

- **Optimization failure of structures containing repeated fixed sequences.**
  If there are more than four copies of a fixed sequence in a structure, Revolvr optimization will often get stuck during sequence-symmetry minimization. This is because the fixed sequence plus one base will fail the test for repeated sequences.

  **Solution:** By extending the fixed sequences by unequal amounts of random sequence, repeats of the exact same sequence in the optimization region can be avoided.

- **No multimeric structures.**
  Revolvr can only perform optimization of one strand at a time.

  **Solution:** Optimization of multi-strand assemblies needs to be done in multiple steps. The user can also replace external KLs with continuous duplexes and optimize the combined structure, replacing the continuous duplex with separate KLs after optimization.

- **Timeout during optimization of large structures.**
  Revolvr can be very slow at optimizing large structures, and due to server constraints the pyFuRNAce webserver has a timeout of 2 hours.

  **Solution:** When running the optimisation through Python scripts, the Revolvr timeout can be set to an arbitrary value. For the 2.5 kb structure shown in Figure 5, we used a 24-hour timeout (see Supplementary Note S5). Designs with particular constraints (e.g. short stems, such as dovetails or long unpaired regions) are typically hard to design, and Revolvr might fail repeatedly. PyFuRNAce also provides a function to run multiple instances of Revolvr in parallel, called ``parallel_road``, which can speed up the optimisation, but requires multiple compute cores. Optimisation through Python scripts can also be run on a high-performance computing cluster, even though the parallelism is shared exclusively on one node.

- **Lack of fine control over KL energies.**
  As RNA structures grow in size and complexity, the number of KL needed to fold the complete structure naturally also increases. By specifying the interaction energy of each KL pair, RNA designers could gain some control over the folding pathway and assembly order of multimeric structures, enabling larger, more complex structures.

  **Solution:** Users can identify orthogonal KLs with the target energies ahead of time and manually insert the desired sequences. The KL sequence can be inserted when the motif is added or edited later in the ``Edit`` tab. We frequently use Table S4 from Geary *et al.* 2021 :cite:`geary2021rna`, which lists KL energies, for this purpose.

We are investigating novel sequence optimization strategies which overcome these limitations. PyFuRNAce is written with this eventuality in mind, and the sequence optimization algorithm is something which can be integrated later with a list of available algorithms.


ViennaRNA
---------

RNAfold is a dependency of Revolvr, used to predict the secondary structure of a given RNA sequence. It is part of the popular ViennaRNA suite :cite:`lorenz2011viennarna`. The folding algorithm is a dynamic programming algorithm :cite:`zuker1984rna`, which uses the Turner nearest neighbor model :cite:`xia1998thermodynamic` to predict the minimum free energy fold of the sequence. There are two major limitations introduced by RNAfold to the optimization of RNA structures:

- **Lack of pseudoknot prediction.**
  Unlike computing nested pairs which can be computed in polynomial time, pseudoknot enumeration is an NP-complete problem :cite:`lyngso2000pseudoknots`. Because of this, RNAfold, like most other secondary structure prediction tools, does not consider pseudoknots. Internal KLs, however, are pseudoknots, which is why Revolvr has to perform KL orthogonalization in a separate step. Unless :math:`P = NP`, this will remain a problem and we can only hope to develop better heuristic algorithms for KL sequence assignment for all but the smallest structures.

- **Lack of stabilizing tertiary interactions.**
  All RNA sequence optimization discussed here is performed at the secondary structure level—only Watson-Crick-Franklin base pairing is considered. This is the most common type of base-pairing in RNA structures; however, RNA can form a huge variety of additional stabilizing hydrogen bonds and π–π stacking interactions, meaning that this approximation is generally an underestimate of the free energy of any given RNA structure. There are two common motifs where this can pose problems in RNA origami design. First, the 180° KLs have an interaction in which two adenines stack across the major groove of the 6-base helix, stabilizing the structure. Second, many fluorescent aptamers, in particular, contain G-quadruplexes where four guanines interact via Hoogsteen-edge interactions. This means that predicted structures from RNAfold involving aptamers are generally inaccurate. We have not seen this cause major problems with designed structures, but it is important to keep in mind as a potential confounding factor. This is not something that we have an immediate solution for, but if cotranscriptional RNA origami becomes a widely adopted technology, quantifying the true energies of 180° KLs would be beneficial.


oxDNA Analysis Tools
--------------------

PyFuRNAce uses oxDNA Analysis Tools (OAT) to convert 3D models from the xoDNA representation (natively supported) to PDB.

- **Inaccurate sub-nucleotide features.**
  As an anisotropic one-bead-per-nucleotide model, oxDNA's file format only tracks the center of mass and orientation of each base, but lower level atomic information is lost. When users download PDB files from pyFuRNAce, the files are generated from the oxDNA representation via OAT's converter, which assembles ideal all-atom representations of nucleotides to the reference frame of the coarse-grained representation. These all-atom representations come from an NMR structure of an RNA helix (PDB ID: 2jxq) :cite:`popenda2008bulged`. This is not an optimized structure and should be considered a low-resolution model.

  **Solution:** Before additional model building or all-atom molecular dynamics simulations, we usually refine the all-atom PDB structures with QRNAS :cite:`stasiewicz2019qrnas` using default settings.

References
----------

.. bibliography::
   :style: unsrt
   :filter: docname in docnames

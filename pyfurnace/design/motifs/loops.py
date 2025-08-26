from . import CONFS_PATH
from ..core.coordinates_3d import Coords
from ..core.strand import Strand
from ..core.motif import Motif


class Loop(Motif):
    """
    Represents a generic unpaired loop motif in an RNA structure.
    Optionally it can be initialized with a sequence or be oriented
    to the left.

    Parameters
    ----------
    open_left : bool, optional
        Whether to flip the loop horizontally and vertically to open to the left.
        Default is False.
    sequence : str, optional
        Nucleotide sequence to include in the loop. If provided, creates a strand
        with the sequence and loop symbols. Default is "".
    **kwargs : dict, optional
        Additional keyword arguments passed to the parent `Motif` class.

    Attributes
    ----------
    strands : list of Strand
        The list of strands that make up the loop motif.
    flipped : bool
        Whether the motif has been flipped horizontally/vertically
        (if `open_left=True`).
    """

    def __init__(self, open_left: bool = False, sequence: str = "", **kwargs) -> None:
        """
        Initialize a loop motif, optionally with a given sequence and orientation.

        If a sequence is provided, a single strand is created for it, flanked by
        loop symbols.

        Parameters
        ----------
        open_left : bool, optional
            Whether to flip the loop horizontally and vertically to open to the
            left (default is False).
        sequence : str, optional
            Nucleotide sequence to insert into the loop. If provided, a strand will
            be created (default is "").
        **kwargs : dict, optional
            Additional keyword arguments passed to the `Motif` superclass.

        Returns
        -------
        None
        """
        # create motif of Cap without basepairs by turing autobasepairing of
        if sequence:
            seq_len = len(sequence)
            ### create the strand
            strand = Strand(
                "─" * seq_len + "╰│╭" + sequence, start=(seq_len, 2), direction=(-1, 0)
            )
            # Add the strand to the list of strands
            kwargs["strands"] = kwargs.get("strands", []) + [strand]

        kwargs["join"] = False
        super().__init__(**kwargs)
        if open_left:
            self.flip(horizontally=True, vertically=True)


class TetraLoop(Loop):
    """
    Represents a specific 4-nucleotide tetraloop RNA motif.
    This class implements a canonical tetraloop structure (e.g., UUCG) with a
    predefined folding pattern and 3D coordinates.

    Parameters
    ----------
    open_left : bool, optional
        Whether to flip the loop to open to the left. Default is False.
    sequence : str, optional
        4-nucleotide RNA sequence for the tetraloop. Default is "UUCG".
    **kwargs : dict, optional
        Additional keyword arguments passed to the `Loop` superclass. You can
        override the default strand using the `strands` argument.

    Attributes
    ----------
    strands : list of Strand
        The single-stranded representation of the tetraloop, including symbolic
        structure markers.
    sequence : str
        The 4-nucleotide sequence used to define the tetraloop.
    _coords : Coords
        3D coordinate data loaded from template structure (e.g., from PDB 2KOC).
    """

    def __init__(
        self, open_left: bool = False, sequence: str = "UUCG", **kwargs
    ) -> None:
        """
        Initialize a tetraloop motif with a specific 4-nucleotide sequence.

        The tetraloop is represented as a single strand folded with a predefined shape.
        Optionally, the loop can be opened to the left.

        Parameters
        ----------
        open_left : bool, optional
            Whether to flip the loop to open to the left (default is False).
        sequence : str, optional
            4-nucleotide sequence to assign to the tetraloop (default is "UUCG").
        **kwargs : dict, optional
            Additional keyword arguments passed to the `Loop` superclass.
            You can override default strand(s) using `strands`.

        Raises
        ------
        ValueError
            If the provided sequence is not exactly 4 nucleotides long.

        Returns
        -------
        None
        """
        """
        Attributes of the class Cap_UUCG, which is a daugther class of the class Motif.
        -------------------------------------------------------------------------------
        UUCG_bool: bool (default= False)
            indicates if a UUCG sequence should be added into the cap
        """
        # create strands deascribing tetraloop
        if len(sequence) != 4:
            raise ValueError(
                "The sequence length doesn't match the length required "
                "for a tetraloop, which is 4."
            )

        # Create new strands if the strand is not provided
        if "strands" in kwargs:
            strands = kwargs.pop("strands")
        else:
            strand = Strand(
                sequence[:2] + "╰│╭" + sequence[2:4], start=(2, 2), direction=(-1, 0)
            )
            ### PDB: 2KOC
            strand._coords = Coords.load_from_file(
                CONFS_PATH / "TetraLoop.dat", dummy_ends=(True, True)
            )
            strands = [strand]

        # create motif of Cap without basepairs by turing autobasepairing of
        kwargs.setdefault("autopairing", False)

        super().__init__(strands=strands, open_left=open_left, **kwargs)

    def set_sequence(self, new_sequence: str) -> None:
        """
        Set a new 4-nucleotide sequence for the tetraloop.

        Parameters
        ----------
        new_sequence : str
            New tetraloop sequence. Must be exactly 4 nucleotides.

        Raises
        ------
        ValueError
            If the provided sequence is not exactly 4 nucleotides long.

        Returns
        -------
        None
        """
        if len(new_sequence) != 4:
            raise ValueError(
                "The sequence length doesn't match the length required "
                "for a tetraloop, which is 4."
            )
        self[0].sequence = new_sequence

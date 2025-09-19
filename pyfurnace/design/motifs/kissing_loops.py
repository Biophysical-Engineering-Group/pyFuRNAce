import os
from typing import Optional, Union, List

# RNA used if a specific sequence is provided to calculate the energy

# pyfurnace imports
from . import CONFS_PATH
from ..core.coordinates_3d import Coords
from ..core.sequence import Sequence
from ..core.strand import Strand
from .loops import Loop

### File Location for the kissing loop energy dictionaries
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


class KissingLoop(Loop):
    """
    Represents a pseudoknotted kissing loop motif in RNA secondary structure.

    Parameters
    ----------
    open_left : bool, optional
        If True, the loop is open on the left side. Default is False.
    sequence : str, optional
        RNA sequence of the internal kissing loop. Default is an empty string.
    seq_len : int, optional
        Expected length of the sequence. Used for validation. Default is 0.
    pk_index : str or int, optional
        Pseudoknot index or symbol identifying the loop. Default is '0'.
    energy : float, optional
        Free energy (in kcal/mol) associated with the loop. Default is -8.5.
    energy_tolerance : float, optional
        Tolerance on the energy for acceptable structural variants. Default is 1.0.
    **kwargs : dict
        Additional arguments passed to the parent `Loop` constructor.

    Attributes
    ----------
    pk_index : str
        Pseudoknot index symbol.
    energy : float
        Free energy of the motif.
    energy_tolerance : float
        Energy tolerance for acceptable variants.
    """

    _KL_coords = Coords()

    def __init__(
        self,
        open_left: bool = False,
        sequence: str = "",
        seq_len: int = 0,
        pk_index: Union[str, int] = "0",
        energy: float = -8.5,
        energy_tolerance: float = 1.0,
        **kwargs,
    ) -> None:
        """
        Initialize a KissingLoop motif representing an internal pseudoknotted
        kissing interaction.

        Parameters
        ----------
        open_left : bool, optional
            If True, flip the loop orientation (default is False).
        sequence : str, optional
            Sequence for the internal strand of the kissing loop (default is "").
        seq_len : int, optional
            Expected sequence length for validation (default is 0).
        pk_index : str or int, optional
            Pseudoknot identifier used to tag the kissing interaction (default "0").
        energy : float, optional
            Free energy associated with the interaction (default is -8.5 kcal/mol).
        energy_tolerance : float, optional
            Energy tolerance for structural variants (default is 1.0 kcal/mol).
        **kwargs : dict
            Additional parameters passed to `Loop`.

        Raises
        ------
        ValueError
            If sequence length mismatches or pk_index has an invalid type.
        """

        pk_index = self._check_pk_index(pk_index)
        self._pk_index = pk_index
        self._seq_len = seq_len
        self._energy = energy
        self._energy_tolerance = energy_tolerance
        if "strands" in kwargs:
            strands = kwargs.pop("strands")
        else:
            strands = self._create_strands(
                sequence=sequence, return_strand=True, pk_index=pk_index
            )

        # create motif with the strands making up the external kissing loop structure
        super().__init__(strands=strands, open_left=open_left, **kwargs)

    ###
    ### Properties
    ###

    @property
    def pk_index(self):
        """Returns the pseudoknot symbol of the kissing loop"""
        return self._pk_index

    @pk_index.setter
    def pk_index(self, new_index):
        """Set the pseudoknot symbol of the kissing loop"""
        new_index = self._check_pk_index(new_index)
        self._create_strands(sequence=self.get_kissing_sequence(), pk_index=new_index)

    @property
    def energy_tolerance(self):
        """Returns the energy tolerance of the internal kissing loop"""
        return self._energy_tolerance

    @energy_tolerance.setter
    def energy_tolerance(self, new_energy_tolerance):
        """Set the energy tolerance of the internal kissing loop"""
        if (
            not isinstance(new_energy_tolerance, (float, int))
            or new_energy_tolerance < 0
        ):
            raise ValueError("The energy tolerance should be a positive number.")
        self._energy_tolerance = new_energy_tolerance
        for strand in self:
            if hasattr(strand, "pk_info"):
                strand.pk_info["dE"] = [new_energy_tolerance]
        self._trigger_callbacks()

    @property
    def energy(self):
        """Returns the energy of the internal kissing loop"""
        return self._energy

    @energy.setter
    def energy(self, new_energy):
        """Set the energy of the internal kissing loop"""
        new_energy = round(float(new_energy), 2)
        self._energy = new_energy
        for strand in self:
            if hasattr(strand, "pk_info"):
                strand.pk_info["E"] = [new_energy]
        self._trigger_callbacks()

    ###
    ### METHODS
    ###

    def complementary_pk_index(self, pk_index: Optional[Union[str, int]] = None) -> str:
        """Returns the complementary pseudoknot symbol of the kissing loop

        Parameters
        ----------
        pk_index : str or int, optional
            The pseudoknot index to complement. If None, uses the current pk_index.

        Returns
        -------
        str
            The complementary pseudoknot index.

        """
        if pk_index is None:
            pk_index = self._pk_index
        elif not isinstance(pk_index, str):
            pk_index = self._check_pk_index(pk_index)

        if "'" in pk_index:
            return pk_index[:-1]
        else:
            return pk_index + "'"

    def get_kissing_sequence(self):
        """Returns the kissing sequence of the kissing loop"""
        return self[0].sequence

    def set_sequence(self, new_seq):
        """Set the sequence of the strand"""
        self._create_strands(
            sequence=new_seq, pk_index=self._pk_index, only_sequence=True
        )

    def _check_pk_index(self, pk_index):
        if pk_index is None:
            pk_index = "0"
        elif not isinstance(pk_index, (int, str)):
            raise ValueError("The pk_index should be an integer or a string.")
        elif isinstance(pk_index, int):
            pk_index = str(abs(pk_index)) + "'" * (pk_index < 0)
        return pk_index

    def _create_strands(
        self,
        sequence: str = "",
        return_strand: bool = False,
        pk_index: Optional[Union[str, int]] = None,
        fold_180kl: bool = False,
        only_sequence: bool = False,
    ) -> Union[None, List[Strand]]:
        """
        Protected class that takes a sequence and a pk_index and creates a strand
        with the kissing loop structure and metadata.

        Parameters
        ----------
        sequence : str, optional
            The sequence for the kissing loop strand (default is "").
        return_strand : bool, default is False
            If True, returns the created strand instead of replacing it in the motif.
        pk_index : str or int, optional
            The pseudoknot index for the kissing loop (default is None).
        fold_180kl : bool, default is False
            If True, folds the kissing loop as a 180-degree loop.
        only_sequence : bool, default is False
            If True, only updates the sequence of the existing strand without
            recreating it (default is False).

        Returns
        -------
        List[Strand] or None
            Returns a list containing the created strand if `return_strand` is True.
            Otherwise, it replaces the existing strands in the motif and returns None.
        """
        self._pk_index = self._check_pk_index(pk_index)

        seq_len = self._seq_len

        if sequence:

            if not isinstance(sequence, Sequence):
                sequence = Sequence(sequence, directionality="53")

            if seq_len and len(sequence) != seq_len:
                raise ValueError(
                    f"The sequence length doesn't match the length "
                    f"for this kissing loop, which is {seq_len}."
                )
            else:
                seq_len = len(sequence)
            if all([s in "ACGU" for s in sequence]):
                from RNA import fold

                if fold_180kl:
                    seq_to_fold = f"A{sequence}A&A{sequence.reverse_complement()}A"
                else:
                    seq_to_fold = f"{sequence}&{sequence.reverse_complement()}"
                self._energy = round(fold(seq_to_fold)[1], 2)
                self._energy_tolerance = 0
        else:
            sequence = Sequence("N" * seq_len, directionality="53")

        # ### if the strands are already created, just update the sequence
        # if hasattr(self, '_strands') and not fold_180kl:
        #     self._strands[0].sequence = sequence
        #     return self._strands

        # if the strands are already created, just update the sequence
        if only_sequence:
            if fold_180kl:
                sequence = "AA" + sequence + "A"
            self._strands[0].sequence = sequence
            self._strands[0].pk_info["id"] = [self._pk_index]
            self._strands[0].pk_info["E"] = [self._energy]
            self._strands[0].pk_info["dE"] = [self._energy_tolerance]
            return self._strands

        ### create the strand
        else:
            strand = Strand(
                f"┼─{sequence}╭╰{'─' * seq_len}─╯┼│╭",
                start=(seq_len + 2, 2),
                direction=(-1, 0),
                directionality=sequence.directionality,
            )
            pk_info = {
                "id": [self._pk_index],
                "ind_fwd": [(0, seq_len - 1)],
                "E": [self._energy],
                "dE": [self._energy_tolerance],
            }
            setattr(strand, "pk_info", pk_info)

        # Add oxDNA coordinates
        strand._coords = self._KL_coords.copy()

        if return_strand:
            return [strand]

        # replace the strands
        self.replace_all_strands([strand], copy=False, join=False)


class KissingLoop120(KissingLoop):
    """
    Kissing loop motif based on structural template from PDB entry 1BJ2.

    Parameters
    ----------
    open_left : bool, optional
        Whether the loop is open on the left side. Default is False.
    sequence : str, optional
        RNA sequence to use. Default is an empty string.
    pk_index : str or int, optional
        Pseudoknot index or label. Default is '0'.
    energy : float, optional
        Free energy (kcal/mol) of the loop. Default is -8.5.
    energy_tolerance : float, optional
        Tolerance on the energy value. Default is 1.0.
    **kwargs : dict
        Additional arguments passed to the parent constructor.
    """

    _KL_coords = Coords.load_from_file(
        CONFS_PATH / "KissingLoop120.dat", dummy_ends=(True, True)
    )

    def __init__(
        self,
        open_left: bool = False,
        sequence: str = "",
        pk_index: str | int = "0",
        energy: float = -8.5,
        energy_tolerance: float = 1.0,
        **kwargs,
    ):
        kwargs["seq_len"] = 7
        super().__init__(
            open_left=open_left,
            sequence=sequence,
            pk_index=pk_index,
            energy=energy,
            energy_tolerance=energy_tolerance,
            **kwargs,
        )


# https://doi.org/10.2210/pdb2D1B/pdb
class KissingLoop180(KissingLoop):
    """
    Kissing loop from HIV with idealized helical geometry for 180-degree interaction.

    Parameters
    ----------
    open_left : bool, optional
        Whether the loop is open on the left. Default is False.
    sequence : str, optional
        RNA sequence for the motif. Default is an empty string.
    pk_index : str or int, optional
        Pseudoknot index for identification. Default is '0'.
    energy : float, optional
        Free energy of the motif. Default is -8.5.
    energy_tolerance : float, optional
        Acceptable energy deviation. Default is 1.0.
    **kwargs : dict
        Additional keyword arguments.
    """

    _KL_coords = Coords.load_from_file(CONFS_PATH / "KissingLoop180.dat")

    def __init__(
        self,
        open_left: bool = False,
        sequence: str = "",
        pk_index: str | int = "0",
        energy: float = -8.5,
        energy_tolerance: float = 1.0,
        **kwargs,
    ):
        kwargs["seq_len"] = 6
        super().__init__(
            open_left=open_left,
            sequence=sequence,
            pk_index=pk_index,
            energy=energy,
            energy_tolerance=energy_tolerance,
            **kwargs,
        )

    def get_kissing_sequence(self):
        """Returns the kissing sequence of the kissing loop"""
        return super().get_kissing_sequence()[2:-1]

    def _create_strands(
        self,
        sequence: str = "",
        return_strand: bool = False,
        pk_index: int = 0,
        only_sequence: bool = False,
    ):
        """
        Protected class that takes a sequence and a pk_index and creates a strand
        with the kissing loop structure and metadata.

        Parameters
        ----------
        sequence : str, optional
            The sequence for the kissing loop strand (default is "").
        return_strand : bool, default is False
            If True, returns the created strand instead of replacing it in the motif.
        pk_index : str or int, optional
            The pseudoknot index for the kissing loop (default is None).
        only_sequence : bool, default is False
            If True, only updates the sequence of the existing strand without
            recreating it (default is False).

        Returns
        -------
        List[Strand] or None
            Returns a list containing the created strand if `return_strand` is True.
            Otherwise, it replaces the existing strands in the motif and returns None.
        """
        strand = super()._create_strands(
            sequence,
            return_strand=True,
            pk_index=pk_index,
            fold_180kl=True,
            only_sequence=only_sequence,
        )[0]

        if not only_sequence:
            # modify the strand start position and pk_info
            strand.start = (10, 2)
            strand.strand = "AA" + strand.strand + "─A"
            strand.pk_info["ind_fwd"] = [(2, 7)]

        # if we don't want to replace the strands, just return the strand
        if return_strand:
            return [strand]

        # replace the strands
        self.replace_all_strands([strand], copy=False, join=False)


class BranchedKissingLoop(KissingLoop):
    """
    Represents a 180° kissing loop motif with a branch connection.
    The motif consists of two strands:
    - 1) A kissing loop strand with a branch
    - 2) A connection strand for the branched kissing loop

    Parameters
    ----------
    open_left : bool, optional
        Whether the motif is open on the left side. Default is False.
    sequence : str, optional
        RNA sequence used for the kissing strand. Default is an empty string.
    pk_index : str or int, optional
        Pseudoknot identifier. Default is '0'.
    energy : float, optional
        Free energy of the structure. Default is -8.5.
    energy_tolerance : float, optional
        Tolerance for energy variations. Default is 1.0.
    **kwargs : dict
        Additional keyword arguments.
    """

    _KL_coords = Coords.load_from_file(
        CONFS_PATH / "BranchedKissingLoop_1.dat", dummy_ends=(True, False)
    )
    _KL_coords2 = Coords.load_from_file(
        CONFS_PATH / "BranchedKissingLoop_2.dat", dummy_ends=(True, True)
    )

    def __init__(
        self,
        open_left=False,
        sequence: str = "",
        pk_index: str | int = "0",
        energy: float = -8.5,
        energy_tolerance: float = 1.0,
        **kwargs,
    ):
        kwargs["seq_len"] = 6
        super().__init__(
            open_left=open_left,
            sequence=sequence,
            pk_index=pk_index,
            energy=energy,
            energy_tolerance=energy_tolerance,
            **kwargs,
        )

    def get_kissing_sequence(self):
        """Returns the kissing sequence of the kissing loop"""
        strand_ind = [i for i, s in enumerate(self) if len(s.sequence) == 7][0]
        return self[strand_ind].sequence[:-1]

    def _create_strands(
        self,
        sequence: str = "",
        return_strand: bool = False,
        pk_index: int = 0,
        only_sequence: bool = False,
    ):
        """
        Protected class that takes a sequence and a pk_index and creates a strand
        with the kissing loop structure and metadata.

        Parameters
        ----------
        sequence : str, optional
            The sequence for the kissing loop strand (default is "").
        return_strand : bool, default is False
            If True, returns the created strand instead of replacing it in the motif.
        pk_index : str or int, optional
            The pseudoknot index for the kissing loop (default is None).
        only_sequence : bool, default is False
            If True, only updates the sequence of the existing strand without
            recreating it (default is False).

        Returns
        -------
        List[Strand] or None
            Returns a list containing the created strand if `return_strand` is True.
            Otherwise, it replaces the existing strands in the motif and returns None.
        """
        strand = super()._create_strands(
            sequence,
            return_strand=True,
            pk_index=pk_index,
            only_sequence=only_sequence,
        )[0]

        if not only_sequence:
            # create the strand
            strand.start = (9, 3)
            strand.direction = (0, -1)
            strand.strand = "│╮" + strand.strand + "─A"
            strand.pk_info["ind_fwd"] = [(0, 5)]

            connect_strand = Strand("╭│", start=(10, 2), direction=(-1, 0))
            connect_strand._coords = self._KL_coords2.copy()
            strands = [strand, connect_strand]
        else:
            self[0].strand = self[0].strand + "A"
            strands = self._strands

        if return_strand:
            return strands

        # replace the strands
        self.replace_all_strands(strands, copy=False, join=False)


class KissingDimer(KissingLoop180):
    """
    A kissing loop dimer composed of two complementary 180-degree loops.

    Parameters
    ----------
    sequence : str, optional
        Sequence of the top strand; bottom is inferred as reverse complement.
    pk_index : str or int, optional
        Pseudoknot identifier for the dimer. Default is '0'.
    energy : float, optional
        Free energy of the dimer. Default is -8.5.
    energy_tolerance : float, optional
        Tolerance for energy deviation. Default is 1.0.
    **kwargs : dict
        Additional keyword arguments.
    """

    _KL_coords2 = Coords.load_from_file(CONFS_PATH / "KissingLoop180_2.dat")

    def __init__(
        self,
        sequence: str = "",
        pk_index: str | int = "0",
        energy: float = -8.5,
        energy_tolerance: float = 1.0,
        **kwargs,
    ):
        super().__init__(
            sequence=sequence,
            pk_index=pk_index,
            energy=energy,
            energy_tolerance=energy_tolerance,
            **kwargs,
        )

    ###
    ### METHODS
    ###

    def _create_strands(
        self,
        sequence: str = "",
        return_strand: bool = False,
        pk_index: int = 0,
        only_sequence: bool = False,
    ):
        """
        Protected class that takes a sequence and a pk_index and creates a strand
        with the kissing loop structure and metadata.

        Parameters
        ----------
        sequence : str, optional
            The sequence for the kissing loop strand (default is "").
        return_strand : bool, default is False
            If True, returns the created strand instead of replacing it in the motif.
        pk_index : str or int, optional
            The pseudoknot index for the kissing loop (default is None).
        only_sequence : bool, default is False
            Dummy parameter for compatibility; not used here.

        Returns
        -------
        List[Strand] or None
            Returns a list containing the created strand if `return_strand` is True.
            Otherwise, it replaces the existing strands in the motif and returns None.
        """
        bottom_pk_index = self.complementary_pk_index(pk_index)
        bottom_strand = super()._create_strands(
            sequence,
            return_strand=True,
            pk_index=bottom_pk_index,
        )[0]
        # add the pk_index to override the pk_index of the bottom strand
        self._pk_index = self._check_pk_index(pk_index)

        seq = bottom_strand.sequence[2:-1]
        rev_comp = seq.reverse_complement()

        ### shift the second strand to make space for the second one
        bottom_strand.start = (13, 3)
        # the second strand is the reverse complement of the first
        bottom_strand.sequence = "AA" + rev_comp + "A"

        ### create the second
        top_strand = KissingLoop180(
            open_left=True,
            sequence=seq,
            pk_index=self._pk_index,
            energy=self._energy,
            energy_tolerance=self._energy_tolerance,
        )[0]

        ## COORDINATES FROM OXVIEW HELIX
        top_strand._coords = self._KL_coords2.copy()

        strands = [top_strand, bottom_strand]

        if return_strand:
            return strands

        # replace the strands
        self.replace_all_strands(strands, copy=False, join=False)

    def set_up_sequence(self, new_seq: Union[str, Sequence]):
        """Set the sequence of the top strand"""
        self.set_sequence(new_seq)

    def set_down_sequence(self, new_seq: Union[str, Sequence]):
        """Set the sequence of the bottom strand"""
        new_seq = Sequence(new_seq, self[1].directionality)
        self.set_sequence(new_seq.reverse_complement())


class KissingDimer120(KissingLoop120):
    """
    Dimer of two complementary kissing loops based on the 120-degree model.

    Parameters
    ----------
    sequence : str, optional
        RNA sequence for the top strand. Default is "".
    pk_index : str or int, optional
        Pseudoknot symbol. Default is '0'.
    energy : float, optional
        Energy of the dimer motif. Default is -8.5.
    energy_tolerance : float, optional
        Allowed energy deviation. Default is 1.0.
    **kwargs : dict
        Additional keyword arguments.
    """

    _KL_coords2 = Coords.load_from_file(
        CONFS_PATH / "KissingLoop120_2.dat", dummy_ends=(True, True)
    )

    def __init__(
        self,
        sequence: str = "",
        pk_index: str | int = "0",
        energy: float = -8.5,
        energy_tolerance: float = 1.0,
        **kwargs,
    ):
        super().__init__(
            sequence=sequence,
            pk_index=pk_index,
            energy=energy,
            energy_tolerance=energy_tolerance,
            **kwargs,
        )

    ###
    ### METHODS
    ###

    def _create_strands(
        self, sequence="", return_strand=False, pk_index=0, only_sequence: bool = False
    ):
        """
        Protected class that takes a sequence and a pk_index and creates a strand
        with the kissing loop structure and metadata.

        Parameters
        ----------
        sequence : str, optional
            The sequence for the kissing loop strand (default is "").
        return_strand : bool, default is False
            If True, returns the created strand instead of replacing it in the motif.
        pk_index : str or int, optional
            The pseudoknot index for the kissing loop (default is None).
        only_sequence : bool, default is False
            Dummy parameter for compatibility; not used here.

        Returns
        -------
        List[Strand] or None
            Returns a list containing the created strand if `return_strand` is True.
            Otherwise, it replaces the existing strands in the motif and returns None.
        """
        # the bottom pk_index is the inverse of the top one
        bottom_pk_index = self.complementary_pk_index(pk_index)
        bottom_strand = super()._create_strands(
            sequence, return_strand=True, pk_index=bottom_pk_index
        )[0]
        # add the pk_index to override the pk_index of the bottom strand
        self._pk_index = self._check_pk_index(pk_index)
        seq = bottom_strand.sequence

        ### shift the second strand to make space for the second one
        bottom_strand.start = (10, 3)
        # the second strand is the reverse complement of the first
        bottom_strand.sequence = seq.reverse_complement()

        ### create the second
        top_strand = KissingLoop120(
            open_left=True,
            sequence=seq,
            pk_index=self._pk_index,
            energy=self._energy,
            energy_tolerance=self._energy_tolerance,
        )[0]
        top_strand._coords = self._KL_coords2.copy()

        strands = [top_strand, bottom_strand]

        if return_strand:
            return strands

        self.replace_all_strands(strands, copy=False, join=False)

    def set_up_sequence(self, new_seq: Union[str, Sequence]):
        """Set the sequence of the top strand"""
        self.set_sequence(new_seq)

    def set_down_sequence(self, new_seq: Union[str, Sequence]):
        """Set the sequence of the bottom strand"""
        new_seq = Sequence(new_seq, self[1].directionality)
        self.set_sequence(new_seq.reverse_complement())


class BranchedDimer(BranchedKissingLoop):
    """
    Branched kissing loop dimer with top and bottom interactions.
    The motif consists of three strands:
    - 1) A 180°KL dimer strand
    - 2) A branched kissing loop strand
    - 3) A connection strand for the branched kissing loop

    Parameters
    ----------
    sequence : str, optional
        Sequence for the top kissing strand. Default is "".
    pk_index : str or int, optional
        Identifier for the pseudoknot. Default is '0'.
    energy : float, optional
        Free energy of the motif. Default is -8.5.
    energy_tolerance : float, optional
        Energy tolerance threshold. Default is 1.0.
    **kwargs : dict
        Additional keyword arguments.
    """

    _KL_coords3 = Coords.load_from_file(CONFS_PATH / "BranchedKissingLoop_3.dat")

    def __init__(
        self,
        sequence: str = "",
        pk_index: str | int = "0",
        energy: float = -8.5,
        energy_tolerance: float = 1.0,
        **kwargs,
    ):
        super().__init__(
            sequence=sequence,
            pk_index=pk_index,
            energy=energy,
            energy_tolerance=energy_tolerance,
            **kwargs,
        )

    ###
    ### METHODS
    ###

    def _create_strands(
        self, sequence="", return_strand=False, pk_index=0, only_sequence: bool = False
    ):
        """
        Protected class that takes a sequence and a pk_index and creates a strand
        with the kissing loop structure and metadata.

        Parameters
        ----------
        sequence : str, optional
            The sequence for the kissing loop strand (default is "").
        return_strand : bool, default is False
            If True, returns the created strand instead of replacing it in the motif.
        pk_index : str or int, optional
            The pseudoknot index for the kissing loop (default is None).
        only_sequence : bool, default is False
            If True, only updates the sequence of the existing strand
            without creating a new one (default is False).

        Returns
        -------
        List[Strand] or None
            Returns a list containing the created strand if `return_strand` is True.
            Otherwise, it replaces the existing strands in the motif and returns None.
        """
        # the bottom pk_index is the inverse of the top one
        bottom_pk_index = self.complementary_pk_index(pk_index)

        if not sequence:
            sequence = Sequence("N" * self._seq_len, directionality="53")
        elif not isinstance(sequence, Sequence):
            sequence = Sequence(sequence, directionality="53")

        rev_comp = sequence.reverse_complement()

        if only_sequence:
            # metti al primo posto il branched KL
            self._strands[0:2] = self._strands[0:2][::-1]

        strands = super()._create_strands(
            rev_comp,
            return_strand=True,
            pk_index=bottom_pk_index,
            only_sequence=only_sequence,
        )
        # add the pk_index to override the pk_index of the bottom strand
        self._pk_index = self._check_pk_index(pk_index)

        if not only_sequence:
            ### shift the bottom branched KL
            strands[0].start = (strands[0].start[0] + 3, strands[0].start[1] + 1)
            strands[1].start = (strands[1].start[0] + 3, strands[1].start[1] + 1)

            ### create the top strand
            top_strand = KissingLoop180(
                open_left=True,
                sequence=sequence,
                pk_index=self._pk_index,
                energy=self._energy,
                energy_tolerance=self._energy_tolerance,
            )[0]
            top_strand._coords = self._KL_coords3.copy()

            strands.insert(0, top_strand)
        else:
            # put back the branched KL at the first position
            self._strands[0:2] = self._strands[0:2][::-1]
            # adjust the energy and energy tolerance of the 180°KL strand
            self[0].sequence = "AA" + sequence + "A"
            self[0].pk_info["E"] = [self._energy]
            self[0].pk_info["dE"] = [self._energy_tolerance]
            strands = self._strands

        if return_strand:
            return strands

        # replace the strands
        self.replace_all_strands(strands, copy=False, join=False)

    def set_up_sequence(self, new_seq: Union[str, Sequence]):
        """Set the sequence of the top strand"""
        self.set_sequence(new_seq)

    def set_down_sequence(self, new_seq: Union[str, Sequence]):
        """Set the sequence of the bottom strand"""
        new_seq = Sequence(new_seq, self[1].directionality)
        self.set_sequence(new_seq.reverse_complement())

import numpy as np
from scipy.spatial.transform import Rotation as R
from typing import Any, Dict, Tuple, List, Union, Literal, TYPE_CHECKING

### OAT IMPORTS
try:
    from oxDNA_analysis_tools.PDB_oxDNA import PDB_oxDNA
    from oxDNA_analysis_tools.UTILS.RyeReader import conf_to_str, get_top_string

    oat_installed = True
except ImportError:
    oat_installed = False

if TYPE_CHECKING:  # only for type checkers / linters, not at runtime
    from .strand import Strand


class ProteinCoords:
    """
    Represents a protein's sequence and its corresponding 3D coordinates.

    Attributes
    ----------
    sequence : str
        The amino acid sequence of the protein.
    coords : np.ndarray
        A NumPy array containing the 3D coordinates associated with the sequence.
    """

    def __init__(self, sequence: str = "", coords: np.ndarray = np.array(())) -> None:
        """
        Initialize the ProteinCoords instance with a sequence and coordinates.

        Parameters
        ----------
        sequence : str, optional
            The amino acid sequence of the protein.
        coords : np.ndarray, optional
            A NumPy array of shape (N, 3) representing the 3D coordinates.
        """
        # Initialize the sequence as an empty string
        self._sequence: str = ""
        # Initialize the coordinates as an empty NumPy array
        self._coords: np.ndarray = np.array(())

        # Set the sequence and coordinates with their setters
        self.sequence = sequence
        self.coords = coords

    def __str__(self) -> str:
        """Return a string representation of the ProteinCoords instance."""
        return f"ProteinCoords({self._sequence})"

    def __repr__(self) -> str:
        """Return a string representation of the ProteinCoords instance."""
        return f"ProteinCoords({self._sequence})"

    def __getitem__(self, key: int) -> np.ndarray:
        """Retrieve coordinates by index."""
        return self.coords[key]

    def __setitem__(self, key: int, value: np.ndarray) -> None:
        """Set coordinates at a specific index."""
        self.coords[key] = value

    def __len__(self) -> int:
        """Return the length of the sequence."""
        return len(self.sequence)

    ###
    ### PROPERTIES
    ###

    @property
    def coords(self) -> np.ndarray:
        """
        Get the 3D coordinates.
        """
        return self._coords

    @coords.setter
    def coords(self, coords: np.ndarray) -> None:
        """
        Set the 3D coordinates and ensure they match the sequence length.

        Parameters
        ----------
        coords : np.ndarray
            A NumPy array of shape (N, 3) representing the coordinates.

        Raises
        ------
        ValueError
            If the coordinates are not a valid NumPy array or do not
            match the sequence length.
        """
        if isinstance(coords, (list, tuple)):
            coords = np.array(coords)
        if not isinstance(coords, np.ndarray):
            raise ValueError("The coords must be a numpy array")
        if (
            self._sequence
            and self._coords.size > 0
            and len(self.sequence) != len(coords)
        ):
            raise ValueError(
                f"The coords and the sequence must have the same length:"
                f" expected len{len(self.sequence)},"
                f"  got {len(coords)} coordinates"
            )
        self._coords = np.array(coords)

    @property
    def sequence(self) -> str:
        """
        Get the protein sequence.
        """
        return self._sequence

    @sequence.setter
    def sequence(self, sequence: str) -> None:
        """
        Set the protein sequence and ensure it matches the coordinate length.

        Parameters
        ----------
        sequence : str
            The amino acid sequence.

        Raises
        ------
        ValueError
            If the sequence length does not match the coordinate length.
        """
        if not isinstance(sequence, str):
            raise ValueError("The sequence must be a string")
        if (
            self._sequence
            and self._coords.size > 0
            and len(sequence) != len(self.coords)
        ):
            raise ValueError("The sequence must have the same length as the coords")
        self._sequence = sequence

    ###
    ### PUBLIC METHODS
    ###

    def copy(self) -> "ProteinCoords":
        """Create a copy of the ProteinCoords instance."""
        return ProteinCoords(self.sequence, self.coords)

    def transform(self, T_matrix: np.ndarray) -> None:
        """
        Apply a transformation matrix to the coordinates.

        Parameters
        ----------
        T_matrix : np.ndarray
            A 4x4 transformation matrix.
        """
        Coords.transform_array(self._coords, T_matrix)


class Coords:
    """
    A class for representing and manipulating 3D molecular coordinates.

    This class provides functionality to store, transform, and manipulate
    3D coordinates used in molecular structures. It supports transformations
    such as rotation and translation, as well as operations for combining
    multiple coordinate sets.
    The Coords class includes functions to build RNA-A helices based on
    the oxRNA model parameters and to handle dovetail crossovers.
    The RNA-A helix parameters and dovetail crossover parameters can be
    adjusted using class methods `set_helix_params` and `set_AE_crossover_params`.

    Parameters
    ----------
    input_array : Union[np.ndarray, List], optional
        A NumPy array or list containing the 3D coordinates, where each element
        represents a nucleotide or molecular unit with associated vectors.
    dummy_ends : Tuple[np.ndarray, np.ndarray], optional
        A tuple containing two NumPy arrays representing dummy end coordinates.
    proteins : List[ProteinCoords], optional
        A list of `ProteinCoords` instances associated with the molecular structure.

    Attributes
    ----------
    array : np.ndarray
        The main coordinate array storing molecular positions.
    dummy_ends : Tuple[np.ndarray, np.ndarray]
        Dummy end coordinates used for alignment or structural constraints.
    proteins : List[ProteinCoords]
        A list of associated `ProteinCoords` objects representing protein components.
    shape : Tuple[int, ...]
        Shape of the coordinates array.
    size : int
        Total number of elements in the coordinates array.

    Class Attributes
    ----------------
    _CACHED_HELICES : dict
        A cache for storing previously computed helical coordinates.
        This is used to optimize performance by avoiding redundant calculations,
        since the Origami are usually built upon repeating helical structures.
    """

    ### --- CACHED DATA ---
    _CACHED_HELICES = dict()
    _CACHED_AE_T = None
    _CACHED_AE_T_INV = None

    ### --- RNA-A HELIX PARAMETERS ---
    # Inclination angle in radians
    inclination = -15.5 * np.pi / 180
    # Distance between base pairs along the backbone
    bp_backbone_distance = 2.0
    # Diameter of the helix
    diameter = 2.35
    # Distance between bases along the helix axis
    base_base_distance = 0.3287
    # Rotation angle between next bases in radians
    rot_angle = 32.73 * np.pi / 180

    ### --- DOVETAIL AE CROSSOVER PARAMETERS ---
    # these values optimize the symmetry of the AE crossover
    # they are obtaiend from a initial guess (manual tuning) and then
    # optimized with a numerical minimization of the asymmetry

    # Shift along the base pair plane
    ae_v_shift = 1.006
    # Rotation around the helix axis in radians
    ae_rot = 40.267 * np.pi / 180
    # Shift in the helix x direction
    ae_hx_shift = 2.293
    # Shift in the helix y direction
    ae_hy_shift = -0.5214

    def __init__(
        self,
        input_array: Union[np.ndarray, List] = np.array(()),
        dummy_ends: Tuple[np.ndarray, np.ndarray] = (np.array(()), np.array(())),
        proteins: List["ProteinCoords"] = None,
    ):
        """
        Initialize Coords object.

        Parameters
        ----------
        input_array : Union[np.ndarray, List]
            A numpy array or list representing the coordinates.
        dummy_ends : Tuple[np.ndarray, np.ndarray]
            A tuple of two numpy arrays representing dummy ends.
        proteins : List[ProteinCoords], optional
            A list of ProteinCoords instances.
        """
        if proteins is None:
            proteins = []
        # initialize the attributes
        self._dummy_ends = ()
        self._array = ()

        # set the attributes
        self.dummy_ends = dummy_ends
        self.proteins = proteins
        self.array = input_array

    def __getitem__(self, key: int) -> np.ndarray:
        """Get item from array."""
        return self.array[key]

    def __setitem__(self, key: int, value: np.ndarray) -> None:
        """Set item in array."""
        self.array[key] = value

    def __str__(self) -> str:
        """Return a string representation of the coordinates."""
        return f"Coords({self.array.tolist()})"

    def __len__(self) -> int:
        """Return the number of elements in the coordinates array."""
        return len(self.array)

    ###
    ### PROPERTIES
    ###

    @property
    def array(self) -> np.ndarray:
        """
        Return the coordinates array.
        """
        return self._array

    @array.setter
    def array(self, new_array: Union[np.ndarray, List, Tuple]) -> None:
        """
        Set the coordinates array.

        Parameters
        ----------
        new_array : Union[np.ndarray, List, Tuple]
            The new array to be set.
        """
        if not isinstance(new_array, (np.ndarray, list, tuple)):
            raise ValueError(
                "The Coordinates array must be a numpy array or "
                "a list of coordinates"
            )
        self._array = np.array(new_array)

    @property
    def dummy_ends(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Return the dummy ends.
        """
        return self._dummy_ends

    @dummy_ends.setter
    def dummy_ends(self, new_dummy: Tuple[np.ndarray, np.ndarray]) -> None:
        """
        Set the dummy ends.

        Parameters
        ----------
        new_dummy : Tuple[np.ndarray, np.ndarray]
            A tuple of numpy arrays representing dummy ends.
        """
        if not isinstance(new_dummy, (list, tuple)) or len(new_dummy) != 2:
            raise ValueError(
                "The dummy_ends argument must be a tuple " "of two numpy arrays."
            )

        self._dummy_ends = (np.array(new_dummy[0]), np.array(new_dummy[1]))

    @property
    def proteins(self) -> List["ProteinCoords"]:
        """
        Return the list of proteins.
        """
        return self._proteins

    @proteins.setter
    def proteins(self, proteins: List["ProteinCoords"]) -> None:
        """
        Set the proteins list.

        Parameters
        ----------
        proteins : List[ProteinCoords]
            A list of ProteinCoords instances.
        """
        if not isinstance(proteins, list) or any(
            not isinstance(p, ProteinCoords) for p in proteins
        ):
            raise ValueError(
                "The proteins argument must be a list of ProteinCoords "
                f"instances. Got: {proteins}"
            )
        self._proteins = proteins

    @property
    def shape(self) -> Tuple[int, ...]:
        """
        Return the shape of the coordinates array.
        """
        return self.array.shape

    @property
    def size(self) -> int:
        """
        Return the size of the coordinates array.
        """
        return self.array.size

    ###
    ### CLASS METHODS
    ###

    @classmethod
    def from_json(cls, data: Dict[str, Any]) -> "Coords":
        """
        Create a Coords instance from a JSON-like dictionary.

        Parameters
        ----------
        data : Dict[str, Any]
            A dictionary containing the coordinate data.

        Returns
        -------
        Coords
            A Coords instance created from the provided data.
        """
        input_array = np.array(data.get("array", []))
        dummy_ends_data = data.get("dummy_ends", ([], []))
        dummy_ends = (np.array(dummy_ends_data[0]), np.array(dummy_ends_data[1]))
        proteins_data = data.get("proteins", [])
        proteins = [
            ProteinCoords(
                sequence=protein_data.get("sequence", ""),
                coords=np.array(protein_data.get("coords", [])),
            )
            for protein_data in proteins_data
        ]
        return cls(input_array=input_array, dummy_ends=dummy_ends, proteins=proteins)

    @classmethod
    def load_from_file(
        cls,
        filename: str,
        dummy_ends: Tuple[bool, bool] = (False, False),
        extend: Tuple[int, int] = (0, 0),
        topology_file: str = None,
        protein: bool = False,
        verbose: bool = False,
    ) -> "Coords":
        """
        Load and cleanup coordinates from an oxDNA configuration or PDB file.

        Parameters
        ----------
        filename : str
            Path to the oxDNA configuration or PDB file.
        dummy_ends : Tuple[bool, bool], default is (False, False)
            A tuple indicating whether the coordinates have dummy ends at the
            beginning or end.
        extend : Tuple[int, int], default is (0, 0)
            Number of nucleotides to extend the coordinates at the beginning and end.
        topology_file : str, optional
            The path to the topology file (only required for proteins).
        protein : bool, default False
            Whether to load protein coordinates.
        verbose : bool, default False
            Whether to print the generated code for the coordinates.

        Returns
        -------
        Coords
            An instance of the Coords class with the cleaned-up coordinates.
        """
        # Check if the file is a pdb file
        pdb = False
        top_text = None  # initialize the topology text

        # check the path
        if not isinstance(filename, str):
            try:
                filename = str(filename)
            except Exception as e:
                raise ValueError(
                    "The filename must be a string or a path-like "
                    f"object, got {type(filename)}. Full error: {e}"
                )

        # check the file extension
        if filename.endswith(".pdb") or filename.endswith(".PDB"):
            pdb = True

        # read the file
        with open(filename, "r") as f:
            conf_text = f.read()

        if not pdb and topology_file:
            with open(topology_file, "r") as f:
                top_text = f.read()

        return cls.load_from_text(
            conf_text,
            dummy_ends,
            extend,
            verbose,
            top_text=top_text,
            pdb_format=pdb,
            protein=protein,
        )

    @classmethod
    def load_from_text(
        cls,
        conf_text: str,
        dummy_ends: Tuple[bool, bool] = (False, False),
        extend: Tuple[int, int] = (0, 0),
        verbose: bool = False,
        top_text: str = None,
        pdb_format: bool = False,
        protein: bool = False,
    ) -> "Coords":
        """
        Load and cleanup coordinates from an oxDNA or PDB configuration text.

        Parameters
        ----------
        conf_text : str
            OxDNA configuration text.
        dummy_ends : Tuple[bool, bool], default is (False, False)
            A tuple indicating whether the coordinates have dummy ends at the
            beginning or end.
        extend : Tuple[int, int], default is (0, 0)
            Number of nucleotides to extend the coordinates at the beginning and end.
        verbose : bool, default False
            Whether to print the generated code for the coordinates.
        top_text : str, optional
            The topology file text (only required for proteins).
        pdb_format : bool, default False
            Whether the configuration text is in PDB format.
        protein : bool, default False
            Whether to load protein coordinates.

        Returns
        -------
        Coords
            A Coords object
        """
        if pdb_format:
            if not oat_installed:
                raise ValueError(
                    "oxDNA_analysis_tools not installed. "
                    "Please install it to use the pdb option"
                )
            confs, systems = PDB_oxDNA(conf_text)
            conf_text = conf_to_str(confs[0])
            top_text = get_top_string(systems[0])

        lines = conf_text.split("\n")

        ### EXTRACT THE COORDINATES
        def extract_coords_from_line(line):
            splitted = line.strip().split()
            return [
                [float(p) for p in splitted[0:3]],
                [float(b) for b in splitted[3:6]],
                [float(n) for n in splitted[6:9]],
            ]

        coords = np.array(
            [extract_coords_from_line(line) for line in lines[3:] if line.strip()]
        )

        ### LOAD PROTEINS
        all_prot_coords = []
        protein_coords = None
        if verbose:
            protein_text = ""
        if protein:

            if top_text is None:
                raise ValueError(
                    "A topology file path must be provided when loading"
                    " a structure with proteins"
                )

            # initialize protein text
            if verbose:
                protein_text = ",\n\tproteins=["
            # Split the topology text into lines
            lines = top_text.split("\n")[1:]
            seq_start_ind = 0
            for line in lines:
                if not line:
                    continue
                # get the sequence
                seq = line.split()[0]
                # Detect a protein sequence
                if "peptide" in line:

                    # load the protein coordinates
                    protein_coords = coords[seq_start_ind : seq_start_ind + len(seq)]
                    if verbose:
                        protein_text += (
                            f'ProteinCoords("{seq}", '
                            f"coords={protein_coords.tolist()}), "
                        )
                    all_prot_coords.append(ProteinCoords(seq, protein_coords))

                    # remove the protein coordinates from the main coordinates
                    coords = np.delete(
                        coords, np.s_[seq_start_ind : seq_start_ind + len(seq)], axis=0
                    )
                    seq_start_ind -= len(seq)

                # update the coordinates index
                seq_start_ind += len(seq)

            if verbose:
                protein_text += "]"

        ### ENTEND THE COORDINATES
        if extend[0] > 0:
            extend_first = Coords.compute_helix_from_nucl(
                coords[0][0],
                coords[0][1],
                coords[0][2],
                extend[0],
                directionality="35",
                double=False,
                use_cached=False,
            )[::-1]
            coords = np.concatenate((extend_first, coords))

        if extend[1] > 0:
            extend_last = Coords.compute_helix_from_nucl(
                coords[-1][0],
                coords[-1][1],
                coords[-1][2],
                extend[1],
                directionality="53",
                double=False,
                use_cached=False,
            )
            coords = np.concatenate((coords, extend_last))

        ### CREATE DUMMIES
        dummy_0 = np.array(())
        dummy_1 = np.array(())
        dummy_text = ""

        if dummy_ends[0]:
            dummy_0 = coords[0]
            coords = coords[1:]

        if dummy_ends[1]:
            dummy_1 = coords[-1]
            coords = coords[:-1]

        if verbose:
            dummy_text = f",\n\tdummy_ends=({dummy_0.tolist()}, {dummy_1.tolist()})"

        if verbose:
            print(f"Coords({coords.tolist()}{dummy_text}{protein_text})")

        return cls(coords, dummy_ends=(dummy_0, dummy_1), proteins=all_prot_coords)

    @classmethod
    def set_helix_params(
        cls,
        inclination: float = -15.5 * np.pi / 180,
        bp_backbone_distance: float = 2.0,
        diameter: float = 2.35,
        base_base_distance: float = 0.3287,
        rot_angle: float = 32.7 * np.pi / 180,
    ) -> None:
        """
        Set the helix parameters for coordinate computations.
        Parameters
        ----------
        inclination : float, default -15.5 * np.pi / 180
            Inclination angle in radians.
        bp_backbone_distance : float, default 2.0
            Distance between base pairs along the backbone.
        diameter : float, default 2.35
            Diameter of the helix.
        base_base_distance : float, default 0.3287
            Distance between bases along the helix axis.
        rot_angle : float, default 32.7 * np.pi / 180
            Rotation angle between next bases in radians.

        Returns
        -------
        None
        """
        cls.inclination = inclination
        cls.bp_backbone_distance = bp_backbone_distance
        cls.diameter = diameter
        cls.base_base_distance = base_base_distance
        cls.rot_angle = rot_angle

        # clear the cached helices since the parameters have changed
        cls._CACHED_HELICES = dict()
        cls._CACHED_AE_T = None
        cls._CACHED_AE_T_INV = None

    @classmethod
    def set_AE_crossover_params(
        cls,
        ae_v_shift: float = 1.006,
        ae_rot: float = 40.267 * np.pi / 180,
        ae_hx_shift: float = 2.293,
        ae_hy_shift: float = -0.5214,
    ) -> None:
        """
        Set the AE crossover parameters for coordinate computations.

        Parameters
        ----------
        ae_v_shift : float, default 0.99
            Shift along the base pair plane.
        ae_rot : float, default 40.52 * np.pi / 180
            Rotation around the helix axis in radians.
        ae_hx_shift : float, default 2.23
            Shift in the helix x direction.
        ae_hy_shift : float, default -0.52
            Shift in the helix y direction.

        Returns
        -------
        None
        """
        cls.ae_v_shift = ae_v_shift
        cls.ae_rot = ae_rot
        cls.ae_hx_shift = ae_hx_shift
        cls.ae_hy_shift = ae_hy_shift
        # clear the cached dovetail transformations
        cls._CACHED_AE_T = None
        cls._CACHED_AE_T_INV = None

    ###
    ### STATIC METHODS
    ###

    @staticmethod
    def apply_transformation(
        T: np.ndarray, p: np.ndarray, b: np.ndarray, n: np.ndarray, local: bool = False
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Apply a transformation matrix to position, base, and normal vectors.

        Parameters
        ----------
        T : np.ndarray
            A 4x4 transformation matrix.
        p : np.ndarray
            Position vector (3D).
        b : np.ndarray
            Base vector (3D).
        n : np.ndarray
            Normal vector (3D).
        local : bool, default False
            Whether to apply transformation in a local reference frame.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray, np.ndarray]
            Transformed position, base, and normal vectors.
        """
        if local:
            p_local, b_local, n_local = p, b, n
            p, b, n = np.array((0, 0, 0)), np.array((1, 0, 0)), np.array((0, 1, 0))

        # Convert position to homogeneous coordinates
        p_homogeneous = np.append(p, 1)

        # Apply the transformation to the position
        p_transformed_homogeneous = T @ p_homogeneous.T
        p_trans = p_transformed_homogeneous.T[:3]

        # Extract the rotation matrix
        rotation_matrix = T[:3, :3]
        rotation = R.from_matrix(rotation_matrix)

        # Apply the rotation to the base and normal vectors
        b_trans = rotation.apply(b)
        n_trans = rotation.apply(n)

        if local:
            # find the original reference system from the local one
            p_trans, b_trans, n_trans = Coords.set_reference(
                p_local, b_local, n_local, p_trans, b_trans, n_trans, reverse=True
            )

        return p_trans, b_trans, n_trans

    @staticmethod
    def combine_coords(
        strand1: "Strand",
        strand2: "Strand",
    ) -> "Coords":
        """
        Combine the coordinates of two strands by applying necessary transformations.
        Use dependency injection (with strand1 and strand2 as arguments) to apply
        the transformation to a whole Strand Block. The second strand is transformed
        to match the first strand's orientation.

        Parameters
        ----------
        strand1 : Strand
            First strand to combine the coordinates.
        strand2 : Strand
            Second strand to combine the coordinates.

        Returns
        -------
        Coords
            Combined coordinate object with the transformed coordinates of the strands.
        """
        coords1 = strand1.coords
        coords2 = strand2.coords

        ### HEDGE CASES
        # there are no coordinates to combine
        if coords1.is_empty() and coords2.is_empty():
            return Coords([])
        # the first strand has no coordinates and no nucleotides
        elif coords1.is_empty():
            return coords2.copy()
        # the second strand has no coordinates and no nucleotides
        elif coords2.is_empty():
            return coords1.copy()

        ### TRANSFORM THE SECOND STRAND
        # transform them only if they are not in the same block os strands
        if strand1.strands_block is not strand2.strands_block:

            ### technically, `strand1.coords will create the coordinates
            ### if they are not present, so this code block is not
            ### needed anymore

            ### -> Start obsolete code
            # # create the coordinates for the first strand
            # if coords1.size == 0 and len(strand1.sequence):
            #     coords1 = Coords.compute_helix_from_nucl(
            #         (0, 0, 0),  # start position
            #         (1, 0, 0),  # base vector
            #         (0, 1, 0),  # normal vector
            #         length=len(strand1.sequence),
            #         directionality=strand1.directionality,
            #     )

            # # create the coordinates for the second strand
            # if coords2.size == 0 and len(strand2.sequence):
            #     coords2 = Coords.compute_helix_from_nucl(
            #         (0, 0, 0),  # start position
            #         (1, 0, 0),  # base vector
            #         (0, 1, 0),  # normal vector
            #         length=len(strand2.sequence),
            #         directionality=strand2.directionality,
            #     )
            ### -> End obsolete code

            ### CHECK THE JOINING POINTS AND TRANSFORM THE ARRAY TWO

            # the first strand has no end dummy
            if coords1.dummy_ends[1].size == 0:
                # the second strand has a start dummy
                if coords2.dummy_ends[0].size > 0:
                    # the dummy start translates to...
                    pos1 = coords2.dummy_ends[0][0]
                    bv1 = coords2.dummy_ends[0][1]
                    nv1 = coords2.dummy_ends[0][2]
                    # ...the last position of the first strand
                    pos2 = coords1[-1][0]
                    bv2 = coords1[-1][1]
                    nv2 = coords1[-1][2]
                else:
                    # the second strand has no dummy: calulcate an helical join from
                    # the last nucleotide of the first strand
                    helix_coords = Coords.compute_helix_from_nucl(
                        coords1[-1][0],
                        coords1[-1][1],
                        coords1[-1][2],
                        length=1,
                        directionality=strand1.directionality,
                    )
                    # the first position of the second strand translates to...
                    pos1 = coords2[0][0]
                    bv1 = coords2[0][1]
                    nv1 = coords2[0][2]
                    # ...the first nucleotide of the helical join
                    pos2 = helix_coords[0][0]
                    bv2 = helix_coords[0][1]
                    nv2 = helix_coords[0][2]

            # the first strand has an end dummy, use it
            else:
                # the array 2 has no coordinate, use the dummy instead
                if coords2.size == 0:
                    # the dummy start translates to...
                    pos1 = coords2.dummy_ends[0][0]
                    bv1 = coords2.dummy_ends[0][1]
                    nv1 = coords2.dummy_ends[0][2]
                else:
                    # the first position of the second strand translates to...
                    pos1 = coords2[0][0]
                    bv1 = coords2[0][1]
                    nv1 = coords2[0][2]
                # ...the dummy end of the first strand
                pos2 = coords1.dummy_ends[1][0]
                bv2 = coords1.dummy_ends[1][1]
                nv2 = coords1.dummy_ends[1][2]
            # Calculate the transformation matrix
            T_matrix = Coords.compute_transformation_matrix(
                pos1, bv1, nv1, pos2, bv2, nv2
            )

            ### Apply the transformation to the second strand motif block
            strand2.strands_block.transform(T_matrix)

        ### COMBINE THE COORDINATES,
        # the coordinates could be empty if there are dummies
        combined_dummy = (coords1.dummy_ends[0], coords2.dummy_ends[1])

        # both strands have no coordinates
        if coords1.size == 0 and coords2.size == 0:
            combined = Coords([], dummy_ends=combined_dummy)

        # the first strand has no coordinates
        elif coords1.size == 0:
            combined = Coords(coords2.array, dummy_ends=combined_dummy)

        # the second strand has no coordinates
        elif coords2.size == 0:
            combined = Coords(coords1.array, dummy_ends=combined_dummy)

        else:
            combined = Coords(
                np.concatenate((coords1.array, coords2.array), axis=0),
                dummy_ends=combined_dummy,
            )

        ### COMBINE THE PROTEINS
        combined.proteins = coords1.proteins + coords2.proteins

        return combined

    @staticmethod
    def compute_helix_from_nucl(
        pos: np.ndarray,
        old_a1: np.ndarray,
        old_a3: np.ndarray,
        length: int,
        directionality: Literal["53", "35"] = "53",
        double: bool = False,
        use_cached: bool = True,
        return_dir_vectors: bool = False,
    ) -> "Coords":
        """
        Compute helical coordinates for nucleotides.
        This function is based on the oxView implementation of the
        helix generation: https://doi.org/10.1093/nar/gkaa417


        Parameters
        ----------
        pos : np.ndarray
            Starting position.
        old_a1 : np.ndarray
            Initial base vector.
        old_a3 : np.ndarray
            Initial normal vector.
        length : int
            Number of nucleotides.
        directionality : str, default "53"
            Direction of the helix ("53" or "35").
        double : bool, default False
            If True, generate double-stranded coordinates.
        use_cached : bool, default False
            If True, use or save cached coordinates for faster computation.
        return_dir_vectors : bool, default False
            If True, return the direction vector and start position also.
            This disables the caching mechanism, because the output is not
            only the coordinates.

        Returns
        -------
        Coords
            Generated helical coordinates.
        """
        if use_cached and not return_dir_vectors:
            info = (
                tuple(pos),
                tuple(old_a1),
                tuple(old_a3),
                directionality,
                double,
                length,
            )
            if info in Coords._CACHED_HELICES:
                return Coords._CACHED_HELICES[info]

        # Adjust arguments: Convert lists or tuples to numpy arrays
        # for consistency in mathematical operations
        direction_35 = directionality == "35"
        if isinstance(pos, (list, tuple)):
            pos = np.array(pos)
        if isinstance(old_a1, (list, tuple)):
            old_a1 = np.array(old_a1)
        if isinstance(old_a3, (list, tuple)):
            old_a3 = np.array(old_a3)

        cord = np.cos(Coords.inclination) * Coords.bp_backbone_distance  # Chord length
        # Distance from center to chord
        center_to_cord = np.sqrt((Coords.diameter / 2) ** 2 - (cord / 2) ** 2)
        fudge = 0.4  # Fudge factor for position adjustment

        # Calculate the axis perpendicular to the old A1 and A3 vectors
        norm_a1_a3 = np.cross(old_a1, old_a3)
        # Create the rotation matrix for the direction vector
        R_dir = R.from_rotvec(norm_a1_a3 * Coords.inclination)
        # Calculate the direction vector,
        # rotating the old normal vector by the inclination angle
        dir_vector = -R_dir.apply(old_a3)

        if direction_35:
            dir_vector *= -1  # Adjust direction if specified

        dir_vector /= np.linalg.norm(dir_vector)  # Normalize direction vector

        # Calculate the coordinates if the helix axis is the Z-axis
        x1, y1, z1 = (
            center_to_cord,
            -cord / 2,
            -(Coords.bp_backbone_distance / 2) * np.sin(Coords.inclination),
        )
        x2, y2, z2 = (
            center_to_cord,
            +cord / 2,
            +(Coords.bp_backbone_distance / 2) * np.sin(Coords.inclination),
        )
        r1 = np.array([x1, y1, z1])
        r2 = np.array([x2, y2, z2])

        ### Set the axis to the correct one
        ref_vector = np.array([0, 0, 1])  # Z-axis is the default axis
        # Calculate the cross product for rotation
        cross_prod = np.cross(ref_vector, dir_vector)
        # Calculate the dot product
        dot_prod = np.dot(ref_vector, dir_vector)
        # Calculate the quaternion components
        # * scalar component to normalize the quaternion
        scalar = np.sqrt((1.0 + dot_prod) * 2.0)

        # The quaternion matrix transforms the ref_vector to the dir_vector
        q_matrix = np.append(cross_prod / scalar, scalar * 0.5)
        q_matrix /= np.linalg.norm(q_matrix)  # Normalize the quaternion
        q1 = R.from_quat(q_matrix)  # Create the rotation from the quaternion
        r1 = q1.apply(r1)  # Apply rotation to r1
        r2 = q1.apply(r2)  # Apply rotation to r2

        # Set a1 to the correct orientation
        r1_to_r2 = r2 - r1
        if direction_35:
            r1_to_r2 = -r1_to_r2  # Adjust direction if specified
        r1_to_r2 /= np.linalg.norm(r1_to_r2)  # Normalize the vector

        ### Rotate the helix axis to the correct orientation
        # Project r1_to_r2 onto dir_vector
        r1_to_r2_proj = r1_to_r2 - np.dot(r1_to_r2, dir_vector) * dir_vector
        # Project old_a1 onto dir_vector
        old_a1_proj = old_a1 - np.dot(old_a1, dir_vector) * dir_vector

        rotAngle2 = np.arccos(
            np.clip(
                np.dot(
                    r1_to_r2_proj / np.linalg.norm(r1_to_r2_proj),
                    old_a1_proj / np.linalg.norm(old_a1_proj),
                ),
                -1.0,
                1.0,
            )
        )

        if np.dot(np.cross(r1_to_r2, old_a1), dir_vector) < 0:
            rotAngle2 *= -1  # Adjust rotation angle based on the cross product
        q2 = R.from_rotvec(dir_vector * rotAngle2)  # Create the rotation matrix
        r1 = q2.apply(r1)  # Apply rotation to r1
        r2 = q2.apply(r2)  # Apply rotation to r2

        # Select the backbone position as reference (r1 for 5'->3', r2 for 3'->5')
        r = r1
        if direction_35:
            r = r2  # Adjust direction if specified
        start_pos = pos - r - old_a1 * fudge  # Calculate the starting position

        # Create per-step rotation matrix
        R_step = R.from_rotvec(
            dir_vector * Coords.rot_angle
        )  # Create rotation for each step

        # Initialize properties of new nucleotide
        # Preallocate result arrays
        out = np.empty((length, 3, 3), dtype=np.float64)
        out_double = np.empty((length, 3, 3), dtype=np.float64) if double else None

        # Precompute the steps
        step_vector = dir_vector * Coords.base_base_distance
        R_step_apply = R_step.apply

        # Generate nucleotide positions and orientations
        for i in range(length):
            # Calculate rotation around central axis and step along axis
            r1 = R_step_apply(r1) + step_vector
            r2 = R_step_apply(r2) + step_vector

            # Calculate a1 orientation
            r1_to_r2 = r2 - r1
            a1 = r1_to_r2 / np.linalg.norm(r1_to_r2)  # Normalize a1

            ### Calculate A3 orientation
            # Project a1 onto dir_vector
            a1proj = a1 - np.dot(a1, dir_vector) * dir_vector
            a1proj /= np.linalg.norm(a1proj)  # Normalize projection
            # Calculate a3
            a3 = (
                -np.cos(Coords.inclination) * dir_vector
                + np.sin(Coords.inclination) * a1proj
            )
            a3 /= np.linalg.norm(a3)  # Normalize a3

            # Calculate position
            r = r1
            if direction_35:
                r = r2  # Adjust direction if specified
                a1 = -a1
                a3 = -a3
            RNA_fudge = a1 * fudge  # Apply fudge factor
            p = r + RNA_fudge + start_pos  # Calculate position

            out[i] = [p, a1, a3]  # Append nucleotide properties

            # Double helix case
            if double:
                # Calculate a1 for double helix
                a1 = -r1_to_r2 / np.linalg.norm(-r1_to_r2)
                # Project a1 onto dir_vector
                a1proj = a1 - np.dot(a1, -dir_vector) * -dir_vector
                # Normalize projection
                a1proj /= np.linalg.norm(a1proj)
                # Calculate a3
                a3 = (
                    np.cos(Coords.inclination) * dir_vector
                    + np.sin(Coords.inclination) * a1proj
                )
                a3 /= np.linalg.norm(a3)  # Normalize a3
                r = r2
                if direction_35:
                    r = r1  # Adjust direction if specified
                    a1 = -a1
                    a3 = -a3
                RNA_fudge = a1 * fudge  # Apply fudge factor
                p = r + RNA_fudge + start_pos  # Calculate position
                # Append nucleotide properties for double helix
                out_double[i] = [p, a1, a3]

        # Combine results
        if double:
            # Combine both strands for double helix, both directions are 5' to 3'
            combined = np.concatenate((out, out_double[::-1]), axis=0)
        else:
            combined = out

        if use_cached:
            # Cache the helical coordinates for future use
            Coords._CACHED_HELICES[info] = Coords(combined)

        if return_dir_vectors:
            return Coords(combined), dir_vector, start_pos
        return Coords(combined)

    @staticmethod
    def compute_AE_crossover(
        helix_length: int = 8,
        return_helices: bool = False,
        use_cached: bool = True,
        compute_all_T: bool = False,
    ) -> Union[Tuple["Coords", "Coords"], Tuple[np.ndarray, np.ndarray]]:
        """
        Compute the coordinates for an A-form AE crossover.
        It generates two parallel helices, then translate and rotate
        a second helix to create the crossover geometry between the middle
        nucleotides of each helix.
        This function is based on the parameters set in `set_AE_crossover_params`.
        The AE parameters are based on manual tuning and scipy refinement so that
        the AE transformation matrix of the two crossover strands is the same.

        Parameters
        ----------
        helix_length : int, default 8
            Length of each helix in nucleotides.
        return_helices : bool, default True
            If True, return the two Coords objects for the helices
            used to create the crossover.
            If False, return the transformation matrices.
        use_cached : bool, default True
            If True, use cached transformation matrices if available.
        compute_all_T : bool, default False
            If True, compute and return additional transformation matrices.

        Returns
        -------
        Union[Tuple[Coords, Coords], Tuple[np.ndarray, np.ndarray]]
            If `return_helices` is True, returns two Coords objects for the helices.
            If `return_helices` is False, returns the transformation matrices.
            If `compute_all_T` is True, returns additional transformation matrices
            for the crossover.

        """
        if (
            use_cached
            and not return_helices
            and not compute_all_T
            and Coords._CACHED_AE_T is not None
            and Coords._CACHED_AE_T_INV is not None
        ):
            return Coords._CACHED_AE_T, Coords._CACHED_AE_T_INV

        pos = np.array((0, 0, 0))
        old_a1 = np.array((1, 0, 0))
        old_a3 = np.array((0, 1, 0))

        h1, dir_vector, start_pos = Coords.compute_helix_from_nucl(
            pos,  # start position
            old_a1,  # base vector
            old_a3,  # normal vector
            length=helix_length,
            double=True,
            # we need the direction vector:
            use_cached=False,
            return_dir_vectors=True,
        )
        h2 = h1.copy()
        # move h2 up along the direction vector
        # equivalent of "step_vector = dir_vector * base_base_distance"
        translate = dir_vector * Coords.ae_v_shift  # 3.55 Amstrongs up
        Transl_mat = np.array(
            [
                [1, 0, 0, translate[0]],
                [0, 1, 0, translate[1]],
                [0, 0, 1, translate[2]],
                [0, 0, 0, 1],
            ]
        )
        h2.transform(Transl_mat)

        ### PREPARE FOR TRANSLATION
        middle = helix_length // 2
        nt = h1[middle][0]  # position of first nucleotide

        # project nt onto the *actual* axis: start_pos + t * dir_vector
        nt_rel = nt - start_pos  # move into axis frame
        # dir_vector is unit, so no denominator needed
        proj_coeff = np.dot(nt_rel, dir_vector)
        proj = start_pos + proj_coeff * dir_vector  # point on the axis closest to nt

        v = nt - proj  # vector from axis to nt
        v /= np.linalg.norm(v)  # normalize

        # get the oterher orthogonal vector
        ortho_vector = np.cross(dir_vector, v)
        ortho_vector /= np.linalg.norm(ortho_vector)
        xy_shift = v * Coords.ae_hx_shift + ortho_vector * Coords.ae_hy_shift

        # PREPARE FOR ROTATION AND ROTATE
        h2_array = h2.array  # shape (N, 3, 3)
        q_rot = R.from_rotvec(dir_vector * Coords.ae_rot)  # rotation quaternion

        # split into position + orientation arrays
        pos = h2_array[:, 0, :]  # (N, 3)
        a1 = h2_array[:, 1, :]  # (N, 3)
        a3 = h2_array[:, 2, :]  # (N, 3)

        # --- rotate positions around the axis ---
        # 1) translate so axis_point is at origin
        pos_shifted = pos - start_pos

        # 2) rotate around axis (through origin)
        pos_rot = q_rot.apply(pos_shifted)

        # 3) translate back
        pos_rot += start_pos

        # --- rotate orientations (NO translation!) ---
        a1_rot = q_rot.apply(a1)
        a3_rot = q_rot.apply(a3)

        # --- recombine ---
        rotated_array = np.stack([pos_rot, a1_rot, a3_rot], axis=1)
        h2 = Coords(rotated_array)

        # TRANSLATE
        Transl_mat2 = np.array(
            [
                [1, 0, 0, xy_shift[0]],
                [0, 1, 0, xy_shift[1]],
                [0, 0, 1, xy_shift[2]],
                [0, 0, 0, 1],
            ]
        )

        h2.transform(Transl_mat2)

        first_nt = h1[middle]
        second_nt = h2[-middle - 1]
        ae_T = Coords.compute_transformation_matrix(
            first_nt[0],
            first_nt[1],
            first_nt[2],
            second_nt[0],
            second_nt[1],
            second_nt[2],
            local=True,
        )
        ae_T_inv = Coords.compute_transformation_matrix(
            second_nt[0],
            second_nt[1],
            second_nt[2],
            first_nt[0],
            first_nt[1],
            first_nt[2],
            local=True,
        )

        Coords._CACHED_AE_T = ae_T
        Coords._CACHED_AE_T_INV = ae_T_inv

        if return_helices:
            return h1, h2

        if compute_all_T:
            first_nt = h2[-middle - 2]
            second_nt = h1[middle + 1]
            ae_T_2 = Coords.compute_transformation_matrix(
                first_nt[0],
                first_nt[1],
                first_nt[2],
                second_nt[0],
                second_nt[1],
                second_nt[2],
                local=True,
            )
            ae_T_inv_2 = Coords.compute_transformation_matrix(
                second_nt[0],
                second_nt[1],
                second_nt[2],
                first_nt[0],
                first_nt[1],
                first_nt[2],
                local=True,
            )

            return ae_T, ae_T_inv, ae_T_2, ae_T_inv_2

        return ae_T, ae_T_inv

    @staticmethod
    def compute_kissing_loop180(
        second_strand: bool = False, branched=False
    ) -> "Coords":
        """
        Compute the coordinates for a 3D kissing loop with 180-degree angle.
        This structure is a key element in RNA nanostructures, so it's important
        to keep it parameterized according to the RNA-A helix parameters.

        Parameters
        ----------
        second_strand : bool, default False
            If True, it generates a complementary kissing loop to the first one.
        branched : bool, default False
            If True, it generates the coordinates with a branched kissing loop,
            (both connecting and kissing strands) with dummy ends for branching.
            The transformation matricess are adapted from the 3D model of branched
            kissing loops from ROAD (https://doi.org/10.1038/s41557-021-00679-1).

        Returns
        -------
        Coords
            Generated coordinates for the 3D kissing loop.
        """

        def middle_nucleotide_coords(n1, n2):
            """
            Compute the coordinates of a middle nucleotide between two nucleotides.

            Parameters
            ----------
            n1 : np.ndarray
                Coordinates of the first nucleotide (position, base, normal).
            n2 : np.ndarray
                Coordinates of the second nucleotide (position, base, normal).

            Returns
            -------
            np.ndarray
                Coordinates of the middle nucleotide (position, base, normal).
            """
            # Alternative more complex method to compute middle nucleotide coords
            # p1, b1, n1 = n1[0], n1[1], n1[2]
            # p2, b2, n2 = n2[0], n2[1], n2[2]

            # r_mid = 0.5 * (p1 + p2)

            # b_mid = (b1 + b2) / np.linalg.norm(b1 + b2)

            # n_raw = (n1 + n2)
            # # Remove component along b_mid
            # n_raw_proj = np.dot(n_raw, b_mid) * b_mid
            # n_ortho = n_raw - n_raw_proj
            # n_norm = np.linalg.norm(n_ortho)
            # n_mid = n_ortho / n_norm if n_norm > 0 else n_raw

            # return np.vstack([r_mid, b_mid, n_mid])
            return (n1 + n2) / 2

        coords = Coords.compute_helix_from_nucl(
            (0, 0, 0),  # start position
            (1, 0, 0),  # base vector
            (0, 1, 0),  # normal vector
            length=9,  # kl180 are equivalent to a stem of 8 nts
            # 1 extra nt will be used for branched loops
            double=True,
        )

        coords = coords.array
        # dummy ends used for branching
        connect_dummy_start = coords[0]
        kl_dummy_end = coords[-1]
        coords = coords[1:-1]  # remove extra nts used for branching dummies

        coords1 = coords[:8]
        coords2 = coords[8:]
        if second_strand:
            coords1 = coords[8:]
            coords2 = coords[:8]

        first_A1 = coords1[0]
        first_kl1_nt = coords2[1]
        second_A1 = (first_A1 + first_kl1_nt) / 2
        second_A1 = middle_nucleotide_coords(first_A1, first_kl1_nt)

        # coords = Coords(np.vstack([[first_A1], [second_A1], coords2[1:7]]))
        if branched:
            connect_T = np.array(
                [
                    [-0.48806033, 0.85282456, 0.18592603, -0.84215779],
                    [-0.71462651, -0.51266626, 0.4758265, -0.42680393],
                    [0.50110485, 0.09931554, 0.8596631, -0.33252832],
                    [0.0, 0.0, 0.0, 1.0],
                ]
            )
            connect_dummy_end = Coords.apply_transformation(
                connect_T,
                connect_dummy_start[0],
                connect_dummy_start[1],
                connect_dummy_start[2],
                local=True,
            )

            kl_T = np.array(
                [
                    [-0.40776493, -0.74831039, 0.52316322, -0.84935859],
                    [-0.73687618, -0.06857565, -0.67253327, -0.13894851],
                    [0.53920428, -0.65979462, -0.52344936, -0.28125316],
                    [0.0, 0.0, 0.0, 1.0],
                ]
            )
            start_KL = coords2[1]
            kl_dummy_start = Coords.apply_transformation(
                kl_T, start_KL[0], start_KL[1], start_KL[2], local=True
            )
            connect_coords = Coords(dummy_ends=(connect_dummy_start, connect_dummy_end))
            kl_coords = Coords(coords2[1:8], dummy_ends=(kl_dummy_start, kl_dummy_end))

            return kl_coords, connect_coords
        else:
            coords = Coords(np.vstack([[first_A1], [second_A1], coords2[1:8]]))
            return coords

    @staticmethod
    def compute_transformation_matrix(
        p1: np.ndarray,
        b1: np.ndarray,
        n1: np.ndarray,
        p2: np.ndarray,
        b2: np.ndarray,
        n2: np.ndarray,
        local: bool = False,
    ) -> np.ndarray:
        """
        Compute a transformation matrix to align one coordinate frame to another.

        Parameters
        ----------
        p1 : np.ndarray
            First position vector.
        b1 : np.ndarray
            First base vector.
        n1 : np.ndarray
            First normal vector.
        p2 : np.ndarray
            Second position vector.
        b2 : np.ndarray
            Second base vector.
        n2 : np.ndarray
            Second normal vector.
        local : bool, default False
            Whether to compute the transformation in a local reference system.

        Returns
        -------
        np.ndarray
            The computed 4x4 transformation matrix.
        """
        # Calculate the third orthogonal vectors
        if local:
            p2, b2, n2 = Coords.set_reference(p1, b1, n1, p2, b2, n2)
            p1, b1, n1 = np.array((0, 0, 0)), np.array((1, 0, 0)), np.array((0, 1, 0))

        v1 = np.cross(b1, n1)
        v2 = np.cross(b2, n2)

        # Construct the rotation matrices
        R1 = np.column_stack((b1, v1, n1))
        R2 = np.column_stack((b2, v2, n2))
        # Calculate the rotation matrix R
        rotation_matrix = R2 @ R1.T
        rotation = R.from_matrix(rotation_matrix)

        # Calculate the translation vector t
        t = p2 - rotation.apply(p1)

        # Construct the 4x4 homogeneous transformation matrix
        T = np.eye(4)
        T[:3, :3] = rotation_matrix
        T[:3, 3] = t
        return T

    @staticmethod
    def set_reference(
        ref_p: np.ndarray,
        ref_b: np.ndarray,
        ref_n: np.ndarray,
        p: np.ndarray,
        b: np.ndarray,
        n: np.ndarray,
        reverse: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Transform coordinates to a new reference system defined by origin, base,
        and normal vectors.

        Parameters
        ----------
        ref_p : np.ndarray
            Reference position vector.
        ref_b : np.ndarray
            Reference base vector.
        ref_n : np.ndarray
            Reference normal vector.
        p : np.ndarray
            Position vector to transform.
        b : np.ndarray
            Base vector to transform.
        n : np.ndarray
            Normal vector to transform.
        reverse : bool, default False
            If True, reverse the transformation.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray, np.ndarray]
            Transformed position, base, and normal vectors.
        """
        # Normalize the reference base and normal vectors
        ref_b /= np.linalg.norm(ref_b)
        ref_n /= np.linalg.norm(ref_n)
        third_axis = np.cross(ref_b, ref_n)

        # Create the rotation matrix to the reference
        rotation_matrix = np.column_stack((ref_b, ref_n, third_axis))

        if reverse:
            rotation = R.from_matrix(rotation_matrix)
            # For reversing, use the transpose of the rotation matrix
            transformed_p = (rotation.apply(p)) + ref_p
        else:
            rotation_matrix = rotation_matrix.T
            rotation = R.from_matrix(rotation_matrix)
            # Apply the transformation to the new reference frame
            transformed_p = rotation.apply(p - ref_p)

        # transform the base and normal vectors
        transformed_b = rotation.apply(b)
        transformed_n = rotation.apply(n)

        return transformed_p, transformed_b, transformed_n

    @staticmethod
    def transform_array(array: np.ndarray, T_matrix: np.ndarray) -> None:
        """
        Apply a transformation matrix to the coordinates inplace.

        Parameters
        ----------
        array : np.ndarray
            A 3D array of coordinates.
        T_matrix : np.ndarray
            A 4x4 transformation matrix.
        """
        if array.size == 0:
            return

        positions = array[:, 0]
        base_vectors = array[:, 1]
        normal_vectors = array[:, 2]

        positions_homogeneous = np.hstack((positions, np.ones((positions.shape[0], 1))))
        transformed_positions = (T_matrix @ positions_homogeneous.T).T[:, :3]

        rotation_matrix = T_matrix[:3, :3]
        rotation = R.from_matrix(rotation_matrix)

        transformed_base_vectors = rotation.apply(base_vectors)
        transformed_normal_vectors = rotation.apply(normal_vectors)

        array[:, 0] = transformed_positions
        array[:, 1] = transformed_base_vectors
        array[:, 2] = transformed_normal_vectors

    ###
    ### PUBLIC METHODS
    ###

    def add_protein(self, protein: "ProteinCoords") -> None:
        """
        Add a ProteinCoords instance to the proteins list.

        Parameters
        ----------
        protein : ProteinCoords
            A ProteinCoords instance to be added.
        """
        if not isinstance(protein, ProteinCoords):
            raise ValueError("The protein argument must be a ProteinCoords instance")
        self._proteins.append(protein)

    def copy(self) -> "Coords":
        """
        Return a copy of the Coords object.

        Returns
        -------
        Coords
            A new Coords object with the same attributes.
        """
        return Coords(
            self.array,
            (self._dummy_ends[0], self._dummy_ends[1]),
            [protein.copy() for protein in self._proteins],
        )

    def is_empty(self) -> bool:
        """
        Check if the Coords object is empty.

        Returns
        -------
        bool
            True if the object is empty, otherwise False.
        """
        return (
            self.size == 0
            and self._dummy_ends[0].size == 0
            and self._dummy_ends[1].size == 0
        )

    def reverse(self) -> "Coords":
        """
        Return a new Coords object with the order of the coordinates reversed.

        Returns
        -------
        Coords
            A new Coords object with reversed coordinates.
        """
        new_dummy = (self._dummy_ends[1], self._dummy_ends[0])
        new_array = np.flip(self._array, axis=0)
        new_proteins = [protein.copy() for protein in self._proteins]
        return Coords(new_array, new_dummy, new_proteins)

    def reverse_in_place(self) -> None:
        """
        Reverse the order of the coordinates in place.
        """
        self._array = np.flip(self._array, axis=0)
        self._dummy_ends = (self._dummy_ends[1], self._dummy_ends[0])

    def to_json(self) -> Dict[str, Any]:
        """
        Serialize the Coords object to a JSON-compatible dictionary.

        Returns
        -------
        Dict[str, Any]
            A dictionary representation of the Coords object.
        """
        return {
            "array": self._array.tolist(),
            "dummy_ends": (
                self._dummy_ends[0].tolist(),
                self._dummy_ends[1].tolist(),
            ),
            "proteins": [
                {"sequence": p.sequence, "coords": p.coords.tolist()}
                for p in self._proteins
            ],
        }

    def transform(self, T: np.ndarray) -> "Coords":
        """
        Apply the transformation matrix T to the coordinates.

        Parameters
        ----------
        T : np.ndarray
            A 4x4 transformation matrix.

        Returns
        -------
        Coords
            The transformed Coords object.
        """
        # Separate position, base vector, and normal vector
        if self._array.size > 0:
            self.transform_array(self._array, T)

        if self._dummy_ends[0].size > 0:
            new_pos, new_bv, new_nv = Coords.apply_transformation(
                T,
                self._dummy_ends[0][0],
                self._dummy_ends[0][1],
                self._dummy_ends[0][2],
            )
            self._dummy_ends[0][0] = new_pos
            self._dummy_ends[0][1] = new_bv
            self._dummy_ends[0][2] = new_nv
        if self._dummy_ends[1].size > 0:
            new_pos, new_bv, new_nv = Coords.apply_transformation(
                T,
                self._dummy_ends[1][0],
                self._dummy_ends[1][1],
                self._dummy_ends[1][2],
            )
            self._dummy_ends[1][0] = new_pos
            self._dummy_ends[1][1] = new_bv
            self._dummy_ends[1][2] = new_nv
        for protein in self._proteins:
            protein.transform(T)
        return self

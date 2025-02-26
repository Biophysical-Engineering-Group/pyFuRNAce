import numpy as np
from scipy.spatial.transform import Rotation as R
from typing import Any, Tuple, List, Union

try:
    from oxDNA_analysis_tools.PDB_oxDNA import PDB_oxDNA
    from oxDNA_analysis_tools.UTILS.RyeReader import conf_to_str, get_top_string
    oat_installed = True
except ImportError:
    oat_installed = False

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
    
    def __init__(self, sequence: str = '', coords: np.ndarray = np.array(())) -> None:
        """
        Initialize the ProteinCoords instance with a sequence and coordinates.
        
        Parameters
        ----------
        sequence : str, optional
            The amino acid sequence of the protein (default is an empty string).
        coords : np.ndarray, optional
            A NumPy array of shape (N, 3) representing the 3D coordinates (default is an empty array).
        """
        self._sequence: str = ''
        self._coords: np.ndarray = np.array(())
        self.sequence = sequence
        self.coords = coords
    
    def __str__(self) -> str:
        """Return a string representation of the ProteinCoords instance."""
        return f'ProteinCoords({self._sequence})'
    
    def __repr__(self) -> str:
        """Return a string representation of the ProteinCoords instance."""
        return f'ProteinCoords({self._sequence})'
    
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
        """Get the 3D coordinates."""
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
            If the coordinates are not a valid NumPy array or do not match the sequence length.
        """
        if isinstance(coords, (list, tuple)):
            coords = np.array(coords)
        if not isinstance(coords, np.ndarray):
            raise ValueError("The coords must be a numpy array")
        if self._sequence and self._coords.size > 0 and len(self.sequence) != len(coords):
            raise ValueError(f"The coords must have the same length as the sequence: expected len {len(self.sequence)}, got {len(coords)} coordinates")
        self._coords = np.array(coords)
    
    @property
    def sequence(self) -> str:
        """Get the protein sequence."""
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
        if self._sequence and self._coords.size > 0 and len(sequence) != len(self.coords):
            raise ValueError("The sequence must have the same length as the coords")
        self._sequence = sequence

    ###
    ### PUBLIC METHODS
    ###
    
    def copy(self) -> 'ProteinCoords':
        """Create a copy of the ProteinCoords instance."""
        return ProteinCoords(self.sequence, np.copy(self.coords))
    
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
    A class to handle 3D coordinate transformations for molecular structures.
    """
    def __init__(self, 
                 input_array: Union[np.ndarray, List] = np.array(()), 
                 dummy_ends: Tuple[np.ndarray, np.ndarray] = (np.array(()), np.array(())), 
                 proteins: List['ProteinCoords'] = None):
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
        self.dummy_ends = dummy_ends
        self.proteins = proteins
        self.array = np.array(input_array)

    def __getitem__(self, key: int) -> np.ndarray:
        """Get item from array."""
        return self.array[key]
    
    def __setitem__(self, key: int, value: np.ndarray) -> None:
        """Set item in array."""
        self.array[key] = value

    def __str__(self) -> str:
        """Return a string representation of the coordinates."""
        return f'Coords({self.array.tolist()})'
    
    def __len__(self) -> int:
        """Return the number of elements in the coordinates array."""
        return len(self.array)
    
    ###
    ### PROPERTIES
    ###

    @property
    def array(self) -> np.ndarray:
        """Return the coordinates array."""
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
            raise ValueError("The Coordinates array must be a numpy array or a list of coordinates")
        self._array = np.array(new_array)

    @property
    def dummy_ends(self) -> Tuple[np.ndarray, np.ndarray]:
        """Return the dummy ends."""
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
            raise ValueError("The dummy_ends argument must be a tuple of two numpy arrays.")
        
        self._dummy_ends = (np.array(new_dummy[0]), np.array(new_dummy[1]))

    @property
    def proteins(self) -> List['ProteinCoords']:
        """Return the list of proteins."""
        return self._proteins
    
    @proteins.setter
    def proteins(self, proteins: List['ProteinCoords']) -> None:
        """
        Set the proteins list.
        
        Parameters
        ----------
        proteins : List[ProteinCoords]
            A list of ProteinCoords instances.
        """
        if not isinstance(proteins, list) or any(not isinstance(p, ProteinCoords) for p in proteins):
            raise ValueError("The proteins argument must be a list of ProteinCoords instances")
        self._proteins = proteins
    
    @property
    def shape(self) -> Tuple[int, ...]:
        """Return the shape of the coordinates array."""
        return self.array.shape
    
    @property
    def size(self) -> int:
        """Return the size of the coordinates array."""
        return self.array.size
    
    ###
    ### STATIC METHODS
    ###

    @staticmethod
    def apply_transformation(T: np.ndarray, 
                             p: np.ndarray, 
                             b: np.ndarray, 
                             n: np.ndarray, 
                             local: bool = False
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
        local : bool, optional
            Whether to apply transformation in a local reference frame, by default False.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray, np.ndarray]
            Transformed position, base, and normal vectors.
        """
        if local:
            p_local, b_local, n_local = p, b, n
            p, b, n = np.array((0,0,0)), np.array((1,0,0)), np.array((0,1,0))

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
            p_trans, b_trans, n_trans = Coords.set_reference(p_local, 
                                                             b_local, 
                                                             n_local, 
                                                             p_trans, 
                                                             b_trans, 
                                                             n_trans, 
                                                             reverse=True)
                                                 
        return p_trans, b_trans, n_trans
    
    @staticmethod
    def combine_coords(strand1: "Strand", 
                       coords1: "Coords", 
                       strand2: "Strand", 
                       coords2: "Coords"
                       ) -> "Coords":
        """
        Combine the coordinates of two strands by applying necessary transformations.
        Use dependency injection (with strand1 and strand2 as arguments) to apply 
        the transformation to a whole Strand Block. The second strand is transformed
        to match the first strand's orientation.

        Parameters
        ----------
        strand1 : Strand
            First strand.
        coords1 : Coords
            Coordinates of the first strand.
        strand2 : Strand
            Second strand.
        coords2 : Coords
            Coordinates of the second strand.

        Returns
        -------
        Coords
            Combined coordinate object containing both transformed strands.
        """
        seq1 = strand1.sequence
        seq2 = strand2.sequence

        ### HEDGE CASES ###
        # there are no coordinates to combine 
        if coords1.is_empty() and not seq1 and coords2.is_empty() and not seq2: 
            return Coords([])
        # the first strand has no coordinates and no nucleotides
        elif coords1.is_empty() and not seq1: 
            return coords2.copy()
        # the second strand has no coordinates and no nucleotides
        elif coords2.is_empty() and not seq2: 
            return coords1.copy()
        
        ### TRANSFORM THE SECOND STRAND ###
        # transform them only if they are not in the same block os strands
        if id(strand1.strands_block) != id(strand2.strands_block): 
            # create the coordinates for the first strand
            if coords1.size == 0 and len(seq1): 
                coords1 = Coords.compute_helix_from_nucl((0,0,0), # start position
                                                         (1,0,0), # base vector
                                                         (0,1,0), # normal vector
                                                         length=len(seq1),
                                                         directionality=seq1.directionality)
            if coords2.size == 0 and len(seq2): # create the coordinates for the second strand
                coords2 = Coords.compute_helix_from_nucl((0,0,0), # start position
                                                         (1,0,0), # base vector
                                                         (0,1,0), # normal vector
                                                         length=len(seq2),
                                                         directionality=seq2.directionality)
                
            ### CHECK THE JOINING POINTS AND TRANSFORM THE ARRAY TWO ###

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
                else: # the second strand has no dummy: calulcate an helical join from the last nucleotide of the first strand
                    helix_coords = Coords.compute_helix_from_nucl(coords1[-1][0], 
                                                                  coords1[-1][1], 
                                                                  coords1[-1][2], 
                                                                  length=1, 
                                                                  directionality=seq1.directionality)
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
            T_matrix = Coords.compute_transformation_matrix(pos1, 
                                                            bv1, 
                                                            nv1, 
                                                            pos2, 
                                                            bv2, 
                                                            nv2)

            ### Apply the transformation to the second strand motif arrays ###
            strand2.strands_block.transform(T_matrix)

        ### COMBINE THE COORDINATES, the coordinates could be empty if there are dummies ### 
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
            combined = Coords(np.concatenate((coords1.array, coords2.array), axis=0), 
                                              dummy_ends=combined_dummy)

        ### COMBINE THE PROTEINS ###
        combined.proteins = coords1.proteins + coords2.proteins

        return combined

    ###
    ### PUBLIC METHODS
    ###

    def copy(self) -> 'Coords':
        """Return a copy of the Coords object."""
        return Coords(np.copy(self.array), 
                      (np.copy(self._dummy_ends[0]), np.copy(self._dummy_ends[1])), 
                      [protein.copy() for protein in self._proteins])

    def add_protein(self, protein: 'ProteinCoords') -> None:
        """Add a ProteinCoords instance to the proteins list."""
        if not isinstance(protein, ProteinCoords):
            raise ValueError("The protein argument must be a ProteinCoords instance")
        self._proteins.append(protein)

    def transform(self, T: np.ndarray) -> 'Coords':
        """Apply the transformation matrix T to the coordinates"""
        # Separate position, base vector, and normal vector
        if self._array.size > 0:
            self.transform_array(self._array, T)

        if self._dummy_ends[0].size > 0:
            new_pos, new_bv, new_nv = Coords.apply_transformation(T, self._dummy_ends[0][0], self._dummy_ends[0][1], self._dummy_ends[0][2])
            self._dummy_ends[0][0] = new_pos
            self._dummy_ends[0][1] = new_bv
            self._dummy_ends[0][2] = new_nv
        if self._dummy_ends[1].size > 0:
            new_pos, new_bv, new_nv = Coords.apply_transformation(T, self._dummy_ends[1][0], self._dummy_ends[1][1], self._dummy_ends[1][2])
            self._dummy_ends[1][0] = new_pos
            self._dummy_ends[1][1] = new_bv
            self._dummy_ends[1][2] = new_nv
        for protein in self._proteins:
            protein.transform(T)
        return self
    
    def reverse_in_place(self):
        """Reverse the order of the coordinates in place"""
        self._array = np.flip(self._array, axis=0)
        self._dummy_ends = (self._dummy_ends[1], self._dummy_ends[0])
    
    def is_empty(self) -> bool:
        """Check if the Coords object is empty."""
        return self.size == 0 and self._dummy_ends[0].size == 0 and self._dummy_ends[1].size == 0

    
    ###
    ### FUNCTIONS FOR MANIPULATING 3D COORDINATES
    ###

    @staticmethod
    def set_reference(ref_p, ref_b, ref_n, p, b, n, reverse=False):
        """
        Transform coordinates to a new reference system defined by origin, base, and normal vectors.
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
    def compute_transformation_matrix(p1, b1, n1, p2, b2, n2, local=False):
        # Calculate the third orthogonal vectors
        if local:
            p2, b2, n2 = Coords.set_reference(p1, b1, n1, p2, b2, n2)
            p1, b1, n1 = np.array((0,0,0)), np.array((1,0,0)), np.array((0,1,0))
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
    def compute_helix_from_nucl(pos, old_a1, old_a3, length, directionality="53", double=False):
        # Adjust arguments: Convert lists or tuples to numpy arrays for consistency in mathematical operations
        direction_35 = directionality == "35"
        if isinstance(pos, (list, tuple)):
            pos = np.array(pos)
        if isinstance(old_a1, (list, tuple)):
            old_a1 = np.array(old_a1)
        if isinstance(old_a3, (list, tuple)):
            old_a3 = np.array(old_a3)

        # Model constants: Define structural parameters for the helix
        inclination = -15.5 * np.pi / 180  # Inclination angle in radians
        bp_backbone_distance = 2  # Distance between base pairs along the backbone
        diameter = 2.35  # Diameter of the helix
        base_base_distance = 0.3287  # Distance between bases along the helix axis
        rot = 32.73 * np.pi / 180  # Rotation angle between successive bases in radians
        cord = np.cos(inclination) * bp_backbone_distance  # Chord length
        center_to_cord = np.sqrt((diameter / 2) ** 2 - (cord / 2) ** 2)  # Distance from center to chord
        fudge = 0.4  # Fudge factor for position adjustment

        norm_a1_a3 = np.cross(old_a1, old_a3)  # Calculate the axis perpendicular to the old A1 and A3 vectors
        R_dir = R.from_rotvec(norm_a1_a3 * inclination)  # Create the rotation matrix for the direction vector
        dir_vector = - R_dir.apply(old_a3) # Calculate the direction vector, rotating the old normal vector by the inclination angle
        if direction_35:
            dir_vector *= -1  # Adjust direction if specified
        dir_vector /= np.linalg.norm(dir_vector)  # Normalize direction vector

        # Calculate the coordinates if the helix axis is the Z-axis
        x1, y1, z1 = center_to_cord, -cord / 2, -(bp_backbone_distance / 2) * np.sin(inclination)
        x2, y2, z2 = center_to_cord, +cord / 2, +(bp_backbone_distance / 2) * np.sin(inclination)
        r1 = np.array([x1, y1, z1])
        r2 = np.array([x2, y2, z2])

        ### Set the axis to the correct one
        ref_vector = np.array([0, 0, 1])  # Z-axis is the default axis
        cross_prod = np.cross(ref_vector, dir_vector)  # Calculate the cross product for rotation
        dot_prod = np.dot(ref_vector, dir_vector)  # Calculate the dot product
        # Calculate the quaternion components
        scalar = np.sqrt((1.0 + dot_prod) * 2.0)  # Scalar component to normalize the quaternion
        # The quaternion matrix transforms the ref_vector to the dir_vector
        q_matrix = np.append(cross_prod / scalar, scalar * 0.5)
        q_matrix /= np.linalg.norm(q_matrix)  # Normalize the quaternion
        q1 = R.from_quat(q_matrix)  # Create the rotation from the quaternion
        r1 = q1.apply(r1)  # Apply rotation to r1
        r2 = q1.apply(r2)  # Apply rotation to r2

        # Set a1 to the correct orientation
        r1_to_r2 = r2 - r1
        if direction_35:
            r1_to_r2 = - r1_to_r2  # Adjust direction if specified
        r1_to_r2 /= np.linalg.norm(r1_to_r2)  # Normalize the vector

        # Rotate the helix axis to the correct orientation
        r1_to_r2_proj = r1_to_r2 - np.dot(r1_to_r2, dir_vector) * dir_vector  # Project r1_to_r2 onto dir_vector
        old_a1_proj = old_a1 - np.dot(old_a1, dir_vector) * dir_vector  # Project old_a1 onto dir_vector
        rotAngle2 = np.arccos(np.clip(np.dot(r1_to_r2_proj / np.linalg.norm(r1_to_r2_proj), old_a1_proj / np.linalg.norm(old_a1_proj)), -1.0, 1.0))
        if np.dot(np.cross(r1_to_r2, old_a1), dir_vector) < 0:
            rotAngle2 *= -1  # Adjust rotation angle based on the cross product
        q2 = R.from_rotvec(dir_vector * rotAngle2)  # Create the rotation matrix
        r1 = q2.apply(r1)  # Apply rotation to r1
        r2 = q2.apply(r2)  # Apply rotation to r2

        # Center point of the helix axis
        r = r1
        if direction_35:
            r = r2  # Adjust direction if specified
        start_pos = pos - r - old_a1 * fudge  # Calculate the starting position

        # Create per-step rotation matrix
        R_step = R.from_rotvec(dir_vector * rot)  # Create rotation for each step

        # Initialize properties of new nucleotide
        out = []
        out_double = []

        # Generate nucleotide positions and orientations
        for _ in range(length):
            # Calculate rotation around central axis and step along axis
            r1 = R_step.apply(r1) + dir_vector * base_base_distance
            r2 = R_step.apply(r2) + dir_vector * base_base_distance

            # Calculate a1 orientation
            r1_to_r2 = r2 - r1
            a1 = r1_to_r2 / np.linalg.norm(r1_to_r2)  # Normalize a1

            # Calculate A3 orientation
            a1proj = a1 - np.dot(a1, dir_vector) * dir_vector  # Project a1 onto dir_vector
            a1proj /= np.linalg.norm(a1proj)  # Normalize projection
            a3 = -np.cos(inclination) * dir_vector + np.sin(inclination) * a1proj  # Calculate a3
            a3 /= np.linalg.norm(a3)  # Normalize a3

            # Calculate position
            r = r1
            if direction_35:
                r = r2  # Adjust direction if specified
                a1 = -a1
                a3 = -a3
            RNA_fudge = a1 * fudge  # Apply fudge factor
            p = r + RNA_fudge + start_pos  # Calculate position

            out.append([p, a1, a3])  # Append nucleotide properties

            # Double helix case
            if double:
                a1 = -r1_to_r2 / np.linalg.norm(-r1_to_r2)  # Calculate a1 for double helix
                a1proj = a1 - np.dot(a1, -dir_vector) * -dir_vector  # Project a1 onto dir_vector
                a1proj /= np.linalg.norm(a1proj)  # Normalize projection
                a3 = np.cos(inclination) * dir_vector + np.sin(inclination) * a1proj  # Calculate a3
                a3 /= np.linalg.norm(a3)  # Normalize a3
                r = r2
                if direction_35:
                    r = r1  # Adjust direction if specified
                    a1 = -a1
                    a3 = -a3
                RNA_fudge = a1 * fudge  # Apply fudge factor
                p = r + RNA_fudge + start_pos  # Calculate position
                out_double.append([p, a1, a3])  # Append nucleotide properties for double helix

        return Coords(out + out_double[::-1])
    
    @staticmethod
    def load_from_file(filename, dummy_ends=(False, False), extend = (0, 0), topology_file= None, protein=False):
        return Coords.cleanup_from_oxdna_file(filename, dummy_ends, extend, return_coords=True, 
                                              topology_file=topology_file,
                                              protein=protein)
    
    @staticmethod
    def load_from_text(conf_text, dummy_ends=(False, False), extend = (0, 0), return_coords=True, top_text=None, protein=False):
        """Return the coordinates array contained in a configuration text"""
        lines = conf_text.split('\n')

        def extract_coords_from_line(line):
            splitted = line.strip().split()
            return [[float(p) for p in splitted[0:3]], [float(b) for b in splitted[3:6]], [float(n) for n in splitted[6:9]]]

        # Extract the coordinates
        coords = np.array([extract_coords_from_line(line) for line in lines[3:] if line.strip()])

        all_prot_coords = []
        protein_coords = None
        protein_text = ''
        if protein: # There is a protein, detect the line and load the coordinates
            if top_text is None:
                raise ValueError("A topology file path must be provided when loading a structure with proteins")
            protein_text = ',\n\tproteins=['
            lines = top_text.split('\n')[1:]
            seq_start_ind = 0
            for line in lines:
                if not line:
                    continue
                seq = line.split()[0]
                if 'peptide' in line:
                    protein_coords = coords[seq_start_ind: seq_start_ind + len(seq)]
                    protein_text += f'ProteinCoords("{seq}", coords={protein_coords.tolist()}), '
                    if return_coords:
                        all_prot_coords.append(ProteinCoords(seq, protein_coords))
                    coords = np.delete(coords, np.s_[seq_start_ind: seq_start_ind + len(seq)], axis=0)
                    seq_start_ind -= len(seq)
                seq_start_ind += len(seq)
            protein_text += ']'

        if extend[0] > 0:
            extend_first = Coords.compute_helix_from_nucl(coords[0][0], coords[0][1], coords[0][2], extend[0], directionality="35", double=False)[::-1]
            coords = np.concatenate((extend_first, coords))

        if extend[1] > 0:
            extend_last = Coords.compute_helix_from_nucl(coords[-1][0], coords[-1][1], coords[-1][2], extend[1], directionality="53", double=False)
            coords = np.concatenate((coords, extend_last))

        dummy_0 = np.array(())
        dummy_1 = np.array(())
        if dummy_ends[0]:
            dummy_0 = coords[0]
            coords = coords[1:]
        if dummy_ends[1]:
            dummy_1 = coords[-1]
            coords = coords[:-1]
        dummy_text = f',\n\tdummy_ends=({dummy_0.tolist()}, {dummy_1.tolist()})'
        if not dummy_ends[0] and not dummy_ends[1]:
            dummy_text = ''

        if return_coords:
            return Coords(coords, dummy_ends=(dummy_0, dummy_1), proteins=all_prot_coords)

        print(f'Coords({coords.tolist()}{dummy_text}{protein_text})')
        
    @staticmethod
    def transform_array(array, T_matrix: np.ndarray) -> None:
        """
        Apply a transformation matrix to the coordinates inplace.
        
        Parameters
        ----------
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

    @staticmethod
    def cleanup_from_oxdna_file(filename, dummy_ends=(False, False), extend = (0, 0), return_coords=False, topology_file=None, protein=False):
        """Return the coordinates array contained in an oxDNA configuration file"""
        # Check if the file is a pdb file
        pdb = False
        top_text = None
        if type(filename) != str:
            try:
                filename = str(filename)
            except Exception as e:
                raise ValueError("The filename must be a string or a path-like object, got ", type(filename))
        if filename.endswith('.pdb') or filename.endswith('.PDB'):
            pdb = True
        if pdb and not oat_installed:
            raise ValueError("oxDNA_analysis_tools not installed. Please install it to use the pdb option")
        
        # read the file
        with open(filename, 'r') as f:
            conf_text = f.read()

        if pdb:
            confs, systems = PDB_oxDNA(conf_text)
            print(confs, systems)
            conf_text = conf_to_str(confs[0])
            top_text = get_top_string(systems[0])
        elif topology_file:
            with open(topology_file, 'r') as f:
                top_text = f.read()

        return Coords.load_from_text(conf_text, dummy_ends, extend, return_coords, top_text=top_text, protein=protein)

    # @staticmethod
    # def compute_dae_crossover4_good(pos, old_a1, old_a3, length, directionality="53", double=False):
    #     # adjust arguments
    #     direction_35 = directionality == "35"
    #     if isinstance(pos, (list, tuple)):
    #         pos = np.array(pos)
    #     if isinstance(old_a1, (list, tuple)):
    #         old_a1 = np.array(old_a1)
    #     if isinstance(old_a3, (list, tuple)):
    #         old_a3 = np.array(old_a3)
        
    #     next_nulc = Coords.compute_helix_from_nucl(pos, old_a1, old_a3, 1, directionality, False)[0]
    #     pos = next_nulc[0]
    #     old_a1 = next_nulc[1]
    #     old_a3 = next_nulc[2]

    #     # Model constants
    #     inclination = -15.5 * np.pi / 180
    #     bp_backbone_distance = 2
    #     diameter = 2.35
    #     base_base_distance = 0.3287
    #     rot = 32.7 * np.pi / 180
    #     cord = np.cos(inclination) * bp_backbone_distance
    #     center_to_cord = np.sqrt((diameter / 2) ** 2 - (cord / 2) ** 2)
    #     fudge = 0.4

    #     # Define the helix axis based on A1 and A3
    #     dir_vector = (old_a3 - old_a1 * np.sin(inclination)) / -(np.cos(inclination) - np.sin(inclination) ** 2)
    #     if direction_35:
    #         dir_vector *= -1
    #     dir_vector /= np.linalg.norm(dir_vector) # normalize

    #     # Calculate the coord if the helix axis is the Z-axis
    #     x1, y1, z1 = center_to_cord, -cord / 2, -(bp_backbone_distance / 2) * np.sin(inclination)
    #     x2, y2, z2 = center_to_cord, +cord / 2, +(bp_backbone_distance / 2) * np.sin(inclination)
    #     r1 = np.array([x1, y1, z1])
    #     r2 = np.array([x2, y2, z2])

    #     ### Set the axis to the correct one
    #     ref_vector = np.array([0, 0, 1]) # Z-axis is the default axis
    #     cross_prod = np.cross(ref_vector, dir_vector) # Calculate the cross product
    #     dot_prod = np.dot(ref_vector, dir_vector) # Calculate the dot product
    #     # Calculate the quaternion components
    #     scalar = np.sqrt((1.0 + dot_prod) * 2.0) # scalar component to normalize the quaternion
    #     # The quaternion matrix transforms the ref_vector to the dir_vector
    #     q_matrix = np.append(cross_prod / scalar, scalar * 0.5)
    #     # Normalize the quaternion
    #     q_matrix /= np.linalg.norm(q_matrix)
    #     # Apply the quaternion to the r1 and r2 vectors
    #     q1 = R.from_quat(q_matrix)
    #     r1 = q1.apply(r1)
    #     r2 = q1.apply(r2)

    #     # Set a1 to the correct orientation
    #     r1_to_r2 = r2 - r1
    #     if direction_35:
    #         r1_to_r2 = - r1_to_r2
    #     r1_to_r2 /= np.linalg.norm(r1_to_r2)

    #     # project the r1_to_r2 vector on the helix axis plane
    #     r1_to_r2_proj = r1_to_r2 - np.dot(r1_to_r2, dir_vector) * dir_vector
    #     # project the base vector on the helix axis plane
    #     old_a1_proj = old_a1 - np.dot(old_a1, dir_vector) * dir_vector
    #     # calculate the rotation angle to orient the helix axis
    #     rotAngle2 = np.arccos(np.clip(np.dot(r1_to_r2_proj / np.linalg.norm(r1_to_r2_proj), old_a1_proj / np.linalg.norm(old_a1_proj)), -1.0, 1.0))
    #     if np.dot(np.cross(r1_to_r2, old_a1), dir_vector) < 0:
    #         rotAngle2 *= -1
    #     # This rotoation matrix rotates the radius 1 and 2 to the point where the base vector is aligned with the old_a1
    #     q2 = R.from_rotvec(dir_vector * rotAngle2)
    #     r1 = q2.apply(r1)
    #     r2 = q2.apply(r2)

    #     # Center point of the helix axis
    #     r = r1
    #     if direction_35:
    #         r = r2
    #     # Calculate the start position
    #     start_pos = pos - r - old_a1 * fudge

    #     # Create per-step rotation matrix
    #     R_step = R.from_rotvec(dir_vector * rot)

    #     # Initialize properties of new nucleotide
    #     out = []
    #     out_double = []

    #     # Generate nucleotide positions and orientations
    #     for i in range(length):
    #         # Calculate rotation around central axis and step along axis
    #         if i > 0:
    #             r1 = R_step.apply(r1) + dir_vector * base_base_distance
    #             r2 = R_step.apply(r2) + dir_vector * base_base_distance

    #         # Calculate a1 orientation
    #         r1_to_r2 = r2 - r1
    #         if i == 0:
    #             r1_to_r2 = - r1_to_r2
    #         a1 = r1_to_r2 / np.linalg.norm(r1_to_r2)

    #         # Calculate A3 orientation
    #         dir_vector = - dir_vector # invert the direction of the helix axis
    #         dir_vector /= np.linalg.norm(dir_vector) # normalize
    #         a1proj = a1 - np.dot(a1, dir_vector) * dir_vector
    #         a1proj /= np.linalg.norm(a1proj)
    #         a3 = -np.cos(inclination) * dir_vector + np.sin(inclination) * a1proj
    #         a3 /= np.linalg.norm(a3)

    #         # Calculate position
    #         r = r1
    #         if direction_35:
    #             r = r2
    #             a1 = -a1
    #             a3 = -a3
    #         RNA_fudge = a1 * fudge
    #         p = r + start_pos + RNA_fudge * 2

    #         out.append([p, a1, a3])

    #         # Double helix case
    #         if double:
    #             a1 = -r1_to_r2 / np.linalg.norm(r1_to_r2)
    #             a1proj = a1 - np.dot(a1, -dir_vector) * -dir_vector
    #             a1proj /= np.linalg.norm(a1proj)
    #             a3 = np.cos(inclination) * dir_vector + np.sin(inclination) * a1proj
    #             a3 /= np.linalg.norm(a3)
    #             r = r2
    #             if direction_35:
    #                 r = r1
    #                 a1 = -a1
    #                 a3 = -a3
    #             RNA_fudge = a1 * fudge
    #             p = r + RNA_fudge + start_pos
    #             out_double.append([p, a1, a3])

    #     return Coords(out + out_double[::-1])
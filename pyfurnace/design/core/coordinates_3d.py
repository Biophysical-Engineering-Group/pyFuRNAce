try:
    from oxDNA_analysis_tools.PDB_oxDNA import PDB_oxDNA
    from oxDNA_analysis_tools.UTILS.RyeReader import conf_to_str, get_top_string
    oat_installed = True
except ImportError:
    oat_installed = False
import numpy as np
from scipy.spatial.transform import Rotation as R

class ProteinCoords():
    def __init__(self, sequence: str='', coords: np.ndarray = np.array(())):
        # initialize the sequence and coordinates
        self._sequence = None
        self._coords = np.array(())
        # assign the values
        self.sequence = sequence
        self.coords = coords

    def __str__(self):
        return f'ProteinCoords({self._sequence})'
    
    def __repr__(self):
        return f'ProteinCoords({self._sequence})'
    
    def __getitem__(self, key):
        return self.coords[key]
    
    def __setitem__(self, key, value):
        self.coords[key] = value

    def __len__(self):
        return len(self.sequence)
    
    # def __deepcopy__(self, memo=None):
    #     new_obj = ProteinCoords(self.sequence, self.coords.copy())
    #     return new_obj

    @property
    def sequence(self):
        return self._sequence
    
    @sequence.setter
    def sequence(self, sequence):
        if not isinstance(sequence, str):
            raise ValueError("The sequence must be a string")
        if self._sequence and self._coords.size > 0 and len(sequence) != len(self.coords):
            raise ValueError("The sequence must have the same length as the coords")
        self._sequence = sequence

    @property
    def coords(self):
        return self._coords
    
    @coords.setter
    def coords(self, coords):
        if isinstance(coords, (list, tuple)):
            coords = np.array(coords)
        if not isinstance(coords, np.ndarray):
            raise ValueError("The coords must be a numpy array")
        if self._sequence and self._coords.size > 0 and len(self.sequence) != len(coords):
            raise ValueError("The coords must have the same length as the sequence: expected len {len(self.sequence)}, got {len(coords)} coordinates")
        self._coords = np.array(coords)

    def copy(self):
        return ProteinCoords(self.sequence, np.copy(self.coords))

    def transform(self, T):
        # Separate position, base vector, and normal vector
        if self._coords.size == 0:
            return self
        
        positions = self._coords[:, 0]
        base_vectors = self._coords[:, 1]
        normal_vectors = self._coords[:, 2]

        # Apply the transformation to all positions
        positions_homogeneous = np.hstack((positions, np.ones((positions.shape[0], 1))))
        transformed_positions = (T @ positions_homogeneous.T).T[:, :3]

        # Extract the rotation matrix
        rotation_matrix = T[:3, :3]

        # Apply the rotation to all base and normal vectors
        transformed_base_vectors = base_vectors @ rotation_matrix.T
        transformed_normal_vectors = normal_vectors @ rotation_matrix.T

        # Update the coordinates in place
        self._coords[:, 0] = transformed_positions
        self._coords[:, 1] = transformed_base_vectors
        self._coords[:, 2] = transformed_normal_vectors

        return self

class Coords():
    def __init__(self, input_array=np.array(()), dummy_ends=(np.array(()), np.array(())), proteins: list = None):
        if proteins is None:
            proteins = []
        self.dummy_ends = dummy_ends
        self.proteins = proteins
        self.array = np.array(input_array)

    def __getitem__(self, key):
        return self.array[key]
    
    def __setitem__(self, key, value):
        self.array[key] = value

    def __str__(self):
        return f'Coords({self.array.tolist()})'
    
    def __len__(self):
        return len(self.array)
    
    @property
    def array(self):
        return self._array
    
    @array.setter
    def array(self, new_array):
        if not isinstance(new_array, (np.ndarray, list, tuple)):
            raise ValueError("The Coordinates array must be a numpy array or a list of coordinates")
        new_array = np.array(new_array)
        shape = new_array.shape
        # if new_array.size > 0 and shape and (len(shape) != 3 or shape[1] != 3 or shape[2] != 3):
            # print(new_array)
            # raise ValueError(f"The Coordinates array must be a 1D array, where each element contains 3 elements for the position, base vector, and normal vector. Got shape {shape}")
        self._array = new_array

    @property
    def size(self):
        return self.array.size
    
    @property
    def shape(self):
        return self.array.shape
    
    @property
    def dummy_ends(self):
        return self._dummy_ends
    
    @dummy_ends.setter
    def dummy_ends(self, new_dummy):
        if not isinstance(new_dummy, (list, tuple, np.ndarray)) or len(new_dummy) != 2:
            raise ValueError(f"The dummy_ends argument must be a list or array of two sets of dummy coordinates.")
        if isinstance(new_dummy, tuple):
            new_dummy = list(new_dummy)
        if not isinstance(new_dummy[0], np.ndarray):
            new_dummy[0] = np.array(new_dummy[0])
        if not isinstance(new_dummy[1], np.ndarray):
            new_dummy[1] = np.array(new_dummy[1])
        self._dummy_ends = new_dummy

    @property
    def proteins(self):
        return self._proteins
    
    @proteins.setter
    def proteins(self, proteins):
        if not isinstance(proteins, list) or proteins and not all([isinstance(protein, ProteinCoords) for protein in proteins]):
            raise ValueError("The proteins argument must be a list of ProteinCoords instances")
        self._proteins = proteins

    def copy(self):
        return Coords(np.copy(self.array), (np.copy(self._dummy_ends[0]), np.copy(self._dummy_ends[1])), [protein.copy() for protein in self._proteins])

    def add_protein(self, protein):
        if not isinstance(protein, ProteinCoords):
            raise ValueError("The protein argument must be a ProteinCoords instance")
        self._proteins.append(protein)

    def transform(self, T):
        """Apply the transformation matrix T to the coordinates"""
        # Separate position, base vector, and normal vector
        if self._array.size > 0:
            positions = self._array[:, 0]
            base_vectors = self._array[:, 1]
            normal_vectors = self._array[:, 2]

            # Apply the transformation to all positions
            positions_homogeneous = np.hstack((positions, np.ones((positions.shape[0], 1))))
            transformed_positions = (T @ positions_homogeneous.T).T[:, :3]

            # Extract the rotation matrix
            R = T[:3, :3]

            # Apply the rotation to all base and normal vectors
            transformed_base_vectors = base_vectors @ R.T
            transformed_normal_vectors = normal_vectors @ R.T

            # Update the coordinates in place
            self._array[:, 0] = transformed_positions
            self._array[:, 1] = transformed_base_vectors
            self._array[:, 2] = transformed_normal_vectors

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

    def is_empty(self):
        return self.size == 0 and self._dummy_ends[0].size == 0 and self._dummy_ends[1].size == 0
    
    @staticmethod # use Dependency injection to apply the transformation to a whole Strand Block
    def combine_coords(strand1, coords1, strand2, coords2):
        """ Add the coordinates of coords2 to coords1 and return the new coordinates"""
        seq1 = strand1.sequence
        seq2 = strand2.sequence

        ### LIMIT CASES ###
        if coords1.is_empty() and not seq1 and coords2.is_empty() and not seq2: # there are no coordinates to combine 
            return Coords([])
        elif coords1.is_empty() and not seq1: # the first strand has no coordinates and no nucleotides
            return coords2.copy()
        elif coords2.is_empty() and not seq2: # the second strand has no coordinates and no nucleotides
            return coords1.copy()
        ### TRANSFORM THE SECOND STRAND ###
        if id(strand1.strands_block) != id(strand2.strands_block): # transform them only if they are not in the same block os strands
            if coords1.size == 0 and len(seq1): # create the coordinates for the first strand
                coords1 = Coords.compute_helix_from_nucl((0,0,0), # start position
                                                            (1,0,0), # base vector
                                                            (0,1,0), # normal vector
                                                            length= len(seq1),
                                                            directionality = seq1.directionality)
            if coords2.size == 0 and len(seq2): # create the coordinates for the second strand
                coords2 = Coords.compute_helix_from_nucl((0,0,0), # start position
                                                            (1,0,0), # base vector
                                                            (0,1,0), # normal vector
                                                            length= len(seq2),
                                                            directionality = seq2.directionality)
                
            ### CHECK THE JOINING POINTS AND TRANSFORM THE ARRAY TWO ###
            if coords1.dummy_ends[1].size == 0: # the first strand has no end dummy
                if coords2.dummy_ends[0].size > 0: # the second strand has a start dummy
                    # the dummy start translates to...
                    pos1 = coords2.dummy_ends[0][0] 
                    bv1 = coords2.dummy_ends[0][1]
                    nv1 = coords2.dummy_ends[0][2]
                    # ...the last position of the first strand
                    pos2 = coords1[-1][0] 
                    bv2 = coords1[-1][1]
                    nv2 = coords1[-1][2]
                else: # the second strand has no dummy: calulcate an helical join from the last nucleotide of the first strand
                    helix_coords = Coords.compute_helix_from_nucl(coords1[-1][0], coords1[-1][1], coords1[-1][2], length=1, directionality=seq1.directionality)
                    # the first position of the second strand translates to...
                    pos1 = coords2[0][0] 
                    bv1 = coords2[0][1]
                    nv1 = coords2[0][2]
                    # ...the first nucleotide of the helical join
                    pos2 = helix_coords[0][0]
                    bv2 = helix_coords[0][1]
                    nv2 = helix_coords[0][2]
            else: # the first strand has an end dummy, use it
                if coords2.size == 0: # the array 2 has no coordinate, use the dummy instead
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
            T_matrix = Coords.compute_transformation_matrix(pos1, bv1, nv1, pos2, bv2, nv2)

            ### Apply the transformation to the second strand motif arrays ###
            strand2.strands_block.transform(T_matrix)

        ### COMBINE THE COORDINATES, the coordinates could be empty if there are dummies ### 
        combined_dummy = (coords1.dummy_ends[0], coords2.dummy_ends[1])
        if coords1.size == 0 and coords2.size == 0: # both strands have no coordinates
            combined = Coords([], dummy_ends=combined_dummy)
        elif coords1.size == 0: # the first strand has no coordinates
            combined = Coords(coords2.array, dummy_ends=combined_dummy)
        elif coords2.size == 0: # the second strand has no coordinates
            combined = Coords(coords1.array, dummy_ends=combined_dummy)
        else:
            combined = Coords(np.concatenate((coords1.array, coords2.array), axis=0), dummy_ends=combined_dummy)

        ### COMBINE THE PROTEINS ###
        combined.proteins = coords1.proteins + coords2.proteins

        return combined
    
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
    def apply_transformation(T, p, b, n, local=False):
        if local:
            p_local, b_local, n_local = p, b, n
            p, b, n = np.array((0,0,0)), np.array((1,0,0)), np.array((0,1,0))
        # Convert position to homogeneous coordinates
        p_homogeneous = np.append(p, 1)
        # Apply the transformation to the position
        p_transformed_homogeneous = T @ p_homogeneous
        p_transformed = p_transformed_homogeneous[:3]
        # Extract the rotation matrix
        rotation_matrix = T[:3, :3]
        rotation = R.from_matrix(rotation_matrix)
        # Apply the rotation to the base and normal vectors
        b_transformed = rotation.apply(b)
        n_transformed = rotation.apply(n)
        if local:
            # # find the original reference system from the local one
            p_transformed, b_transformed, n_transformed = Coords.set_reference(p_local, b_local, n_local, p_transformed, b_transformed, n_transformed, reverse=True)
        return p_transformed, b_transformed, n_transformed

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
    def compute_dae_crossover4_good(pos, old_a1, old_a3, length, directionality="53", double=False):
        # adjust arguments
        direction_35 = directionality == "35"
        if isinstance(pos, (list, tuple)):
            pos = np.array(pos)
        if isinstance(old_a1, (list, tuple)):
            old_a1 = np.array(old_a1)
        if isinstance(old_a3, (list, tuple)):
            old_a3 = np.array(old_a3)
        
        next_nulc = Coords.compute_helix_from_nucl(pos, old_a1, old_a3, 1, directionality, False)[0]
        pos = next_nulc[0]
        old_a1 = next_nulc[1]
        old_a3 = next_nulc[2]

        # Model constants
        inclination = -15.5 * np.pi / 180
        bp_backbone_distance = 2
        diameter = 2.35
        base_base_distance = 0.3287
        rot = 32.7 * np.pi / 180
        cord = np.cos(inclination) * bp_backbone_distance
        center_to_cord = np.sqrt((diameter / 2) ** 2 - (cord / 2) ** 2)
        fudge = 0.4

        # Define the helix axis based on A1 and A3
        dir_vector = (old_a3 - old_a1 * np.sin(inclination)) / -(np.cos(inclination) - np.sin(inclination) ** 2)
        if direction_35:
            dir_vector *= -1
        dir_vector /= np.linalg.norm(dir_vector) # normalize

        # Calculate the coord if the helix axis is the Z-axis
        x1, y1, z1 = center_to_cord, -cord / 2, -(bp_backbone_distance / 2) * np.sin(inclination)
        x2, y2, z2 = center_to_cord, +cord / 2, +(bp_backbone_distance / 2) * np.sin(inclination)
        r1 = np.array([x1, y1, z1])
        r2 = np.array([x2, y2, z2])

        ### Set the axis to the correct one
        ref_vector = np.array([0, 0, 1]) # Z-axis is the default axis
        cross_prod = np.cross(ref_vector, dir_vector) # Calculate the cross product
        dot_prod = np.dot(ref_vector, dir_vector) # Calculate the dot product
        # Calculate the quaternion components
        scalar = np.sqrt((1.0 + dot_prod) * 2.0) # scalar component to normalize the quaternion
        # The quaternion matrix transforms the ref_vector to the dir_vector
        q_matrix = np.append(cross_prod / scalar, scalar * 0.5)
        # Normalize the quaternion
        q_matrix /= np.linalg.norm(q_matrix)
        # Apply the quaternion to the r1 and r2 vectors
        q1 = R.from_quat(q_matrix)
        r1 = q1.apply(r1)
        r2 = q1.apply(r2)

        # Set a1 to the correct orientation
        r1_to_r2 = r2 - r1
        if direction_35:
            r1_to_r2 = - r1_to_r2
        r1_to_r2 /= np.linalg.norm(r1_to_r2)

        # project the r1_to_r2 vector on the helix axis plane
        r1_to_r2_proj = r1_to_r2 - np.dot(r1_to_r2, dir_vector) * dir_vector
        # project the base vector on the helix axis plane
        old_a1_proj = old_a1 - np.dot(old_a1, dir_vector) * dir_vector
        # calculate the rotation angle to orient the helix axis
        rotAngle2 = np.arccos(np.clip(np.dot(r1_to_r2_proj / np.linalg.norm(r1_to_r2_proj), old_a1_proj / np.linalg.norm(old_a1_proj)), -1.0, 1.0))
        if np.dot(np.cross(r1_to_r2, old_a1), dir_vector) < 0:
            rotAngle2 *= -1
        # This rotoation matrix rotates the radius 1 and 2 to the point where the base vector is aligned with the old_a1
        q2 = R.from_rotvec(dir_vector * rotAngle2)
        r1 = q2.apply(r1)
        r2 = q2.apply(r2)

        # Center point of the helix axis
        r = r1
        if direction_35:
            r = r2
        # Calculate the start position
        start_pos = pos - r - old_a1 * fudge

        # Create per-step rotation matrix
        R_step = R.from_rotvec(dir_vector * rot)

        # Initialize properties of new nucleotide
        out = []
        out_double = []

        # Generate nucleotide positions and orientations
        for i in range(length):
            # Calculate rotation around central axis and step along axis
            if i > 0:
                r1 = R_step.apply(r1) + dir_vector * base_base_distance
                r2 = R_step.apply(r2) + dir_vector * base_base_distance

            # Calculate a1 orientation
            r1_to_r2 = r2 - r1
            if i == 0:
                r1_to_r2 = - r1_to_r2
            a1 = r1_to_r2 / np.linalg.norm(r1_to_r2)

            # Calculate A3 orientation
            dir_vector = - dir_vector # invert the direction of the helix axis
            dir_vector /= np.linalg.norm(dir_vector) # normalize
            a1proj = a1 - np.dot(a1, dir_vector) * dir_vector
            a1proj /= np.linalg.norm(a1proj)
            a3 = -np.cos(inclination) * dir_vector + np.sin(inclination) * a1proj
            a3 /= np.linalg.norm(a3)

            # Calculate position
            r = r1
            if direction_35:
                r = r2
                a1 = -a1
                a3 = -a3
            RNA_fudge = a1 * fudge
            p = r + start_pos + RNA_fudge * 2

            out.append([p, a1, a3])

            # Double helix case
            if double:
                a1 = -r1_to_r2 / np.linalg.norm(r1_to_r2)
                a1proj = a1 - np.dot(a1, -dir_vector) * -dir_vector
                a1proj /= np.linalg.norm(a1proj)
                a3 = np.cos(inclination) * dir_vector + np.sin(inclination) * a1proj
                a3 /= np.linalg.norm(a3)
                r = r2
                if direction_35:
                    r = r1
                    a1 = -a1
                    a3 = -a3
                RNA_fudge = a1 * fudge
                p = r + RNA_fudge + start_pos
                out_double.append([p, a1, a3])

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

import pyfurnace as pf

# MAKE A SIMPLE TILE AND CLEANUP
origami = pf.simple_origami(
    dt_list=[120],
    kl_columns=1,
    main_stem=[22],
    add_terminal_helix=True,
    align="first",
    use_angles=True,
)

# Adjust the length of the stems with Tetraloops
origami[(1, 7)].length = 5
origami[(1, 1)].length = 5

# Adjust the length of the stems with external KL
origami[(0, 7)].length = 7
origami[(2, 7)].length = 8
origami[(2, 1)].length = 6
origami[(0, 1)].length = 7

# Replace the tetraloops with kissing loops in the first and last helix
origami[(0, 0)] = pf.KissingLoop180(open_left=False, pk_index="1")
origami[(2, 0)] = pf.KissingLoop180(open_left=False, pk_index="2")
origami[(0, 8)] = pf.KissingLoop180(open_left=True, pk_index="2'")
origami[(2, 8)] = pf.KissingLoop180(open_left=True, pk_index="1'")

# Make three tiles
tile = origami.copy()
n_tiles = 3

for i in range(n_tiles - 1):
    # add and empty line, it will be filled when
    # the columns are repeated
    origami.append([])

    # add the tile to the origami
    for line in tile:
        origami.append(line, copy=True)


# REPEATE THE COLUMNS
n_columns = 2

# This is gonna be the basic column
column = origami.copy()
for col_ind in range(1, n_columns):

    # add two empty lines, it will be filled when
    # the columns are repeated
    origami.append([])
    origami.append([])

    # # add the column to the origami
    for line_ind, line in enumerate(column):
        line_index = line_ind + 2 * col_ind
        insert_at = (line_index, len(origami[line_index]))
        origami.insert(insert_at, line)


# Join adjacent 180 KL to kissing dimers
for y, line in enumerate(origami):

    # go through the line in reverse order
    for i in range(len(line) - 1, 0, -1):

        # if the two adjacent motifs are kissing loops
        # pop the one and replace the other with a kissing dimer
        if isinstance(line[i], pf.KL180) and isinstance(line[i - 1], pf.KL180):
            origami.pop((y, i))
            origami[y, i - 1] = pf.KissingDimer()

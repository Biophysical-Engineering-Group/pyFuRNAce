from itertools import product
import pyfurnace as pf

### ORIGAMI PARAMETERS
dovetails_list = [-3, -3]
kl_columns = 3
kl_delay = 150  # nucleotides delay to claculate the folding barrier

# collection of best configurations
# a configuration is a tuple of the form (stem_pos, best_start_position)
best_confs = []
total_min = float("inf")
n_helix = len(dovetails_list) + 2

starts = list(product(range(0, n_helix), repeat=kl_columns))[::-1]

### Find the best start position
for stem_pos in starts:

    ### Usually the central helices are not a good start position (heuristic)
    # Avoiding central helices screens out a lot of configurations, so
    # uncomment the next lines if you want to test them

    ### FILTERING OUT CENTRAL HELICES
    # skip = False
    # for h in range(3, n_helix - 3):
    #     if h in stem_pos:
    #         skip = True
    #         break
    # if skip:
    #     continue
    ### END OF FILTERING

    print("Stem pos:", stem_pos)

    # Create the origami with the given stem positions
    # and no start/end stem
    origami = pf.simple_origami(
        dovetails_list,
        kl_columns=kl_columns,
        main_stem=33,
        stem_pos=stem_pos,
        align="first",
        add_start_end=False,
    )

    start_barrier = origami.assembled.folding_barriers(kl_delay=kl_delay)[1]
    print("Start barrier:", start_barrier)

    db, stacks = pf.dot_bracket_to_stacks(origami.structure)
    min = start_barrier
    best_middle = 0
    for db, (start, end) in zip(db, stacks):
        if db not in "()":
            continue
        pos = (start + end) // 2
        new_strucutre = pf.rotate_dot_bracket(origami.structure, pos)
        new_bar = pf.folding_barriers(structure=new_strucutre, kl_delay=kl_delay)[1]
        if new_bar < min:
            min = new_bar
            best_middle = pos
            best_structure = new_strucutre
        if new_bar < total_min:
            total_min = new_bar
            best_confs = [(stem_pos, pos)]
        elif new_bar == total_min:
            best_confs.append((stem_pos, pos))

    print("Best start:", best_middle, "Barrier:", min)

    print()

print("Best min penalty:", total_min)
print("Best confs:")
for conf in best_confs:
    print("Main stem positions:", conf[0], "Best start position:", conf[1])

# Make the origami with the best configuration
best_stem_pos = best_confs[0][0]
best_start_pos = best_confs[0][1]
origami = pf.simple_origami(
    dovetails_list,
    kl_columns=kl_columns,
    main_stem=33,
    stem_pos=best_stem_pos,
    align="first",
    add_start_end=False,
)

# Add the start/end stem, but check the orientation
for flip in range(2):

    # create a copy of the origami, get the slice and the motif
    ori_copy = origami.copy()
    start_slice = origami.get_slice_at_seq_index(best_start_pos)
    m = origami[start_slice]

    ### IMPORTANT REMARK:
    ### THIS ASSUMENTS THE MOTIF HAS A LENGTH PROPERTY
    stem_1 = pf.Stem(m.length // 2)
    start_end = pf.start_end_stem()
    if flip:
        start_end.flip()
    stem2 = pf.Stem(m.length - stem_1.length)
    ori_copy[start_slice] = [stem_1, start_end, stem2]

    # This is the good origami, save it
    if ori_copy.folding_barriers(kl_delay=kl_delay)[1] == total_min:
        origami = ori_copy
        break

print("Example origami:")
print(origami)

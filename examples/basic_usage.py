import pyfurnace as pf
from pyfurnace.generate import generate_road

line1 = [
    pf.TetraLoop(),
    pf.Stem(7),
    pf.Dovetail(-2, up_cross=False),
    pf.Stem(6),
    pf.KissingDimer(),
    pf.Stem(6),
    pf.Dovetail(-2, up_cross=False),
    pf.Stem(4),
    pf.Broccoli(),
    pf.Stem(4),
    pf.TetraLoop(True),
]

line2 = [
    pf.TetraLoop(),
    pf.Stem(7),
    pf.Dovetail(-2, down_cross=False),
    pf.Stem(10),
    pf.start_end_stem(),
    pf.Stem(10),
    pf.Dovetail(-2, down_cross=False),
    pf.Stem(7),
    pf.TetraLoop(True),
]

origami = pf.Origami(line1, line2, aling="center")

print(origami)
print(origami.structure)
print(origami.sequence)

### Optimization
constraints = origami.sequence
constraints = "GGGA" + constraints[4:]

sequence = generate_road(
    origami.structure, sequence=constraints, pseudoknots=origami.pseudoknots
)

origami.sequence = sequence
print(origami)
origami.save_3d_model("basic_origami")
origami.save_fasta("basic_origami")

from pyfurnace.design import template_2_helix, pseudo_to_dot
from pyfurnace.generate import generate_road, parallel_road, fold


def test_road():
    origami = template_2_helix()
    sequence = generate_road(
        structure=origami.structure,
        sequence=origami.sequence,
        pseudoknots=origami.pseudoknots,
    )
    seq_fold = fold(sequence)
    assert len(sequence) == len(origami.sequence)
    assert seq_fold == origami.structure.translate(pseudo_to_dot)


def test_parallel_road():
    origami = template_2_helix()
    sequence = parallel_road(
        structure=origami.structure,
        sequence=origami.sequence,
        pseudoknots=origami.pseudoknots,
        zip_directory=False,
    )
    seq_fold = fold(sequence)
    assert len(sequence) == len(origami.sequence)
    assert seq_fold == origami.structure.translate(pseudo_to_dot)

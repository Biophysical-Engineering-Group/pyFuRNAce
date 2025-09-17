import os
import re
import pytest
import tempfile
from pathlib import Path

from pyfurnace import simple_origami
from pyfurnace.design import (
    template_2_helix,
    template_rectangle_10H_3X,
    pseudo_to_dot,
    KissingLoop180,
)
from pyfurnace.generate import generate_road, parallel_road, fold


def test_road():
    origami = template_2_helix()
    origami[0, -1] = KissingLoop180(open_left=1, pk_index=1)
    origami[1, -1] = KissingLoop180(open_left=1, pk_index=-1)
    with tempfile.TemporaryDirectory() as tmpdir:
        dir_path = Path(tmpdir) / "road_test"
        pattern = rf"{re.escape(str(dir_path))} does not exist\. Creating it\."
        with pytest.warns(UserWarning, match=pattern):
            sequence = generate_road(
                structure=origami.structure,
                sequence=origami.sequence,
                pseudoknots=origami.pseudoknots,
                initial_sequence=origami.sequence.get_random_sequence(),
                verbose=True,
                directory=str(dir_path),
            )
    seq_fold = fold(sequence)
    assert len(sequence) == len(origami.sequence)
    assert seq_fold == origami.structure.translate(pseudo_to_dot)
    assert (
        origami[0, -1].get_kissing_sequence()
        == origami[1, -1].get_kissing_sequence().reverse_complement()
    )


def test_generate_road_raises_on_multistrand_in_structure():
    structure = "..&..(((...)))"  # contains '&'
    sequence = "AUAUAGCGCGCGC"  # any length is fine; we fail before length matters
    with pytest.raises(ValueError, match=r"does not support multistranded"):
        generate_road(structure, sequence)


def test_generate_road_raises_on_multistrand_in_sequence():
    structure = "....(((...)))"  # no '&' here
    sequence = "AU&AUGCGCGCGC"  # '&' in sequence
    with pytest.raises(ValueError, match=r"does not support multistranded"):
        generate_road(structure, sequence)


def test_generate_road_raises_on_length_mismatch():
    structure = "....(((...)))"  # length 12
    sequence = "A" * 11  # mismatch
    with pytest.raises(
        ValueError, match=r"length of the sequence and structure must match"
    ):
        generate_road(structure, sequence)


def test_timeout_road():
    origami = template_rectangle_10H_3X()
    with pytest.warns(
        UserWarning,
        match=r"Optimization failed. Please check the ROAD algorithm output.",
    ):
        generate_road(
            structure=origami.structure,
            sequence=origami.sequence,
            pseudoknots=origami.pseudoknots,
            initial_sequence=origami.sequence.get_random_sequence(),
            timeout=0.001,
        )


def test_road_zipping():
    # also don't add pseudoknots here
    origami = simple_origami([-2])
    sequence, zip_path = generate_road(
        structure=origami.structure,
        sequence=origami.sequence,
        callback=print,
        origami_code="import pyfurnace as pf\npf.simple_origami([-2])",
        zip_directory=True,
    )
    origami[1, 2].structure = "..&.."
    origami[1, 6].structure = "..&.."
    seq_fold = fold(sequence)
    assert len(sequence) == len(origami.sequence)
    assert seq_fold == origami.structure.translate(pseudo_to_dot)
    assert zip_path.endswith(".zip")
    assert os.path.exists(zip_path)


def test_timeout_parallel_road():
    origami = template_rectangle_10H_3X()
    with pytest.raises(
        RuntimeError, match=r"No successful ROAD design found in any trial."
    ):
        parallel_road(
            structure=origami.structure,
            sequence=origami.sequence,
            pseudoknots=origami.pseudoknots,
            initial_sequence=origami.sequence.get_random_sequence(),
            timeout=0.001,
        )


def test_parallel_road():
    origami = template_2_helix()
    sequence = parallel_road(
        structure=origami.structure,
        sequence=origami.sequence,
        pseudoknots=origami.pseudoknots,
        zip_directory=False,
        wait_for_all=False,
    )
    seq_fold = fold(sequence)
    assert len(sequence) == len(origami.sequence)
    assert seq_fold == origami.structure.translate(pseudo_to_dot)


def test_parallel_road_wait():
    n_trials = 2
    origami = template_2_helix()
    sequences = parallel_road(
        structure=origami.structure,
        sequence=origami.sequence,
        pseudoknots=origami.pseudoknots,
        n_trials=n_trials,
        zip_directory=False,
        wait_for_all=True,
    )
    seq_fold = fold(sequences[0])
    assert len(sequences) == 2
    assert len(sequences[0]) == len(origami.sequence)
    assert seq_fold == origami.structure.translate(pseudo_to_dot)


def test_parallel_road_wait_zip():
    n_trials = 2
    origami = template_2_helix()
    with tempfile.TemporaryDirectory() as tmpdir:
        dir_path = Path(tmpdir) / "parallel_road_test"
        pattern = rf"{re.escape(str(dir_path))} did not exist\. Created it\."
        with pytest.warns(UserWarning, match=pattern):
            sequences = parallel_road(
                structure=origami.structure,
                sequence=origami.sequence,
                pseudoknots=origami.pseudoknots,
                n_trials=n_trials,
                zip_directory=True,
                wait_for_all=True,
                save_to=str(dir_path),
            )
    seq_fold = fold(sequences[0])
    assert len(sequences) == n_trials
    assert len(sequences[0]) == len(origami.sequence)
    assert seq_fold == origami.structure.translate(pseudo_to_dot)

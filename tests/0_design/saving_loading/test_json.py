import pytest
import json
from pyfurnace import Strand, Motif, Stem, Dovetail, L7Ae, Origami, simple_origami


def test_strand_to_json_includes_all_attributes():
    strand = Strand("--AUGCUAGCUAGC-", start=(1, 2), direction=(1, 0))
    strand.custom_attribute = "custom_value"
    strand_json = strand.to_json()
    assert "custom_attribute" in strand_json
    assert strand_json["custom_attribute"] == "custom_value"
    assert strand_json["start"] == "1,2"
    assert strand_json["direction"] == "1,0"
    assert strand_json["strand"] in ("──AUGCUAGCUAGC-", "--AUGCUAGCUAGC-")
    assert (
        "pyfurnace_version" not in strand_json
    )  # version is added only when saving to file

    strand_json_save = json.loads(strand.save_json(return_data=True))
    assert (
        "pyfurnace_version" in strand_json_save
    )  # version is added when saving to file
    strand_json_save.pop("pyfurnace_version")
    # doing dumps comparison to avoid issues with ordering or tuples/lists in json
    assert json.loads(json.dumps(strand_json)) == strand_json_save


def test_strand_save_json(tmp_path):
    strand = Strand("--AUGCUAGCUAGC-", start=(1, 2), direction=(1, 0))
    strand.custom_attribute = "custom_value"
    strand.save_json(tmp_path / "test_strand.json")

    with open(tmp_path / "test_strand.json", "r") as f:
        strand_json = json.load(f)

    assert "custom_attribute" in strand_json
    assert strand_json["custom_attribute"] == "custom_value"
    assert strand_json["start"] == "1,2"
    assert strand_json["direction"] == "1,0"
    assert strand_json["strand"] in ("──AUGCUAGCUAGC-", "--AUGCUAGCUAGC-")


def test_strand_from_json(tmp_path):
    strand = Strand("--AUGCUAGCUAGC-", start=(1, 2), direction=(1, 0))
    strand.custom_attribute = "custom_value"
    strand.save_json(tmp_path / "test_strand.json")
    with open(tmp_path / "test_strand.json", "r") as f:
        strand_string = f.read()

    with pytest.raises(ValueError):
        Strand.from_json(123)  # invalid input type

    loaded_strand = Strand.from_json(strand_string)
    loaded_from_file = Strand.from_json_file(tmp_path / "test_strand.json")
    assert loaded_strand.start == strand.start
    assert loaded_strand.direction == strand.direction
    assert loaded_strand.strand == strand.strand
    assert loaded_strand.custom_attribute == strand.custom_attribute
    assert loaded_strand == loaded_from_file

    # check if you can load it back correctly
    new_strand = Strand.from_json_file(tmp_path / "test_strand.json")
    assert new_strand.start == strand.start
    assert new_strand.direction == strand.direction
    assert new_strand.strand == strand.strand
    assert new_strand.custom_attribute == strand.custom_attribute


def test_motif_to_json_includes_all_attributes():
    motif = Motif(Strand("--AUGCUAGCUAGC-", start=(1, 2), direction=(1, 0)))
    motif.custom_attribute = "custom_value"
    motif_json = motif.to_json()

    assert "custom_attribute" in motif_json
    assert motif_json["custom_attribute"] == "custom_value"
    assert motif_json["motif_type"] == "Motif"
    assert motif_json["basepair"] == dict()
    assert motif_json["strands"][0]["start"] == "1,2"
    assert motif_json["strands"][0]["direction"] == "1,0"
    assert motif_json["strands"][0]["strand"] in ("──AUGCUAGCUAGC-", "--AUGCUAGCUAGC-")
    assert (
        "pyfurnace_version" not in motif_json
    )  # version is added only when saving to file

    motif_json_save = json.loads(motif.save_json(return_data=True))
    assert (
        "pyfurnace_version" in motif_json_save
    )  # version is added when saving to file
    motif_json_save.pop("pyfurnace_version")
    assert json.loads(json.dumps(motif_json)) == motif_json_save


def test_stem_save_json_includes_all_attributes(tmp_path):
    stem = Stem(10, wobble_interval=0, wobble_tolerance=0)
    stem.custom_attribute = "custom_value"
    stem.save_json(tmp_path / "test_stem.json")

    with open(tmp_path / "test_stem.json", "r") as f:
        stem_dict = json.load(f)

    assert "custom_attribute" in stem_dict
    assert stem_dict["custom_attribute"] == "custom_value"
    assert stem_dict["motif_type"] == "Stem"
    assert stem_dict["_length"] == 10
    assert stem_dict["strands"][0]["strand"] == "N" * 10
    assert stem_dict["strands"][0]["start"] == "0,0"


def test_strands_block_dummy_ends_in_json():
    dt = Dovetail(-2)

    # check reloading of dummy ends
    strand1 = dt.strands[0]
    s1_json = strand1.to_json()
    assert "dummy_ends" in s1_json["coords"]
    s1_reload = Strand.from_json(json.dumps(s1_json))
    assert s1_reload.coords.dummy_ends[0].size == strand1.coords.dummy_ends[0].size
    assert s1_reload.coords.dummy_ends[1].size == strand1.coords.dummy_ends[1].size

    # check reloading of dovetail motif
    dt_reload = Motif.from_json(json.dumps(dt.to_json()))
    # it's only a unique strands block for dovetail
    assert len(set(id(s.strands_block) for s in dt_reload.strands)) == 1
    for s in dt_reload:
        assert not s.coords.is_empty()


def test_motif_save_load_with_proteins(tmp_path):
    mot = L7Ae()
    mot.custom_attribute = "custom_value"
    mot.save_json(tmp_path / "test_mot.json")
    with open(tmp_path / "test_mot.json", "r") as f:
        mot_string = f.read()

    mot_dict = json.loads(mot_string)
    assert "custom_attribute" in mot_dict
    assert mot_dict["custom_attribute"] == "custom_value"
    assert mot_dict["motif_type"] == "Aptamer"
    assert mot_dict["strands"][0]["start"] == "0,0"
    assert mot_dict["strands"][1]["coords"]["proteins"]

    with pytest.raises(ValueError):
        Motif.from_json(123)  # invalid input type

    new_mot = Motif.from_json(mot_string)
    assert new_mot.__class__.__name__ == mot.__class__.__name__
    assert new_mot.custom_attribute == mot.custom_attribute
    assert new_mot.strands[0].strand == mot.strands[0].strand
    assert new_mot.strands[1].strand == mot.strands[1].strand

    # load directly from file
    new_mot_file = Motif.from_json_file(tmp_path / "test_mot.json")
    assert new_mot_file.__class__.__name__ == mot.__class__.__name__
    assert new_mot_file.custom_attribute == mot.custom_attribute
    assert new_mot_file.strands[0].strand == mot.strands[0].strand
    assert new_mot_file.strands[1].strand == mot.strands[1].strand
    assert (
        new_mot_file[1].coords.proteins[0].sequence
        == mot.strands[1].coords.proteins[0].sequence
    )
    assert (
        new_mot_file[1].coords.proteins[0].coords
        == mot.strands[1].coords.proteins[0].coords
    ).all()


def test_origami_to_json(tmp_path):
    ori = simple_origami(dt_list=[-1], kl_columns=1)
    ori_from_json = Origami.from_json(
        ori.to_json()
    )  # the json is consumed to create the new origami

    ori_from_json.save_json(tmp_path / "test_origami.json")
    ori_from_file2 = Origami.from_json_file(tmp_path / "test_origami.json")
    for i, s in enumerate(ori.assembled):
        assert s.strand == ori_from_file2.assembled[i].strand
    assert ori.assembled.basepair == ori_from_file2.assembled.basepair
    assert ori.assembled == ori_from_file2.assembled

    assert (
        "pyfurnace_version" not in ori.to_json()
    )  # version is added only when saving to file
    ori_json_save = json.loads(ori.save_json(return_data=True))
    assert "pyfurnace_version" in ori_json_save  # version is added when saving to file
    ori_json_save.pop("pyfurnace_version")
    assert json.loads(json.dumps(ori.to_json())) == ori_json_save

    with open("ori1.json", "w") as f:
        json.dump(ori.to_json(), f, indent=4)
    with open("ori2.json", "w") as f:
        json.dump(ori_json_save, f, indent=4)


if __name__ == "__main__":
    # from pathlib import Path
    # test_origami_to_json(Path('./'))
    test_strands_block_dummy_ends_in_json()

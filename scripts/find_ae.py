"""
Find AEs in Origami
This script allows users to interactively find and visualize the parameters for
antiparallel-even turns (AE) crossovers in RNA origami structures. Users can adjust
the helix parameters and AE crossover parameters using sliders, and visualize the
resulting 3D structure with highlighted crossover points.

Usage:
    streamlit run scripts/find_ae.py
"""

import pyfurnace as pf
import numpy as np
import streamlit as st
from st_oxview import oxview_from_text

if __name__ == "__main__":

    st.set_page_config(layout="wide")
    st.title("Find AEs in Origami")

    st.write("RNA Helix parameters:")
    with st.container(horizontal=1, horizontal_alignment="distribute"):
        diameter = st.slider("Diameter", min_value=0.0, max_value=10.0, value=2.35)
        inclination_deg = st.slider(
            "Inclination (degrees)", min_value=-45.0, max_value=45.0, value=-15.5
        )
        inclination = inclination_deg * np.pi / 180  # convert to radians
        bp_backbone_distance = st.slider(
            "Base pair backbone distance",
            min_value=0.0,
            max_value=5.0,
            value=2.0,
        )
        base_base_distance = st.slider(
            "Base to base distance",
            min_value=0.0,
            max_value=1.0,
            value=0.3287,
        )

    pf.Coords.set_helix_params(
        diameter=diameter,
        inclination=inclination,
        bp_backbone_distance=bp_backbone_distance,
        base_base_distance=base_base_distance,
    )

    st.write("AE crossover parameters:")
    with st.container(horizontal=1, horizontal_alignment="distribute"):
        ae_v_shift = st.slider(
            "AE vertical shift", min_value=-3.0, max_value=3.0, value=1.006
        )
        ae_rot = st.slider(
            "AE rotation (degrees)", min_value=-45.0, max_value=90.0, value=40.267
        )
        # convert to radians
        ae_rot = ae_rot * np.pi / 180
        ae_hx_shift = st.slider(
            "AE horizontal x shift", min_value=-3.0, max_value=3.0, value=2.293
        )
        ae_hy_shift = st.slider(
            "AE horizontal y shift", min_value=-3.0, max_value=3.0, value=-0.5214
        )

    pf.Coords.set_AE_crossover_params(
        ae_v_shift=ae_v_shift,
        ae_rot=ae_rot,
        ae_hx_shift=ae_hx_shift,
        ae_hy_shift=ae_hy_shift,
    )

    with st.container(
        horizontal=1, horizontal_alignment="distribute", vertical_alignment="center"
    ):
        helix_length = st.number_input(
            "Helix length (number of base pairs)", min_value=1, max_value=100, value=8
        )
        dt_list = st.text_input("Dovetail list (comma-separated integers)", value="-3")
        dt_list = [int(x.strip()) for x in dt_list.split(",")]
        ss_assembly = st.toggle(
            "Single-stranded assembly", value=False, key="ss_assembly"
        )

    h1, h2 = pf.Coords.compute_AE_crossover(
        helix_length, return_helices=True, use_cached=False
    )

    stem1 = pf.Stem(length=len(h1) // 2)
    stem2 = pf.Stem(length=len(h2) // 2)

    # load coordinates into stems
    stem1[0].coords = h1[: len(h1) // 2]
    stem1[1].coords = h1[len(h1) // 2 :]
    stem2[0].coords = h2[: len(h2) // 2]
    stem2[1].coords = h2[len(h2) // 2 :]
    origami = pf.Origami([[stem1], [stem2]])

    # clear cached AE matrices to force recomputation with new parameters
    conf, top = origami.save_3d_model(return_text=True)

    simp_ori = pf.simple_origami(dt_list, main_stem=33)

    simp_ori.ss_assembly = ss_assembly
    conf2, top2 = simp_ori.save_3d_model(return_text=True)

    tot_len = len(h1) + len(h2)
    middle = helix_length // 2
    with st.container(horizontal=1, horizontal_alignment="distribute"):
        oxview_from_text(
            configuration=conf,
            topology=top,
            frame_id=1,
            colormap="Reds",
            index_colors=[
                (
                    0
                    if i
                    not in (
                        middle,
                        middle + 1,
                        tot_len - middle - 1,
                        tot_len - middle - 2,
                    )
                    else 1
                )
                for i in range(tot_len)
            ],
            key="display_nano",
        )
        oxview_from_text(
            configuration=conf2,
            topology=top2,
            frame_id=2,
            key="display_ss",
        )

"""
Plot oxDNA-style configuration files with Plotly.
The plot shows nucleotide positions as points, base vectors as solid lines,
and normal vectors as dotted lines.
The plot is interactive and allows zooming and rotation.

Requires pyFuRNAce: pip install pyfurnace
Requires plotly: pip install plotly
Requires `structure.dat` file in the same directory.
Author: Luca Monari
Date: December 2025

Usage:
    python scripts/plot_oxDNA_conf.py


"""

import numpy as np
import pyfurnace as pf
import plotly.graph_objects as go


def load_oxdna_file(filename):
    """
    Load oxDNA-style file where each row has:
    x y z  bx by bz  nx ny nz

    Lines starting with '#' are treated as comments.
    """
    data = np.loadtxt(filename, comments="#")
    pos = data[:, 0:3]  # x, y, z
    base = data[:, 3:6]  # bx, by, bz
    normal = data[:, 6:9]  # nx, ny, nz
    return pos, base, normal


def _make_vector_lines(pos, vec, scale=1.0):
    """
    Build x, y, z arrays for line segments representing vectors.

    Each vector is drawn as a line from pos[i] to pos[i] + scale * vec[i].
    `None` is inserted between segments for Plotly.
    """
    xs, ys, zs = [], [], []
    for p, v in zip(pos, vec):
        xs.extend([p[0], p[0] + scale * v[0], None])
        ys.extend([p[1], p[1] + scale * v[1], None])
        zs.extend([p[2], p[2] + scale * v[2], None])
    return xs, ys, zs


def plot_oxdna_vectors_plotly(
    coords,
    base_scale=1.0,
    normal_scale=1.0,
    line_length_factor=0.4,  # <--- fraction of vector to draw as a line
    base_line_width=6,  # <--- thickness of base vector lines
    normal_line_width=4,  # <--- thickness of normal vector lines
):
    array = coords.array
    pos = array[:, 0]
    base = array[:, 1]
    normal = array[:, 2]

    # Scatter points for nucleotide positions
    pos_trace = go.Scatter3d(
        x=pos[:, 0],
        y=pos[:, 1],
        z=pos[:, 2],
        mode="markers",
        marker=dict(size=3),
        name="Positions",
    )

    # Lines for base vectors (shorter & thicker)
    bx, by, bz = _make_vector_lines(pos, base, scale=base_scale * line_length_factor)
    base_trace = go.Scatter3d(
        x=bx,
        y=by,
        z=bz,
        mode="lines",
        line=dict(width=base_line_width),  # thicker
        name="Base vectors",
    )

    # Lines for normal vectors (shorter & thicker)
    nx, ny, nz = _make_vector_lines(
        pos, normal, scale=normal_scale * line_length_factor
    )
    normal_trace = go.Scatter3d(
        x=nx,
        y=ny,
        z=nz,
        mode="lines",
        line=dict(width=normal_line_width, dash="dot"),  # thicker
        name="Normal vectors",
    )

    fig = go.Figure(data=[pos_trace, base_trace, normal_trace])

    fig.update_layout(
        title="oxDNA positions, base vectors, and normal vectors (interactive)",
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z",
            aspectmode="data",  # equal aspect ratio
        ),
        showlegend=True,
    )

    fig.show()


if __name__ == "__main__":
    filename = "structure.dat"
    coords = pf.Coords.load_from_file(filename)
    plot_oxdna_vectors_plotly(
        coords,
        base_scale=1.0,
        normal_scale=1.0,
        line_length_factor=0.4,  # try 0.2–0.5
        base_line_width=6,
        normal_line_width=4,
    )

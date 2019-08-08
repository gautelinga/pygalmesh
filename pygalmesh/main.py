# -*- coding: utf-8 -*-
#
import os
import tempfile

import meshio

from _pygalmesh import (
    _generate_mesh,
    _generate_from_inr,
    _generate_from_off,
    _generate_from_off_with_features,
    _generate_periodic_mesh,
    _generate_surface_mesh,
    _remesh_surface,
    _orient_surface_mesh,
)


def generate_mesh(
    domain,
    feature_edges=None,
    bounding_sphere_radius=0.0,
    lloyd=False,
    odt=False,
    perturb=True,
    exude=True,
    edge_size=0.0,
    facet_angle=0.0,
    facet_size=0.0,
    facet_distance=0.0,
    cell_radius_edge_ratio=0.0,
    cell_size=0.0,
    verbose=True,
):
    feature_edges = [] if feature_edges is None else feature_edges

    fh, outfile = tempfile.mkstemp(suffix=".mesh")
    os.close(fh)

    _generate_mesh(
        domain,
        outfile,
        feature_edges=feature_edges,
        bounding_sphere_radius=bounding_sphere_radius,
        lloyd=lloyd,
        odt=odt,
        perturb=perturb,
        exude=exude,
        edge_size=edge_size,
        facet_angle=facet_angle,
        facet_size=facet_size,
        facet_distance=facet_distance,
        cell_radius_edge_ratio=cell_radius_edge_ratio,
        cell_size=cell_size,
        verbose=verbose,
    )

    mesh = meshio.read(outfile)
    os.remove(outfile)
    return mesh


def generate_periodic_mesh(
    domain,
    bounding_cuboid,
    lloyd=False,
    odt=False,
    perturb=True,
    exude=True,
    edge_size=0.0,
    facet_angle=0.0,
    facet_size=0.0,
    facet_distance=0.0,
    cell_radius_edge_ratio=0.0,
    cell_size=0.0,
    number_of_copies_in_output=1,
    verbose=True,
):
    fh, outfile = tempfile.mkstemp(suffix=".mesh")
    os.close(fh)

    assert number_of_copies_in_output in [1, 2, 4, 8]

    _generate_periodic_mesh(
        domain,
        outfile,
        bounding_cuboid,
        lloyd=lloyd,
        odt=odt,
        perturb=perturb,
        exude=exude,
        edge_size=edge_size,
        facet_angle=facet_angle,
        facet_size=facet_size,
        facet_distance=facet_distance,
        cell_radius_edge_ratio=cell_radius_edge_ratio,
        cell_size=cell_size,
        number_of_copies_in_output=number_of_copies_in_output,
        verbose=verbose,
    )

    mesh = meshio.read(outfile)
    os.remove(outfile)
    return mesh


def generate_surface_mesh(
    domain,
    bounding_sphere_radius=0.0,
    angle_bound=0.0,
    radius_bound=0.0,
    distance_bound=0.0,
    verbose=True,
):
    fh, outfile = tempfile.mkstemp(suffix=".off")
    os.close(fh)

    _generate_surface_mesh(
        domain,
        outfile,
        bounding_sphere_radius=bounding_sphere_radius,
        angle_bound=angle_bound,
        radius_bound=radius_bound,
        distance_bound=distance_bound,
        verbose=verbose,
    )

    mesh = meshio.read(outfile)
    os.remove(outfile)
    return mesh


def remesh_surface(
    mesh,
    edge_size=0.0,
    facet_angle=0.0,
    facet_size=0.0,
    facet_distance=0.0,
    verbose=True
):
    fh_in, infile = tempfile.mkstemp(suffix=".off")
    os.close(fh_in)
    meshio.write(infile, mesh)

    fh_out, outfile = tempfile.mkstemp(suffix=".off")
    os.close(fh_out)

    _remesh_surface(
        infile,
        outfile,
        edge_size=edge_size,
        facet_angle=facet_angle,
        facet_size=facet_size,
        facet_distance=facet_distance,
        verbose=verbose,
    )
    mesh = meshio.read(outfile)
    os.remove(infile)
    os.remove(outfile)
    return mesh


def orient_surface_mesh(
    mesh,
    verbose=True
):
    fh_in, infile = tempfile.mkstemp(suffix=".off")
    os.close(fh_in)
    meshio.write(infile, mesh)

    fh_out, outfile = tempfile.mkstemp(suffix=".off")
    os.close(fh_out)

    _orient_surface_mesh(
        infile,
        outfile,
        verbose=verbose,
    )
    mesh_out = meshio.read(outfile)
    os.remove(infile)
    os.remove(outfile)
    return mesh_out


def generate_volume_mesh_from_surface_mesh(
    mesh,
    lloyd=False,
    odt=False,
    perturb=True,
    exude=True,
    detect_features=False,
    edge_size=0.0,
    facet_angle=0.0,
    facet_size=0.0,
    facet_distance=0.0,
    cell_radius_edge_ratio=0.0,
    cell_size=0.0,
    verbose=True,
):
    if isinstance(mesh, str):
        mesh = meshio.read(mesh)

    fh, off_file = tempfile.mkstemp(suffix=".off")
    os.close(fh)
    meshio.write(off_file, mesh)

    fh, outfile = tempfile.mkstemp(suffix=".mesh")
    os.close(fh)

    if detect_features:
        _generate_from_off_with_features(
            off_file,
            outfile,
            lloyd=lloyd,
            odt=odt,
            perturb=perturb,
            exude=exude,
            edge_size=edge_size,
            facet_angle=facet_angle,
            facet_size=facet_size,
            facet_distance=facet_distance,
            cell_radius_edge_ratio=cell_radius_edge_ratio,
            cell_size=cell_size,
            verbose=verbose,
        )
    else:
        _generate_from_off(
            off_file,
            outfile,
            lloyd=lloyd,
            odt=odt,
            perturb=perturb,
            exude=exude,
            edge_size=edge_size,
            facet_angle=facet_angle,
            facet_size=facet_size,
            facet_distance=facet_distance,
            cell_radius_edge_ratio=cell_radius_edge_ratio,
            cell_size=cell_size,
            verbose=verbose,
        )

    mesh = meshio.read(outfile)
    os.remove(off_file)
    os.remove(outfile)
    return mesh


def generate_from_inr(
    inr_filename,
    lloyd=False,
    odt=False,
    perturb=True,
    exude=True,
    edge_size=0.0,
    facet_angle=0.0,
    facet_size=0.0,
    facet_distance=0.0,
    cell_radius_edge_ratio=0.0,
    cell_size=0.0,
    verbose=True,
):
    fh, outfile = tempfile.mkstemp(suffix=".mesh")
    os.close(fh)

    _generate_from_inr(
        inr_filename,
        outfile,
        lloyd=lloyd,
        odt=odt,
        perturb=perturb,
        exude=exude,
        edge_size=edge_size,
        facet_angle=facet_angle,
        facet_size=facet_size,
        facet_distance=facet_distance,
        cell_radius_edge_ratio=cell_radius_edge_ratio,
        cell_size=cell_size,
        verbose=verbose,
    )

    mesh = meshio.read(outfile)
    os.remove(outfile)
    return mesh

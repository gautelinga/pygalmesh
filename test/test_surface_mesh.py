# -*- coding: utf-8 -*-
#
import numpy
import meshio

import pygalmesh

import helpers


def test_sphere():
    radius = 1.0
    s = pygalmesh.Ball([0.0, 0.0, 0.0], radius)
    mesh = pygalmesh.generate_surface_mesh(
        s, angle_bound=30, radius_bound=0.1, distance_bound=0.1, verbose=False
    )

    tol = 1.0e-2
    assert abs(max(mesh.points[:, 0]) - radius) < tol
    assert abs(min(mesh.points[:, 0]) + radius) < tol
    assert abs(max(mesh.points[:, 1]) - radius) < tol
    assert abs(min(mesh.points[:, 1]) + radius) < tol
    assert abs(max(mesh.points[:, 2]) - radius) < tol
    assert abs(min(mesh.points[:, 2]) + radius) < tol

    areas = helpers.compute_triangle_areas(mesh.points, mesh.cells["triangle"])
    surface_area = sum(areas)
    assert abs(surface_area - 4 * numpy.pi * radius ** 2) < 0.1
    return


def test_remesh_lion():
    helpers.download("lion-head.off", "e4b00c47fee69a37a0d95ab71a63cc74",
                     url="https://raw.githubusercontent.com/CGAL/cgal/master/Mesh_3/examples/Mesh_3/data/")
    mesh_in = meshio.read("/tmp/lion-head.off")
    mesh_out = pygalmesh.remesh_surface(
        mesh_in,
        edge_size=0.025,
        facet_angle=25,
        facet_size=0.1,
        facet_distance=0.001)

    area_in = sum(helpers.compute_triangle_areas(
        mesh_in.points,
        mesh_in.cells["triangle"]))
    area_out = sum(helpers.compute_triangle_areas(
        mesh_out.points,
        mesh_out.cells["triangle"]))
    err_rel = abs(area_in - area_out) / area_in

    assert err_rel < 0.005
    return


def test_orient_lion():
    helpers.download("lion-head.off", "e4b00c47fee69a37a0d95ab71a63cc74",
                     url="https://raw.githubusercontent.com/CGAL/cgal/master/Mesh_3/examples/Mesh_3/data/")
    mesh = meshio.read("/tmp/lion-head.off")
    faces = mesh.cells["triangle"]
    area = sum(helpers.compute_triangle_areas(mesh.points, faces))

    # Shuffle some normals
    for i in [10, 20, 30, 40, 50]:
        faces[i, [0, 1]] = faces[i, [1, 0]]

    mesh.cells["triangle"] = faces
    mesh_oriented = pygalmesh.orient_surface_mesh(mesh)
    area_oriented = sum(helpers.compute_triangle_areas(mesh_oriented.points, mesh_oriented.cells["triangle"]))

    assert area == area_oriented

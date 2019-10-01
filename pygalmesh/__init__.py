# https://github.com/pybind/pybind11/issues/1004
from _pygalmesh import (
    Ball,
    Cone,
    Cuboid,
    Cylinder,
    Difference,
    DomainBase,
    Ellipsoid,
    Extrude,
    HalfSpace,
    Intersection,
    Polygon2D,
    RingExtrude,
    Rotate,
    Scale,
    SizingFieldBase,
    Stretch,
    Tetrahedron,
    Torus,
    Translate,
    Union,
)

from .main import (
    generate_mesh,
    generate_periodic_mesh,
    generate_surface_mesh,
    remesh_surface,
    orient_surface_mesh,
    generate_volume_mesh_from_surface_mesh,
    generate_from_inr,
    generate_with_sizing_field,
)

from .__about__ import (
    __author__,
    __author_email__,
    __copyright__,
    __license__,
    __maintainer__,
    __status__,
    __version__,
)

__all__ = [
    "__author__",
    "__author_email__",
    "__copyright__",
    "__license__",
    "__version__",
    "__maintainer__",
    "__status__",
    #
    "DomainBase",
    "SizingFieldBase",
    "Translate",
    "Rotate",
    "Scale",
    "Stretch",
    "Intersection",
    "Union",
    "Difference",
    "Extrude",
    "Ball",
    "Cuboid",
    "Ellipsoid",
    "Tetrahedron",
    "Cone",
    "Cylinder",
    "Torus",
    "HalfSpace",
    "Polygon2D",
    "RingExtrude",
    #
    "generate_mesh",
    "generate_with_sizing_field",
    "generate_periodic_mesh",
    "generate_surface_mesh",
    "remesh_surface",
    "orient_surface_mesh",
    "generate_volume_mesh_from_surface_mesh",
    "generate_from_inr",
]

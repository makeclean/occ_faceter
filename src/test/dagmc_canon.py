"""Generate some canonical representation of a DAGMC file such that it can be
compared via the Unix utility diff(1).  The intention is that this tool can
be run on a file processed by occ_faceter and another file that was faceted by
another tool and get comparable output.  In practice, this only works for a
limited set of simple shapes as curves will be discreteized differently by
each faceter.
"""

from argparse import ArgumentParser
from collections import Counter
from dataclasses import dataclass
from typing import Any, Dict, List, Union

import h5py


@dataclass(frozen=True)
class Coords:
    coords: List[float]


@dataclass(frozen=True)
class Entity:
    etype: str
    data: object


@dataclass(frozen=True)
class MeshSet:
    contents: List[int]
    children: List[int]
    parents: List[int]
    flags: int


@dataclass(frozen=True)
class MoabEnts:
    coords: Dict[int, Coords]
    edges: Dict[int, Entity]
    tris: Dict[int, Entity]
    meshsets: Dict[int, MeshSet]

    # repeat of the above
    ents: Dict[int, Union[Coords, Entity, MeshSet]]


def _grouper(data):
    "splits data up into groups as lambda is called specifying the end of the group"
    prev = 0

    def next_group(pos: int):
        nonlocal prev
        pos = pos + 1
        assert pos >= prev
        result = data[prev:pos]
        prev = pos
        return result

    return next_group


def _expand_compressed(result):
    "expand a compressed group"
    it = iter(result)
    result = []
    for i in it:
        # this will raise if group doesn't contain an even number of elements
        n = next(it)
        result.extend(range(i, i + n))
    return result


def test_grouper():
    g = _grouper([1, 2, 101, 2])
    assert g(-1) == []
    assert g(1) == [1, 2]
    assert g(1) == []
    assert g(3) == [101, 2]
    # and with compressed groups
    assert _expand_compressed([]) == []
    assert _expand_compressed([101, 2]) == [101, 102]


def _enum_start(ds: h5py.Dataset) -> Any:
    return enumerate(ds[:].tolist(), ds.attrs["start_id"])


def load_moab_entities(tstt: h5py.Group) -> MoabEnts:
    coords = {
        idx: Coords(coords=row) for idx, row in _enum_start(tstt["nodes/coordinates"])
    }

    edges: Dict[int, Entity] = {}
    tris: Dict[int, Entity] = {}
    others: Dict[int, Entity] = {}
    entlookup = {
        "Edge": edges,
        "Tri": tris,
    }

    for grp in tstt["elements"].values():
        if not isinstance(grp, h5py.Group):
            break

        # get name of group's element type
        et = grp.attrs.get_id("element_type").get_type()
        typename = et.enum_nameof(grp.attrs["element_type"]).decode()
        d = entlookup.get(typename, others)

        for idx, row in _enum_start(grp["connectivity"]):
            d[idx] = Entity(typename, row)

    next_contents = _grouper(tstt["sets/contents"][:].tolist())
    next_children = _grouper(tstt["sets/children"][:].tolist())
    next_parents = _grouper(tstt["sets/parents"][:].tolist())

    meshsets: Dict[int, MeshSet] = {}

    for idx, [content, child, parent, flags] in _enum_start(tstt["sets/list"]):

        # pull out contents, expanding groups as necessary
        contents = next_contents(content)
        if flags & 0x8:
            contents = _expand_compressed(contents)

        meshsets[idx] = MeshSet(
            contents=contents,
            children=next_children(child),
            parents=next_parents(parent),
            flags=flags,
        )

    return MoabEnts(
        coords=coords,
        edges=edges,
        tris=tris,
        meshsets=meshsets,
        ents={**coords, **edges, **tris, **others, **meshsets},
    )


def main():
    parser = ArgumentParser(
        description=(
            "Generate a canonical representation of a DAGMC file such that "
            "two facetings can be compared via diff(1)."
        )
    )
    parser.add_argument("file", help="DAGMC file to open")
    args = parser.parse_args()

    with h5py.File(args.file) as h5m:
        ents = load_moab_entities(h5m["/tstt"])

    coords = {i: c.coords for i, c in ents.coords.items()}
    for coord in sorted(coords.values()):
        print("coord", coord)

    edge_coords = [sorted(coords[i] for i in edge.data) for edge in ents.edges.values()]
    edge_coords.sort()
    for c1, c2 in edge_coords:
        print("edge", c1, c2)

    # no canonical way of turning faces into triangles, only thing that
    # seems to make sense is their count
    print("tris", len(ents.tris))

    sets = Counter(
        (len(s.contents), len(s.children), len(s.parents))
        for s in ents.meshsets.values()
    )

    for ccp, n in sorted(sets.items()):
        print("meshset", ccp, n)


if __name__ == "__main__":
    main()

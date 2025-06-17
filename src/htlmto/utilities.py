import numpy as np
from collections import defaultdict
from cifkit import Cif
import datetime


from cif_reader import read_cif


def print_to_console(func):
    """Print current step and its staus at the end to console."""

    def wrapper(**kwargs):
        func_name = func.__name__[4:]
        if func_name == "lm":
            print("\tRunning optimization ", end=" ")
        else:
            print(f"\tRunning {func_name}", end=" ")
        res = func(**kwargs)

        if func_name == "lm":
            print("error" if res[0] else " ok")
        else:
            print("error" if res[0] else "ok")
        return res

    return wrapper


def get_d_by_dmin_CN(v):

    points_wd = [[p[3], p[1], p[0]] for p in v]

    # sort
    points_wd = sorted(points_wd, key=lambda x: x[1])[:20]
    distances = np.array([p[1] for p in points_wd])
    distances /= distances.min()

    gaps = np.array(
        [distances[i] - distances[i - 1] for i in range(1, len(distances))]
    )
    ind_gaps = np.argsort(gaps)

    inner_gap_index = np.array(ind_gaps[::-1])
    outer_CN = int(inner_gap_index[0]) + 1

    return outer_CN


def get_distances_from_cifkit(cifpath):
    cif = Cif(cifpath)
    cif.compute_connections()

    max_distances = defaultdict(dict)
    # {site: {site1: d1, site2: d2, ...}, site2: {}, ...}
    conns = cif.connections
    for k, v in conns.items():
        cn = get_d_by_dmin_CN(v)

        v = sorted([p[:2] for p in v], key=lambda x: x[1])[:cn]
        for _site in set([_v[0] for _v in v]):
            neigh_d_w_site_label = [_v[1] for _v in v if _v[0] == _site]
            # print(k, _site, max(neigh_d_w_site_label), neigh_d_w_site_label)
            max_distances[k][_site] = max(neigh_d_w_site_label)

    return dict(max_distances)


def extract_data_from_cif(cif_path):
    """Parse CIF and extract the information needed for LMTO calculation."""

    cif = read_cif(cif_path)
    cell = cif.cell

    if cif.id:
        name = f"{cif.formula}-{cif.id}"
    else:
        name = f"{cif.formula}-{str(datetime.datetime.now())[:-7]}"

    data = {
        "_cell_length_a": cell[0],
        "_cell_length_b": cell[1],
        "_cell_length_c": cell[2],
        "_cell_angle_alpha": cell[3],
        "_cell_angle_beta": cell[4],
        "_cell_angle_gamma": cell[5],
        "_space_group_IT_number": cif.space_group_number,
        "origin": 2 if cif.has_origin_choice_2 else 1,
        "atom_site_data": cif.site_data,
        "name": name,
    }

    return data

"""Certified OUTER enclosures for shared-raw-error EMA audit equations.

Standard library only for the certificate; numpy only ranks candidate bases.
All input decimals are interpreted as exact rationals. No target nodes, true
beta, execution roots, or true report poles are accepted by the auditor.
The algebraic conditions are necessary, not sufficient for positive-model
realizability. A nonempty enclosure is never called a feasible-model proof.
"""
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
import argparse, json, math


def rat(x):
    return x if isinstance(x, Q) else Q(str(x))


def inv(a):
    n = len(a)
    t = [list(row) + [Q(int(i == j)) for j in range(n)] for i, row in enumerate(a)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if t[i][j]), None)
        if pivot is None:
            return None
        t[j], t[pivot] = t[pivot], t[j]
        v = t[j][j]
        t[j] = [x / v for x in t[j]]
        for i in range(n):
            if i != j:
                v = t[i][j]
                t[i] = [x - v * y for x, y in zip(t[i], t[j])]
    return [r[n:] for r in t]


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), Q(0))


def mv(a, b):
    return [dot(r, b) for r in a]


def equations(costs, rho, lengths):
    """Each row [A_1,A_2,A_3,b] consists of linear raw-data forms."""
    offsets = [sum(lengths[:i]) for i in range(len(lengths))]
    size = sum(lengths)
    zero = [Q(0)] * size
    rows, labels = [], []
    for a, b in combinations(range(len(costs)), 2):
        n = min(lengths[a], lengths[b])
        C, D = [], []
        for j in range(n):
            c, d = zero.copy(), zero.copy()
            c[offsets[b] + j], c[offsets[a] + j] = costs[b], -costs[a]
            d[offsets[b] + j], d[offsets[a] + j] = Q(1), Q(-1)
            C.append(c); D.append(d)
        W = [[x - y for x, y in zip(C[j], C[j - 1] if j else zero)] for j in range(n)]
        H = [[C[j + 1][l] - C[j][l] - rho * D[j][l] for l in range(size)] for j in range(n - 1)]
        for j in range(n - 1):
            rows.append([H[j], W[j - 1] if j else zero.copy(),
                         [-x for x in (H[j - 1] if j else C[0])], W[j]])
            labels.append([a, b, j])
    return rows, labels


def select_bases(A, labels, limit):
    import numpy as np
    m = np.array([[float(v) for v in r] for r in A])
    candidates = []
    for I in combinations(range(len(A)), 3):
        sub = m[list(I)]
        scales = np.linalg.norm(sub, axis=1)
        if np.any(scales == 0):
            continue
        # This is a heuristic ordering only. Exact rational inversion follows.
        s = np.linalg.svd(sub / scales[:, None], compute_uv=False)
        score = s[-1] / s[0] if s[0] else 0
        candidates.append((score, I))
    candidates.sort(reverse=True)
    chosen = [I for _, I in candidates[:limit]]
    # Always retain every original two-fee basis, so adding channels cannot
    # weaken the interval when the same original data and bounds are retained.
    original = [i for i, lab in enumerate(labels) if lab[:2] == [0, 1]]
    for I in combinations(original, 3):
        if I not in chosen:
            chosen.append(I)
    return chosen


def enclosure_for_basis(rows, A, b, eps, I):
    T = inv([A[i] for i in I])
    if T is None:
        return None
    u0 = mv(T, [b[i] for i in I])
    G = [[Q(0) for _ in range(3)] for _ in range(3)]
    d = [Q(0)] * 3
    for j, radius in enumerate(eps):
        if not radius:
            continue
        Aj = [[rows[i][k][j] for k in range(3)] for i in I]
        bj = [rows[i][3][j] for i in I]
        qj = mv(T, [x - y for x, y in zip(bj, mv(Aj, u0))])
        # The SAME e_j enters every row, and b_j-A_j*u0 is combined
        # before taking absolute values. There are no formula-specific errors.
        Gj = [[dot(T[k], [Aj[l][v] for l in range(3)]) for v in range(3)] for k in range(3)]
        for k in range(3):
            d[k] += radius * abs(qj[k])
            for v in range(3):
                G[k][v] += radius * abs(Gj[k][v])
    q = max(sum(row) for row in G)
    if q >= 1:
        return None
    M = [[Q(int(i == j)) - G[i][j] for j in range(3)] for i in range(3)]
    w = mv(inv(M), d)
    assert all(x >= 0 for x in w)
    return [(u0[k] - w[k], u0[k] + w[k]) for k in range(3)], q


def audit(payload):
    costs = list(map(rat, payload['costs']))
    rho = rat(payload['rho'])
    reports = [list(map(rat, row)) for row in payload['reports']]
    lengths = list(map(len, reports))
    if len(costs) != len(reports) or len(costs) < 2 or min(lengths) < 2:
        raise ValueError('At least two cost channels and two readings per channel are required.')
    if rho <= 0 or any(k <= 0 for k in costs) or len(set(costs)) != len(costs):
        raise ValueError('Known rho and distinct positive internal costs are required.')
    raw = sum(reports, [])
    error = payload['error_bounds']
    if not isinstance(error, list):
        eps = [rat(error)] * len(raw)
    else:
        eps = [rat(v) for v in sum(error, [])] if isinstance(error[0], list) else list(map(rat, error))
    if len(eps) != len(raw) or any(e < 0 for e in eps):
        raise ValueError('One nonnegative simultaneous error bound per raw reading is required.')
    rows, labels = equations(costs, rho, lengths)
    A = [[dot(form, raw) for form in row[:3]] for row in rows]
    b = [dot(row[3], raw) for row in rows]
    import numpy as np
    singular = np.linalg.svd(np.array([[float(v) for v in r] for r in A]), compute_uv=False)
    bounds = [(Q(0), Q(1)) for _ in range(3)]
    used, qs = [], []
    reason = None
    for I in select_bases(A, labels, int(payload.get('basis_limit', 8))):
        result = enclosure_for_basis(rows, A, b, eps, I)
        if result is None:
            continue
        intervals, q = result
        used.append(list(I)); qs.append(q)
        bounds = [(max(l, a), min(h, c)) for (l, h), (a, c) in zip(bounds, intervals)]
    if any(l > h for l, h in bounds):
        reason = 'Disjoint certified enclosures or physical parameter bounds.'
    # Propagate t=beta*r through the parameter box. This is inclusion-preserving.
    if reason is None:
        for _ in range(4):
            (bl, bh), (rl, rh), (tl, th) = bounds
            tl, th = max(tl, bl * rl), min(th, bh * rh)
            if rl > 0: bh = min(bh, th / rl)
            if rh > 0: bl = max(bl, tl / rh)
            if bl > 0: rh = min(rh, th / bl)
            if bh > 0: rl = max(rl, tl / bh)
            bounds = [(bl, bh), (rl, rh), (tl, th)]
            if any(l > h for l, h in bounds):
                reason = 'Product constraint t=beta*r has empty outer enclosure.'
                break
    rejected_rows = []
    if reason is None:
        center = [(l + h) / 2 for l, h in bounds]
        radii = [(h - l) / 2 for l, h in bounds]
        for i, row in enumerate(rows):
            residual = dot(A[i], center) - b[i]
            rad = dot([abs(v) for v in A[i]], radii)
            for j, e in enumerate(eps):
                shared_coeff = dot([row[k][j] for k in range(3)], center) - row[3][j]
                rad += e * (abs(shared_coeff) + dot([abs(row[k][j]) for k in range(3)], radii))
            if abs(residual) > rad:
                rejected_rows.append(labels[i])
        if rejected_rows:
            reason = 'A full-prefix cross-channel necessary equation excludes zero.'
        elif any(v + e < 0 for v, e in zip(raw, eps)):
            reason = 'Report positivity is incompatible with the simultaneous error box.'
    goal = rat(payload.get('max_beta_width', '.02'))
    width = bounds[0][1] - bounds[0][0]
    status = ('MODEL_INCOMPATIBLE' if reason else
              'ABSTAIN_NO_CONTRACTION' if not used else
              'ABSTAIN_WIDE_INTERVAL' if width > goal else 'CONDITIONAL_INTERVAL')
    threshold = rat(payload.get('threshold', '.7'))
    decision = 'ABSTAIN'
    if status == 'CONDITIONAL_INTERVAL':
        if bounds[0][0] > threshold: decision = 'ABOVE_THRESHOLD_CONDITIONAL'
        elif bounds[0][1] < threshold: decision = 'BELOW_THRESHOLD_CONDITIONAL'
    def serial_interval(pair):
        l, h = pair
        return {'exact': [str(l), str(h)],
                'outward_float': [math.nextafter(float(l), -math.inf), math.nextafter(float(h), math.inf)]}
    return {
        'status': status, 'decision': decision, 'reason': reason,
        'beta_outer': None if reason else serial_interval(bounds[0]),
        'report_pole_outer': None if reason else serial_interval(bounds[1]),
        'beta_width': None if reason else float(width),
        'contract_validated': False, 'positive_model_feasibility_proved': False,
        'enclosure_type': 'certified outer enclosure of legal-model compatibility projection; not exact feasible set',
        'coverage': 'conditional on ALL raw errors lying simultaneously in the supplied box',
        'contracting_bases': len(used), 'minimum_contraction_bound': float(min(qs)) if qs else None,
        'basis_rows': used, 'rejected_equations': rejected_rows,
        'readings': lengths, 'rho_and_costs_assumed_known': True,
        'nominal_smallest_singular_value': float(singular[-1]),
        'nominal_condition_number': float(singular[0] / singular[-1]) if singular[-1] else None,
        'arithmetic': 'exact rational through all certificate decisions; float basis ranking only',
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = audit(json.loads(args.input.read_text(encoding='utf-8')))
    output = json.dumps(result, ensure_ascii=False, indent=2)
    if args.output: args.output.write_text(output, encoding='utf-8')
    else: print(output)

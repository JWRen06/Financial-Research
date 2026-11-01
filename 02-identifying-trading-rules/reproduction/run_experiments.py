from pathlib import Path
from decimal import Decimal as D, localcontext
from fractions import Fraction as Q
import json, csv, hashlib, time
from audit_interval import audit, equations, dot
from hidden_system import HiddenSystem, make_request

OUT = Path(__file__).resolve().parent


def contains(result, beta):
    if result['beta_outer'] is None: return False
    lo, hi = map(Q, result['beta_outer']['exact'])
    return lo <= Q(beta) <= hi


def contains_pole(result, pole):
    if result['report_pole_outer'] is None: return False
    lo, hi = map(Q, result['report_pole_outer']['exact'])
    return lo <= Q(pole) <= hi


def reduced(result):
    out = {k: result[k] for k in ('status', 'decision', 'beta_width', 'contracting_bases',
                                  'minimum_contraction_bound', 'readings',
                                  'nominal_smallest_singular_value', 'nominal_condition_number')}
    out['beta_interval'] = result['beta_outer']['outward_float'] if result['beta_outer'] else None
    out['report_pole_interval'] = result['report_pole_outer']['outward_float'] if result['report_pole_outer'] else None
    return out


def run():
    start = time.time(); rows = []; requests = []; coverage_count = 0
    def evaluate(group, label, system, costs, length, eps, beta='.8', **meta):
        nonlocal coverage_count
        req, covered = make_request(system, costs, length, eps)
        assert covered
        result = audit(req)
        if beta is not None:
            assert result['status'] != 'MODEL_INCOMPATIBLE', (group, label, result)
            assert contains(result, beta), (group, label, result)
            assert contains_pole(result, meta.get('pole', '.6')), (group, label, result)
            coverage_count += 1
        entry = dict(group=group, label=label, epsilon=eps, costs=costs,
                     **meta, **reduced(result))
        rows.append(entry); requests.append((entry, req, result))
        print(group, label, result['status'], result['beta_width'], flush=True)
        return result

    # Paired sweeps vary one contract feature at a time.
    for eps in ('1e-12', '1e-10', '1e-8', '1e-6', '1e-4'):
        evaluate('noise', eps, HiddenSystem(), ['1', '2'], 6, eps)
    for gap in ('1', '.1', '.01', '.001'):
        for eps in ('1e-10', '1e-8', '1e-6'):
            costs = ['1', str(D(1) + D(gap))]
            evaluate('cost_gap', gap + ':' + eps, HiddenSystem(), costs, 6, eps, gap=gap)
    for s in ('.8', '.95', '.99', '.999', '.9999'):
        for eps in ('1e-10', '1e-8', '1e-6'):
            evaluate('slow_target', s + ':' + eps, HiddenSystem(nodes=(s,), weights=('1',)),
                     ['1', '2'], 6, eps, persistence=s)
    for r in ('.1', '.6', '.9', '.99'):
        evaluate('report_memory', r, HiddenSystem(pole=r), ['1', '2'], 6, '1e-8', pole=r)
    # Fixed total twelve readings, followed by nested eighteen-reading design.
    for eps in ('1e-10', '1e-8', '1e-6'):
        system = HiddenSystem()
        base = evaluate('third_cost', '2x6:' + eps, system, ['1', '2'], 6, eps)
        evaluate('third_cost', '3x4:' + eps, system, ['1', '2', '4'], 4, eps)
        nested = evaluate('third_cost', '3x6:' + eps, system, ['1', '2', '4'], 6, eps)
        assert nested['beta_width'] <= base['beta_width'] + 1e-30
    # Report ALL candidates in this exploratory scan, including worse designs.
    for third_cost in ('.05', '.1', '.25', '.5', '1.5', '3', '4', '8', '20'):
        for n in (4, 6):
            evaluate('third_cost_scan', third_cost + ':' + str(n), HiddenSystem(),
                     ['1', '2', third_cost], n, '1e-8', third_cost=third_cost)
    for fault in ('different_poles', 'two_poles'):
        for eps in ('1e-12', '1e-8', '1e-4'):
            result = evaluate('misspecification', fault + ':' + eps,
                              HiddenSystem(misspec=fault), ['1', '2'], 6, eps,
                              beta=None, fault=fault,
                              counted_as_recovery_success=False)
            if eps == '1e-12':
                assert result['status'] == 'MODEL_INCOMPATIBLE', (fault, result)

    # Two DISTINCT legal worlds share exactly the same uncertain input.
    with localcontext() as ctx:
        ctx.prec = 65
        worlds = [HiddenSystem(beta=b, nodes=('.999999',), weights=('1',)) for b in ('.6', '.85')]
        values = [[w.reset_and_report(k, 6) for k in ['1', '2']] for w in worlds]
        centers = [[(a + b) / 2 for a, b in zip(x, y)] for x, y in zip(*values)]
        eps = max(abs(a - b) / 2 for x, y in zip(*values) for a, b in zip(x, y)) + D('1e-48')
        req = {'rho': '1', 'costs': ['1', '2'], 'reports': [[str(v) for v in x] for x in centers],
               'error_bounds': str(eps), 'max_beta_width': '.02', 'threshold': '.7'}
        result = audit(req)
        assert contains(result, '.6') and contains(result, '.85')
        assert result['decision'] == 'ABSTAIN' and result['status'].startswith('ABSTAIN')
        entry = dict(group='indistinguishable', label='midpoint_two_legal_worlds',
                     epsilon=str(eps), costs=['1', '2'], true_betas=['.6', '.85'], **reduced(result))
        rows.append(entry); requests.append((entry, req, result))

    # Exact rational legal world tests the zero-error and shared-dependency logic.
    beta, ea, eb, rho = Q(4, 5), Q(1, 4), Q(2, 5), Q(3)
    costs = [rho * beta * e / ((1 - e) * (1 - beta * e)) for e in (ea, eb)]
    def exact_channel(e):
        p, z = Q(0), Q(0); seq = []
        for n in range(6):
            forcing = sum(w * (1 - s) * s ** n / (1 - beta * e * s)
                          for s, w in [(Q(1, 5), Q(2, 5)), (Q(3, 5), Q(3, 5))])
            p = e * p + (1 - e) * (1 - beta * e) * forcing
            z = Q(3, 5) * z + Q(2, 5) * p
            seq.append(str(z))
        return seq
    req = {'rho': str(rho), 'costs': list(map(str, costs)),
           'reports': [exact_channel(ea), exact_channel(eb)], 'error_bounds': '0'}
    exact = audit(req)
    assert exact['status'] == 'CONDITIONAL_INTERVAL'
    assert exact['beta_outer']['exact'] == ['4/5', '4/5']
    assert exact['report_pole_outer']['exact'] == ['3/5', '3/5']

    # An exact three-cost/four-prefix witness, with no floating rank decision.
    ec = Q(3, 5)
    kc = rho * beta * ec / ((1 - ec) * (1 - beta * ec))
    req3 = {'rho': str(rho), 'costs': list(map(str, costs + [kc])),
            'reports': [exact_channel(e)[:4] for e in (ea, eb, ec)], 'error_bounds': '0'}
    result3 = audit(req3)
    assert result3['beta_outer']['exact'] == ['4/5', '4/5']
    assert result3['report_pole_outer']['exact'] == ['3/5', '3/5']
    forms3, labels3 = equations(costs + [kc], rho, [4, 4, 4])
    raw3 = list(map(Q, sum(req3['reports'], [])))
    A3 = [[dot(form, raw3) for form in row[:3]] for row in forms3]
    I = result3['basis_rows'][0]
    M = [A3[i] for i in I]
    a,b,c = M[0]; d,e,f = M[1]; g,h,i = M[2]
    determinant = a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
    assert determinant != 0
    rank_witness = {'beta': str(beta), 'rho': str(rho), 'costs': req3['costs'],
                    'prefix_length': 4, 'selected_equations': [labels3[i] for i in I],
                    'exact_nonzero_determinant': str(determinant),
                    'identification': 'this exact dataset, not a uniform minimum theorem',
                    'request': req3, 'result': result3}
    (OUT/'three_cost_exact_witness.json').write_text(json.dumps(rank_witness, indent=2), encoding='utf-8')

    # Adversarial one-coordinate error and correlated all-coordinate errors.
    adversarial = 0
    for mode in range(13):
        radius = Q(1, 10 ** 10)
        test = json.loads(json.dumps(req)); test['error_bounds'] = str(radius)
        flat = sum(test['reports'], [])
        for i in range(12):
            e = radius * (1 if i % 2 else -1) if mode == 12 or mode == i else Q(0)
            flat[i] = str(Q(flat[i]) + e)
        test['reports'] = [flat[:6], flat[6:]]
        res = audit(test)
        assert contains(res, '4/5')
        adversarial += 1

    varied = 0
    for bb in ('.2', '.55', '.8', '.97'):
        for rr in ('.2', '1', '5'):
            system = HiddenSystem(beta=bb, rho=rr, nodes=('0', '.4', '.95'),
                                  weights=('.2', '.5', '.3'), pole='.85', gain='.15')
            request, covered = make_request(system, ['.1', '3'], 6, '1e-10', known_rho=rr)
            answer = audit(request)
            assert covered and contains(answer, bb)
            assert contains_pole(answer, '.85')
            varied += 1

    for i, (entry, req, result) in enumerate(requests):
        target = OUT / 'audit_inputs'; target.mkdir(exist_ok=True)
        (target / f'case_{i:03d}.json').write_text(json.dumps(req, indent=2), encoding='utf-8')
    (OUT / 'results.json').write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding='utf-8')
    keys = sorted(set().union(*(r.keys() for r in rows)))
    with (OUT / 'results.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys); writer.writeheader(); writer.writerows(rows)
    summary = {'scenarios': len(rows), 'legal_single_world_coverage_checks': coverage_count,
               'indistinguishable_worlds_covered': 2, 'exact_zero_error_check': True,
               'three_cost_four_prefix_exact_rank_witness': True,
               'adversarial_shared_error_checks': adversarial,
               'additional_beta_rho_grid_checks': varied,
               'misspecification_cases': sum(r['group'] == 'misspecification' for r in rows),
               'misspecification_rejections': sum(r['group'] == 'misspecification' and r['status'] == 'MODEL_INCOMPATIBLE' for r in rows),
               'synthetic_contract_test_only': True, 'seconds': round(time.time() - start, 2),
               'audit_source_sha256': hashlib.sha256((OUT / 'audit_interval.py').read_bytes()).hexdigest()}
    (OUT / 'verification.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False), flush=True)


if __name__ == '__main__': run()

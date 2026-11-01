"""Controlled hidden optimizer/report service. Audit inputs contain reports only.

This is a synthetic contract test, not an institution or public-strategy result.
High precision Decimal simulation is independent of the rational auditor.
"""
from decimal import Decimal as D, localcontext
import hashlib


class HiddenSystem:
    def __init__(self, beta='.8', rho='1', nodes=('.1', '.45', '.8'),
                 weights=('.3', '.4', '.3'), pole='.6', gain='.4',
                 misspec=None):
        self.__beta, self.__rho = D(beta), D(rho)
        self.__nodes = list(map(D, nodes)); self.__weights = list(map(D, weights))
        self.__pole, self.__gain = D(pole), D(gain)
        self.__misspec = misspec

    def reset_and_report(self, internal_cost, prefix_length):
        """Only authorized control: reset to zero and set internal cost."""
        with localcontext() as ctx:
            ctx.prec = 70
            k = D(str(internal_cost)); beta, rho = self.__beta, self.__rho
            b = 1 + beta + rho * beta / k
            E = 2 / (b + (b * b - 4 * beta).sqrt())
            # Analytic optimal target look-ahead enters a state recursion.
            # The auditor receives neither this recursion nor these states.
            inventory, report1, report2 = D(0), D(0), D(0)
            values = []
            r = self.__pole
            if self.__misspec == 'different_poles':
                r += D('.02') * (k - 1)
            for n in range(prefix_length):
                forcing = sum(w * (1 - s) * (s ** n if n else D(1)) /
                              (1 - beta * E * s) for s, w in zip(self.__nodes, self.__weights))
                inventory = E * inventory + (1 - E) * (1 - beta * E) * forcing
                report1 = r * report1 + self.__gain * inventory
                if self.__misspec == 'two_poles':
                    report2 = D('.15') * report2 + self.__gain * inventory
                    value = D('.98') * report1 + D('.02') * report2
                else:
                    value = report1
                values.append(value)
            return values


def make_request(system, costs, length, epsilon, max_width='.02', known_rho='1'):
    """The serialized boundary has no hidden beta, poles, targets or roots."""
    with localcontext() as ctx:
        ctx.prec = 65
        eps = D(str(epsilon)); rounding = D('1e-48')
        reports = []
        covered = True
        for cost in costs:
            exact = system.reset_and_report(cost, length)
            channel = []
            for n, v in enumerate(exact):
                # Noise is indexed by cost and lag; adding a channel does not
                # change the already observed data. No stochastic coverage claim.
                code = int(hashlib.sha256(f'{cost}:{n}'.encode()).hexdigest()[:8], 16)
                sign = D(1 if code % 2 else -1)
                amplitude = D(2 + code % 5) / 10
                y = (v + sign * amplitude * eps).quantize(rounding)
                covered &= abs(y - v) <= eps + rounding
                channel.append(str(y))
            reports.append(channel)
        return {'rho': known_rho, 'costs': list(map(str, costs)),
                'reports': reports, 'error_bounds': str(eps + rounding),
                'max_beta_width': max_width, 'threshold': '.7'}, covered

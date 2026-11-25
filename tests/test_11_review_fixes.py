"""Regression tests for the corrections made after the v0.2.0 review."""

from __future__ import annotations

import numpy as np
import pytest

from combusta import (
    DRY_AIR,
    DeadStateConditions,
    compute_gci,
    compute_thermomechanical_exergy,
    hardenberg_hase_delay,
)


@pytest.mark.unit
class TestGCIFormalOrder:
    def test_formal_order_caps_an_optimistic_apparent_order(self) -> None:
        # phi(h) = 10 + 2 h^3 has apparent order exactly 3
        h = (1.0, 1.5, 2.25)
        phi = [10.0 + 2.0 * x**3 for x in h]
        free = compute_gci(phi[0], phi[1], phi[2], h1=h[0], h2=h[1], h3=h[2])
        capped = compute_gci(phi[0], phi[1], phi[2], h1=h[0], h2=h[1], h3=h[2], formal_order=2.0)

        assert free.apparent_order_p == pytest.approx(3.0, abs=1e-10)
        assert not free.order_limited and free.order_used_p == pytest.approx(3.0, abs=1e-10)
        assert capped.apparent_order_p == pytest.approx(3.0, abs=1e-10)
        assert capped.order_limited and capped.order_used_p == 2.0

        ea = abs(phi[0] - phi[1]) / phi[0]
        assert capped.gci_fine_21 == pytest.approx(1.25 * ea / (1.5**2 - 1.0), rel=1e-12)
        assert capped.gci_fine_21 > free.gci_fine_21  # conservative

    def test_formal_order_not_applied_below_it(self) -> None:
        h = (1.0, 1.5, 2.25)
        phi = [10.0 + 2.0 * x**2 for x in h]
        res = compute_gci(phi[0], phi[1], phi[2], h1=h[0], h2=h[1], h3=h[2], formal_order=2.0 + 1e-6)
        assert not res.order_limited

    def test_divergent_triplet_warns(self) -> None:
        with pytest.warns(UserWarning, match="Divergent grid triplet"):
            res = compute_gci(10.0, 12.0, 12.5, r21=1.5, r32=1.5)
        assert res.convergence_regime == "DIVERGENT"

    def test_invalid_formal_order(self) -> None:
        with pytest.raises(Exception, match="formal_order"):
            compute_gci(12.0, 14.5, 20.125, r21=1.5, r32=1.5, formal_order=-1.0)


@pytest.mark.unit
class TestHardenbergHaseForm:
    def test_default_is_published_product_form(self) -> None:
        t, p, sp, rpm, cn = 900.0, 6.0e6, 7.2, 2400.0, 50.0
        e_a = 618840.0 / (cn + 25.0)
        expected_deg = (0.36 + 0.22 * sp) * np.exp(
            e_a * (1.0 / (8.314462618 * t) - 1.0 / 17190.0) * (21.2 / (p * 1e-5 - 12.4)) ** 0.63
        )
        res = hardenberg_hase_delay(t, p, sp, rpm, cn)
        assert float(res.ignition_delay_deg) == pytest.approx(expected_deg, rel=1e-12)

    def test_additive_form_is_deprecated(self) -> None:
        with pytest.warns(DeprecationWarning):
            hardenberg_hase_delay(900.0, 6.0e6, 7.2, 2400.0, 50.0, pressure_coupling="additive")

    def test_low_pressure_message_names_a_validity_limit(self) -> None:
        with pytest.raises(ValueError, match="not a physical auto-ignition limit"):
            hardenberg_hase_delay(900.0, 1.0e6, 7.2, 2400.0, 50.0)


@pytest.mark.unit
class TestThermomechanicalExergyUnclipped:
    def test_exergy_is_nonnegative_without_clipping(self) -> None:
        ds = DeadStateConditions.from_altitude(0.0)
        t, p = np.meshgrid(np.linspace(220.0, 2500.0, 40), np.geomspace(2.0e4, 2.0e7, 40))
        e = compute_thermomechanical_exergy(t.ravel(), p.ravel(), ds, 1.0e-3, gas_mixture=DRY_AIR)
        assert np.all(e >= -1e-9)  # Gibbs inequality, checked on the raw values
        assert float(compute_thermomechanical_exergy(np.array([ds.t0_k]), np.array([ds.p0_pa]), ds, 1.0e-3)[0]) == pytest.approx(0.0, abs=1e-9)

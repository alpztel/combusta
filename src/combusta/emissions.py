"""Exhaust emission post-processing: concentrations to masses, specific emissions, and
combustion efficiency.

Conventions follow common engine-testing practice (Heywood, 2018, Sec. 2.12 and 4.9):

* Specific emission $\\mathrm{is}_i = m_i / W$ [g/kWh] relative to indicated work $W$.
* Emission index $\\mathrm{EI}_i = m_i / m_f$ [g/kg fuel].
* NO$_x$ is reported as NO$_2$-equivalent mass, the convention of exhaust-emission
  regulations (e.g. ISO 8178, 40 CFR 1065), so a NO$_x$ mole fraction is converted
  with the molar mass of NO$_2$.
* Combustion efficiency from the chemical energy left in the exhaust
  (Heywood, 2018, Eq. 4.69): $\\eta_c = 1 - \\sum_i m_i \\mathrm{LHV}_i / (m_f \\mathrm{LHV}_f)$.

Mole fractions are converted to masses with the exhaust molar amount of complete
combustion at the global air-excess ratio (``product_moles``); trace species change the
total by far less than measurement uncertainty. Dry-basis concentrations are converted
to wet basis with the complete-combustion water fraction,
$x_\\text{wet} = x_\\text{dry} (1 - x_{\\text{H}_2\\text{O}})$, a simplification of the
dry-to-wet corrections in ISO 8178-1 that neglects intake humidity.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Final, Literal, TypeAlias

import numpy as np

from combusta.exceptions import PhysicalBoundsError
from combusta.thermo.fuels import ATOMIC_MASS_C, ATOMIC_MASS_H, ATOMIC_MASS_O, Fuel, product_moles
from combusta.thermo.nasa_polynomials import SPECIES_DB

ATOMIC_MASS_N: Final[float] = 14.0067e-3

#: Molar masses of common exhaust species [kg/mol] (IUPAC atomic weights).
SPECIES_MOLAR_MASS: Final[Mapping[str, float]] = {
    "NO": ATOMIC_MASS_N + ATOMIC_MASS_O,
    "NO2": ATOMIC_MASS_N + 2.0 * ATOMIC_MASS_O,
    "NOx": ATOMIC_MASS_N + 2.0 * ATOMIC_MASS_O,  # reported as NO2 equivalent
    "N2O": 2.0 * ATOMIC_MASS_N + ATOMIC_MASS_O,
    "NH3": ATOMIC_MASS_N + 3.0 * ATOMIC_MASS_H,
    "CO": ATOMIC_MASS_C + ATOMIC_MASS_O,
    "CO2": ATOMIC_MASS_C + 2.0 * ATOMIC_MASS_O,
    "H2": 2.0 * ATOMIC_MASS_H,
    "H2O": 2.0 * ATOMIC_MASS_H + ATOMIC_MASS_O,
    "CH4": ATOMIC_MASS_C + 4.0 * ATOMIC_MASS_H,
    "O2": 2.0 * ATOMIC_MASS_O,
}

#: Lower heating values of combustible exhaust species [J/kg] (Heywood, 2018, App. D).
#: Unburned liquid-fuel hydrocarbons are passed as ``"HC"`` and use the fuel LHV.
UNBURNED_SPECIES_LHV: Final[Mapping[str, float]] = {
    "H2": 119.96e6,
    "CO": 10.10e6,
    "CH4": 50.03e6,
}

ConcentrationBasis: TypeAlias = Literal["wet", "dry"]


@dataclass(frozen=True, slots=True)
class ExhaustComposition:
    """Complete-combustion exhaust amount per cycle.

    Attributes
    ----------
    total_moles : float
        Exhaust amount per cycle [mol].
    water_mole_fraction : float
        Wet-basis water mole fraction $x_{\\text{H}_2\\text{O}}$ [-].
    mean_molar_mass_kg_per_mol : float
        Mean exhaust molar mass [kg/mol].
    """

    total_moles: float
    water_mole_fraction: float
    mean_molar_mass_kg_per_mol: float


def exhaust_composition(fuel: Fuel, lambda_air: float, fuel_mass_kg: float) -> ExhaustComposition:
    """Exhaust amount and water fraction of complete combustion of ``fuel_mass_kg``."""
    if not np.isfinite(fuel_mass_kg) or fuel_mass_kg <= 0.0:
        raise PhysicalBoundsError("fuel_mass_kg must be finite and > 0.")
    moles = product_moles(fuel, lambda_air)
    per_fuel = sum(moles.values())
    n_fuel = fuel_mass_kg / fuel.molar_mass_kg_per_mol
    mass_per_fuel = sum(n * SPECIES_DB[sp].molar_mass_kg_per_mol for sp, n in moles.items())
    return ExhaustComposition(
        total_moles=n_fuel * per_fuel,
        water_mole_fraction=moles.get("H2O", 0.0) / per_fuel,
        mean_molar_mass_kg_per_mol=mass_per_fuel / per_fuel,
    )


def species_mass_from_mole_fraction(
    mole_fraction: float,
    species: str,
    *,
    fuel: Fuel,
    lambda_air: float,
    fuel_mass_kg: float,
    basis: ConcentrationBasis = "wet",
) -> float:
    """Convert an exhaust mole fraction to species mass per cycle [kg].

    Parameters
    ----------
    mole_fraction : float
        Mole fraction [-] (``ppm * 1e-6``). For ``species="NOx"`` the result is the
        NO2-equivalent mass.
    species : str
        Key of ``SPECIES_MOLAR_MASS``.
    fuel, lambda_air, fuel_mass_kg
        Fuel, global air-excess ratio, and fuel mass per cycle [kg].
    basis : {"wet", "dry"}, default="wet"
        Measurement basis; dry values are corrected with the complete-combustion water
        fraction.
    """
    if not np.isfinite(mole_fraction) or not 0.0 <= mole_fraction < 1.0:
        raise PhysicalBoundsError("mole_fraction must lie in [0, 1).")
    if species not in SPECIES_MOLAR_MASS:
        raise KeyError(f"Unknown species '{species}'. Known: {sorted(SPECIES_MOLAR_MASS)}.")
    exhaust = exhaust_composition(fuel, lambda_air, fuel_mass_kg)
    if basis == "dry":
        x_wet = mole_fraction * (1.0 - exhaust.water_mole_fraction)
    elif basis == "wet":
        x_wet = mole_fraction
    else:
        raise ValueError("basis must be 'wet' or 'dry'.")
    return x_wet * exhaust.total_moles * SPECIES_MOLAR_MASS[species]


def specific_emission(species_mass_kg: float, work_j: float) -> float:
    """Specific emission $m_i / W$ [kg/J]; multiply by 3.6e9 for g/kWh."""
    if work_j <= 0.0:
        raise PhysicalBoundsError("Work must be > 0 J for a specific emission.")
    if species_mass_kg < 0.0:
        raise PhysicalBoundsError("Species mass must be >= 0 kg.")
    return species_mass_kg / work_j


def emission_index(species_mass_kg: float, fuel_mass_kg: float) -> float:
    """Emission index $m_i / m_f$ [kg/kg]; multiply by 1e3 for g/kg fuel."""
    if fuel_mass_kg <= 0.0:
        raise PhysicalBoundsError("Fuel mass must be > 0 kg for an emission index.")
    if species_mass_kg < 0.0:
        raise PhysicalBoundsError("Species mass must be >= 0 kg.")
    return species_mass_kg / fuel_mass_kg


def combustion_efficiency(
    fuel: Fuel,
    fuel_mass_kg: float,
    unburned_mass_kg: Mapping[str, float],
    *,
    lhv_table: Mapping[str, float] = UNBURNED_SPECIES_LHV,
) -> float:
    """Combustion efficiency from chemical energy left in the exhaust (Heywood, 2018, Eq. 4.69).

    $$\\eta_c = 1 - \\frac{\\sum_i m_i\\,\\mathrm{LHV}_i}{m_f\\,\\mathrm{LHV}_f}$$

    ``unburned_mass_kg`` maps species to mass per cycle [kg]. Species in ``lhv_table``
    use their own heating value; ``"HC"`` (unburned fuel-like hydrocarbons) uses the fuel
    LHV. Species with no heating value (e.g. NOx) are ignored.
    """
    if fuel_mass_kg <= 0.0:
        raise PhysicalBoundsError("Fuel mass must be > 0 kg.")
    lost = 0.0
    for species, mass in unburned_mass_kg.items():
        if mass < 0.0:
            raise PhysicalBoundsError(f"Unburned mass of '{species}' must be >= 0 kg.")
        if species == "HC":
            lost += mass * fuel.lhv_j_per_kg
        elif species in lhv_table:
            lost += mass * lhv_table[species]
    eta = 1.0 - lost / (fuel_mass_kg * fuel.lhv_j_per_kg)
    if eta < 0.0:
        raise PhysicalBoundsError("Unburned chemical energy exceeds the fuel energy.")
    return eta


def has_heating_value(species: str) -> bool:
    """``True`` if ``species`` contributes to the combustion-efficiency loss."""
    return species == "HC" or species in UNBURNED_SPECIES_LHV

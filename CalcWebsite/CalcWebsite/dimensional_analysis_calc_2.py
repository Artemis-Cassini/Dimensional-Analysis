# Density-based quarts ↔ pounds conversions (density in kg/quart)
def quarts_to_pounds_using_density(qts, density_kg_per_qt):
    """
    Convert quarts to pounds using density (kg/quart).
    """
    kg = qts * density_kg_per_qt
    lbs = kg * 2.20462
    steps = [
        f"Step 1: Given {qts} quarts and density = {density_kg_per_qt} kg/quart",
        f"Step 2: Mass in kg = {qts} × {density_kg_per_qt} = {kg:.4f} kg",
        f"Step 3: Convert to pounds: {kg:.4f} × 2.20462 = {lbs:.4f} lb",
        f"Final Answer: {lbs:.4f} lb"
    ]
    return steps

def pounds_to_quarts_using_density(lbs, density_kg_per_qt):
    """
    Convert pounds to quarts using density (kg/quart).
    """
    kg = lbs / 2.20462
    qts = kg / density_kg_per_qt
    steps = [
        f"Step 1: Given {lbs} lb and density = {density_kg_per_qt} kg/quart",
        f"Step 2: Convert to kg: {lbs} ÷ 2.20462 = {kg:.4f} kg",
        f"Step 3: Volume in quarts = {kg:.4f} ÷ {density_kg_per_qt} = {qts:.4f} qts",
        f"Final Answer: {qts:.4f} quarts"
    ]
    return steps

# Time ↔ Months (assume 1 month = 30 days)
SECONDS_PER_MONTH = 30 * 24 * 3600
def months_to_seconds(months):
    sec = months * SECONDS_PER_MONTH
    steps = [
        f"Step 1: Given {months} months",
        f"Step 2: 1 month = {SECONDS_PER_MONTH} seconds",
        f"{months} × {SECONDS_PER_MONTH} = {sec:.2f} s",
        f"Final Answer: {sec:.2f} s"
    ]
    return steps

def seconds_to_months(sec):
    months = sec / SECONDS_PER_MONTH
    steps = [
        f"Step 1: Given {sec} seconds",
        f"Step 2: Divide by {SECONDS_PER_MONTH} to convert s → months",
        f"{sec} ÷ {SECONDS_PER_MONTH} = {months:.6f} months",
        f"Final Answer: {months:.6f} months"
    ]
    return steps

# Fuel consumption: miles ↔ gallons using rate (gal per km)
def miles_to_gallons_using_rate(miles, fuel_gal, dist_km):
    km = miles * 1.60934
    gal = fuel_gal * (km / dist_km)
    steps = [
        f"Step 1: Given rate = {fuel_gal} gal per {dist_km} km",
        f"Step 2: Convert {miles} miles to km: ×1.60934 = {km:.4f} km",
        f"Step 3: Gallons needed = {fuel_gal} × ({km} ÷ {dist_km}) = {gal:.4f} gal",
        f"Final Answer: {gal:.4f} gallons"
    ]
    return steps

def gallons_to_miles_using_rate(gal, fuel_gal, dist_km):
    # inverse: given gallons, compute miles
    km = dist_km * (gal / fuel_gal)
    miles = km / 1.60934
    steps = [
        f"Step 1: Given rate = {fuel_gal} gal per {dist_km} km",
        f"Step 2: km traveled = {dist_km} × ({gal} ÷ {fuel_gal}) = {km:.4f} km",
        f"Step 3: Convert to miles: {km:.4f} ÷ 1.60934 = {miles:.4f} miles",
        f"Final Answer: {miles:.4f} miles"
    ]
    return steps

# Weeks ↔ Minutes
def weeks_to_minutes(weeks):
    minutes = weeks * 7 * 24 * 60
    steps = [
        f"Step 1: Given {weeks} weeks",
        f"Step 2: 1 week = 7×24×60 = {7*24*60} minutes",
        f"{weeks} × {7*24*60} = {minutes:.2f} minutes",
        f"Final Answer: {minutes:.2f} minutes"
    ]
    return steps

def minutes_to_weeks(minutes):
    weeks = minutes / (7 * 24 * 60)
    steps = [
        f"Step 1: Given {minutes} minutes",
        f"Step 2: Divide by {7*24*60} to convert minutes → weeks",
        f"{minutes} ÷ {7*24*60} = {weeks:.6f} weeks",
        f"Final Answer: {weeks:.6f} weeks"
    ]
    return steps

# Distance-time ↔ ft/s
def miles_time_to_fps(miles, seconds):
    feet = miles * 5280
    fps = feet / seconds
    steps = [
        f"Step 1: Convert {miles} miles to feet: ×5280 = {feet:.2f} ft",
        f"Step 2: Divide by {seconds} s → {fps:.4f} ft/s",
        f"Final Answer: {fps:.4f} ft/s"
    ]
    return steps

def fps_to_miles_time(fps, miles=1):
    # inverse: time to travel one mile
    seconds = (miles * 5280) / fps
    steps = [
        f"Step 1: Given speed = {fps} ft/s",
        f"Step 2: Distance 1 mile = 5280 ft, time = 5280 ÷ {fps} = {seconds:.4f} s",
        f"Final Answer: {seconds:.4f} seconds"
    ]
    return steps
# Speed conversions: miles/hour ↔ feet/second
def mph_to_fps(mph):
    fps = mph * 5280 / 3600
    return [
        f"Step 1: Given {mph} miles/hour",
        f"Step 2: Convert miles to feet: × 5280 ft/mile",
        f"Step 3: Convert hours to seconds: ÷ 3600 s/hour",
        f"{mph} × 5280 ÷ 3600 = {fps:.4f} ft/s",
        f"Final Answer: {fps:.4f} ft/s"
    ]

def fps_to_mph(fps):
    mph = fps * 3600 / 5280
    return [
        f"Step 1: Given {fps} feet/second",
        f"Step 2: Convert seconds to hours: × 3600 s/hour",
        f"Step 3: Convert feet to miles: ÷ 5280 ft/mile",
        f"{fps} × 3600 ÷ 5280 = {mph:.4f} mi/h",
        f"Final Answer: {mph:.4f} mi/h"
    ]

# Weight conversions: pounds & ounces ↔ milligrams
def lb_oz_to_mg(lb, oz):
    mg = lb * 453592.37 + oz * 28349.523
    return [
        f"Step 1: Convert {lb} lb to mg: {lb} × 453592.37 = {lb*453592.37:.0f} mg",
        f"Step 2: Convert {oz} oz to mg: {oz} × 28349.523 = {oz*28349.523:.0f} mg",
        f"Step 3: Sum = {lb*453592.37 + oz*28349.523:.0f} mg",
        f"Final Answer: {mg:.0f} mg"
    ]

def mg_to_lb_oz(mg):
    lb = int(mg // 453592.37)
    rem_mg = mg - lb * 453592.37
    oz = rem_mg / 28349.523
    return [
        f"Step 1: Given {mg:.0f} mg",
        f"Step 2: Pounds = {mg:.0f} ÷ 453592.37 = {lb} lb",
        f"Step 3: Remainder = {rem_mg:.0f} mg → {oz:.2f} oz",
        f"Final Answer: {lb} lb and {oz:.2f} oz"
    ]

# Distance conversions: miles ↔ centimeters
def miles_to_cm(miles):
    cm = miles * 160934.4
    return [
        f"Step 1: Given {miles} miles",
        f"Step 2: Multiply by 160934.4 cm/mile",
        f"{miles} × 160934.4 = {cm:.0f} cm",
        f"Final Answer: {cm:.0f} cm"
    ]

def cm_to_miles(cm):
    miles = cm / 160934.4
    return [
        f"Step 1: Given {cm:.0f} cm",
        f"Step 2: Divide by 160934.4 cm/mile",
        f"{cm} ÷ 160934.4 = {miles:.6f} miles",
        f"Final Answer: {miles:.6f} miles"
    ]
# Polyprotic acid titration concentration calculation
def titration_concentration(V_titrant_L, M_titrant, titrant_comp, V_analyte_L, analyte_comp):
    """
    Calculate analyte molarity from titrant usage.
    Assumes analyte releases H+ equal to count of H in formula.
    """
    # 1) Moles of titrant used
    n_titrant = M_titrant * V_titrant_L
    # 2) Determine number of protons per analyte molecule
    counts = parse_formula(analyte_comp)
    if 'H' not in counts:
        return [f"Error: no ionizable H found in {analyte_comp}."]
    H_count = counts['H']
    # 3) Moles of analyte neutralized
    n_analyte = n_titrant / H_count
    # 4) Molarity of analyte
    M_analyte = n_analyte / V_analyte_L
    steps = []
    steps.append(f"Step 1: Moles of {titrant_comp} = {M_titrant} M × {V_titrant_L:.4f} L = {n_titrant:.4f} mol")
    steps.append(f"Step 2: Stoichiometry: {H_count} H⁺ per {analyte_comp}, so moles of {analyte_comp} = {n_titrant:.4f} ÷ {H_count} = {n_analyte:.4f} mol")
    steps.append(f"Step 3: Molarity of {analyte_comp} = {n_analyte:.4f} mol ÷ {V_analyte_L:.4f} L = {M_analyte:.4f} M")
    steps.append(f"Final Answer: {M_analyte:.4f} M")
    return steps
# Acid-base titration volume calculation (1:1 stoichiometry)
def titration_volume(M_base, base_compound, V_acid_L, M_acid, acid_compound):
    """
    Calculate the volume of base needed to titrate acid (1:1 stoichiometry).
    M_base, M_acid in mol/L; V_acid_L in L.
    """
    # 1) Calculate moles of acid
    acid_moles = M_acid * V_acid_L
    # 2) Moles of base needed (1:1)
    base_moles = acid_moles
    # 3) Calculate volume of base in liters
    V_base_L = base_moles / M_base
    # 4) Convert to mL
    V_base_mL = V_base_L * 1000
    steps = []
    steps.append(f"Step 1: Moles of {acid_compound} = {M_acid} M × {V_acid_L:.4f} L = {acid_moles:.4f} mol")
    steps.append(f"Step 2: For 1:1 reaction, moles of {base_compound} needed = {base_moles:.4f} mol")
    steps.append(f"Step 3: Volume of {base_compound} = moles ÷ M = {base_moles:.4f} mol ÷ {M_base} M = {V_base_L:.4f} L")
    steps.append(f"Step 4: Convert to mL: {V_base_L:.4f} L × 1000 = {V_base_mL:.1f} mL")
    steps.append(f"Final Answer: {V_base_L:.4f} L ({V_base_mL:.1f} mL)")
    return steps
import re
from collections import defaultdict

# Constants
AVOGADRO = 6.022e23  # molecules per mole

# Periodic table with molar masses (partial, but expandable)
PERIODIC_TABLE = {
    "H": 1.008, 
    "He": 4.0026,
    "Li": 6.94,
    "Be": 9.0122,
    "B": 10.81,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "Ne": 20.180,
    "Na": 22.990,
    "Mg": 24.305,
    "Al": 26.982,
    "Si": 28.085,
    "P": 30.974,
    "S": 32.06,
    "Cl": 35.45,
    "Ar": 39.948,
    "K": 39.098,
    "Ca": 40.078,
    "Sc": 44.956,
    "Ti": 47.867,
    "V": 50.942,
    "Cr": 51.996,
    "Mn": 54.938,
    "Fe": 55.845,
    "Co": 58.933,
    "Ni": 58.693,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ga": 69.723,
    "Ge": 72.630,
    "As": 74.922,
    "Se": 78.971,
    "Br": 79.904,
    "Kr": 83.798,
    "Rb": 85.469,
    "Sr": 87.62,
    "Y": 88.906,
    "Zr": 91.224,
    "Nb": 92.906,
    "Mo": 95.95,
    "Tc": 98,
    "Ru": 101.07,
    "Rh": 102.91,
    "Pd": 106.42,
    "Ag": 107.87,
    "Cd": 112.41,
    "In": 114.82,
    "Sn": 118.71,
    "Sb": 121.76,
    "Te": 127.60,
    "I": 126.90,
    "Xe": 131.29,
    "Cs": 132.91,
    "Ba": 137.33,
    "La": 138.91,
    "Ce": 140.12,
    "Pr": 140.91,
    "Nd": 144.24,
    "Pm": 145,
    "Sm": 150.36,
    "Eu": 151.96,
    "Gd": 157.25,
    "Tb": 158.93,
    "Dy": 162.50,
    "Ho": 164.93,
    "Er": 167.26,
    "Tm": 168.93,
    "Yb": 173.05,
    "Lu": 174.97,
    "Hf": 178.49,
    "Ta": 180.95,
    "W": 183.84,
    "Re": 186.21,
    "Os": 190.23,
    "Ir": 192.22,
    "Pt": 195.08,
    "Au": 196.97,
    "Hg": 200.59,
    "Tl": 204.38,
    "Pb": 207.2,
    "Bi": 208.98,
    "Po": 209,
    "At": 210,
    "Rn": 222,
    "Fr": 223,
    "Ra": 226,
    "Ac": 227,
    "Th": 232.04,
    "Pa": 231.04,
    "U": 238.03,
    "Np": 237,
    "Pu": 244,
    "Am": 243,
    "Cm": 247,
    "Bk": 247,
    "Cf": 251,
    "Es": 252,
    "Fm": 257, 
    "Md": 258,
    "No": 259,
    "Lr": 266,
    "Rf": 267,
    "Db": 268,
    "Sg": 269,
    "Bh": 270,
    "Hs": 277,
    "Mt": 278,
    "Ds": 281,
    "Rg": 282,
    "Cn": 285,
    "Nh": 286,
    "Fl": 289,
    "Mc": 290,
    "Lv": 293,
    "Ts": 294,
    "Og": 294,
}

def parse_formula(formula):
    formula = fix_case(formula)
    tokens = re.findall(r'([A-Z][a-z]?|\d+|\(|\)|\*|\.)', formula)
    stack = [defaultdict(int)]
    i = 0

    while i < len(tokens):
        token = tokens[i]

        if token == '(':
            stack.append(defaultdict(int))
        elif token == ')':
            temp = stack.pop()
            i += 1
            multiplier = int(tokens[i]) if i < len(tokens) and tokens[i].isdigit() else 1
            for elem, count in temp.items():
                stack[-1][elem] += count * multiplier
        elif re.match(r'[A-Z][a-z]?', token):
            elem = token
            i += 1
            count = int(tokens[i]) if i < len(tokens) and tokens[i].isdigit() else 1
            if not elem in PERIODIC_TABLE:
                raise ValueError(f"Unknown element: {elem}")
            stack[-1][elem] += count
            if i < len(tokens) and str(count) == tokens[i]:
             i += 1
            continue
        elif token in ('*', '.'):
            # hydrate separator with multiplier support
            i += 1
            # read leading multiplier if present
            if i < len(tokens) and tokens[i].isdigit():
                hydrate_mult = int(tokens[i])
                i += 1
            else:
                hydrate_mult = 1
            # parse the rest of the hydrate formula
            hydrate_formula = ''.join(tokens[i:])
            hydrate_counts = parse_formula(hydrate_formula)
            # apply the multiplier to each element
            for elem, count in hydrate_counts.items():
                stack[-1][elem] += count * hydrate_mult
            break
        i += 1

    return dict(stack[-1])

def fix_case(compound):
    i = 0
    result = ''
    while i < len(compound):
        if compound[i].isalpha():
            # Try two-letter element
            if i + 1 < len(compound) and compound[i:i+2] in PERIODIC_TABLE:
                result += compound[i:i+2]
                i += 2
            elif compound[i] in PERIODIC_TABLE:
                result += compound[i]
                i += 1
            else:
                # Try capitalizing properly to find match
                two_letter = compound[i:i+2].capitalize()
                one_letter = compound[i].upper()
                if two_letter in PERIODIC_TABLE:
                    result += two_letter
                    i += 2
                elif one_letter in PERIODIC_TABLE:
                    result += one_letter
                    i += 1
                else:
                    result += compound[i]
                    i += 1
        else:
            result += compound[i]
            i += 1
    return result


def molar_mass(formula):
    elements = parse_formula(formula)
    mass = 0
    for elem, count in elements.items():
        if elem not in PERIODIC_TABLE:
            raise ValueError(f"Unknown element: {elem}")
        mass += PERIODIC_TABLE[elem] * count
    return mass

# Grams ↔ Moles
def grams_to_moles(grams, molar_mass):
    moles = grams / molar_mass
    return [
        f"Step 1: Given {grams} g",
        f"Step 2: Molar mass = {molar_mass:.3f} g/mol",
        f"Step 3: Divide grams by molar mass",
        f"{grams} g ÷ {molar_mass:.3f} g/mol = {moles:.4f} mol",
        f"Final Answer: {moles:.4f} mol"
    ]

def moles_to_grams(moles, molar_mass):
    grams = moles * molar_mass
    return [
        f"Step 1: Given {moles} mol",
        f"Step 2: Molar mass = {molar_mass:.3f} g/mol",
        f"Step 3: Multiply moles by molar mass",
        f"{moles} mol × {molar_mass:.3f} g/mol = {grams:.2f} g",
        f"Final Answer: {grams:.2f} g"
    ]

# Moles ↔ Particles
def moles_to_particles(moles):
    particles = moles * AVOGADRO
    return [
        f"Step 1: Given {moles} mol",
        f"Step 2: Avogadro's number = {AVOGADRO:.3e} particles/mol",
        f"Step 3: Multiply moles by Avogadro’s number",
        f"{moles} mol × {AVOGADRO:.3e} = {particles:.3e} particles",
        f"Final Answer: {particles:.3e} particles"
    ]

def particles_to_moles(particles):
    moles = particles / AVOGADRO
    return [
        f"Step 1: Given {particles:.3e} particles",
        f"Step 2: Divide by Avogadro's number",
        f"{particles:.3e} ÷ {AVOGADRO:.3e} = {moles:.4f} mol",
        f"Final Answer: {moles:.4f} mol"
    ]

# Moles ↔ Liters (STP)
STP_VOLUME = 22.4  # L/mol

def moles_to_liters_stp(moles):
    liters = moles * STP_VOLUME
    return [
        f"Step 1: Given {moles} mol",
        f"Step 2: 1 mol = 22.4 L at STP",
        f"Step 3: Multiply moles by 22.4",
        f"{moles} mol × 22.4 L/mol = {liters:.2f} L",
        f"Final Answer: {liters:.2f} L"
    ]

def liters_to_moles_stp(liters):
    moles = liters / STP_VOLUME
    return [
        f"Step 1: Given {liters} L at STP",
        f"Step 2: 1 mol = 22.4 L at STP",
        f"Step 3: Divide liters by 22.4",
        f"{liters} L ÷ 22.4 L/mol = {moles:.4f} mol",
        f"Final Answer: {moles:.4f} mol"
    ]

# Molarity (mol ↔ L)
def molarity_moles_to_liters(moles, molarity):
    liters = moles / molarity
    return [
        f"Step 1: Given {moles} mol and {molarity} M",
        f"Step 2: Molarity = mol/L",
        f"Step 3: Volume = mol ÷ M",
        f"{moles} mol ÷ {molarity} mol/L = {liters:.3f} L",
        f"Final Answer: {liters:.3f} L"
    ]

def molarity_liters_to_moles(liters, molarity):
    moles = liters * molarity
    return [
        f"Step 1: Given {liters} L and {molarity} M",
        f"Step 2: Molarity = mol/L",
        f"Step 3: Moles = M × L",
        f"{molarity} mol/L × {liters} L = {moles:.4f} mol",
        f"Final Answer: {moles:.4f} mol"
    ]

# Molarity (Explicit Input)
def molarity_from_moles_and_volume(moles, liters):
    molarity = moles / liters
    return [
        f"Step 1: Given {moles} mol and {liters} L",
        f"Step 2: Molarity = mol ÷ L",
        f"{moles} ÷ {liters} = {molarity:.3f} M",
        f"Final Answer: {molarity:.3f} M"
    ]

# Density conversions
def grams_to_volume(grams, density):
    volume = grams / density
    return [
        f"Step 1: Given {grams} g and density = {density} g/mL",
        f"Step 2: Volume = grams ÷ density",
        f"{grams} g ÷ {density} g/mL = {volume:.2f} mL",
        f"Final Answer: {volume:.2f} mL"
    ]

def volume_to_grams(volume, density):
    grams = volume * density
    return [
        f"Step 1: Given {volume} mL and density = {density} g/mL",
        f"Step 2: Mass = volume × density",
        f"{volume} mL × {density} g/mL = {grams:.2f} g",
        f"Final Answer: {grams:.2f} g"
    ]


# Percent Yield
def percent_yield(actual, theoretical):
    percent = (actual / theoretical) * 100
    return [
        f"Step 1: Given actual = {actual} g and theoretical = {theoretical} g",
        f"Step 2: Percent yield = (actual ÷ theoretical) × 100",
        f"({actual} ÷ {theoretical}) × 100 = {percent:.2f}%",
        f"Final Answer: {percent:.2f}%"
    ]

# Limiting reagent and theoretical yield (full support)
def limiting_reagent_and_theoretical_yield(val1, comp1, val2, comp2, reaction_str):
    """
    Determine limiting reagent and theoretical yield of the first product.
    val1, val2 in moles; comp1, comp2 as formula strings; reaction_str like '2A+3B->4C+2D'.
    """
    import re
    # Parse reaction into reactants and products
    sides = reaction_str.replace(' ', '').split('->')
    reactants_list = sides[0].split('+')
    products_list = sides[1].split('+') if len(sides) > 1 else []

    # Parse stoichiometric coefficients for reactants
    stoich = {}
    for r in reactants_list:
        m = re.match(r'(\d*)([A-Za-z][A-Za-z0-9().*]+)', r)
        coeff = int(m.group(1)) if m.group(1) else 1
        form = m.group(2)
        stoich[form] = coeff

    # Map provided amounts
    amounts = {comp1: val1, comp2: val2}

    # Calculate available ratio for each reactant
    ratios = {}
    for form, coeff in stoich.items():
        if form in amounts:
            ratios[form] = amounts[form] / coeff

    # Identify limiting reagent
    limiting = min(ratios, key=ratios.get)
    min_ratio = ratios[limiting]

    # Build step-by-step output
    steps = ["Step 1: Compute moles available ÷ stoichiometric coefficient for each reactant"]
    for form, coeff in stoich.items():
        if form in amounts:
            steps.append(f"{amounts[form]} mol ÷ {coeff} = {ratios[form]:.4f}")
    steps.append(f"Limiting reagent is {limiting}")

    # Theoretical yield for the first product (if any)
    if products_list:
        prod = products_list[0]
        m = re.match(r'(\d*)([A-Za-z][A-Za-z0-9().*]+)', prod)
        prod_coeff = int(m.group(1)) if m.group(1) else 1
        prod_form = m.group(2)
        prod_moles = min_ratio * prod_coeff
        steps.append(f"Step 2: Theoretical yield of {prod_form} (in moles) = {min_ratio:.4f} × {prod_coeff} = {prod_moles:.4f} mol")
        # Convert to grams if possible
        try:
            mm = molar_mass(prod_form)
            grams = prod_moles * mm
            steps.append(f"Step 3: Convert to grams: {prod_moles:.4f} mol × {mm:.3f} g/mol = {grams:.2f} g")
        except Exception:
            pass

    return steps

# Stoichiometric mol-to-mol conversion helper
def stoichiometry_moles(val, from_comp, to_comp, reaction_str):
    """
    Convert moles of one species to moles of another via balanced reaction.
    val: moles of from_comp; reaction_str like '3H2+N2->2NH3'
    """
    import re
    sides = reaction_str.replace(' ', '').split('->')
    terms = sides[0].split('+') + (sides[1].split('+') if len(sides)>1 else [])
    stoich = {}
    for term in terms:
        m = re.match(r'(\d*)([A-Za-z][A-Za-z0-9().*]+)', term)
        coeff = int(m.group(1)) if m.group(1) else 1
        frm = m.group(2)
        stoich[frm] = coeff
    if from_comp not in stoich or to_comp not in stoich:
        return [f"Error: '{from_comp}' or '{to_comp}' not in reaction."]
    ratio = stoich[to_comp] / stoich[from_comp]
    result = val * ratio
    steps = [
        f"Step 1: Stoichiometric coefficients → {from_comp} = {stoich[from_comp]}, {to_comp} = {stoich[to_comp]}",
        f"Step 2: {val:.4f} mol {from_comp} × ({stoich[to_comp]}/{stoich[from_comp]}) = {result:.4f} mol {to_comp}",
        f"Final Answer: {result:.4f} mol {to_comp}"
    ]
    return steps

# Compound-to-element mole conversion
def moles_compound_to_element(val, compound_formula, element):
    """
    Convert moles of a compound to moles of a specific element within that compound.
    """
    # Parse formula to get element counts
    counts = parse_formula(compound_formula)
    if element not in counts:
        return [f"Error: element '{element}' not found in compound {compound_formula}."]
    atom_count = counts[element]
    result = val * atom_count
    steps = [
        f"Step 1: Compound formula {compound_formula} contains {atom_count} {element} atom(s) per formula unit",
        f"Step 2: {val:.4f} mol {compound_formula} × {atom_count} = {result:.4f} mol {element}",
        f"Final Answer: {result:.4f} mol {element}"
    ]
    return steps

# Unit conversions: mL↔L and g↔kg

# Unit conversions: yards ↔ centimeters, yards ↔ meters
def yards_to_centimeters(yards):
    cm = yards * 91.44
    return [
        f"Step 1: Given {yards} yards",
        "Step 2: Multiply by 91.44 cm/yard",
        f"{yards} × 91.44 = {cm:.2f} cm",
        f"Final Answer: {cm:.2f} cm"
    ]

def centimeters_to_yards(cm):
    yards = cm / 91.44
    return [
        f"Step 1: Given {cm} cm",
        "Step 2: Divide by 91.44 cm/yard",
        f"{cm} ÷ 91.44 = {yards:.4f} yards",
        f"Final Answer: {yards:.4f} yards"
    ]

def yards_to_meters(yards):
    m = yards * 0.9144
    return [
        f"Step 1: Given {yards} yards",
        "Step 2: Multiply by 0.9144 m/yard",
        f"{yards} × 0.9144 = {m:.4f} m",
        f"Final Answer: {m:.4f} m"
    ]

def meters_to_yards(m):
    yards = m / 0.9144
    return [
        f"Step 1: Given {m} meters",
        "Step 2: Divide by 0.9144 m/yard",
        f"{m} ÷ 0.9144 = {yards:.4f} yards",
        f"Final Answer: {yards:.4f} yards"
    ]

def ounces_to_grams(oz):
    g = oz * 28.3495
    return [
        f"Step 1: Given {oz} ounces",
        "Step 2: Multiply by 28.3495 g/oz",
        f"{oz} oz × 28.3495 = {g:.2f} g",
        f"Final Answer: {g:.2f} g"
    ]

def grams_to_ounces(g):
    oz = g / 28.3495
    return [
        f"Step 1: Given {g} grams",
        "Step 2: Divide by 28.3495 g/oz",
        f"{g} g ÷ 28.3495 = {oz:.4f} oz",
        f"Final Answer: {oz:.4f} oz"
    ]

# Length conversions: inches ↔ millimeters
def inches_to_mm(inches):
    mm = inches * 25.4
    return [
        f"Step 1: Given {inches} inches",
        "Step 2: Multiply by 25.4 mm/inch",
        f"{inches} in × 25.4 = {mm:.4f} mm",
        f"Final Answer: {mm:.4f} mm"
    ]

def mm_to_inches(mm):
    inches = mm / 25.4
    return [
        f"Step 1: Given {mm} mm",
        "Step 2: Divide by 25.4 mm/inch",
        f"{mm} mm ÷ 25.4 = {inches:.4f} in",
        f"Final Answer: {inches:.4f} in"
    ]

# Volume conversions: mL ↔ dm³
def ml_to_dm3(ml):
    dm3 = ml / 1000
    return [
        f"Step 1: Given {ml} mL",
        "Step 2: Divide by 1000 to convert mL → dm³",
        f"{ml} mL ÷ 1000 = {dm3:.4f} dm³",
        f"Final Answer: {dm3:.4f} dm³"
    ]

def dm3_to_ml(dm3):
    ml = dm3 * 1000
    return [
        f"Step 1: Given {dm3} dm³",
        "Step 2: Multiply by 1000 to convert dm³ → mL",
        f"{dm3} dm³ × 1000 = {ml:.1f} mL",
        f"Final Answer: {ml:.1f} mL"
    ]

# Mass conversions: grams ↔ metric tons
def grams_to_tons(g):
    tons = g / 1e6
    return [
        f"Step 1: Given {g} grams",
        "Step 2: Divide by 1e6 to convert g → metric tons",
        f"{g} g ÷ 1e6 = {tons:.6f} tons",
        f"Final Answer: {tons:.6f} tons"
    ]

def tons_to_grams(tons):
    g = tons * 1e6
    return [
        f"Step 1: Given {tons} tons",
        "Step 2: Multiply by 1e6 to convert tons → g",
        f"{tons} tons × 1e6 = {g:.0f} g",
        f"Final Answer: {g:.0f} g"
    ]

# Time conversions: minutes ↔ years
def minutes_to_years(mins):
    years = mins / 525600
    return [
        f"Step 1: Given {mins} minutes",
        "Step 2: Divide by 525600 to convert minutes → years",
        f"{mins} min ÷ 525600 = {years:.6f} years",
        f"Final Answer: {years:.6f} years"
    ]

def years_to_minutes(years):
    mins = years * 525600
    return [
        f"Step 1: Given {years} years",
        "Step 2: Multiply by 525600 to convert years → minutes",
        f"{years} years × 525600 = {mins:.0f} minutes",
        f"Final Answer: {mins:.0f} minutes"
    ]

# Volume to mass: quarts ↔ pounds (water)
def quarts_to_pounds(qts):
    lbs = qts * 2.08645
    return [
        f"Step 1: Given {qts} quarts",
        "Step 2: Multiply by 2.08645 lb/quart",
        f"{qts} qt × 2.08645 = {lbs:.2f} lb",
        f"Final Answer: {lbs:.2f} lb"
    ]

def pounds_to_quarts(lbs):
    qts = lbs / 2.08645
    return [
        f"Step 1: Given {lbs} pounds",
        "Step 2: Divide by 2.08645 lb/quart",
        f"{lbs} lb ÷ 2.08645 = {qts:.4f} qt",
        f"Final Answer: {qts:.4f} qt"
    ]
def ml_to_liter(ml):
    liters = ml / 1000
    steps = [
        f"Step 1: Given {ml} mL",
        "Step 2: Divide by 1000 to convert mL → L",
        f"{ml} mL ÷ 1000 = {liters:.3f} L",
        f"Final Answer: {liters:.3f} L"
    ]
    return steps

def g_to_kg(grams):
    kilograms = grams / 1000
    steps = [
        f"Step 1: Given {grams} g",
        "Step 2: Divide by 1000 to convert g → kg",
        f"{grams} g ÷ 1000 = {kilograms:.3f} kg",
        f"Final Answer: {kilograms:.3f} kg"
    ]
    return steps

# Distance unit conversions
def miles_to_kilometers(miles):
    km = miles * 1.60934
    return [
        f"Step 1: Given {miles} miles",
        f"Step 2: Multiply by 1.60934 km/mile",
        f"{miles} miles × 1.60934 = {km:.4f} km",
        f"Final Answer: {km:.4f} km"
    ]

def kilometers_to_miles(km):
    miles = km / 1.60934
    return [
        f"Step 1: Given {km} kilometers",
        f"Step 2: Divide by 1.60934 km/mile",
        f"{km} km ÷ 1.60934 = {miles:.4f} miles",
        f"Final Answer: {miles:.4f} miles"
    ]

def miles_to_feet(miles):
    feet = miles * 5280
    return [
        f"Step 1: Given {miles} miles",
        f"Step 2: Multiply by 5280 ft/mile",
        f"{miles} miles × 5280 = {feet:.2f} ft",
        f"Final Answer: {feet:.2f} ft"
    ]

def feet_to_miles(feet):
    miles = feet / 5280
    return [
        f"Step 1: Given {feet} feet",
        f"Step 2: Divide by 5280 ft/mile",
        f"{feet} ft ÷ 5280 = {miles:.4f} miles",
        f"Final Answer: {miles:.4f} miles"
    ]

# Volume unit conversions: gallons ↔ cubic meters
def gallons_to_cubic_meters(gallons):
    m3 = gallons * 0.00378541
    return [
        f"Step 1: Given {gallons} gallons",
        "Step 2: Multiply by 0.00378541 m³/gallon",
        f"{gallons} gal × 0.00378541 = {m3:.6f} m³",
        f"Final Answer: {m3:.6f} m³"
    ]

def cubic_meters_to_gallons(m3):
    gal = m3 / 0.00378541
    return [
        f"Step 1: Given {m3} cubic meters",
        "Step 2: Divide by 0.00378541 m³/gallon",
        f"{m3} m³ ÷ 0.00378541 = {gal:.4f} gal",
        f"Final Answer: {gal:.4f} gal"
    ]

# Time unit conversions: seconds↔minutes, minutes↔hours, hours↔days
def seconds_to_minutes(sec):
    minutes = sec / 60
    return [
        f"Step 1: Given {sec} seconds",
        "Step 2: Divide by 60 to convert seconds → minutes",
        f"{sec} s ÷ 60 = {minutes:.4f} min",
        f"Final Answer: {minutes:.4f} min"
    ]

def minutes_to_seconds(mins):
    sec = mins * 60
    return [
        f"Step 1: Given {mins} minutes",
        "Step 2: Multiply by 60 to convert minutes → seconds",
        f"{mins} min × 60 = {sec:.2f} s",
        f"Final Answer: {sec:.2f} s"
    ]

def minutes_to_hours(mins):
    hrs = mins / 60
    return [
        f"Step 1: Given {mins} minutes",
        "Step 2: Divide by 60 to convert minutes → hours",
        f"{mins} min ÷ 60 = {hrs:.4f} hr",
        f"Final Answer: {hrs:.4f} hr"
    ]

def hours_to_minutes(hrs):
    mins = hrs * 60
    return [
        f"Step 1: Given {hrs} hours",
        "Step 2: Multiply by 60 to convert hours → minutes",
        f"{hrs} hr × 60 = {mins:.2f} min",
        f"Final Answer: {mins:.2f} min"
    ]

def hours_to_days(hrs):
    days = hrs / 24
    return [
        f"Step 1: Given {hrs} hours",
        "Step 2: Divide by 24 to convert hours → days",
        f"{hrs} hr ÷ 24 = {days:.4f} day",
        f"Final Answer: {days:.4f} day"
    ]

def days_to_hours(days):
    hrs = days * 24
    return [
        f"Step 1: Given {days} days",
        "Step 2: Multiply by 24 to convert days → hours",
        f"{days} day × 24 = {hrs:.2f} hr",
        f"Final Answer: {hrs:.2f} hr"
    ]

# Rate conversion: minutes per day → hours per day
def min_per_day_to_hr_per_day(min_per_day):
    hr_per_day = min_per_day / 60
    steps = [
        f"Step 1: Given {min_per_day} minutes/day",
        "Step 2: Divide by 60 to convert minutes to hours",
        f"{min_per_day} ÷ 60 = {hr_per_day:.4f} hours/day",
        f"Final Answer: {hr_per_day:.4f} hours/day"
    ]
    return steps

# Time unit conversions: seconds ↔ years
SECONDS_PER_YEAR = 365 * 24 * 3600

def seconds_to_years(sec):
    years = sec / SECONDS_PER_YEAR
    steps = [
        f"Step 1: Given {sec} seconds",
        f"Step 2: Divide by {SECONDS_PER_YEAR} to convert seconds → years",
        f"{sec} s ÷ {SECONDS_PER_YEAR} = {years:.6f} years",
        f"Final Answer: {years:.6f} years"
    ]
    return steps

def years_to_seconds(years):
    sec = years * SECONDS_PER_YEAR
    steps = [
        f"Step 1: Given {years} years",
        f"Step 2: Multiply by {SECONDS_PER_YEAR} to convert years → seconds",
        f"{years} years × {SECONDS_PER_YEAR} = {sec:.2f} seconds",
        f"Final Answer: {sec:.2f} seconds"
    ]
    return steps

# Distance from time & speed: seconds at m/s to miles
def seconds_speed_to_miles(seconds, speed_m_s):
    meters = speed_m_s * seconds
    km = meters / 1000
    miles = km / 1.60934
    steps = [
        f"Step 1: Given {seconds} s at {speed_m_s} m/s",
        f"Step 2: Distance = speed × time = {speed_m_s} m/s × {seconds} s = {meters:.2f} m",
        f"Step 3: Convert to kilometers: {meters:.2f} m ÷ 1000 = {km:.4f} km",
        f"Step 4: Convert to miles: {km:.4f} km ÷ 1.60934 = {miles:.4f} miles",
        f"Final Answer: {miles:.4f} miles"
    ]
    return steps

# Distance-speed to time: miles at m/s to seconds
def miles_speed_to_seconds(distance_miles, speed_m_s):
    """
    Convert distance in miles at a given speed (m/s) to time in seconds.
    """
    # 1) Convert miles to meters
    meters = distance_miles * 1609.34
    # 2) Time = distance ÷ speed
    sec = meters / speed_m_s
    steps = [
        f"Step 1: Given {distance_miles} miles at {speed_m_s} m/s",
        f"Step 2: Convert miles to meters: {distance_miles} × 1609.34 = {meters:.2f} m",
        f"Step 3: Time = distance ÷ speed = {meters:.2f} m ÷ {speed_m_s:.2e} m/s = {sec:.2e} s",
        f"Final Answer: {sec:.2e} s"
    ]
    return steps

# Unit conversions: feet ↔ kilometers
def feet_to_kilometers(feet):
    km = feet * 0.0003048
    return [
        f"Step 1: Given {feet} feet",
        f"Step 2: Multiply by 0.0003048 km/ft",
        f"{feet} ft × 0.0003048 = {km:.4f} km",
        f"Final Answer: {km:.4f} km"
    ]

def kilometers_to_feet(km):
    feet = km / 0.0003048
    return [
        f"Step 1: Given {km} kilometers",
        f"Step 2: Divide by 0.0003048 km/ft",
        f"{km} km ÷ 0.0003048 = {feet:.2f} ft",
        f"Final Answer: {feet:.2f} ft"
    ]

# Time conversions: minutes ↔ days
def minutes_to_days(minutes):
    days = minutes / 1440
    return [
        f"Step 1: Given {minutes} minutes",
        "Step 2: Divide by 1440 to convert minutes → days",
        f"{minutes} min ÷ 1440 = {days:.4f} days",
        f"Final Answer: {days:.4f} days"
    ]

def days_to_minutes(days):
    minutes = days * 1440
    return [
        f"Step 1: Given {days} days",
        "Step 2: Multiply by 1440 to convert days → minutes",
        f"{days} days × 1440 = {minutes:.2f} minutes",
        f"Final Answer: {minutes:.2f} minutes"
    ]

# Time conversions: weeks ↔ years
def weeks_to_years(weeks):
    years = weeks * 7 / 365
    return [
        f"Step 1: Given {weeks} weeks",
        "Step 2: Convert weeks to days: × 7",
        "Step 3: Divide by 365 to convert days → years",
        f"{weeks} × 7 ÷ 365 = {years:.4f} years",
        f"Final Answer: {years:.4f} years"
    ]

def years_to_weeks(years):
    weeks = years * 365 / 7
    return [
        f"Step 1: Given {years} years",
        "Step 2: Multiply by 365 to convert years → days",
        "Step 3: Divide by 7 to convert days → weeks",
        f"{years} × 365 ÷ 7 = {weeks:.2f} weeks",
        f"Final Answer: {weeks:.2f} weeks"
    ]

# Volume conversions: gallons → liters
def gallons_to_liters(gallons):
    liters = gallons * 3.78541
    return [
        f"Step 1: Given {gallons} gallons",
        "Step 2: Multiply by 3.78541 L/gallon",
        f"{gallons} gal × 3.78541 = {liters:.4f} L",
        f"Final Answer: {liters:.4f} L"
    ]

# Drops conversion factor
DROPS_PER_ML = 20

def ml_to_drops(ml):
    drops = ml * DROPS_PER_ML
    return [
        f"Step 1: Given {ml} mL",
        f"Step 2: 1 mL = {DROPS_PER_ML} drops",
        f"{ml} mL × {DROPS_PER_ML} = {drops:.0f} drops",
        f"Final Answer: {drops:.0f} drops"
    ]

def gallons_to_gallons_per_second(gallons_per_day):
    # Convert flow rate from gallons/day to gallons/second
    gal_per_s = gallons_per_day / 86400
    return [
        f"Step 1: Given {gallons_per_day} gallons/day",
        f"Step 2: Divide by number of seconds in a day (86400) to convert to gal/s",
        f"{gallons_per_day} gal/day ÷ 86400 s/day = {gal_per_s:.6f} gal/s",
        f"Final Answer: {gal_per_s:.6f} gal/s"
    ]


def cubic_meters_to_drops(m3):
    ml = m3 * 1e6
    drops = ml * DROPS_PER_ML
    return [
        f"Step 1: Given {m3} m³",
        f"Step 2: Convert to mL: {m3} × 10^6 = {ml:.0f} mL",
        f"Step 3: 1 mL = {DROPS_PER_ML} drops, so {ml:.0f} mL × {DROPS_PER_ML} = {drops:.0f} drops",
        f"Final Answer: {drops:.0f} drops"
    ]

# Drops → Cubic meters helper
def drops_to_cubic_meters(drops):
    """
    Convert drops to cubic meters.
    """
    # Convert drops to mL
    ml = drops / DROPS_PER_ML
    # Convert mL to cubic meters
    m3 = ml / 1e6
    steps = [
        f"Step 1: Given {drops:.0f} drops",
        f"Step 2: Convert drops to mL: {drops:.0f} drops ÷ {DROPS_PER_ML} drops/mL = {ml:.4f} mL",
        f"Step 3: Convert mL to m³: {ml:.4f} mL ÷ 1e6 = {m3:.10f} m³",
        f"Final Answer: {m3:.10f} m³"
    ]
    return steps

# Supply duration calculation (total supply and daily rate)
def days_supply_from_rate(total_gallons, daily_rate_gpd):
    days = total_gallons / daily_rate_gpd
    steps = [
        f"Step 1: Total supply = {total_gallons} gallons",
        f"Step 2: Daily usage rate = {daily_rate_gpd} gallons/day",
        f"Step 3: Days supply = {total_gallons} ÷ {daily_rate_gpd} = {days:.2f} days",
        f"Final Answer: {days:.2f} days"
    ]
    return steps

# Inverse supply calculation: gallons/day × days → total gallons
def total_supply_from_rate_and_days(daily_rate_gpd, days):
    total = daily_rate_gpd * days
    steps = [
        f"Step 1: Given daily usage rate = {daily_rate_gpd} gallons/day",
        f"Step 2: Given duration = {days} days",
        f"Step 3: Total supply = {daily_rate_gpd} × {days} = {total:.2f} gallons",
        f"Final Answer: {total:.2f} gallons"
    ]
    return steps

# Inverse unit conversions

def liters_to_gallons(liters):
    gal = liters / 3.78541
    return [
        f"Step 1: Given {liters} liters",
        "Step 2: Divide by 3.78541 L/gallon",
        f"{liters} L ÷ 3.78541 = {gal:.4f} gallons",
        f"Final Answer: {gal:.4f} gallons"
    ]

def drops_to_milliliters(drops):
    ml = drops / DROPS_PER_ML
    return [
        f"Step 1: Given {drops} drops",
        f"Step 2: Convert drops to mL: {drops} ÷ {DROPS_PER_ML} = {ml:.4f} mL",
        f"Final Answer: {ml:.4f} mL"
    ]

def gallons_per_second_to_gallons_per_day(gal_per_s):
    galpd = gal_per_s * 86400
    return [
        f"Step 1: Given {gal_per_s} gallons/second",
        "Step 2: Multiply by 86400 s/day to convert to gallons/day",
        f"{gal_per_s} gal/s × 86400 = {galpd:.2f} gallons/day",
        f"Final Answer: {galpd:.2f} gallons/day"
    ]

def hours_per_day_to_minutes_per_day(hrpd):
    minpd = hrpd * 60
    return [
        f"Step 1: Given {hrpd} hours/day",
        "Step 2: Multiply by 60 to convert hours → minutes",
        f"{hrpd} hr/day × 60 = {minpd:.2f} minutes/day",
        f"Final Answer: {minpd:.2f} minutes/day"
    ]

 # Consumption time calculation (e.g., dog food rate)
def consumption_time_to_target(target_qty, target_unit, rate_qty, rate_unit, rate_time_s):
    """
    Calculate time to consume target_qty of unit, given rate_qty per rate_time_s (seconds).
    Returns time in seconds and minutes.
    """
    # seconds per unit
    sec_per_unit = rate_time_s / rate_qty
    total_sec = sec_per_unit * target_qty
    minutes = total_sec / 60
    steps = [
        f"Step 1: Consumption rate = {rate_qty} {rate_unit} / {rate_time_s} s",
        f"Step 2: Time per {rate_unit} = {rate_time_s} s ÷ {rate_qty} = {sec_per_unit:.4f} s/{rate_unit}",
        f"Step 3: Total time = {sec_per_unit:.4f} s/{rate_unit} × {target_qty} {target_unit} = {total_sec:.4f} s",
        f"Step 4: Convert to minutes: {total_sec:.4f} s ÷ 60 = {minutes:.4f} min",
        f"Final Answer: {minutes:.4f} min"
    ]
    return steps

# Volume to mass using density for liters
def volume_to_mass(volume_L, density_g_per_mL):
    """
    Convert volume in liters to mass in grams using density (g/mL).
    """
    mL = volume_L * 1000
    mass_g = mL * density_g_per_mL
    steps = [
        f"Step 1: Given density = {density_g_per_mL} g/mL",
        f"Step 2: Convert volume to mL: {volume_L} L × 1000 = {mL:.1f} mL",
        f"Step 3: Mass = {mL:.1f} mL × {density_g_per_mL} g/mL = {mass_g:.2f} g",
        f"Final Answer: {mass_g:.2f} g"
    ]
    return steps

# Mass → volume in liters using density (g/mL)
def mass_to_liters(grams, density_g_per_mL):
    """
    Convert mass in grams to volume in liters using density (g/mL).
    """
    mL = grams / density_g_per_mL
    L = mL / 1000
    steps = [
        f"Step 1: Given {grams} g and density = {density_g_per_mL} g/mL",
        f"Step 2: Volume in mL = {grams} g ÷ {density_g_per_mL} g/mL = {mL:.2f} mL",
        f"Step 3: Convert to liters: {mL:.2f} mL ÷ 1000 = {L:.4f} L",
        f"Final Answer: {L:.4f} L"
    ]
    return steps

def dimensional_analysis(query):
    import re
    query = query.strip()
    print(f"DEBUG: Received query: '{query}'")

    # Yards ↔ Centimeters
    match_yd_to_cm = re.match(r'convert\s+([\deE.+-]+)\s*yards?\s+to\s+centimeters?', query, re.IGNORECASE)
    if match_yd_to_cm:
        return '\n'.join(yards_to_centimeters(float(match_yd_to_cm.group(1))))

    match_cm_to_yd = re.match(r'convert\s+([\deE.+-]+)\s*centimeters?\s+to\s+yards?', query, re.IGNORECASE)
    if match_cm_to_yd:
        return '\n'.join(centimeters_to_yards(float(match_cm_to_yd.group(1))))

    # Yards ↔ Meters
    match_yd_to_m = re.match(r'convert\s+([\deE.+-]+)\s*yards?\s+to\s+meters?', query, re.IGNORECASE)
    if match_yd_to_m:
        return '\n'.join(yards_to_meters(float(match_yd_to_m.group(1))))

    match_m_to_yd = re.match(r'convert\s+([\deE.+-]+)\s*meters?\s+to\s+yards?', query, re.IGNORECASE)
    if match_m_to_yd:
        return '\n'.join(meters_to_yards(float(match_m_to_yd.group(1))))

    # Inches ↔ Millimeters
    match_in_to_mm = re.match(
        r'convert\s+([\deE.+-]+)\s*inches?\s+to\s+millimeters?',
        query, re.IGNORECASE
    )
    if match_in_to_mm:
        return '\n'.join(inches_to_mm(float(match_in_to_mm.group(1))))

    match_mm_to_in = re.match(
        r'convert\s+([\deE.+-]+)\s*millimeters?\s+to\s+inches?',
        query, re.IGNORECASE
    )
    if match_mm_to_in:
        return '\n'.join(mm_to_inches(float(match_mm_to_in.group(1))))

    # mL ↔ dm³
    match_ml_to_dm3 = re.match(
        r'convert\s+([\deE.+-]+)\s*(?:mL|ml|milliliters?)\s+to\s+dm3|dm³',
        query, re.IGNORECASE
    )
    if match_ml_to_dm3:
        return '\n'.join(ml_to_dm3(float(match_ml_to_dm3.group(1))))

    match_dm3_to_ml = re.match(
        r'convert\s+([\deE.+-]+)\s*(?:dm3|dm³)\s+to\s+(?:mL|ml|milliliters?)',
        query, re.IGNORECASE
    )
    if match_dm3_to_ml:
        return '\n'.join(dm3_to_ml(float(match_dm3_to_ml.group(1))))

    # Ounces ↔ Grams
    match_oz_to_g = re.match(r'convert\s+([\deE.+-]+)\s*ounces?\s+to\s+grams?', query, re.IGNORECASE)
    if match_oz_to_g:
        return '\n'.join(ounces_to_grams(float(match_oz_to_g.group(1))))

    match_g_to_oz = re.match(r'convert\s+([\deE.+-]+)\s*grams?\s+to\s+ounces?', query, re.IGNORECASE)
    if match_g_to_oz:
        return '\n'.join(grams_to_ounces(float(match_g_to_oz.group(1))))

    # Grams ↔ Tons
    match_g_to_ton = re.match(r'convert\s+([\deE.+-]+)\s*grams?\s+to\s+tons?', query, re.IGNORECASE)
    if match_g_to_ton:
        return '\n'.join(grams_to_tons(float(match_g_to_ton.group(1))))

    match_ton_to_g = re.match(r'convert\s+([\deE.+-]+)\s*tons?\s+to\s+grams?', query, re.IGNORECASE)
    if match_ton_to_g:
        return '\n'.join(tons_to_grams(float(match_ton_to_g.group(1))))

    # Minutes ↔ Years
    match_min_to_yr = re.match(r'convert\s+([\deE.+-]+)\s*minutes?\s+to\s+years?', query, re.IGNORECASE)
    if match_min_to_yr:
        return '\n'.join(minutes_to_years(float(match_min_to_yr.group(1))))

    match_yr_to_min = re.match(r'convert\s+([\deE.+-]+)\s*years?\s+to\s+minutes?', query, re.IGNORECASE)
    if match_yr_to_min:
        return '\n'.join(years_to_minutes(float(match_yr_to_min.group(1))))

    # Quarts ↔ Pounds
    match_qt_to_lb = re.match(r'convert\s+([\deE.+-]+)\s*quarts?\s+to\s+pounds?', query, re.IGNORECASE)
    if match_qt_to_lb:
        return '\n'.join(quarts_to_pounds(float(match_qt_to_lb.group(1))))

    match_lb_to_qt = re.match(r'convert\s+([\deE.+-]+)\s*pounds?\s+to\s+quarts?', query, re.IGNORECASE)
    if match_lb_to_qt:
        return '\n'.join(pounds_to_quarts(float(match_lb_to_qt.group(1))))

    # Quarts ↔ Pounds via density
    match_qt_to_lb_den = re.match(
        r'convert\s+([\deE.+-]+)\s*quarts?\s+to\s+pounds?\s+using\s+density\s+([\deE.+-]+)\s*kg/quart',
        query, re.IGNORECASE
    )
    if match_qt_to_lb_den:
        qts, dens = match_qt_to_lb_den.groups()
        return '\n'.join(quarts_to_pounds_using_density(float(qts), float(dens)))

    match_lb_to_qt_den = re.match(
        r'convert\s+([\deE.+-]+)\s*pounds?\s+to\s+quarts?\s+using\s+density\s+([\deE.+-]+)\s*kg/quart',
        query, re.IGNORECASE
    )
    if match_lb_to_qt_den:
        lbs, dens = match_lb_to_qt_den.groups()
        return '\n'.join(pounds_to_quarts_using_density(float(lbs), float(dens)))

    # Time ↔ Months
    match_mon_to_s = re.match(
        r'convert\s+([\deE.+-]+)\s*months?\s+to\s+seconds?',
        query, re.IGNORECASE
    )
    if match_mon_to_s:
        return '\n'.join(months_to_seconds(float(match_mon_to_s.group(1))))

    match_s_to_mon = re.match(
        r'convert\s+([\deE.+-]+)\s*seconds?\s+to\s+months?',
        query, re.IGNORECASE
    )
    if match_s_to_mon:
        return '\n'.join(seconds_to_months(float(match_s_to_mon.group(1))))

    # Fuel consumption: miles ↔ gallons using rate
    match_mi_to_gal_rate = re.match(
        r'convert\s+([\deE.+-]+)\s*miles?\s+to\s+gallons?\s+using\s+([\deE.+-]+)\s*gal\s+per\s+([\deE.+-]+)\s*km',
        query, re.IGNORECASE
    )
    if match_mi_to_gal_rate:
        mi, gal, km = map(float, match_mi_to_gal_rate.groups())
        return '\n'.join(miles_to_gallons_using_rate(mi, gal, km))

    match_gal_to_mi_rate = re.match(
        r'convert\s+([\deE.+-]+)\s*gallons?\s+to\s+miles?\s+using\s+([\deE.+-]+)\s*gal\s+per\s+([\deE.+-]+)\s*km',
        query, re.IGNORECASE
    )
    if match_gal_to_mi_rate:
        gal, gal1, km = map(float, match_gal_to_mi_rate.groups())
        return '\n'.join(gallons_to_miles_using_rate(gal, gal1, km))

    # Weeks ↔ Minutes
    match_w_to_min = re.match(
        r'convert\s+([\deE.+-]+)\s*weeks?\s+to\s+minutes?',
        query, re.IGNORECASE
    )
    if match_w_to_min:
        return '\n'.join(weeks_to_minutes(float(match_w_to_min.group(1))))

    match_min_to_w = re.match(
        r'convert\s+([\deE.+-]+)\s*minutes?\s+to\s+weeks?',
        query, re.IGNORECASE
    )
    if match_min_to_w:
        return '\n'.join(minutes_to_weeks(float(match_min_to_w.group(1))))

    # Distance-time ↔ ft/s
    match_mi_time_to_fps = re.match(
        r'convert\s+([\deE.+-]+)\s*miles?\s+in\s+([\deE.+-]+)\s*seconds?\s+to\s+feet/second',
        query, re.IGNORECASE
    )
    if match_mi_time_to_fps:
        mi, sec = map(float, match_mi_time_to_fps.groups())
        return '\n'.join(miles_time_to_fps(mi, sec))

    match_fps_to_time = re.match(
        r'convert\s+([\deE.+-]+)\s*ft/s\s+to\s+seconds?',
        query, re.IGNORECASE
    )
    if match_fps_to_time:
        fps = float(match_fps_to_time.group(1))
        return '\n'.join(fps_to_miles_time(fps))

    # Convert mL to L
    match_ml_to_L = re.match(
        r'convert\s+([\deE.+-]+)\s*(?:mL|ml|milliliters?)\s+to\s*(?:L|l|liters?)$',
        query, re.IGNORECASE
    )
    if match_ml_to_L:
        ml = float(match_ml_to_L.group(1))
        return '\n'.join(ml_to_liter(ml))

    # Convert g to kg
    match_g_to_kg = re.match(
        r'convert\s+([\deE.+-]+)\s*(?:g|grams?)\s+to\s*(?:kg|kilograms?)$',
        query, re.IGNORECASE
    )
    if match_g_to_kg:
        grams = float(match_g_to_kg.group(1))
        return '\n'.join(g_to_kg(grams))

    # Distance conversions
    match_mi_to_km = re.match(
        r'convert\s+([\deE.+-]+)\s*miles?\s+to\s+kilometers?',
        query, re.IGNORECASE
    )
    if match_mi_to_km:
        mi = float(match_mi_to_km.group(1))
        return '\n'.join(miles_to_kilometers(mi))

    match_km_to_mi = re.match(
        r'convert\s+([\deE.+-]+)\s*kilometers?\s+to\s+miles?',
        query, re.IGNORECASE
    )
    if match_km_to_mi:
        km = float(match_km_to_mi.group(1))
        return '\n'.join(kilometers_to_miles(km))

    match_mi_to_ft = re.match(
        r'convert\s+([\deE.+-]+)\s*miles?\s+to\s+feet?',
        query, re.IGNORECASE
    )
    if match_mi_to_ft:
        mi = float(match_mi_to_ft.group(1))
        return '\n'.join(miles_to_feet(mi))

    match_ft_to_mi = re.match(
        r'convert\s+([\deE.+-]+)\s*feet?\s+to\s+miles?',
        query, re.IGNORECASE
    )
    if match_ft_to_mi:
        ft = float(match_ft_to_mi.group(1))
        return '\n'.join(feet_to_miles(ft))

    # Volume conversions: gallons ↔ cubic meters
    match_gal_to_m3 = re.match(
        r'convert\s+([\deE.+-]+)\s*gallons?\s+to\s+cubic\s+meters?',
        query, re.IGNORECASE
    )
    if match_gal_to_m3:
        gal = float(match_gal_to_m3.group(1))
        return '\n'.join(gallons_to_cubic_meters(gal))

    match_m3_to_gal = re.match(
        r'convert\s+([\deE.+-]+)\s*cubic\s+meters?\s+to\s+gallons?',
        query, re.IGNORECASE
    )
    if match_m3_to_gal:
        m3 = float(match_m3_to_gal.group(1))
        return '\n'.join(cubic_meters_to_gallons(m3))

    # Additional liquid/unit conversions
    match_gal_to_Lit = re.match(
        r'convert\s+([\deE.+-]+)\s*gallons?\s+to\s+liters?',
        query, re.IGNORECASE
    )
    if match_gal_to_Lit:
        gal = float(match_gal_to_Lit.group(1))
        return '\n'.join(gallons_to_liters(gal))

    match_ml_to_drops = re.match(
        r'convert\s+([\deE.+-]+)\s*(?:mL|ml|milliliters?)\s+to\s+drops?',
        query, re.IGNORECASE
    )
    if match_ml_to_drops:
        ml = float(match_ml_to_drops.group(1))
        return '\n'.join(ml_to_drops(ml))

    match_galpd_to_galps = re.match(
        r'convert\s+([\deE.+-]+)\s*gallons?/day\s+to\s+gallons?/second',
        query, re.IGNORECASE
    )
    if match_galpd_to_galps:
        galpd = float(match_galpd_to_galps.group(1))
        return '\n'.join(gallons_to_gallons_per_second(galpd))


    # Handle drops → cubic meters
    match_drops_to_m3 = re.match(
        r'convert\s+([\deE.+-]+)\s*drops?\s+to\s+cubic\s+meters?',
        query, re.IGNORECASE
    )
    if match_drops_to_m3:
        drops = float(match_drops_to_m3.group(1))
        return '\n'.join(drops_to_cubic_meters(drops))

    match_m3_to_drops = re.match(
        r'convert\s+([\deE.+-]+)\s*cubic\s+meters?\s+to\s+drops?',
        query, re.IGNORECASE
    )
    if match_m3_to_drops:
        m3 = float(match_m3_to_drops.group(1))
        return '\n'.join(cubic_meters_to_drops(m3))

    # Time conversions
    match_s_to_min = re.match(
        r'convert\s+([\deE.+-]+)\s*seconds?\s+to\s+minutes?',
        query, re.IGNORECASE
    )
    if match_s_to_min:
        sec = float(match_s_to_min.group(1))
        return '\n'.join(seconds_to_minutes(sec))

    match_min_to_s = re.match(
        r'convert\s+([\deE.+-]+)\s*minutes?\s+to\s+seconds?',
        query, re.IGNORECASE
    )
    if match_min_to_s:
        mins = float(match_min_to_s.group(1))
        return '\n'.join(minutes_to_seconds(mins))

    match_min_to_hr = re.match(
        r'convert\s+([\deE.+-]+)\s*minutes?\s+to\s+hours?',
        query, re.IGNORECASE
    )
    if match_min_to_hr:
        mins = float(match_min_to_hr.group(1))
        return '\n'.join(minutes_to_hours(mins))

    match_hr_to_min = re.match(
        r'convert\s+([\deE.+-]+)\s*hours?\s+to\s+minutes?',
        query, re.IGNORECASE
    )
    if match_hr_to_min:
        hrs = float(match_hr_to_min.group(1))
        return '\n'.join(hours_to_minutes(hrs))

    match_hr_to_day = re.match(
        r'convert\s+([\deE.+-]+)\s*hours?\s+to\s+days?',
        query, re.IGNORECASE
    )
    if match_hr_to_day:
        hrs = float(match_hr_to_day.group(1))
        return '\n'.join(hours_to_days(hrs))

    match_day_to_hr = re.match(
        r'convert\s+([\deE.+-]+)\s*days?\s+to\s+hours?',
        query, re.IGNORECASE
    )
    if match_day_to_hr:
        days = float(match_day_to_hr.group(1))
        return '\n'.join(days_to_hours(days))

    # Rate conversion: minutes/day → hours/day
    match_minpd_to_hrpd = re.match(
        r'convert\s+([\deE.+-]+)\s*minutes?/day\s+to\s*hours?/day',
        query, re.IGNORECASE
    )
    if match_minpd_to_hrpd:
        rate = float(match_minpd_to_hrpd.group(1))
        return '\n'.join(min_per_day_to_hr_per_day(rate))

    # Time conversions: seconds ↔ years
    match_s_to_yr = re.match(
        r'convert\s+([\deE.+-]+)\s*seconds?\s+to\s*years?',
        query, re.IGNORECASE
    )
    if match_s_to_yr:
        sec = float(match_s_to_yr.group(1))
        return '\n'.join(seconds_to_years(sec))

    match_yr_to_s = re.match(
        r'convert\s+([\deE.+-]+)\s*years?\s+to\s*seconds?',
        query, re.IGNORECASE
    )
    if match_yr_to_s:
        yrs = float(match_yr_to_s.group(1))
        return '\n'.join(years_to_seconds(yrs))

    # Hours/day → Minutes/day
    match_hrpd_to_minpd = re.match(
        r'convert\s+([\deE.+-]+)\s*hours?/day\s+to\s*minutes?/day',
        query, re.IGNORECASE
    )
    if match_hrpd_to_minpd:
        hrpd = float(match_hrpd_to_minpd.group(1))
        return '\n'.join(hours_per_day_to_minutes_per_day(hrpd))

    # Liters → Gallons
    match_L_to_gal = re.match(
        r'convert\s+([\deE.+-]+)\s*(?:L|l|liters?)\s+to\s+gallons?',
        query, re.IGNORECASE
    )
    if match_L_to_gal:
        liters = float(match_L_to_gal.group(1))
        return '\n'.join(liters_to_gallons(liters))

    # Drops → mL
    match_drops_to_ml = re.match(
        r'convert\s+([\deE.+-]+)\s*drops?\s+to\s*(?:mL|ml|milliliters?)',
        query, re.IGNORECASE
    )
    if match_drops_to_ml:
        drops = float(match_drops_to_ml.group(1))
        return '\n'.join(drops_to_milliliters(drops))

    # Gallons/second → Gallons/day
    match_galps_to_galpd = re.match(
        r'convert\s+([\deE.+-]+)\s*gallons?/second\s+to\s+gallons?/day',
        query, re.IGNORECASE
    )
    if match_galps_to_galpd:
        galps = float(match_galps_to_galpd.group(1))
        return '\n'.join(gallons_per_second_to_gallons_per_day(galps))

    # ft/sec → mi/hr (allow ft/s, ft/sec, mi/h, mi/hr)
    match_ftps_to_mihr = re.match(
        r'convert\s+([\deE.+-]+)\s*ft/?s(?:ec(onds?)?)?\s+to\s*mi/?h(?:r|ours?)?',
        query, re.IGNORECASE
    )
    if match_ftps_to_mihr:
        fps_val = float(match_ftps_to_mihr.group(1))
        return '\n'.join(fps_to_mph(fps_val))

    # Inverse supply: daily usage and duration → total gallons
    match_supply_inv = re.search(
        r'uses\s+([\deE.+-]+)\s*gallons?/day\s+for\s+([\deE.+-]+)\s*days',
        query, re.IGNORECASE
    )
    if match_supply_inv:
        rate_str, days_str = match_supply_inv.groups()
        rate = float(rate_str)
        days = float(days_str)
        return '\n'.join(total_supply_from_rate_and_days(rate, days))

    # Supply duration: total gallons and daily usage rate
    match_supply = re.search(
        r'holds\s+([\deE.+-]+)\s*gallons?.*?uses\s+([\deE.+-]+)\s*gallons?/day',
        query, re.IGNORECASE | re.DOTALL
    )
    if match_supply:
        total_str, rate_str = match_supply.groups()
        total = float(total_str)
        rate = float(rate_str)
        return '\n'.join(days_supply_from_rate(total, rate))

    # Generic consumption time: rate_qty rate_unit in rate_time seconds → target_qty target_unit
    match_cons = re.search(
        r'(\d+(?:\.\d+)?)\s*([A-Za-z]+)\s+in\s+(\d+(?:\.\d+)?)\s*seconds.*?(\d+(?:\.\d+)?)\s*([A-Za-z]+)',
        query, re.IGNORECASE
    )
    if match_cons:
        rate_qty, rate_unit, rate_time_s, target_qty, target_unit = match_cons.groups()
        rate_qty = float(rate_qty)
        rate_time_s = float(rate_time_s)
        target_qty = float(target_qty)
        return '\n'.join(consumption_time_to_target(
            target_qty, target_unit, rate_qty, rate_unit, rate_time_s
        ))

    # Zinc block mass (density g/mL, volume in L)
    match_zinc = re.search(
        r'density of zinc is approximately\s*([\deE.+-]+)\s*g/mL.*?mass of a\s*([\deE.+-]+)\s*L block of zinc',
        query, re.IGNORECASE | re.DOTALL
    )
    if match_zinc:
        density, vol_L = match_zinc.groups()
        density = float(density)
        vol_L = float(vol_L)
        return '\n'.join(volume_to_mass(vol_L, density))

    # Distance-speed to time: miles at m/s to seconds
    match_mi_at_speed_to_s = re.match(
        r'convert\s+([\deE.+-]+)\s*miles?\s+at\s+([\deE.+-]+)\s*m/s\s+to\s+seconds?',
        query, re.IGNORECASE
    )
    if match_mi_at_speed_to_s:
        dist, speed = map(float, match_mi_at_speed_to_s.groups())
        return '\n'.join(miles_speed_to_seconds(dist, speed))

    # Time-speed to distance: seconds at m/s to miles
    match_ts_to_mi = re.match(
        r'convert\s+([\deE.+-]+)\s*seconds?\s+at\s+([\deE.+-]+)\s*m/s\s+to\s+miles?',
        query, re.IGNORECASE
    )
    if match_ts_to_mi:
        sec, speed = map(float, match_ts_to_mi.groups())
        return '\n'.join(seconds_speed_to_miles(sec, speed))

    # Feet ↔ Kilometers
    match_ft_to_km = re.match(
        r'convert\s+([\deE.+-]+)\s*feet?\s+to\s+kilometers?',
        query, re.IGNORECASE
    )
    if match_ft_to_km:
        ft = float(match_ft_to_km.group(1))
        return '\n'.join(feet_to_kilometers(ft))

    match_km_to_ft = re.match(
        r'convert\s+([\deE.+-]+)\s*kilometers?\s+to\s+feet?',
        query, re.IGNORECASE
    )
    if match_km_to_ft:
        km = float(match_km_to_ft.group(1))
        return '\n'.join(kilometers_to_feet(km))

    # Minutes ↔ Days
    match_min_to_day = re.match(
        r'convert\s+([\deE.+-]+)\s*minutes?\s+to\s+days?',
        query, re.IGNORECASE
    )
    if match_min_to_day:
        mins = float(match_min_to_day.group(1))
        return '\n'.join(minutes_to_days(mins))

    match_day_to_min = re.match(
        r'convert\s+([\deE.+-]+)\s*days?\s+to\s+minutes?',
        query, re.IGNORECASE
    )
    if match_day_to_min:
        days = float(match_day_to_min.group(1))
        return '\n'.join(days_to_minutes(days))

    # Weeks ↔ Years
    match_w_to_yr = re.match(
        r'convert\s+([\deE.+-]+)\s*weeks?\s+to\s+years?',
        query, re.IGNORECASE
    )
    if match_w_to_yr:
        w = float(match_w_to_yr.group(1))
        return '\n'.join(weeks_to_years(w))

    match_yr_to_w = re.match(
        r'convert\s+([\deE.+-]+)\s*years?\s+to\s+weeks?',
        query, re.IGNORECASE
    )
    if match_yr_to_w:
        y = float(match_yr_to_w.group(1))
        return '\n'.join(years_to_weeks(y))

    # Speed conversions: mph ↔ fps
    match_mph_to_fps = re.match(
        r'convert\s+([\deE.+-]+)\s*miles/hour\s+to\s+feet/second',
        query, re.IGNORECASE
    )
    if match_mph_to_fps:
        return '\n'.join(mph_to_fps(float(match_mph_to_fps.group(1))))

    match_fps_to_mph = re.match(
        r'convert\s+([\deE.+-]+)\s*feet/second\s+to\s+miles/hour',
        query, re.IGNORECASE
    )
    if match_fps_to_mph:
        return '\n'.join(fps_to_mph(float(match_fps_to_mph.group(1))))

    # Weight conversions: lb & oz ↔ mg
    match_lb_oz_to_mg = re.match(
        r'convert\s+([\deE.+-]+)\s*lb(?:s)?\s+and\s+([\deE.+-]+)\s*oz\s+to\s+mg',
        query, re.IGNORECASE
    )
    if match_lb_oz_to_mg:
        lb, oz = map(float, match_lb_oz_to_mg.groups())
        return '\n'.join(lb_oz_to_mg(lb, oz))

    match_mg_to_lb_oz = re.match(
        r'convert\s+([\deE.+-]+)\s*mg\s+to\s+lb(?:s)?\s+and\s+oz',
        query, re.IGNORECASE
    )
    if match_mg_to_lb_oz:
        mg = float(match_mg_to_lb_oz.group(1))
        return '\n'.join(mg_to_lb_oz(mg))

    # Distance conversions: miles ↔ cm
    match_mi_to_cm = re.match(
        r'convert\s+([\deE.+-]+)\s*miles?\s+to\s+centimeters?',
        query, re.IGNORECASE
    )
    if match_mi_to_cm:
        return '\n'.join(miles_to_cm(float(match_mi_to_cm.group(1))))

    match_cm_to_mi = re.match(
        r'convert\s+([\deE.+-]+)\s*centimeters?\s+to\s+miles?',
        query, re.IGNORECASE
    )
    if match_cm_to_mi:
        return '\n'.join(cm_to_miles(float(match_cm_to_mi.group(1))))

    # Pattern for general conversions: convert X unit1 of compound to unit2
    conv_pattern = re.compile(
    r'convert\s+([\deE.+-]+)\s*'
    r'(molecules?|particles?|atoms?|moles?|grams?|liters?|mol|g|l)\s+of\s+'
    r'([A-Za-z0-9().*]+)\s+to\s+'
    r'(molecules?|particles?|atoms?|moles?|grams?|liters?|mol|g|l)',
    re.IGNORECASE
)
    

    m = conv_pattern.fullmatch(query)   # requires that the entire string conforms
    print(f"DEBUG: conv_pattern match: {bool(m)}")
    if m:
        value_str, raw_from, raw_compound, raw_to = m.groups()
        try:
            value = float(value_str)
        except ValueError:
            return "Invalid numeric value."
        # Normalize units to singular lowercase
        unit_map = {
            'mole': 'mole', 'moles': 'mole',
            'mol': 'mole',
            'gram': 'gram', 'grams': 'gram',
            'g': 'gram',
            'molecule': 'particle', 'molecules': 'particle',
            'particle': 'particle', 'particles': 'particle',
            'atom': 'particle', 'atoms': 'particle',
            'liter': 'liter', 'liters': 'liter',
            'l': 'liter',
        }
        from_unit = unit_map.get(raw_from.lower(), None)
        to_unit = unit_map.get(raw_to.lower(), None)
        print(f"DEBUG: Normalized units: from='{from_unit}', to='{to_unit}'")
        if not from_unit or not to_unit:
            return "Unsupported unit."

        # Normalize compound formatting
        compound = fix_case(raw_compound)

        steps = [f"Step 1: Starting with {value} {from_unit}(s) of {compound}"]

        # Branch for each conversion
        if from_unit == 'mole' and to_unit == 'gram':
            try:
                mm = molar_mass(compound)
            except ValueError as e:
                return str(e)
            steps.extend(moles_to_grams(value, mm))
            return '\n'.join(steps)

        elif from_unit == 'gram' and to_unit == 'mole':
            try:
                mm = molar_mass(compound)
            except ValueError as e:
                return str(e)
            steps.extend(grams_to_moles(value, mm))
            return '\n'.join(steps)

        elif from_unit == 'mole' and to_unit == 'particle':
            steps.extend(moles_to_particles(value))
            return '\n'.join(steps)

        elif from_unit == 'particle' and to_unit == 'mole':
            steps.extend(particles_to_moles(value))
            return '\n'.join(steps)

        elif from_unit == 'gram' and to_unit == 'particle':
            try:
                mm = molar_mass(compound)
            except ValueError as e:
                return str(e)
            # 1) grams → moles
            steps.extend(grams_to_moles(value, mm))
            # 2) moles → particles
            moles_value = value / mm
            steps.extend(moles_to_particles(moles_value))
            return '\n'.join(steps)

        elif from_unit == 'particle' and to_unit == 'gram':
            # 1) particles → moles
            moles_value = value / AVOGADRO
            steps.extend(particles_to_moles(value))
            # 2) moles → grams
            try:
                mm = molar_mass(compound)
            except ValueError as e:
                return str(e)
            steps.extend(moles_to_grams(moles_value, mm))
            return '\n'.join(steps)

        elif from_unit == 'mole' and to_unit == 'liter':
            steps.extend(moles_to_liters_stp(value))
            return '\n'.join(steps)

        elif from_unit == 'liter' and to_unit == 'mole':
            steps.extend(liters_to_moles_stp(value))
            return '\n'.join(steps)

        else:
            return "Unsupported conversion type."

    # Molarity from moles and volume
    molar_pattern = re.compile(
        r'(?:molarity\s+from|convert)\s+([\deE.+-]+)\s*mol\s+and\s+([\deE.+-]+)\s*(?:L|l|liters?)',
        re.IGNORECASE
    )
    m2 = molar_pattern.match(query)
    if m2:
        mol, vol = map(float, m2.groups())
        return '\n'.join(molarity_from_moles_and_volume(mol, vol))

    # Percent yield
    yield_pattern = re.compile(
        r'percent\s+yield\s+from\s+([\deE.+-]+)\s*g\s+actual\s+and\s+([\deE.+-]+)\s*g\s+theoretical',
        re.IGNORECASE
    )
    m3 = yield_pattern.match(query)
    if m3:
        actual, theoretical = map(float, m3.groups())
        return '\n'.join(percent_yield(actual, theoretical))

    # Handle density: volume → mass
    match_density_vol_to_mass = re.match(
        r'convert\s+([\deE.+-]+)\s*(?:mL|ml|milliliters?)\s+of\s+[A-Za-z0-9().*]+\s+to\s+grams\s+using\s+density\s+([\deE.+-]+)\s*g/mL',
        query, re.IGNORECASE
    )
    if match_density_vol_to_mass:
        volume_str, density_str = match_density_vol_to_mass.groups()
        volume = float(volume_str)
        density = float(density_str)
        return '\n'.join(volume_to_grams(volume, density))

    # Handle density: mass → volume
    match_density_mass_to_vol = re.match(
        r'convert\s+([\deE.+-]+)\s*grams?\s+of\s+[A-Za-z0-9().*]+\s+to\s+volume\s+using\s+density\s+([\deE.+-]+)\s*g/mL',
        query, re.IGNORECASE
    )
    if match_density_mass_to_vol:
        grams_str, density_str = match_density_mass_to_vol.groups()
        grams = float(grams_str)
        density = float(density_str)
        return '\n'.join(grams_to_volume(grams, density))

    # Mass → volume in liters using density
    match_mass_to_L = re.match(
        r'convert\s+([\deE.+-]+)\s*grams?\s+of\s+[A-Za-z0-9().*]+\s+to\s+liters?\s+using\s+density\s+([\deE.+-]+)\s*g/mL',
        query, re.IGNORECASE
    )
    if match_mass_to_L:
        grams_str, density_str = match_mass_to_L.groups()
        grams = float(grams_str)
        density = float(density_str)
        return '\n'.join(mass_to_liters(grams, density))

    # Acid-base titration volume calculation
    match_titr = re.match(
        r'what volume of\s+([\deE.+-]+)\s*M\s+([A-Za-z0-9().*]+)\s+is needed to titrate\s+([\deE.+-]+)\s*(mL|L)\s+of\s+([\deE.+-]+)\s*M\s+([A-Za-z0-9().*]+)',
        query, re.IGNORECASE
    )
    if match_titr:
        base_M_str, base_comp, acid_V_str, acid_V_unit, acid_M_str, acid_comp = match_titr.groups()
        # Parse numbers
        base_M = float(base_M_str)
        acid_V = float(acid_V_str)
        # Convert acid volume to liters
        if acid_V_unit.lower() == 'ml':
            acid_V = acid_V / 1000
        acid_M = float(acid_M_str)
        # Normalize compounds
        base_comp = fix_case(base_comp)
        acid_comp = fix_case(acid_comp)
        return '\n'.join(titration_volume(base_M, base_comp, acid_V, acid_M, acid_comp))

    # Polyprotic acid titration concentration
    match_titr2 = re.match(
        r'it takes\s+([\deE.+-]+)\s*(?:mL|L)\s+of\s+([\deE.+-]+)\s*M\s+([A-Za-z0-9().*]+)\s+solution\s+to\s+completely neutralize\s+([\deE.+-]+)\s*(?:mL|L)\s+of\s+([A-Za-z0-9().*]+)\s+solution',
        query, re.IGNORECASE
    )
    if match_titr2:
        Vt_str, Mt_str, titrant_comp, Va_str, analyte_comp = match_titr2.groups()
        Vt = float(Vt_str) / 1000
        Mt = float(Mt_str)
        Va = float(Va_str) / 1000
        titrant_comp = fix_case(titrant_comp)
        analyte_comp = fix_case(analyte_comp)
        return '\n'.join(titration_concentration(Vt, Mt, titrant_comp, Va, analyte_comp))

    # Compound-to-element conversion
    match_elem = re.match(
        r'convert\s+([\deE.+-]+)\s*moles?\s+of\s+([A-Za-z0-9().*]+)\s+to\s+moles?\s+of\s+([A-Za-z][a-z]?)$',
        query, re.IGNORECASE
    )
    if match_elem:
        val_str, comp_formula, elem_symbol = match_elem.groups()
        val = float(val_str)
        comp_formula = fix_case(comp_formula)
        elem_symbol = elem_symbol.capitalize()
        return '\n'.join(moles_compound_to_element(val, comp_formula, elem_symbol))

    # Stoichiometric mol-to-mol conversion
    match_stoich = re.match(
        r'convert\s+([\deE.+-]+)\s*moles?\s+of\s+([A-Za-z0-9().*]+)\s+to\s+moles?\s+of\s+([A-Za-z0-9().*]+)\s+in reaction\s+(.+)',
        query, re.IGNORECASE
    )
    if match_stoich:
        val_str, comp_from, comp_to, reaction = match_stoich.groups()
        val = float(val_str)
        comp_from = fix_case(comp_from)
        comp_to = fix_case(comp_to)
        return '\n'.join(stoichiometry_moles(val, comp_from, comp_to, reaction))

    # Limiting reagent and theoretical yield
    match_lr = re.match(
        r'limiting reagent and theoretical yield from\s+([\deE.+-]+)\s*mol\s+of\s+([A-Za-z0-9().*]+)\s+and\s+([\deE.+-]+)\s*mol\s+of\s+([A-Za-z0-9().*]+)\s+in reaction\s+(.+)',
        query, re.IGNORECASE
    )
    if match_lr:
        val1, comp1, val2, comp2, reaction_str = match_lr.groups()
        comp1 = fix_case(comp1)
        comp2 = fix_case(comp2)
        return '\n'.join(
            limiting_reagent_and_theoretical_yield(
                float(val1), comp1, float(val2), comp2, reaction_str
            )
        )

    return "Sorry, I couldn't understand that conversion request."

# Example usage
# … all your function definitions above …

if __name__ == "__main__":
    # Interactive REPL for local testing only
    while True:
        q = input("Enter your dimensional analysis query (or type 'exit' to quit): ")
        if q.strip().lower() == "exit":
            break
        print(dimensional_analysis(q))   # or however you print results
    
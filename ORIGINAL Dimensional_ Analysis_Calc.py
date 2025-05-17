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

def dimensional_analysis(query):
    import re
    query = query.strip()
    print(f"DEBUG: Received query: '{query}'")

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

    # Pattern for general conversions: convert X unit1 of compound to unit2
    conv_pattern = re.compile(
    r'convert\s+([\deE.+-]+)\s*'
    r'(molecules?|particles?|moles?|grams?|liters?|mol|g|l)\s+of\s+'
    r'([A-Za-z0-9().*]+)\s+to\s+'
    r'(molecules?|particles?|moles?|grams?|liters?|mol|g|l)',
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
while True:
    query = input("Enter your dimensional analysis query (or type 'exit' to quit): ")
    if query.lower() == 'exit':
        break
    print(dimensional_analysis(query))
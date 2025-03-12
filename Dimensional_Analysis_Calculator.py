import re

# Constants
AVOGADRO = 6.022e23  # molecules per mole
MOLAR_MASSES = {
    "H2O": 18.015,
    "CO2": 44.01,
    "NaCl": 58.44,
    "O2": 32.00,
    "H2": 2.016
}

def dimensional_analysis(query):
    query = query.lower().strip()

    # Parsing logic
    mole_match = re.search(r'(\d*\.?\d+)\s*moles?\s*of\s*(\w+)', query)
    gram_match = re.search(r'(\d*\.?\d+)\s*grams?\s*of\s*(\w+)', query)
    mole_to_molecule_match = re.search(r'(\d*\.?\d+)\s*moles?\s*to\s*molecules?', query)
    molecule_to_mole_match = re.search(r'(\d*\.?\d+)\s*molecules?\s*to\s*moles?', query)

    if mole_match:
        moles, compound = float(mole_match[1]), mole_match[2].upper()
        return moles_to_grams(moles, compound)
    elif gram_match:
        grams, compound = float(gram_match[1]), gram_match[2].upper()
        return grams_to_moles(grams, compound)
    elif mole_to_molecule_match:
        moles = float(mole_to_molecule_match[1])
        return moles_to_molecules(moles)
    elif molecule_to_mole_match:
        molecules = float(molecule_to_mole_match[1])
        return molecules_to_moles(molecules)
    else:
        return "Sorry, I couldn't understand that conversion request."

def moles_to_grams(moles, compound):
    if compound not in MOLAR_MASSES:
        return f"Unknown compound: {compound}"
    molar_mass = MOLAR_MASSES[compound]
    grams = moles * molar_mass

    explanation = (
        f"Step 1: Start with {moles} moles of {compound}.\n"
        f"Step 2: Multiply by the molar mass of {compound} ({molar_mass} g/mol).\n"
        f"Step 3: {moles} mol × {molar_mass} g/mol = {grams:.2f} g\n"
    )
    return explanation + f"Final Answer: {grams:.2f} g"

def grams_to_moles(grams, compound):
    if compound not in MOLAR_MASSES:
        return f"Unknown compound: {compound}"
    molar_mass = MOLAR_MASSES[compound]
    moles = grams / molar_mass

    explanation = (
        f"Step 1: Start with {grams} grams of {compound}.\n"
        f"Step 2: Divide by the molar mass of {compound} ({molar_mass} g/mol).\n"
        f"Step 3: {grams} g ÷ {molar_mass} g/mol = {moles:.4f} mol\n"
    )
    return explanation + f"Final Answer: {moles:.4f} mol"

def moles_to_molecules(moles):
    molecules = moles * AVOGADRO
    explanation = (
        f"Step 1: Start with {moles} moles.\n"
        f"Step 2: Multiply by Avogadro's number ({AVOGADRO:.3e} molecules/mol).\n"
        f"Step 3: {moles} mol × {AVOGADRO:.3e} = {molecules:.3e} molecules\n"
    )
    return explanation + f"Final Answer: {molecules:.3e} molecules"

def molecules_to_moles(molecules):
    moles = molecules / AVOGADRO
    explanation = (
        f"Step 1: Start with {molecules:.3e} molecules.\n"
        f"Step 2: Divide by Avogadro's number ({AVOGADRO:.3e} molecules/mol).\n"
        f"Step 3: {molecules:.3e} ÷ {AVOGADRO:.3e} = {moles:.4f} mol\n"
    )
    return explanation + f"Final Answer: {moles:.4f} mol"

# Example usage
while True:
    query = input("Enter your dimensional analysis query (or type 'exit' to quit): ")
    if query.lower() == 'exit':
        break
    print(dimensional_analysis(query))

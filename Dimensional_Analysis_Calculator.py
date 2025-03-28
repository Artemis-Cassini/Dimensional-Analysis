import re

# Constants
AVOGADRO = 6.022e23  # molecules per mole, molar masses of regular elements, or atomic mass
MOLAR_MASSES = {
    "H": 1.008, 
    "HE": 4.0026,
    "LI": 6.94,
    "BE": 9.0122,
    "B": 10.81,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "NE": 20.180,
    "NA": 22.990,
    "MG": 24.305,
    "AL": 26.982,
    "SI": 28.085,
    "P": 30.974,
    "S": 32.06,
    "CL": 35.45,
    "AR": 39.948,
    "K": 39.098,
    "CA": 40.078,
    "SC": 44.956,
    "TI": 47.867,
    "V": 50.942,
    "CR": 51.996,
    "MN": 54.938,
    "FE": 55.845,
    "CO": 58.933,
    "NI": 58.693,
    "CU": 63.546,
    "ZN": 65.38,
    "GA": 69.723,
    "GE": 72.630,
    "AS": 74.922,
    "SE": 78.971,
    "BR": 79.904,
    "KR": 83.798,
    "RB": 85.469,
    "SR": 87.62,
    "Y": 88.906,
    "ZR": 91.224,
    "NB": 92.906,
    "MO": 95.95,
    "TC": 98,
    "RU": 101.07,
    "RH": 102.91,
    "PD": 106.42,
    "AG": 107.87,
    "CD": 112.41,
    "IN": 114.82,
    "SN": 118.71,
    "SB": 121.76,
    "TE": 127.60,
    "I": 126.90,
    "XE": 131.29,
    "CS": 132.91,
    "BA": 137.33,
    "LA": 138.91,
    "CE": 140.12,
    "PR": 140.91,
    "ND": 144.24,
    "PM": 145,
    "SM": 150.36,
    "EU": 151.96,
    "GD": 157.25,
    "TB": 158.93,
    "DY": 162.50,
    "HO": 164.93,
    "ER": 167.26,
    "TM": 168.93,
    "YB": 173.05,
    "LU": 174.97,
    "HF": 178.49,
    "TA": 180.95,
    "W": 183.84,
    "RE": 186.21,
    "OS": 190.23,
    "IR": 192.22,
    "PT": 195.08,
    "AU": 196.97,
    "HG": 200.59,
    "TL": 204.38,
    "PB": 207.2,
    "BI": 208.98,
    "PO": 209,
    "AT": 210,
    "RN": 222,
    "FR": 223,
    "RA": 226,
    "AC": 227,
    "TH": 232.04,
    "PA": 231.04,
    "U": 238.03,
    "NP": 237,
    "PU": 244,
    "AM": 243,
    "CM": 247,
    "BK": 247,
    "CF": 251,
    "ES": 252,
    "FM": 257, 
    "MD": 258,
    "NO": 259,
    "LR": 266,
    "RF": 267,
    "DB": 268,
    "SG": 269,
    "BH": 270,
    "HS": 277,
    "MT": 278,
    "DS": 281,
    "RG": 282,
    "CN": 285,
    "NH": 286,
    "FL": 289,
    "MC": 290,
    "LV": 293,
    "TS": 294,
    "OG": 294,
#compounds starting with A, this is the only Actinium
    "AC2O3": 502.053, 
#Silver starts here ----------------------------------
    "AGBF4": 194.673,
    "AGBR": 187.77,
    "AGBRO": 203.771,
    "AGBRO2": 219.772, #NO PAGE
    "AGBRO3": 235.770,
    "AGBRO4": 251.77, #NO PAGE
    "AGCL": 143.32,
    "AGCL3CU2": 341.312, #NO PAGE
    "AGCLO3": 191.319,
    "AGCLO4": 207.319,
    "AGCN": 133.8856,
    "AGCNO": 149.885,
    "AGF": 126.8666,
    "AGF2": 145.865,
    "AGI": 234.77,
    "AGIO": 250.772,
    "AGIO2": 266.768, #NO PAGE
    "AGIO3": 282.77,
    "AGIO4": 298.766, #NO PAGE
    "AGMNO4": 226.804,
    "AGN3": 149.888,
    "AGNO3": 169.872,
    "AGO": 123.87,
    "AGONC": 149.885,
    "AGPF6": 252.83,
    "AGSNC": 165.948,
    "AG2C2": 239.758,
    "AG2CO3": 275.75,
    "AG2C2O4": 303.755,
    "AG2CL2": 286.64, #NO PAGE
    "AG2CRO4": 331.73,
    "AG2CR2O7": 431.76,
    "AG2F": 234.734,
    "AG2MOO4": 375.67,
    "AG2O": 231.735,
    "AG2S": 247.80,
    "AG2SO4": 311.79,
    "AG2SE": 294.7,
    "AG2SEO3": 342.69,
    "AG2SEO4": 358.707, #NO PAGE
    "AG2TE": 343.3364,
    "AG3BR2": 483.418, #NO PAGE
    "AG3BR3": 563.322, #NO PAGE
    "AG3CL3": 429.96, #NO PAGE
    "AG3I3": 704.31, #NO PAGE
    "AG3PO4": 418.574,
#start of aluminium ----------------------------------
    "ALBO3": 85.789, #NO PAGE
    "ALBR": 106.886,
    "ALBR3": 266.694,
    "ALCL": 62.43,
    "ALCLF": 81.43, #NO PAGE
    "ALCL2F": 116.88, #NO PAGE
    "ALCLO": 78.431, #NO PAGE
    "ALCL2H": 98.89, #NO PAGE
    "ALCL3": 133.341,
    "ALCL4CS": 301.692, #NO PAGE
    "ALCL4K": 207.88, #NO PAGE
    "ALCL4NA": 191.78331,
    "ALCL4RB": 254.25, #NO PAGE
    "ALCL6K3": 356.976, #NO PAGE
    "ALCL6NA3": 308.652, #NO PAGE
    "ALCO": 85.915, #NO PAGE
    "ALF": 45.98,
    "ALFO": 61.979, #NO PAGE
    "ALF2": 64.978, #NO PAGE
    "ALF2O": 80.977, #NO PAGE
    "ALF3": 83.977,
    "ALF4K": 142,
    "ALF4LI": 109.914, #NO PAGE
    "ALF6K3": 258.274, #NO PAGE
    "ALF6LI3": 161.79,
    "ALF6NA3": 209.9,
    "ALGAINP": 242.499, #NO PAGE
    "ALH3": 30.006,
    "ALH4": 31.014, #NO PAGE
    "AL(OH)3": 78.003,
    "ALI": 153.886,
    "ALI3": 407.695,
    "ALLIO2": 65.92,
    "ALN": 40.989,
    "AL(NO2)3": 164.997, #NO PAGE
    "AL(NO3)3": 212.996,
    "ALNAO2": 81.97,
    "ALO": 42.98,
    "ALOSI": 71.066, #NO PAGE
    "ALO2": 58.98, #NO PAGE
    "ALP": 57.9552,
    "ALPO4": 121.9529,
    "ALTE": 154.582, #NO PAGE
    "ALTE2": 282.182, #NO PAGE
    "AL2BEO4": 126.9722,
    "AL2BR6": 533.388, #NO PAGE
    "AL2(CO3)3": 233.988,
    "AL2CL9K3": 490.308, #NO PAGE
    "AL2COO4": 176.892,
    "AL2F6": 167.952,
    "AL2I6": 407.695,
    "AL2MGO4": 142.265,
    "AL2O": 69.963, #NO PAGE
    "AL2O2": 85.962, #NO PAGE
    "AL2O3": 101.960,
    "AL2O5SI": 162.0456,
    "AL2SIO5": 162.044,
    "AL2O7SI2": 222.127, #NO PAGE
    "AL2S": 86.024, #NO PAGE
    "AL2S3": 150.158,
    "AL2(SO4)3": 342.15,
    "AL2SE": 132.935, #NO PAGE
    "AL2SI2O5(OH)4": 258.157,
    "AL2TE": 181.564, #NO PAGE
    "AL3F14NA5": 461.8711070,
    "AL4C3": 143.95853,
    "AL6BEO10": 330.8942,
    "AL6O13SI2": 426.049,
# start of Argon here YAY ----------------------------------
    "ARCLF": 94.396, #NO PAGE
    "ARCLH": 76.406, #NO PAGE
    "ARFH": 59.954, #NO PAGE
# start of Arsenic here YIPPEE ----------------------------------
    "ASBRO": 170.825,
    "ASBR3": 314.634,
    "ASCLO": 126.371,
    "ASCL3": 181.272,
    "ASCL3O": 197.271,
    "ASCL4F": 235.72,
    "ASF3": 131.916,
    "ASF5": 169.912,
    "ASH3": 77.946,
    "ASI3": 455.622,
    "ASO": 90.921,
    "ASO2": 106.92,
    "ASP": 105.896,
    "ASP3": 167.844,
    "ASTL": 279.302,
    "AS2I4": 657.444,
    "AS2O3": 197.841,
    "AS2P2": 211.792,
    "AS2O5": 229.839,
    "AS2S4": 278.084,
    "AS2S5": 310.144,
    "AS2SE": 228.815,
    "AS2SE3": 386.757,
    "AS2SE5": 544.699,
    "AS3O4": 288.762,
    "AS3P": 255.74,
    "AS4O3": 347.685,
    "AS4O5": 379.683,
    "AS4S3": 395.868,
    "AS4S4": 427.928,
# we start GOLD GOLD GOLD  ----------------------------------
    "AUB": 207.78,
    "AUBO": 223.779,
    "AUBR": 276.874,
    "AUBR3": 436.682,
    "AUCN": 222.988,
    "AUCL": 232.42,
    "AUCL3": 303.32,
    "AUF3": 253.964,
    "AUI": 323.87,
    "AUI3": 577.67,
    "AU(OH)3": 247.991,
    "AUTE": 324.57,
    "AU2O3": 441.937,
    "AU2S": 426,
    "AU2S3": 490.12,
    "AU2(SEO4)3": 822.841,
    "AU2SE3": 630.853,
# B starts here yippee we start with Boron  ----------------------------------
    "BAS": 85.732,
    "BBR3": 250.522,
    "BCL3": 117.16,
    "BF3": 67.804,
    "BH3": 13.834,
    "BH4": 14.842,
    "BH6N": 30.865,
    "BI3": 391.57,
    "BN": 24.817,
    "B(NO3)3": 100.828,
    "B(NO3)4": 258.826,
    "B(OH)3": 61.831,
    "BO3": 58.807,
    "BP": 41.784,
    "BPO4": 105.78,
    "B2CL4": 163.42,
    "B2F4": 97.612,
    "B2H6": 27.668,
    "B2MG": 45.925,
    "B2H2SE3": 260.549,
    "B2O3": 69.617,
    "B2S3": 117.8,
    "B3N3H6": 80.499,
    "B4C": 55.251,
    "B4CAO7": 195.311,
    "B4NA2O710H2O": 381.36,
# BARIUM BARIUM GO GO GO GO GO GO GO GO GO  ----------------------------------
    "BA(ALO2)2": 255.29,
    "BA(ASO3)2": 383.168,
    "BA(ASO4)2": 415.166,
    "BAB6": 202.19,
    "BABR2": 297.138,
    "BA(BRO)2": 329.136,
    "BA(BRO2)2": 361.134,
    "BA(BRO3)2": 393.132, 
    "BA(BRO3)2H2O": 411.147,
    "BA(BRO3)22H2O": 429.162,
    "BA(BRO4)2": 425.13,
    "BA(BRO4)23H2O": 479.19,
    "BA(BRO4)24H2O": 497.19,
    "BA(CHO2)2": 227.364,
    "BA(C2H3O2)2": 255.418,
    "BA(CN)2": 189.366,
    "BAHFO3": 363.817,
    "BAHGI4": 845.52,
    "BA(HS)2": 203.466,
    "BAI2": 391.13, 
    "BA(IO)2": 423.128,
    


}

def dimensional_analysis(query):
    query = query.lower().strip()

    # Parsing logic
    mole_match = re.search(r'(convert\s*)?(\d*\.?\d+)\s*moles?\s*(of\s*)?(\w+)', query)
    gram_match = re.search(r'(convert\s*)?(\d*\.?\d+)\s*grams?\s*(of\s*)?(\w+)', query)
    kg_match = re.search(r'(convert\s*)?(\d*\.?\d+)\s*(kilograms?|kg)\s*(of\s*)?(\w+)', query)
    mg_match = re.search(r'(convert\s*)?(\d*\.?\d+)\s*(milligrams?|mg)\s*(of\s*)?(\w+)', query)
    liter_match = re.search(r'(convert\s*)?(\d*\.?\d+)\s*(liters?|l)\s*(of\s*)?(\w+)', query)
    mole_to_molecule_match = re.search(r'(convert\s*)?(\d*\.?\d+)\s*moles?\s*(to\s*)?molecules?', query)
    molecule_to_mole_match = re.search(r'(convert\s*)?(\d*\.?\d+)\s*molecules?\s*(to\s*)?moles?', query)

    if mole_match:
        moles, compound = float(mole_match[2]), mole_match[4].upper()
        return moles_to_grams(moles, compound)
    elif gram_match:
        grams, compound = float(gram_match[2]), gram_match[4].upper()
        return grams_to_moles(grams, compound)
    elif kg_match:
        kg, compound = float(kg_match[2]), kg_match[5].upper()
        grams = mass_conversion(kg, "kg")
        return grams_to_moles(grams, compound)
    elif mg_match:
        mg, compound = float(mg_match[2]), mg_match[4].upper()
        grams = mass_conversion(mg, "mg")
        return grams_to_moles(grams, compound)
    elif liter_match:
        volume_liters = float(liter_match[2])
        return liters_to_moles(volume_liters)
    elif mole_to_molecule_match:
        moles = float(mole_to_molecule_match[2])
        return moles_to_molecules(moles)
    elif molecule_to_mole_match:
        molecules = float(molecule_to_mole_match[2])
        return molecules_to_moles(molecules)
    else:
        return "Sorry, I couldn't understand that conversion request. Try phrasing it like: 'Convert 5 grams of CO2 to moles' or 'How many molecules are in 2 moles?'"

def moles_to_grams(moles, compound):
    if compound not in MOLAR_MASSES:
        return f"Error: Unknown compound '{compound}'. Please check spelling or add it to the database."
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
        return f"Error: Unknown compound '{compound}'. Please check spelling or add it to the database."
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

def liters_to_moles(volume_liters):
    mol = volume_liters / 22.4
    explanation = (
        f"Step 1: Start with {volume_liters} L.\n"
        f"Step 2: Divide by molar volume at STP (22.4 L/mol).\n"
        f"Step 3: {volume_liters} L ÷ 22.4 L/mol = {mol:.4f} mol\n"
    )
    return explanation + f"Final Answer: {mol:.4f} mol"

def mass_conversion(mass, from_unit):
    if from_unit == "kg":
        return mass * 1000
    elif from_unit == "mg":
        return mass / 1000
    elif from_unit == "g":
        return mass
    return mass

# Example usage
while True:
    query = input("Enter your dimensional analysis query (or type 'exit' to quit): ")
    if query.lower() == 'exit':
        break
    print(dimensional_analysis(query))

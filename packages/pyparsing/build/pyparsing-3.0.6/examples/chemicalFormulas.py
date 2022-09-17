#
# chemicalFormulas.py
#
# Copyright (c) 2003,2019 Paul McGuire
#

import pyparsing as pp

atomicWeight = {
    "O": 15.9994,
    "H": 1.00794,
    "Na": 22.9897,
    "Cl": 35.4527,
    "C": 12.0107,
}

digits = "0123456789"

# Version 1
element = pp.Word(pp.alphas.upper(), pp.alphas.lower(), max=2).set_name("element")
# for stricter matching, use this Regex instead
# element = Regex("A[cglmrstu]|B[aehikr]?|C[adeflmorsu]?|D[bsy]|"
#                 "E[rsu]|F[emr]?|G[ade]|H[efgos]?|I[nr]?|Kr?|L[airu]|"
#                 "M[dgnot]|N[abdeiop]?|Os?|P[abdmortu]?|R[abefghnu]|"
#                 "S[bcegimnr]?|T[abcehilm]|U(u[bhopqst])?|V|W|Xe|Yb?|Z[nr]")
elementRef = pp.Group(element + pp.Optional(pp.Word(digits), default="1"))
formula = elementRef[...]


def sum_atomic_weights(element_list):
    return sum(atomicWeight[elem] * int(qty) for elem, qty in element_list)


formula.runTests(
    """\
    H2O
    C6H5OH
    NaCl
    """,
    fullDump=False,
    postParse=lambda _, tokens: "Molecular weight: {}".format(
        sum_atomic_weights(tokens)
    ),
)
print()

# Version 2 - access parsed items by results name
elementRef = pp.Group(
    element("symbol") + pp.Optional(pp.Word(digits), default="1")("qty")
)
formula = elementRef[...]


def sum_atomic_weights_by_results_name(element_list):
    return sum(atomicWeight[elem.symbol] * int(elem.qty) for elem in element_list)


formula.runTests(
    """\
    H2O
    C6H5OH
    NaCl
    """,
    fullDump=False,
    postParse=lambda _, tokens: "Molecular weight: {}".format(
        sum_atomic_weights_by_results_name(tokens)
    ),
)
print()

# Version 3 - convert integers during parsing process
integer = pp.Word(digits).setParseAction(lambda t: int(t[0])).setName("integer")
elementRef = pp.Group(element("symbol") + pp.Optional(integer, default=1)("qty"))
formula = elementRef[...].setName("chemical_formula")


def sum_atomic_weights_by_results_name_with_converted_ints(element_list):
    return sum(atomicWeight[elem.symbol] * int(elem.qty) for elem in element_list)


formula.runTests(
    """\
    H2O
    C6H5OH
    NaCl
    """,
    fullDump=False,
    postParse=lambda _, tokens: "Molecular weight: {}".format(
        sum_atomic_weights_by_results_name_with_converted_ints(tokens)
    ),
)
print()

# Version 4 - parse and convert integers as subscript digits
subscript_digits = "₀₁₂₃₄₅₆₇₈₉"
subscript_int_map = {e[1]: e[0] for e in enumerate(subscript_digits)}


def cvt_subscript_int(s):
    ret = 0
    for c in s[0]:
        ret = ret * 10 + subscript_int_map[c]
    return ret


subscript_int = pp.Word(subscript_digits).addParseAction(cvt_subscript_int).set_name("subscript")

elementRef = pp.Group(element("symbol") + pp.Optional(subscript_int, default=1)("qty"))
formula = elementRef[1, ...].setName("chemical_formula")
formula.runTests(
    """\
    H₂O
    C₆H₅OH
    NaCl
    """,
    fullDump=False,
    postParse=lambda _, tokens: "Molecular weight: {}".format(
        sum_atomic_weights_by_results_name_with_converted_ints(tokens)
    ),
)
print()

#!/usr/bin/env python
"""
Parses the KEGG module databases. Entries in this database have certain format,
that includes information of different combinations of genes to carry out the full
module function. This script downloads the module 'definition' in memory and parses
out the possible combinations of genes.
"""
__author__ = "Joris van Steenbrugge"
__email__ = "joris.vansteenbrugge@wur.nl"

import requests
import re
from sympy.parsing.sympy_parser import parse_expr
from sympy.logic.boolalg import to_dnf
from sympy import symbols
from bs4 import BeautifulSoup


def pprint(l):
    for val in l:
        print(val)


def retrieve_definition(module_name):
    """Retrieves the html page based on the database entry. The module definition
    is then retrieved from the html text.

        Keyword arguments:
            module_name -- str, name of the module, should be matchable with
                            r'M[0-9]{5}'.

        Returns:
            definition -- str, indicating AND and OR relationship between gene
                            annotations (e.g. (A,B) C -> (A OR B) AND C).
    """
    base_url = "https://www.genome.jp/dbget-bin/www_bget?md:{}"

    # Retrieve the html using web requests
    html_page = requests.get(base_url.format(module_name))

    # Parse the base HTML using BeautifulSoup
    soup = BeautifulSoup(html_page.text, 'html.parser')
    all_lines = soup.text.split('\n')

    # Extract the line containing the module definition
    definition = ""
    read = False
    for line in all_lines:
        if line.startswith("Definition"):
            read = True
            continue
        if read:
            definition = line
            break
    return definition


def clean_definition(definition):
    """Give the definition without any optional terms.
       Optional terms are indicated using '-' (e.g. (A-B) C <- A optional B AND C)

       Keyword arguments:
           definition -- str, indicating AND and OR relationship between gene
                           annotations (e.g. (A,B) C -> (A OR B) AND C).

    """
    if '-' not in definition:  # then we don't need cleaning
        return definition

    # double -- quick and dirty fix
    definition = definition.replace(" -- --", "").replace(" -- ", "-").replace("-- ", '-').replace(" --", "")

    # Keep doing this until all optional terms are gone
    while '-' in definition:
        min_pos = definition.index('-')

        # If there is only one KO term afterwards we skip ahead the length of the
        # KO term to remove itself.
        if definition[min_pos + 1] == 'K':
            definition = definition.replace(definition[min_pos:min_pos + 7], "")  # remove the -K
            continue

        # In this case, there are multiple KO terms afterwards, so we want to
        # remove them all.
        elif definition[min_pos + 1] == '(':

            opens = ['(']
            end_pos = None
            # We go find the corresponding closing parenthesis
            for idx in range(min_pos + 2, len(definition)):
                char = definition[idx]
                if char == ')':
                    opens.pop()
                elif char == '(':
                    opens.append('(')

                # We found the closing parenthesis
                if len(opens) == 0:
                    end_pos = idx
                    break

            if end_pos:
                # Remove all characters between the optional parentheses
                definition = definition.replace(definition[min_pos:end_pos + 1], "")
        else:
            raise ValueError("Unsupported optional term found", definition[min_pos + 1])

    definition = definition.replace("  ", " ")
    if definition.startswith(' '):
        definition = definition[1:]

    while definition.endswith(" "):
        definition = definition[0:len(definition) - 1]

    return definition.replace("  ", " ")


def get_symbol_dict(definitions):
    """Create a symbol dictionary, as required for the sympy_parserself.

        Keyword arguments:
            definition -- str, indicating AND and OR relationship between gene
                            annotations (e.g. (A,B) C -> (A OR B) AND C).

    """
    ko_symbols = re.findall(r'K[0-9]{5}', definitions)
    symbol_dict = {ko: symbols(ko) for ko in ko_symbols}

    return symbol_dict


def parse_definition(definitions):
    """Convert a definitions set to its disjunctive normal form

        Keyword arguments:
            definition -- str, indicating AND and OR relationship between gene
                            annotations (e.g. (A,B) C -> (A OR B) AND C).
    """
    symbol_dict = get_symbol_dict(definitions)

    expr_text = definitions.replace(" ", " & ").replace(",", " | ").replace('+', ' & ')
    if expr_text.endswith('& ') or expr_text.endswith('| '):
        expr_text = expr_text[:-3]
        print(expr_text)
    try:
        expr = parse_expr(expr_text, local_dict=symbol_dict)
    except SyntaxError:
        print(definitions)
        print(expr_text)
        return

    dnf = to_dnf(expr)
    combinations = [d.strip("() ").replace(" ", "") for d in str(dnf).split("|")]

    combinations = [x.split("&") for x in combinations]

    return combinations


def parse_module(module_name):
    definition = retrieve_definition(module_name)

    definition = clean_definition(definition)

    try:
        definition_parsed = parse_definition(definition)
    except ValueError as err:
        print(err.args)
        definition_parsed = []

    return definition_parsed


if __name__ == '__main__':
    mod_list = []
    with open('../data/kegg_modules_2019_07_23.tsv') as in_file:
        in_file.readline()
        for line in in_file:
            line = line.rstrip().split()
            mod_list.append(line[0])

    mod_list = list(set(mod_list))

    with open("../data/kegg_modules_2019_07_23_DNF_flast.tsv", 'w') as dnf_file:
        for module in mod_list:
            print(module)
            options = parse_module(module)

            opt_line = ";".join([','.join(option) for option in options])
            dnf_file.write(f"{module};{opt_line}\n")

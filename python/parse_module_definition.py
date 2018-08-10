#!/usr/bin/env python

import requests
import re
from sympy.parsing.sympy_parser import parse_expr
from sympy.logic.boolalg import to_dnf
from sympy import symbols
from bs4 import BeautifulSoup

def pprint(l):
    for val in l:
        print val

def retrieve_definition(module_name):
    base_url = "https://www.genome.jp/dbget-bin/www_bget?md:{}"
    html_page = requests.get(base_url.format(module_name))

    soup = BeautifulSoup(html_page.text, 'html.parser')
    all_lines = soup.text.split('\n')


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
    """Give the definition without any optional terms
    """
    if '-' not in definition: #then we don't need cleaning
        return definition

    while '-' in definition:
        min_pos = definition.index('-')
        if definition[min_pos+1] == 'K':
            definition = definition.replace(definition[min_pos:min_pos+7], "") #remove the -K


            continue
        elif definition[min_pos+1] == '(':

            opens = ['(']
            end_pos = None
            # We go find corresponding closing parenthesis
            for idx in range(min_pos+2, len(definition)):
                char = definition[idx]
                if char == ')':
                    opens.pop()
                elif char == '(':
                    opens.append('(')

                if len(opens) == 0:
                    end_pos = idx
                    break

            if end_pos:
                definition = definition.replace(definition[min_pos:end_pos+1], "")
        else:
            print 'wtf'

    return definition


def get_symbol_dict(definitions):
    ko_symbols = re.findall(r'K[0-9]{5}', definitions)
    symbol_dict = {ko:symbols(ko) for ko in ko_symbols}

    return symbol_dict

def parse_definition(definitions):
    symbol_dict = get_symbol_dict(definitions)

    expr_text = definitions.replace(" ", " & ").replace(",", " | ").replace('+', ' & ')
    expr = parse_expr(expr_text, local_dict = symbol_dict)

    dnf = to_dnf(expr)
    combinations = [d.translate(None,"() ") for d in str(dnf).split("|")]

    combinations = [x.split("&") for x in combinations]

    return combinations


def parse_module(module_name):
    definition = retrieve_definition(module_name)
    definition = clean_definition(definition)


    return parse_definition(definition)

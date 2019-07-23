#!/usr/bin/env python2


import requests
import re
from datetime import date



def get_module_names(file):
    modules = []
    with open(file) as kegg_brite:
        for line in kegg_brite:
            if not line.startswith("D"):
                continue

            modules.append(line.strip().split()[1])

    return modules

def parseModule(module):
    module_specific = "http://www.genome.jp/kegg-bin/show_module?{}"
    r = requests.get(module_specific.format(module))
    r = r.text

    matches = list(set([ match.group(1) for match in re.finditer(r'>(K[0-9]{5})',r)]))

    return matches


def writeModule(name, definition, file):
    for ko in definition:
        file.write(f"{name}\t{ko}\n")



if __name__ == "__main__":
    date_today = str(date.today()).replace('-',"_")
    outfile = open(f'kegg_modules_{date_today}.tsv','w')
    outfile.write("Module\tAnnotation\n")


    modules = get_module_names('tmp_ko00002_2019_7_23.keg')

    total = len(modules)

    c = 1
    for module in modules:
        definition = parseModule(module)
        writeModule(module, definition, outfile)

        print(f'done {c}/{total}')
        c += 1

    outfile.close()

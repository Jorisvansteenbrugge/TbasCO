import json

import os

def Load_json(json_file):
    with open(json_file) as f:
        return json.load(f)




def Parse_json(data, levels = 3):

    if levels == 3:
        
        for mod_type in data['children']:
            module_type = mod_type['name']

            for category in mod_type['children']:
                current_cat = category['name']

                

                if 'children' not in category['children'][0]:
                    current_trait = category['name']
                    

                    module_name    = current_trait
                    module_description = module_name

                    print(f"{module_name}\t{module_type}\t{current_cat}\t{module_description}")
                
                else:
                    
                    for trait in category['children']:
                        current_trait = trait['name']

                        
                        module_name    = current_trait
                        module_description = module_name

                        print(f"{module_name}\t{module_type}\t{current_cat}\t{module_description}")

    if levels == 2:

        for mod_type in data['children']:
            module_type = mod_type['name']

            for trait in mod_type['children']:
                current_trait = trait['name']

                current_trait = trait['name']
                module_name    = current_trait
                module_description = " ".join(current_trait)

                print(f"{module_name}\t\t{module_type}\t{module_description}")

if __name__ == '__main__':
    Parse_json(Load_json('data/brite/2020_05_kegg_brite_transporters.json'), levels = 3)
    Parse_json(Load_json('data/brite/2020_05_kegg_brite_two_component_systems.json'), levels = 2)
    Parse_json(Load_json('data/brite/2020_05_kegg_brite_secretion_systems.json'), levels = 2)

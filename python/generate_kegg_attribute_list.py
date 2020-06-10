import json

def Load_json(json_file):
    with open(json_file) as f:
        return json.load(f)

def Parse_json(data):
    header = "Module_name\tModule_type\tGeneral_category\tSpecific_category\tModule_description"
    print(header)

    for mod_type in data['children']:
        module_type = mod_type['name']

        for category in mod_type['children']:
            current_cat = category['name']

            for specific_category in category['children']:
                current_specific_category = specific_category['name']

                for module in specific_category['children']:
                    current_module = module['name'].split(" ")
                    module_name    = current_module[0]
                    module_description = " ".join(current_module[1:])

                    print(f"{module_name}\t{module_type}\t{current_cat}\t{current_specific_category}\t{module_description}")


if __name__ == '__main__':
    data = Load_json('data/kegg_categories_2020_05.json')
    Parse_json(data)

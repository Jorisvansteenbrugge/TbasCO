

def import_att_list(file_name):
    node_map = {}

    with open(file_name) as in_file:
        in_file.readline()

        for line in in_file:
            line = line.strip().split(';')
            node_map[line[0]] = line[1:]

    return node_map

def parse_node_name(node_name):
    return node_name.split(".")[0]
    

def lookup_node(node_map, node_name, outfile = 'expanded_attribute_table.csv'):
    with open(outfile, 'a') as out:
        node_name_clean = parse_node_name(node_name)

        try:
            outline = [node_name] + node_map[node_name_clean]
            out.write(";".join(outline) + "\n")
        except:
            pass


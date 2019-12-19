import re

def get_match_objs(BRITE_instance):

  BRITE_content = BRITE_instance.read()
  if BRITE_content.startswith('+E'):
    print('+e')
    re_pattern = r"C {4}(.+[^\n])\n(D      K[0-9]{5}[^\n]+\n)+"
  elif BRITE_content.startswith('+C'):
    print("+c")
    re_pattern = r"B  (.+)\n(C    K[0-9]{5}[^\n]+\n)+"

  for match in re.finditer(re_pattern, BRITE_content):
    yield match


def get_ko_ids(full_match_group, re_pattern = r"K[0-9]{5}"):
  return re.findall(re_pattern, full_match_group)


def Parse_brite_db(BRITE_DIR = "../data/brite"):
  from glob import glob

  database = {}

  for BRITE_FILE in glob(f"{BRITE_DIR}/*.keg"):
    print(BRITE_FILE)
    with open(BRITE_FILE) as BRITE_instance:
      for match_group in get_match_objs(BRITE_instance):

        group_name = match_group[1]
        ko_terms = get_ko_ids(match_group[0])
        database[group_name] = ko_terms

  print(len(database.keys()))
  return database



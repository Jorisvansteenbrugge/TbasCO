def get_module_categories(categories):
  cat_names = []
  output = []
  current = []
  with open(categories) as in_file:
    for line in in_file:
      line = line.strip()
      if line.startswith('C'):
        line = line.split()
        cat_names.append(" ".join(line[1:]))
        if len(current) != 0:
          output.append(current)
      elif line.startswith('D'):
        line = line.split()
        current.append(line[1])
    output.append(current)

  return (output, cat_names)



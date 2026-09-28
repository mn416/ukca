#!/usr/bin/python3 

import sys
import re

asad_vars = [
  "cdt",
  "cdt_diag",
  "co3",
  "deriv",
  "dpd",
  "dpw",
  "ej",
  "emr",
  "f",
  "fdot",
  "fj",
  "fpsc1",
  "fpsc2",
  "ftilde",
  "interval",
  "ipa",
  "jsubs",
  "lati",
  "linfam",
  "ltrig",
  "modified_map",
  "ncsteps",
  "ncsteps_factor",
  "p",
  "pd",
  "pmintnd",
  "prk",
  "prod",
  "qa",
  "ratio",
  "rk",
  "sh2o",
  "shno3",
  "slos",
  "spfj",
  "sph2o",
  "sphno3",
  "t",
  "t300",
  "tnd",
  "wp",
  "co2",
  "y",
  "ydot",
  "za",
  "firstcall",
]

if len(sys.argv) != 2:
  print("Usage: asad_rename.py FILE.F90")
  sys.exit(-1)

with open(sys.argv[1], "r") as f:
  contents = f.read()

routines = []
rest = contents
while rest:
  result = rest.split("SUBROUTINE", 1)
  routines.append(result[0])
  if len(result) == 1: break
  rest = "SUBROUTINE" + result[1]
  result = rest.split("END SUBROUTINE", 1)
  routines.append(result[0])
  if len(result) == 1: break
  rest = "END SUBROUTINE" + result[1]

new_routines = []
for r in routines:
  use_pattern = r"(USE asad_mod, *ONLY: )(.*?[^&])\n"
  match = re.search(use_pattern, r, re.DOTALL)
  if match:
    state_vars = []
    non_state_vars = []
    for item in match.group(2).split(","):
        item = item.replace(r"&", "")
        item = item.replace("\n", "")
        item = item.strip()
        if item == "&": continue
        if item in asad_vars:
            state_vars.append(item)
        else:
            non_state_vars.append(item)
    if state_vars:
      non_state_vars.append("s=>asad_state")
      # Remove state vars from use statement
      def chunk(xs, n):
        return [xs[i:i + n] for i in range(0, len(xs), n)]
      lines = []
      chunks = chunk(non_state_vars, 5)
      for (i, c) in enumerate(chunks):
        prefix = match.group(1)
        if lines: prefix = " " * len(prefix)
        line = prefix + ", ".join(c)
        if i < len(chunks)-1:
          line += ","
          line += " " * (79 - len(line)) + "&"
        lines.append(line)
      new_use_stmt = "\n".join(lines)
      r = re.sub(use_pattern, new_use_stmt + "\n", r, flags=re.DOTALL)
      # Rename all state vars
      non_fortran_id_char = r"(^|$|\n|[^a-zA-Z0-9_])"
      asad_var_pattern = non_fortran_id_char + \
                         r"(" + "|".join(state_vars) + r")" + \
                         non_fortran_id_char
      r = re.sub(asad_var_pattern, r"\1s%\2\3", r)
  new_routines.append(r)

lines = "".join(new_routines).splitlines()
new_lines = []
for line in lines:
  if line.endswith("&") and len(line) > 80:
    space_to_remove = len(line) - 80
    trim_str = " " * space_to_remove + "&"
    line = line.replace(trim_str, "&")
  new_lines.append(line)

new_contents = "\n".join(new_lines) + "\n"
with open(sys.argv[1], "w", encoding="utf-8") as f:
  f.write(new_contents)

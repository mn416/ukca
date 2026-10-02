#!/usr/bin/python3

import re
import sys
import copy
from pathlib import Path

if len(sys.argv) <= 1:
  print("Usage: top_sort.py root_1 root_2 ... root_n")
  sys.exit(-1)

def all_calls(s):
  calls = {}
  call_pattern = r"CALL\s+([A-Za-z0-9_]+)\s*\("
  for match in re.finditer(call_pattern,s):
    sub_loc = s.rfind("SUBROUTINE ", 0, match.start())
    if sub_loc < 0: continue
    call_name = match.group(1)
    sub_str = s[sub_loc:]
    sub_match = re.search(r"SUBROUTINE\s+([A-Za-z0-9_]+)", sub_str)
    if sub_match:
      sub_name = sub_match.group(1)
      if sub_name not in calls:
        calls[sub_name] = set()
      calls[sub_name].add(call_name)
  return calls


files = [p.name for p in Path('.').glob('*.F90') if p.is_file()]
calls = {}
for file_name in files:
  with open(file_name, "r") as f:
    contents = f.read()
  for (k, v) in all_calls(contents).items():
    if k not in calls: calls[k] = set()
    calls[k].update(v)

# Prune subs not reachable from roots
roots = sys.argv[1:]
pruned = {}
frontier = set(roots)
visited = set()
while frontier:
    sub_name = frontier.pop()
    if sub_name in visited: continue
    visited.add(sub_name)
    if sub_name not in calls: continue
    cs = calls[sub_name]
    pruned[sub_name] = cs
    frontier.update(cs - visited)
calls = pruned

# Add called subs that don't exist
exclude = {"dr_hook", "umPrint", "ereport"}
for cs in list(calls.values()):
  for c in cs:
    if c not in exclude and c not in calls: calls[c] = set()

# Topological sort
sorted_calls = []
while calls:
   no_dep = None
   for (sub, cs) in calls.items():
     if not cs:
       no_dep = sub
       break
   if not no_dep:
     no_dep = list(calls.keys())[0]
   sorted_calls.append(no_dep)
   del calls[no_dep]
   for uses in calls.values():
       uses.discard(no_dep)

for call in sorted_calls:
  print(call)

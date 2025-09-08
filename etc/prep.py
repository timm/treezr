import re, sys
old = r'(?m)^(def[^\n]*:\s*\n)([ \t]+)[ru]?["\']([^"\']+?)["\']\s*\n'
new = r'# \3.\n\1'
with open(sys.argv[1]) as f: 
  print(re.sub(old, new, f.read(), flags=re.MULTILINE))

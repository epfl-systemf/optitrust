import os

transfo = "lib/framework/transfo/"
out = "out.ml"
cnt = 30
gap = 2

with open(out, 'w') as out_f:
  for f_path in os.listdir(transfo):
    with open(transfo + f_path) as f:
      um = (cnt - len(f_path))//2
      rem = (cnt - len(f_path))%2
      docstring = f"""
(* {'/'*cnt}
   {'/'*(um-gap)}{' '*gap}{f_path.capitalize()}{' '*gap}{'/'*(um-gap+rem)}
   {'/'*cnt} *)\n
      """
      incomment = False # waiting for line
      for line in f.readlines():
        if not incomment and line.startswith("(** ["):
          docstring += "(*========================*)\n"
          incomment = True
        if incomment:
          docstring += line
          if line.endswith("*)\n"):
            incomment = False
      out_f.write(f"{docstring}")

import os

print("sample-id,forward-absolute-filepath,reverse-absolute-filepath")

with open("data/manifest.txt", 'r') as fin:
   next(fin)
   for line in fin:
       llist = line.strip().split('\t')
       fwd = llist[1]
       if os.path.isfile("/".join(fwd.split("/")[1:])):
           print(",".join(llist))
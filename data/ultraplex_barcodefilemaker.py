"""parse brcodes into weird ultraplex barcode format
"""


revdict = []
with open("barcodes_rev.txt") as rhandle:
    for line in rhandle:
       revdict.append(line.strip())


with open("barcodes_fwd.txt") as fhandle:
    for line in fhandle:
       flist = [line.strip()]
       flist.extend(revdict)
       print(",".join(flist))


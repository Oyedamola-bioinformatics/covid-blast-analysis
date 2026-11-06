import difflib, itertools

def read(p):
    return "".join(l.strip() for l in open(p) if not l.startswith(">")).upper()

ref, smp = read("data/NC_045512.2.fasta"), read("data/sample_omicron.fasta")
print("lens", len(ref), len(smp), "N in sample:", smp.count("N"))

bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
tab = {a+b+c: aas[i] for i, (a, b, c) in enumerate(itertools.product(bases, repeat=3))}

def tr(s):
    return "".join(tab.get(s[i:i+3], "X") for i in range(0, len(s) - len(s) % 3, 3))

rs = ref[21562:25384]
rp = tr(rs).rstrip("*")

seed = rs[:20]
start = smp.find(seed)
print("seed start (0-based) in sample:", start)

end = smp.find(rs[-20:], start) + 20
print("end:", end, "len", end - start)

sp = tr(smp[start:end])
print("sample spike aa len", len(sp), "ends with stop:", sp.endswith("*"))
sp = sp.rstrip("*")

sm = difflib.SequenceMatcher(None, rp, sp, autojunk=False)
out = []
for tag, i1, i2, j1, j2 in sm.get_opcodes():
    if tag == "equal":
        continue
    if tag == "replace" and (i2 - i1) == (j2 - j1):
        for k in range(i2 - i1):
            out.append(f"{rp[i1+k]}{i1+k+1}{sp[j1+k]}")
    elif tag == "delete":
        out.append(f"del{i1+1}" if i2 - i1 == 1 else f"del{i1+1}-{i2}")
    elif tag == "insert":
        out.append(f"ins{i1}_{i1+1}{sp[j1:j2]}")
    else:
        out.append(f"{tag}:{rp[i1:i2]}{i1+1}-{i2}->{sp[j1:j2]}")

print(len(out))
print(out)

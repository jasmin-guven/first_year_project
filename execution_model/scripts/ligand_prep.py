import BioSimSpace as bss

print (f"{sys.argv[0]} {sys.argv[1]}")
idx = int( sys.argv[1] )
stream = open("ligands.dat","r")
lines = stream.readlines()
lig_name = lines[idx].rstrip()
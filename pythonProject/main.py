import requests
from bs4 import BeautifulSoup
from PIL import Image, ImageDraw
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP

# 636x637

aminos_posibles = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
colores_posibles = [(128,0,0), (128,128,0), (0,128,0), (128,0,128), (0,128,128), (0,0,128)]

proteinas_MW = []
proteinas_PI = []
input_proteinas = []
while True:
    entrada = input("Ingresa entrada UniProt de una proteina ('X' para terminar):")
    if entrada == "X":
        break
    else:
        input_proteinas.append(entrada)

for proteina in input_proteinas:
    uniprot_entry = proteina

    r = requests.get('https://rest.uniprot.org/uniprotkb/{0}.txt'.format(uniprot_entry))

    soup = BeautifulSoup(r.content, 'html.parser').prettify()

    datos_uniprot = soup.split("SQ   SEQUENCE")
    datos_lista = datos_uniprot[-1].split(";")

    peso_molecular = datos_lista[1].lstrip()
    peso_molecular = peso_molecular.rstrip(" MW")
    peso_molecular = float(peso_molecular)/1000
    proteinas_MW.append(peso_molecular)

    seq_amino = datos_lista[3].lstrip()
    seq_amino = seq_amino.rstrip()
    seq_amino = seq_amino.replace("\n", " ")
    seq_amino = seq_amino.replace("/", " ")
    seq_amino = seq_amino.replace(" ", "")
    lista_seq_amino = list(seq_amino)

    # calcular punto isoelectrico
    protein = IP(seq_amino)
    punto_iso = protein.pi()
    proteinas_PI.append(punto_iso)

#####
img = Image.open('grid2.PNG')
draw = ImageDraw.Draw(img)

# RANGO ISO 0-14
# RANGO MW kDa 0-20
# 604x321
# ISO x MW
for x in range(len(input_proteinas)):
    punto_iso_grafica = (604*proteinas_PI[x])/14
    peso_molecular_grafica = (321*proteinas_MW[x])/20

    draw.ellipse(xy=(punto_iso_grafica-5, peso_molecular_grafica-5,
                     punto_iso_grafica+5, peso_molecular_grafica+5),
                fill=colores_posibles[x],
                width=5)
    draw.text((punto_iso_grafica, peso_molecular_grafica-20),
              input_proteinas[x], fill=colores_posibles[x])

img.show()

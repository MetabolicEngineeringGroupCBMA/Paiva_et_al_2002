#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: py:percent,ipynb
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python (pydna312_bp185)
#     language: python
#     name: pydna312_bp185
# ---

# %% [markdown]
# ![image.png](attachment:359a6097-a0a9-45aa-abd0-125156555dbf.png)
#
# [link to paper](https://portlandpress.com/biochemj/article-abstract/363/3/737/40440/Utilization-of-green-fluorescent-protein-as-a?redirectedFrom=fulltext)

# %%
# if pydna_utils is not istalled, change the import below to:
# from pydna.genbank import genbank 
from pydna_utils.genbank import genbank

# %% [markdown]
# # Figure 1 (page 739) 
# Strategies followed for the construction of the JEN1–GFP fusion (A) and for the verification of gene fusion by analytical PCR (B)
#
# ![xxx.png](attachment:cb4cbe49-62dd-438e-ba74-98d478664ea9.png)

# %%
import math
import re
import regex
from pydna.parsers import parse_primers
from pydna.amplify import pcr
from pydna.assembly2 import homologous_recombination_integration

# %% [markdown]
# # Table 2 (page 738)
#
# ![image.png](attachment:585ccbcb-91a2-4ef2-b1be-2fdf03e96e1b.png)

# %% [markdown]
# Primers below are from Table 2 above. [Free online OCR](https://www.newocr.com/) was used to extract the primer sequences, since the table seems to be an image in the PDF document.
#
# The underlined part of the S1 primer above is different from the one recommended in Wach et al. 1997 (The source of the (pFA6a-GFPS65T-kanMX6 plasmid). 
#
#     GAT-GCT-TAT-TGT-TTT-GCC-AAA-AAT-GAT-CTT-TTT-TTC-AAT-AAC-
#                                                             -GGTCGACGGATCCCCGGG
#
#     CT-GAT-GTA-CAG-AGA-TAT-ATT-ATA-CAT-TAC-AAT-TGT-ACA-AAC-
#                                                            -ATCGATGAATTCGAGCTCG
#
#
# This change removed a linker formed by the MCS between the two proteins:
#
#     Upstream protein - GRRIPGLIN - SKGEELF...(GFP)...GMDELYK
#                        --linker--
#                        
# ![image.png](attachment:32e69498-48b9-4ea3-b046-c0ea16cf8b12.png)
#
# ***
#
# The reverse primer sits slightly upstream of the recommended primer. This avoids adding the PmeI restriction site into the PCR product.
#
# ![image.png](attachment:79e9b0f6-9a27-456c-ba03-431b615a2442.png)
#
#
#     ATT GAT TCG AAC GTC TCA AAG ACA TAT GAG GAG CAT ATT GAG ACC GTT TAA    <== JEN1 open reading frame
#         GAT TCG AAC GTC TCA AAG ACA TAT GAG GAG CAT ATT GAG ACC GTT-
#                                                                     |
#                                                                      -AGTAAAGGAGAAGAACTTTTC    s1 primer
#
#

# %%
s1, s2, k2, k3, a1, a2 = parse_primers("""
>S1
GAT TCG AAC GTC TCA AAG ACA TAT GAG GAG CAT ATT GAG ACC GTT  AGTAAAGGAGAAGAACTTTTC
>S2
GTTACATAGAGAAGCGAACACGCCCTAGAGAGCAATGAAAAGTGA  GGATGGCGGCGTTAGTATC
>K2
CGATAGATTGTCGCACCTG
>K3
CCATCCTATGGAACTGCCTC
>Al
GGCCTATCCAAGGATGCTGTC
>A2
GGCCCATTCAGTGCAAGAACC
""")

# %%
# import os
# os.environ["pydna_email"] = "bjornjobb@gmail.com"

# %%
template = genbank("AJ002682.1")
assert len(template) == 4882
assert template.seguid() == 'ldseguid=zH5VDEw2UM9jjLVdEco_oZxzwss'

# %% [markdown]
# The [pFA6a-GFPS65T-kanMX6](https://www.ncbi.nlm.nih.gov/nuccore/AJ002682) plasmid has the Genbank accession number `AJ002682.1`. 
#
# The Genbank entry has only this feature key/value related to GFP:
#
#      source          63..776
#                      /organism="Aequorea victoria"
#                      /mol_type="other DNA"
#                      /db_xref="taxon:6100"
#                      /note="GFP gene (green fluorescent protein)"
#
# The slice 62 .. 776 contain the protein below, a good candidate for GFP. We use 62 instead of 63 
# since Python sequences are 0-indexed, while Genbank files start at 1. 

# %%
GFP_aa = template[62:776].seq.translate()
print(GFP_aa)

# %% [markdown]
# The protein sequence above was analyzed by protein [BLAST](https://www.uniprot.org/blast) at uniprot:
#
# ![image.png](attachment:4bc6ee58-4af5-451a-a4ed-2d427c111529.png)
#
# The protein sequence is identified as "GFP", which is what we want.
#
# The cassette is amplified using primers s1 and s2 as described in the "A" panel above. 

# %%
cassette = pcr(s1, s2, template)

# %% [markdown]
# JEN1 sequence [SGD](https://www.yeastgenome.org/locus/JEN1)

# %%
locus = genbank("NC_001143.9", 22234-1000, 24084+1000)
assert locus.seguid() == "ldseguid=JztP6Y8QAjTIWOHIWS_dFz1JVHE"

# %%
locus

# %% [markdown]
# The locus contain JEN1 and 1000 bp up & downstream.

# %%
JEN1_aa = locus[1000:-1000].seq.translate()

# %%
JEN1_aa

# %%
len(JEN1_aa)

# %% [markdown]
# The protein sequence above was analyzed by protein [BLAST](https://www.uniprot.org/blast) at uniprot:
#
# ![image.png](attachment:4c06b6cd-7e91-4eb3-b22a-be8c91e17846.png)

# %%
JEN1_GFP_fusion_protein = (JEN1_aa+GFP_aa).replace("*", "")

# %%
JEN1_GFP_fusion_protein

# %%
JEN1_GFP_locus, *_ = homologous_recombination_integration(locus, [cassette])
assert JEN1_GFP_locus.seguid() == 'ldseguid=4fA-Rwp79zn53w-mlVBM89WodVg'

# %%
JEN1_GFP_locus.name = "JEN1_GFP_locus"
JEN1_GFP_locus.stamp()

# %%
JEN1_GFP_locus.comment("""\
The S. cerevisiae JEN1 orf c-terminally tagged with GFP from pFA6a-GFPS65T-KanMX6.
""")

# %%
print(JEN1_GFP_locus.format())

# %% [markdown]
# Diagnostic PCR from the B panel above. The published results are given as comments:

# %%
pcr(a1, a2, locus) # 1000 bp

# %%
pcr(a1, k2, JEN1_GFP_locus) # 2127 bp

# %%
pcr(k3, a2, JEN1_GFP_locus) # 770 bp

# %%
pcr(a1, a2, JEN1_GFP_locus) # 3040 bp

# %%
seq = JEN1_GFP_locus.seq

# %%
frame0 = seq + "n" * (-len(seq) % 3) 
frame1 = "n" + frame0[:-1]
frame2 = "n" + frame1[:-1]

# %%
assert len(frame0) == len(frame1) == len(frame2)

# %%
translated_frames = tuple(str(f.translate()) for f in (frame0, frame1, frame2))

# %%
singleframelength = len(translated_frames[0])

# %%
target = "msssi.+delyk".upper()

# %%
mobj = re.search(target, "".join(translated_frames))

# %%
mobj

# %%
frame = math.ceil(mobj.start()/singleframelength) - 1
frame

# %%
start = mobj.start() * 3 - (3 * frame * singleframelength) - frame
start

# %%
end = start + (mobj.end() - mobj.start())*3 + 3

# %%
end

# %%
seq[start:end]

# %%
extracted_fusion_protein = seq[start:end].translate(stop_symbol="")

# %%
JEN1_GFP_fusion_protein

# %%
extracted_fusion_protein

# %%
assert extracted_fusion_protein == JEN1_GFP_fusion_protein

# %%
from pydna.align import align

# %%
aln, editlist = align(JEN1_GFP_fusion_protein, extracted_fusion_protein)

# %%
aln.counts()

# %%
print(aln)

# %%
assert JEN1_GFP_fusion_protein.seguid() == 'lsseguid=XPZS5iMC990Pyp1k9GjdCOmecbM'

# %%
UTR5 = seq[start-150:start+3]

# %%
list(re.finditer("ATG", str(UTR5)))

# %%
(len(UTR5) - UTR5.find("ATG"))%3 # 0 = in frame, 1 or 2 = out of frame

# %%
print(UTR5[12:])

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
#       jupytext_version: 1.18.1
#   kernelspec:
#     display_name: Python (pydna312_bp185)
#     language: python
#     name: pydna312_bp185
# ---

# %% [markdown]
# ![image.png](attachment:359a6097-a0a9-45aa-abd0-125156555dbf.png)

# %% [markdown]
# [link to paper](https://portlandpress.com/biochemj/article-abstract/363/3/737/40440/Utilization-of-green-fluorescent-protein-as-a?redirectedFrom=fulltext)

# %% [markdown]
# ![image.png](attachment:585ccbcb-91a2-4ef2-b1be-2fdf03e96e1b.png)

# %% [markdown]
# ![xxx.png](attachment:cb4cbe49-62dd-438e-ba74-98d478664ea9.png)

# %%
from pydna_utils.genbank import genbank
from pydna.parsers import parse_primers
from pydna.amplify import pcr
from pydna.assembly2 import homologous_recombination_integration

# %%
s1, s2, k2, k3, a1, a2 = parse_primers("""
>S1
GATTCGAACGTCTCAAAGACATATGAGGAGCATATTGAGACCGTTAGTAAAGGAGAAGAACTTTTC
>S2
GTTACATAGAGAAGCGAACACGCCCTAGAGAGCAATGAAAAGTGAGGATGGCGGCGTTAGTATC
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
template = genbank("AJ002682.1")
assert template.seguid() == 'ldseguid=zH5VDEw2UM9jjLVdEco_oZxzwss'

# %%
cassette = pcr(s1, s2, template)

# %%
locus = genbank("BK006944.2", 22234, 26159)
assert locus.seguid() == 'ldseguid=ZlCe_azWtYuHbKGWof0KX1LBlFs'

# %%
JEN1_GFP_locus, = homologous_recombination_integration(locus, [cassette])
assert JEN1_GFP_locus.seguid() == 'ldseguid=yNNY448oU9oK0PGj3lAJMMQgtLs'

# %%
JEN1_GFP_locus.name = "JEN1_GFP_locus"
JEN1_GFP_locus.stamp()
JEN1_GFP_locus.comment("""\
The S. cerevisiae JEN1 orf c-terminally tagged with GFP from pFA6a-GFPS65T-KanMX6.
""")

# %%
print(JEN1_GFP_locus.format())

# %%
pcr(a1, a2, locus) # 1000 bp

# %%
pcr(a1, k2, JEN1_GFP_locus) # 2127 bp

# %%
pcr(k3, a2, JEN1_GFP_locus) # 770 bp

# %%
pcr(a1, a2, JEN1_GFP_locus) # 3040 bp

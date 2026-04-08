from pydna.primer import Primer
from pydna_utils.genbank import genbank
#from pydna.amplify import pcr
#from pydna.primer_screen import reverse_primers
from pydna.design import primer_design
from pydna.readers import read
from Bio.SeqFeature import SeqFeature, FeatureLocation
template = genbank("AJ002682.1")
assert len(template) == 4882
assert template.seguid() == 'ldseguid=zH5VDEw2UM9jjLVdEco_oZxzwss'

GFP_aa = template[62:776].seq.translate()
start = 62



recommended_reverse_primer = Primer("ggatggcggcgttagtatc")
# reverse_primers(template, [recommended_reverse_primer])

end = template.find(recommended_reverse_primer.rc()) + len(recommended_reverse_primer)

template[start:end].seq

amp = primer_design(template[start:end], limit=16)

amp.figure()

def feature_span(s, start="start", end="end", where = ('locus_tag','label')):
    for f in s.features:
        if f.qualifiers["label"][0] == start:
            x = int(f.location.start)
            y = int(f.location.end)
            break
    for f in s.features:
        if f.qualifiers["label"][0] == end:
            y = int(f.location.end)
            break
    return slice(x, y)

s = read("pFA6a-GFPS65T-kanMX6.gb")

cassette_feature = SeqFeature(FeatureLocation(*feature_span(s)))

locus = genbank("NC_001143.9", 22234-1000, 24084+1000)
assert locus.seguid() == "ldseguid=JztP6Y8QAjTIWOHIWS_dFz1JVHE"
JEN1_aa = locus[1000:-1000].seq.translate()
JEN1_aa
locus.find_aminoacids(JEN1_aa)

jen1_orf_feature = SeqFeature(FeatureLocation(1000, 2851))

assert jen1_orf_feature.extract(locus).isorf()

primer_design(jen1_orf_feature)
